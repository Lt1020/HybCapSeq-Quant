"""Accession-resolved reference chemistry and declared theoretical mass models."""

from __future__ import annotations

import concurrent.futures
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from .core import NA, MISSING, STANDARD_PATH, number, read_tsv, sha256, write_tsv

CATALOG_NAME = "virus_pathogen_database.molecular_annotations.tsv"
MANIFEST_NAME = "virus_pathogen_database.molecular_annotations.manifest.json"
METADATA_NAME = "virus_pathogen_database.all_metadata.tsv"
SCHEMA_VERSION = "1.0"
NCBI_ENDPOINT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ICTV_URL = "https://ictv.global/sites/default/files/VMR/VMR_MSL41.v1.20260729.xlsx"
ICTV_VERSION = "VMR_MSL41.v1.20260729"
ICTV_SHA256 = "b79b5d82a1b3b8e9dd5e19afe8fe1a8f441267474918a7cefa8ae4913adf45bb"
MASS_SOURCE = (
    "https://www.promega.com/-/media/files/resources/technical-references/"
    "nucleotides-and-nucleic-acids.pdf"
)
MASS_FACTORS = {"ssRNA": 340.0, "dsRNA": 680.0, "ssDNA": 330.0, "dsDNA": 660.0}
MASS_MODEL = "average_reference_position_v1"
PACKAGE_REFERENCE = Path(__file__).parent / "data" / "reference"
ACCESSION_PATTERN = re.compile(r"\b[A-Z][A-Z0-9_]*\d+(?:\.\d+)?\b")


def normalized_taxon(value):
    """Return a lookup key without changing any original taxonomic annotation."""
    if value is None or str(value).strip().upper() in MISSING:
        return None
    return re.sub(r"^[a-z]__", "", str(value).strip())


def genome_chemistry(value):
    """Resolve only unambiguous RNA/DNA and single/double strand annotations."""
    text = str(value or "").replace(" ", "")
    found = set(re.findall(r"(?:ss|ds)(?:RNA|DNA)", text, flags=re.I))
    found = {v[:2].lower() + v[2:].upper() for v in found}
    return next(iter(found)) if len(found) == 1 else None


class NcbiSummaryClient:
    """Cache public summaries with a shared rate limit below three calls/second."""

    def __init__(self, cache_dir, interval=0.4):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.lock = threading.Lock()
        self.last_start = 0.0
        self.interval = max(0.4, interval)

    def fetch(self, accessions):
        key = hashlib.sha256("\n".join(accessions).encode()).hexdigest()
        path = self.cache_dir / (key + ".json")
        if path.is_file():
            return json.loads(path.read_text()), sha256(path)
        query = urllib.parse.urlencode(
            {
                "db": "nuccore",
                "id": ",".join(accessions),
                "retmode": "json",
                "retmax": 1000,
                "tool": "hybcapseq_quant",
            }
        )
        for attempt in range(5):
            with self.lock:
                pause = self.interval - (time.monotonic() - self.last_start)
                if pause > 0:
                    time.sleep(pause)
                self.last_start = time.monotonic()
            request = urllib.request.Request(
                NCBI_ENDPOINT + "?" + query,
                headers={"User-Agent": "HybCapSeq-Quant-reference-annotation"},
            )
            try:
                with urllib.request.urlopen(request, timeout=50) as response:
                    raw = response.read()
                data = json.loads(raw)
                if "result" not in data:
                    raise ValueError("NCBI summary lacks a result object")
                path.write_bytes(raw)
                return data, hashlib.sha256(raw).hexdigest()
            except (OSError, ValueError) as exc:
                if attempt == 4:
                    raise RuntimeError(
                        f"NCBI summary retrieval failed for batch {key}: {exc}"
                    ) from exc
                time.sleep(min(2**attempt, 8))
        raise AssertionError("Unreachable retry state")


def collect_ncbi_summaries(metadata, cache_dir, workers=3, batch_size=150):
    """Query every viral accession; never substitute a different sequence version."""
    accessions = sorted(
        {r["Accession"] for r in metadata if not r["Accession"].startswith("ERCC-")}
    )
    client = NcbiSummaryClient(cache_dir)
    batches = [
        accessions[i : i + batch_size] for i in range(0, len(accessions), batch_size)
    ]
    records = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(3, workers)) as pool:
        futures = {pool.submit(client.fetch, batch): batch for batch in batches}
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            batch = futures[future]
            data, digest = future.result()
            result = data["result"]
            returned = {
                str(result[uid].get("accessionversion")): result[uid]
                for uid in result.get("uids", [])
                if isinstance(result.get(uid), dict)
            }
            for accession in batch:
                record = returned.get(accession)
                if record is None:
                    records[accession] = {
                        "accession": accession,
                        "retrieval_status": "exact_accession_version_not_returned",
                        "response_sha256": digest,
                    }
                else:
                    records[accession] = {
                        "accession": accession,
                        "retrieval_status": "exact_accession_version",
                        "response_sha256": digest,
                        **{
                            key: record.get(key)
                            for key in (
                                "accessionversion",
                                "uid",
                                "taxid",
                                "slen",
                                "moltype",
                                "strand",
                                "biomol",
                                "topology",
                                "organism",
                                "title",
                                "updatedate",
                            )
                        },
                    }
            if done % 10 == 0 or done == len(batches):
                print(
                    f"NCBI summaries: {done}/{len(batches)} batches; "
                    f"{len(records)} requested accessions",
                    flush=True,
                )
    return records


def read_ictv_vmr(path):
    """Read an official VMR workbook; conflicting taxon chemistries stay unresolved."""
    try:
        import openpyxl
    except ImportError as exc:
        raise ValueError(
            "Building annotations from ICTV VMR requires: pip install '.[references]'"
        ) from exc
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = next(s for s in workbook if s.title.startswith("VMR MSL"))
    source = sheet.iter_rows(values_only=True)
    header = next(source)
    keys = [str(v or "") for v in header]
    index = defaultdict(list)
    for values in source:
        row = dict(zip(keys, values))
        if not row.get("Genome"):
            continue
        evidence = {
            "genome": str(row["Genome"]),
            "chemistry": genome_chemistry(row["Genome"]),
            "isolate_id": str(row.get("Isolate ID") or "NA"),
        }
        for rank in ("Species", "Genus", "Family", "Order"):
            if row.get(rank):
                index[(rank.lower(), str(row[rank]))].append(evidence)
        for accession in ACCESSION_PATTERN.findall(
            str(row.get("Virus GENBANK accession") or "")
        ):
            index[("accession", accession)].append(evidence)
    workbook.close()
    return index


def ictv_assignment(row, index):
    accession = row["Accession"]
    candidates = [("accession", accession), ("accession", accession.split(".")[0])]
    candidates.extend(
        (rank, normalized_taxon(row.get(rank)))
        for rank in ("species", "genus", "family", "order")
    )
    for rank, name in candidates:
        if name is None or str(name).startswith("unclassified_"):
            continue
        evidence = index.get((rank, name), [])
        if not evidence:
            continue
        chemistries = {e["chemistry"] for e in evidence}
        if len(chemistries) == 1 and None not in chemistries:
            return {
                "chemistry": next(iter(chemistries)),
                "genome": ";".join(sorted({e["genome"] for e in evidence})),
                "rank": rank,
                "taxon": name,
                "isolate_ids": (
                    ",".join(sorted({e["isolate_id"] for e in evidence}))
                    if len(evidence) <= 8
                    else None
                ),
                "supporting_records": len(evidence),
                "evidence_sha256": hashlib.sha256(
                    json.dumps(evidence, sort_keys=True).encode()
                ).hexdigest(),
            }
        return {
            "chemistry": None,
            "genome": ";".join(sorted({e["genome"] for e in evidence})),
            "rank": rank,
            "taxon": name,
            "isolate_ids": (
                ",".join(sorted({e["isolate_id"] for e in evidence}))
                if len(evidence) <= 8
                else None
            ),
            "supporting_records": len(evidence),
            "evidence_sha256": hashlib.sha256(
                json.dumps(evidence, sort_keys=True).encode()
            ).hexdigest(),
        }
    return {}


def annotate_reference(row, record, ictv, commercial_ercc):
    accession = row["Accession"]
    length = float(row["Length"])
    assigned = ictv_assignment(row, ictv)
    ncbi_mol = str(record.get("moltype") or "").upper()
    ncbi_strand = str(record.get("strand") or "").lower()
    ncbi_type = (
        ncbi_strand + ncbi_mol
        if ncbi_mol in {"RNA", "DNA"} and ncbi_strand in {"ss", "ds"}
        else None
    )
    molecular_type = None
    status = "unresolved_molecular_type"
    evidence = None
    notes = []
    is_ercc = accession.startswith("ERCC-")
    exact = record.get("retrieval_status") == "exact_accession_version"
    length_matches = exact and int(record.get("slen") or -1) == int(length)
    if is_ercc:
        molecular_type, status = "ssRNA", "assigned_ERCC_transcript_reference"
        evidence = "ERCC_RNA_transcript_reference_definition"
        notes.append("Mass model describes an RNA transcript, not a DNA plasmid")
        if accession not in commercial_ercc:
            notes.append("Not in commercial Mix 1; no concentration assignment")
    elif not exact:
        status = "exact_accession_version_unavailable"
    elif not length_matches:
        status = "ncbi_database_length_mismatch"
    elif assigned.get("chemistry"):
        molecular_type, status = assigned["chemistry"], "resolved_ictv_genome"
        evidence = "ICTV_genome_type_with_NCBI_accession_and_length_verification"
        if ncbi_mol in {"RNA", "DNA"} and not molecular_type.endswith(ncbi_mol):
            status = "resolved_ictv_genome_record_chemistry_differs"
            notes.append("Genome-equivalent mass; deposited record chemistry differs")
        elif ncbi_type and ncbi_type != molecular_type:
            status = "resolved_ictv_genome_record_strandedness_differs"
            notes.append(
                "Genome-equivalent mass; deposited record strandedness differs"
            )
    elif ncbi_type:
        molecular_type, status = ncbi_type, "resolved_ncbi_accession"
        evidence = "NCBI_accession_moltype_and_strand"
    if molecular_type == "dsDNA" and "-RT" in str(assigned.get("genome", "")):
        molecular_type = "partially_dsDNA"
        status = "partially_double_stranded_DNA_requires_material_assignment"
        notes.append("Unspecified strand completion prevents a unique genome mass")
    factor = MASS_FACTORS.get(molecular_type)
    mass = factor * length if factor is not None else None
    return {
        "accession": accession,
        "assembly": row.get("Assembly", "NA"),
        "reference_name": row.get("Name", "NA"),
        "length_nt": int(length),
        "assembly_length_nt": int(float(row["Asm_length"])),
        "family_raw": row.get("family", "NA"),
        "family_normalized": normalized_taxon(row.get("family")),
        "ncbi_taxid": record.get("taxid"),
        "ncbi_accession_version": record.get("accessionversion"),
        "ncbi_length_nt": record.get("slen"),
        "ncbi_molecule_type": record.get("moltype"),
        "ncbi_strandedness": record.get("strand"),
        "ncbi_biomol": record.get("biomol"),
        "ncbi_record_update": record.get("updatedate"),
        "ncbi_retrieval_status": record.get("retrieval_status"),
        "ncbi_source_url": (
            f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}" if not is_ercc else None
        ),
        "ncbi_response_sha256": record.get("response_sha256"),
        "ictv_genome_type": assigned.get("genome"),
        "ictv_match_rank": assigned.get("rank"),
        "ictv_match_taxon": assigned.get("taxon"),
        "ictv_isolate_ids": assigned.get("isolate_ids"),
        "ictv_supporting_record_count": assigned.get("supporting_records"),
        "ictv_evidence_sha256": assigned.get("evidence_sha256"),
        "ictv_version": ICTV_VERSION if assigned else None,
        "ictv_source_url": ICTV_URL if assigned else None,
        "genome_nucleic_acid_type": (
            assigned["chemistry"][-3:]
            if assigned.get("chemistry")
            else (
                ncbi_mol if ncbi_mol in {"RNA", "DNA"} else ("RNA" if is_ercc else None)
            )
        ),
        "molecule_type": molecular_type,
        "strandedness": (
            "partially_ds"
            if molecular_type == "partially_dsDNA"
            else molecular_type[:2] if molecular_type else None
        ),
        "molecular_mass_per_nt_g_mol": factor,
        "reference_molecular_mass_g_mol": mass,
        "reference_mass_g_per_molecule": mass / NA if mass else None,
        "reference_mass_ng_per_molecule": mass / NA * 1e9 if mass else None,
        "mass_model": MASS_MODEL if mass else None,
        "mass_model_source": MASS_SOURCE if mass else None,
        "classification_evidence": evidence,
        "annotation_status": status,
        "notes": "; ".join(notes) or None,
    }


def build_catalog(
    metadata_path, vmr_path, cache_dir, output_dir, workers=3, ncbi_snapshot=None
):
    """Create a complete accession sidecar and a checksum-bound provenance manifest."""
    if sha256(vmr_path) != ICTV_SHA256:
        raise ValueError(f"ICTV workbook checksum does not match {ICTV_VERSION}")
    _, metadata = read_tsv(metadata_path)
    if len({r["Accession"] for r in metadata}) != len(metadata):
        raise ValueError("Database metadata contains duplicate accessions")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    catalog = output_dir / CATALOG_NAME
    if catalog.exists():
        raise ValueError(f"Refusing to overwrite reference annotations: {catalog}")
    if ncbi_snapshot:
        with gzip.open(ncbi_snapshot, "rt", encoding="utf-8") as source:
            snapshot = [json.loads(line) for line in source]
        records = {row["accession"]: row for row in snapshot}
        expected = {
            r["Accession"] for r in metadata if not r["Accession"].startswith("ERCC-")
        }
        if len(records) != len(snapshot) or set(records) != expected:
            raise ValueError(
                "NCBI snapshot accession set does not match database metadata"
            )
    else:
        records = collect_ncbi_summaries(metadata, cache_dir, workers=workers)
    index = read_ictv_vmr(vmr_path)
    _, ercc = read_tsv(STANDARD_PATH)
    commercial = {r["ercc_id"] for r in ercc}
    rows = [
        annotate_reference(r, records.get(r["Accession"], {}), index, commercial)
        for r in metadata
    ]
    assemblies = defaultdict(list)
    for row in rows:
        assemblies[row["assembly"]].append(row)
    for group in assemblies.values():
        lengths = {r["assembly_length_nt"] for r in group}
        factors = {r["molecular_mass_per_nt_g_mol"] for r in group}
        group_mass = None
        if len(lengths) != 1 or sum(r["length_nt"] for r in group) != next(
            iter(lengths)
        ):
            group_status = "inconsistent_assembly_reference_length"
        elif None in factors:
            group_status = "unresolved_assembly_member_mass"
        elif len(factors) != 1:
            group_status = "inconsistent_assembly_molecular_chemistry"
        else:
            group_mass = sum(r["reference_molecular_mass_g_mol"] for r in group)
            group_status = "resolved_complete_assembly"
        for row in group:
            row["assembly_reference_count"] = len(group)
            row["assembly_molecular_mass_g_mol"] = group_mass
            row["assembly_mass_status"] = group_status
    write_tsv(catalog, rows)
    source = output_dir / "ncbi_accession_summary.jsonl.gz"
    with gzip.GzipFile(filename=str(source), mode="wb", mtime=0) as binary:
        dest = io.TextIOWrapper(binary, encoding="utf-8")
        for accession in sorted(records):
            dest.write(json.dumps(records[accession], sort_keys=True) + "\n")
        dest.flush()
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "catalog_filename": CATALOG_NAME,
        "catalog_sha256": sha256(catalog),
        "database_metadata_filename": METADATA_NAME,
        "database_metadata_sha256": sha256(metadata_path),
        "reference_count": len(rows),
        "resolved_mass_count": sum(
            r["reference_molecular_mass_g_mol"] is not None for r in rows
        ),
        "annotation_status_counts": dict(Counter(r["annotation_status"] for r in rows)),
        "molecule_type_counts": dict(Counter(r["molecule_type"] or "NA" for r in rows)),
        "ncbi_summary_endpoint": NCBI_ENDPOINT,
        "ncbi_summary_filename": source.name,
        "ncbi_summary_sha256": sha256(source),
        "ncbi_snapshot_input_sha256": sha256(ncbi_snapshot) if ncbi_snapshot else None,
        "ictv_vmr_url": ICTV_URL,
        "ictv_vmr_version": ICTV_VERSION,
        "ictv_vmr_sha256": sha256(vmr_path),
        "classification_policy": (
            "ICTV genome chemistry has priority over deposited NCBI record chemistry; "
            "both annotations are retained. Mass represents one viral genome equivalent, "
            "not an inferred physical chemical form in an individual sample."
        ),
        "mass_model": MASS_MODEL,
        "mass_factors_g_mol_per_reference_position": MASS_FACTORS,
        "mass_interpretation": (
            "Calculated approximate molar mass of the annotated reference molecule; "
            "not an experimentally measured molecular mass or a virion mass. "
            "Reference length represents one strand; duplex factors include both strands."
        ),
        "avogadro_constant_mol_inverse": NA,
    }
    (output_dir / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def load_catalog(path, metadata_path=None):
    """Validate chemistry, mass arithmetic, exact identifiers, and provenance."""
    path = Path(path)
    manifest_path = path.with_name(MANIFEST_NAME)
    if not manifest_path.is_file():
        raise ValueError(f"Molecular annotation manifest is required: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported molecular annotation schema")
    if sha256(path) != manifest.get("catalog_sha256"):
        raise ValueError("Molecular annotation catalog checksum mismatch")
    _, rows = read_tsv(path)
    records = {}
    for row in rows:
        acc = row["accession"]
        if acc in records:
            raise ValueError(f"Duplicate molecular annotation accession: {acc}")
        length = number(row.get("length_nt"), positive=True)
        factor = number(row.get("molecular_mass_per_nt_g_mol"), positive=True)
        mass = number(row.get("reference_molecular_mass_g_mol"), positive=True)
        kind = row.get("molecule_type")
        if factor is not None:
            if kind not in MASS_FACTORS or factor != MASS_FACTORS[kind]:
                raise ValueError(f"Inconsistent molecular chemistry or factor: {acc}")
            if mass is None or not math.isclose(mass, factor * length, rel_tol=1e-10):
                raise ValueError(f"Inconsistent reference molecular mass: {acc}")
        elif mass is not None:
            raise ValueError(
                f"Reference mass supplied without a molecular factor: {acc}"
            )
        if mass is not None:
            for key, expected in (
                ("reference_mass_g_per_molecule", mass / NA),
                ("reference_mass_ng_per_molecule", mass / NA * 1e9),
            ):
                value = number(row.get(key), positive=True)
                if value is None or not math.isclose(value, expected, rel_tol=1e-10):
                    raise ValueError(f"Inconsistent {key}: {acc}")
        records[acc] = row
    if len(records) != manifest.get("reference_count"):
        raise ValueError("Molecular annotation record count mismatch")
    if metadata_path:
        metadata_path = Path(metadata_path)
        if sha256(metadata_path) != manifest.get("database_metadata_sha256"):
            raise ValueError(
                "Molecular annotations belong to a different database metadata file"
            )
        _, metadata = read_tsv(metadata_path)
        if {r["Accession"] for r in metadata} != set(records):
            raise ValueError("Database and molecular annotation accession sets differ")
        for row in metadata:
            annotation = records[row["Accession"]]
            if float(row["Length"]) != float(annotation["length_nt"]):
                raise ValueError(f"Reference length mismatch: {row['Accession']}")
            if row["Assembly"] != annotation["assembly"]:
                raise ValueError(f"Reference assembly mismatch: {row['Accession']}")
            if float(row["Asm_length"]) != float(annotation["assembly_length_nt"]):
                raise ValueError(f"Assembly length mismatch: {row['Accession']}")
    return records, {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": sha256(manifest_path),
        "reference_count": len(records),
        "mass_model": manifest.get("mass_model"),
        "classification_policy": manifest.get("classification_policy"),
        "database_metadata_sha256": manifest.get("database_metadata_sha256"),
    }


def database_catalog(directory):
    """Require a checksum-matched sidecar in the selected database directory."""
    directory = Path(directory)
    return load_catalog(directory / CATALOG_NAME, directory / METADATA_NAME)


def install_catalog(directory, replace=False):
    """Install bundled annotations without changing alignment database files."""
    import shutil
    import tempfile

    directory = Path(directory)
    metadata = directory / METADATA_NAME
    if not metadata.is_file():
        raise ValueError(f"Database metadata is required: {metadata}")
    manifest = json.loads((PACKAGE_REFERENCE / MANIFEST_NAME).read_text())
    if sha256(metadata) != manifest["database_metadata_sha256"]:
        raise ValueError(
            "Bundled annotations require the documented v3.3 metadata checksum. "
            "Build a new catalog for a different database."
        )
    destinations = [directory / CATALOG_NAME, directory / MANIFEST_NAME]
    if any(p.exists() for p in destinations) and not replace:
        if all(p.is_file() for p in destinations):
            records, context = database_catalog(directory)
            if context["sha256"] == manifest["catalog_sha256"]:
                return records, context
        raise ValueError(
            "Existing annotations differ; review them before using --replace"
        )
    with tempfile.TemporaryDirectory(
        prefix=".molecular_annotations_", dir=directory
    ) as tmp:
        staging = Path(tmp)
        with gzip.open(PACKAGE_REFERENCE / (CATALOG_NAME + ".gz"), "rb") as source:
            with (staging / CATALOG_NAME).open("wb") as dest:
                shutil.copyfileobj(source, dest)
        shutil.copyfile(PACKAGE_REFERENCE / MANIFEST_NAME, staging / MANIFEST_NAME)
        load_catalog(staging / CATALOG_NAME, metadata)
        for dest in destinations:
            (staging / dest.name).replace(dest)
    return database_catalog(directory)


def annotation_summary(records, context):
    """Return auditable counts without assigning quantities to unknown references."""
    return {
        **context,
        "resolved_mass_count": sum(
            number(r.get("reference_molecular_mass_g_mol"), positive=True) is not None
            for r in records.values()
        ),
        "molecule_type_counts": dict(
            Counter(r["molecule_type"] for r in records.values())
        ),
        "annotation_status_counts": dict(
            Counter(r["annotation_status"] for r in records.values())
        ),
    }
