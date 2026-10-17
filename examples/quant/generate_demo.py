"""Generate explicitly simulated evidence with analytically known RNA copies.

Signal model: aligned bases = 1e9 * input RNA mass in ng.
This is a software check, not experimental calibration or validation data.
"""

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "src/hybcapseq_quant/data"
OUT = Path(__file__).resolve().parent
AVOGADRO = 6.02214076e23


def write(path, rows, fields=None):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def create(out=OUT, sid="SIMULATED_01"):
    out.mkdir(parents=True, exist_ok=True)
    with (DATA / "ercc92_mix1.tsv").open() as f:
        standards = list(csv.DictReader(f, delimiter="\t"))
    rows = []

    def detection(acc, length, n, assembly=None, genome_length=None):
        return {
            "sample_ID": sid,
            "Name": acc,
            "description": "SIMULATED SOFTWARE TEST ONLY",
            "Length": length,
            "Segment": "NA",
            "Accession": acc,
            "Assembly": assembly or acc,
            "Asm_length": genome_length or length,
            "kingdom": "Orthornavirae",
            "phylum": "NA",
            "tclass": "NA",
            "order": "NA",
            "family": "NA",
            "genus": "NA",
            "species": "NA",
            "subspecies": "NA",
            "RPKMF": max(1, int(n / 150)) / (length / 1000) / 20,
            "read_count": max(1, int(n / 150)),
            "covered_bases": min(length, int(n)),
            "mean_coverage": n / length,
            "avg_read_identity": 1,
            "Pi": 0,
            "filtered_reads_in_sample": 20000000,
        }

    for s in standards:
        length = int(s["length_nt"])
        molecules = float(s["stock_attomol_ul"]) * 1e-18 * AVOGADRO / 1000 * 2
        mass = molecules * length * 340 / AVOGADRO * 1e9
        rows.append(detection(s["ercc_id"], length, mass * 1e9))
    # 10,000 molecules in 14uL, originating from 2uL of extract.
    total = 10000 * 10000 * 340 / AVOGADRO * 1e9 * 1e9
    rows.append(detection("SIM_RNA", 10000, total))
    rows.append(detection("SIM_SEG_A", 4000, total * 0.4, "SIM_SEGMENTED", 10000))
    rows.append(detection("SIM_SEG_B", 6000, total * 0.6, "SIM_SEGMENTED", 10000))
    write(out / (sid + ".detected_virus.info.tsv"), rows)
    sample = {
        "sample_id": sid,
        "sample_type": "wastewater",
        "batch_id": "SIMULATED",
        "replicate_group": "NA",
        "library_input_volume_ul": 2,
        "library_input_mass_ng": 120,
        "extract_nucleic_acid_concentration_ng_ul": 60,
        "ercc_dilution_factor": 1000,
        "ercc_spike_volume_ul": 2,
        "ercc_mix": "Mix1",
        "reaction_volume_ul": 14,
        "total_effective_bases": 3000000000,
        "extract_elution_volume_ul": 80,
        "extraction_input_volume_ul": 300,
        "concentrate_volume_ml": 0.460,
        "wastewater_volume_ml": 40,
        "pmmov_recovery_fraction": 0.25,
        "alignment_qc_verified": 1,
        "detection_qc_verified": 1,
        "calibration_transfer_verified": 1,
        "loq_reaction_copies_ul": 1,
    }
    write(out / "samples.tsv", [sample])
    return sample, rows


if __name__ == "__main__":
    create()
