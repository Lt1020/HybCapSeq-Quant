"""Integrity checks for synthetic ERCC reference files, without sequence mapping."""
import hashlib
from pathlib import Path
from .core import InputError, numeric

IUPAC = frozenset('ACGTURYSWKMBDHVN')

def audit_ercc_reference(path, standard_rows):
    expected = {}
    standard_keys = set()
    for row in standard_rows:
        rid = row.get('ercc_id', '').strip()
        lot = row.get('lot_id', '').strip()
        if not rid.startswith('ERCC-') or not lot:
            raise InputError('standards require ERCC- identifiers and nonempty lot_id')
        key = lot, rid
        if key in standard_keys:
            raise InputError(f'duplicate standard {key}')
        standard_keys.add(key)
        length = numeric(row.get('length_nt'), f'{rid}.length_nt', positive=True)
        if length != int(length):
            raise InputError(f'{rid}: length_nt must be an integer')
        if rid in expected and expected[rid] != int(length):
            raise InputError(f'{rid}: conflicting reference lengths between lots')
        expected[rid] = int(length)
    if not expected:
        raise InputError('standards table is empty')
    sequences = {}
    current = None
    parts = []

    def flush():
        if current is None:
            return
        sequence = ''.join(parts)
        if not sequence:
            raise InputError(f'{current}: empty sequence')
        invalid = set(sequence) - IUPAC
        if invalid:
            raise InputError(f'{current}: invalid sequence characters {sorted(invalid)}')
        sequences[current] = sequence

    with Path(path).open(encoding='utf-8-sig') as stream:
        for line_no, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                flush()
                words = line[1:].split()
                if not words:
                    raise InputError(f'{path}:{line_no}: empty FASTA identifier')
                current = words[0]
                if current not in expected:
                    raise InputError(f'{path}:{line_no}: unexpected reference {current!r}')
                if current in sequences:
                    raise InputError(f'{current}: duplicate FASTA identifier')
                parts = []
            else:
                if current is None:
                    raise InputError(f'{path}:{line_no}: sequence before first FASTA header')
                parts.append(''.join(line.split()).upper())
    flush()
    missing = sorted(set(expected)-set(sequences))
    if missing:
        raise InputError(f'missing ERCC references: {missing}')
    result = []
    for rid, sequence in sequences.items():
        if len(sequence) != expected[rid]:
            raise InputError(f'{rid}: FASTA length {len(sequence)} != declared length {expected[rid]}')
        canonical = sum(sequence.count(x) for x in 'ACGTU')
        result.append(dict(ercc_id=rid, length_nt=len(sequence),
            gc_percent=100*(sequence.count('G')+sequence.count('C'))/canonical if canonical else None,
            ambiguous_nt=len(sequence)-canonical,
            contains_u=int('U' in sequence),
            sequence_sha256=hashlib.sha256(sequence.encode('ascii')).hexdigest()))
    return result

