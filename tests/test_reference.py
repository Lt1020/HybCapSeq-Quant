from pathlib import Path
import json
import tempfile
import unittest
from ercc_audit.core import InputError
from ercc_audit.reference import audit_ercc_reference
from ercc_audit.cli import main, write_tsv, STANDARD_FIELDS


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.path=Path(self.tmp.name)/'synthetic.fasta'
        self.standard=[dict(lot_id='SIMULATED',ercc_id='ERCC-SIM-001',
                            length_nt='8',stock_copies_per_ul='1000')]

    def tearDown(self):
        self.tmp.cleanup()

    def test_multiline_case_and_sequence_checksum(self):
        self.path.write_text('>ERCC-SIM-001 simulation only\nacgt\nACGU\n')
        row=audit_ercc_reference(self.path,self.standard)[0]
        self.assertEqual(row['length_nt'],8)
        self.assertEqual(row['gc_percent'],50)
        self.assertEqual(row['contains_u'],1)
        self.assertEqual(len(row['sequence_sha256']),64)

    def test_missing_and_unexpected_reference(self):
        for text in ['', '>OTHER\nACGTACGT\n']:
            self.path.write_text(text)
            with self.subTest(text=text), self.assertRaises(InputError):
                audit_ercc_reference(self.path,self.standard)

    def test_duplicate_and_empty_reference(self):
        for text in ['>ERCC-SIM-001\n', '>ERCC-SIM-001\nACGTACGT\n>ERCC-SIM-001\nACGTACGT\n']:
            self.path.write_text(text)
            with self.subTest(text=text), self.assertRaises(InputError):
                audit_ercc_reference(self.path,self.standard)

    def test_bad_characters_and_length(self):
        for seq in ['ACGTXCGT','ACGT']:
            self.path.write_text('>ERCC-SIM-001\n'+seq+'\n')
            with self.subTest(seq=seq), self.assertRaises(InputError):
                audit_ercc_reference(self.path,self.standard)

    def test_conflicting_lot_lengths(self):
        self.path.write_text('>ERCC-SIM-001\nACGTACGT\n')
        other=dict(self.standard[0],lot_id='SECOND',length_nt='9')
        with self.assertRaisesRegex(InputError,'conflicting reference lengths'):
            audit_ercc_reference(self.path,self.standard+[other])

    def test_cli_writes_audit_and_manifest(self):
        self.path.write_text('>ERCC-SIM-001\nACGTACGT\n')
        standards=Path(self.tmp.name)/'standards.tsv'
        write_tsv(standards,self.standard,STANDARD_FIELDS)
        out=Path(self.tmp.name)/'out'
        status=main(['audit-reference','--standards',str(standards),
                     '--fasta',str(self.path),'--outdir',str(out)])
        self.assertEqual(status,0)
        self.assertTrue((out/'ercc_reference_audit.tsv').is_file())
        self.assertEqual(len(json.loads((out/'run_manifest.json').read_text())['inputs']),2)


if __name__=='__main__':
    unittest.main()
