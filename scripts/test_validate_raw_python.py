"""Analytic checks for raw aligned/unaligned measurements."""
from pathlib import Path
import tempfile
import unittest

import validate_raw_python as raw


class RawMeasurements(unittest.TestCase):
    def test_unaligned_length_composition_and_duplicates(self):
        result = raw.unaligned_metrics({'T1': 'AC', 'T2': 'ACGT', 'T3': 'ATGT', 'T4': 'ATGT'}, 4)
        self.assertEqual(result['mean_length'], 3.5)
        self.assertEqual(result['min_length'], 2)
        self.assertEqual(result['max_length'], 4)
        self.assertEqual(result['total_residues'], 14)
        self.assertEqual(result['duplicate'], 1)
        self.assertEqual(result['unique_fraction'], 0.75)
        self.assertEqual(result['freq_A'], 4 / 14)

    def test_unaligned_rejects_gaps(self):
        with self.assertRaises(AssertionError):
            raw.unaligned_metrics({'T1': 'A-', 'T2': 'CC'}, 2)

    def test_tip_matching_and_all_gap_columns_are_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            tree = Path(temp) / 'tree.nwk'
            tree.write_text('(T1:1,T2:1,(T3:1,T4:1):1);')
            task = dict(bank='u4', sample=0, tree=str(tree), model='LG+GC{0.5}', seed=42, alpha=0.5)
            aligned = dict(T1='A--', T2='C--', T3='A--', T4='G--', Node1='AAA')
            unaligned = dict(T1='A', T2='C', T3='A', T4='G', Node1='AAA')
            result = raw.measure(task, aligned, unaligned)
            self.assertEqual(result['aligned']['records'], 4)
            self.assertEqual(result['aligned']['aligned_length'], 3)
            self.assertEqual(result['aligned']['all_gap_columns'], 2)
            self.assertEqual(result['structure']['source_aligned_records'], 5)
            self.assertEqual(result['unaligned']['mean_length'], 1)
            unaligned['T1'] = 'C'
            with self.assertRaises(AssertionError):
                raw.measure(task, aligned, unaligned)


if __name__ == '__main__':
    unittest.main()
