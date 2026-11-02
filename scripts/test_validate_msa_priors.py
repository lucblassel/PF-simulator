"""Analytic checks for the measurements used in the standalone comparison."""
import math
import unittest

from validate_msa_priors import canonical_hash, compare_groups, ks, metrics


class MeasurementTests(unittest.TestCase):
    def test_known_alignment(self):
        rows = {'T1': 'AC--', 'T2': 'ACGT', 'T3': 'ATGT', 'T4': 'ATGT'}
        result = metrics(rows, 4)
        self.assertEqual(result['aligned_length'], 4)
        self.assertEqual(result['ungapped_mean'], 3.5)
        self.assertEqual(result['gap_fraction'], 1 / 8)
        self.assertEqual(result['gapped_columns'], 1 / 2)
        self.assertEqual(result['invariant_fraction'], 3 / 4)
        self.assertEqual(result['complete_invariant_fraction'], 1 / 2)
        self.assertAlmostEqual(result['entropy_mean'], math.log(2) / 4)
        self.assertEqual(result['first_pair_distance'], 0)
        self.assertEqual(result['mean_gap_run'], 2)
        self.assertEqual(result['max_gap_run'], 2)
        self.assertEqual(result['duplicate'], 1)
        self.assertEqual(result['unique_fraction'], 3 / 4)
        self.assertEqual(result['freq_A'], 4 / 14)
        self.assertEqual(result['freq_T'], 5 / 14)
        self.assertAlmostEqual(sum(v for k, v in result.items() if k.startswith('freq_')), 1)

    def test_gap_removal_can_reveal_duplicate(self):
        result = metrics({'T1': 'A-C', 'T2': 'AC-', 'T3': 'AAC', 'T4': 'ACC'}, 4)
        self.assertEqual(result['duplicate'], 0)
        self.assertEqual(result['ungapped_duplicate'], 1)

    def test_discrete_ks(self):
        self.assertEqual(ks([0, 0, 1, 1], [0, 1, 1, 1]), 0.25)
        self.assertEqual(ks([0, 0], [1, 1]), 1)
        self.assertEqual(ks([1, 0], [0, 1]), 0)

    def test_hash_is_independent_of_record_order(self):
        self.assertEqual(canonical_hash({'T1': 'AC', 'T2': 'CC'}),
                         canonical_hash({'T2': 'CC', 'T1': 'AC'}))

    def test_control_and_constants(self):
        rows, constants = [], []
        compare_groups([{'metric': 0, 'constant': 3}] * 458,
                       [{'metric': 1, 'constant': 3}] * 458,
                       'analytic', 'control', rows, constants)
        self.assertEqual(rows[0]['ks_d'], 1)
        self.assertLess(rows[0]['p_upper'], 1e-100)
        self.assertEqual(len(constants), 1)

    def test_small_retained_ensemble_uses_safe_bound(self):
        rows, constants = [], []
        compare_groups([{'metric': 0}] * 100, [{'metric': 1}] * 100,
                       'small', 'control', rows, constants)
        self.assertEqual(rows[0]['bound'], 'union_of_one_sample_DKW')
        self.assertAlmostEqual(rows[0]['p_upper'], 4 * math.exp(-50))


if __name__ == '__main__':
    unittest.main()
