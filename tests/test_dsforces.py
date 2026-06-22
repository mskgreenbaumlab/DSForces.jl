import math
import unittest

from dsforces import (
    ComputeAlpha,
    ComputeDSForce,
    FindLongestComplementarySegment,
)


class DSForcesTests(unittest.TestCase):
    def test_compute_alpha(self):
        self.assertAlmostEqual(ComputeAlpha("AACCGGTT"), 0.375)
        self.assertAlmostEqual(ComputeAlpha("AACCCCAA"), 0.0)
        self.assertAlmostEqual(ComputeAlpha("GGCCCCGG"), 0.5)
        self.assertAlmostEqual(ComputeAlpha("GTTTG"), 0.48)
        self.assertAlmostEqual(ComputeAlpha("GTCTCTGA"), 0.40625)
        self.assertEqual(ComputeAlpha("AACCNGGTT"), 0.0)

    def test_find_longest_complementary_segment(self):
        self.assertEqual(FindLongestComplementarySegment("AACCGGTT"), 4)
        self.assertEqual(
            FindLongestComplementarySegment("AACCGGTT", return_coords=True),
            ((1, 4), (5, 8)),
        )

        self.assertEqual(FindLongestComplementarySegment("CCAAAGG"), 2)
        self.assertEqual(
            FindLongestComplementarySegment("CCAAAGG", return_coords=True),
            ((1, 2), (6, 7)),
        )

        self.assertEqual(FindLongestComplementarySegment("CAACCCAGAAGGGG"), 3)
        self.assertEqual(
            FindLongestComplementarySegment("CAACCCAGAAGGGG", return_coords=True),
            ((4, 6), (11, 13)),
        )

        self.assertEqual(FindLongestComplementarySegment("AAAAACCCCC"), 2)
        self.assertEqual(
            FindLongestComplementarySegment("AAAAACCCCC", return_coords=True),
            ((10, 10), (10, 10)),
        )

    def test_compute_ds_force(self):
        self.assertIsNone(ComputeDSForce("AGATCGAGACCATCCTGGCCAACATGGTGAAATN"))

        seq = "AACCGGTT"
        ds_force, seq_a, seq_b = ComputeDSForce(seq, return_lcs_positions=True)
        self.assertEqual(seq_a, (1, 4))
        self.assertEqual(seq_b, (5, 8))
        self.assertFalse(math.isnan(ds_force))

    def test_constrained_range_accepts_python_range(self):
        seq = "ACAACGTAACGGTCGAGTCG"
        ds_force, seq_a, _ = ComputeDSForce(
            seq,
            return_lcs_positions=True,
            seq_a_range=range(1, 11),
        )
        self.assertFalse(math.isnan(ds_force))
        self.assertLessEqual(1, seq_a[0])
        self.assertLessEqual(seq_a[1], 10)

    def test_sliding_matches_single_window(self):
        seq = "AACCGGTTACGT"
        sliding = ComputeDSForce(seq, sliding_window_length=8)
        sequential = [ComputeDSForce(seq[i : i + 8]) for i in range(len(seq) - 8 + 1)]
        for x, y in zip(sliding, sequential, strict=True):
            if math.isnan(x) or math.isnan(y):
                self.assertTrue(math.isnan(x) and math.isnan(y))
            else:
                self.assertAlmostEqual(x, y, places=10)


if __name__ == "__main__":
    unittest.main()
