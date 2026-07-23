import unittest

import numpy as np

from run_external_benchmark import deterministic_batch_dp_means


class DeterministicBatchDPMeansTests(unittest.TestCase):
    def test_separates_two_compact_groups(self):
        data = np.array(
            [
                [0.0, 0.0],
                [0.1, 0.0],
                [5.0, 5.0],
                [5.1, 5.0],
            ]
        )

        labels, iterations, reason = deterministic_batch_dp_means(
            data, penalty=0.1, max_iter=20
        )

        self.assertEqual(reason, "fixed-point")
        self.assertLessEqual(iterations, 3)
        self.assertEqual(labels[0], labels[1])
        self.assertEqual(labels[2], labels[3])
        self.assertNotEqual(labels[0], labels[2])

    def test_partition_is_invariant_to_row_permutation(self):
        data = np.array(
            [
                [0.0, 0.0],
                [0.1, 0.0],
                [2.5, 2.5],
                [2.6, 2.5],
                [5.0, 5.0],
                [5.1, 5.0],
            ]
        )
        permutation = np.array([4, 0, 3, 5, 1, 2])

        labels, _, _ = deterministic_batch_dp_means(
            data, penalty=0.1, max_iter=20
        )
        permuted_labels, _, _ = deterministic_batch_dp_means(
            data[permutation], penalty=0.1, max_iter=20
        )
        restored = np.empty_like(permuted_labels)
        restored[permutation] = permuted_labels

        expected_pairs = labels[:, None] == labels[None, :]
        restored_pairs = restored[:, None] == restored[None, :]
        np.testing.assert_array_equal(expected_pairs, restored_pairs)


if __name__ == "__main__":
    unittest.main()
