import sys
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp


PYTHON_DIR = Path(__file__).resolve().parents[1]
REPO_DIR = PYTHON_DIR.parent
sys.path.insert(0, str(PYTHON_DIR))

from fasc_core import run_fasc  # noqa: E402


def dual_fixture():
    group_a = np.array(
        [
            [1.00, 0.02, 0.01, 0.98, 0.03, 0.01],
            [0.97, 0.05, 0.02, 1.00, 0.02, 0.01],
            [0.99, 0.03, 0.02, 0.97, 0.05, 0.02],
            [1.00, 0.01, 0.03, 0.99, 0.03, 0.02],
            [0.98, 0.04, 0.01, 0.98, 0.04, 0.01],
            [0.96, 0.06, 0.02, 0.97, 0.04, 0.03],
        ]
    )
    group_b = group_a[:, [1, 0, 2, 4, 3, 5]]
    group_c = group_a[:, [2, 1, 0, 5, 4, 3]]
    return np.vstack([group_a, group_b, group_c])


def co_clustering(labels):
    labels = np.asarray(labels)
    assigned = labels > 0
    return (
        labels[:, None] == labels[None, :]
    ) & assigned[:, None] & assigned[None, :]


def run_dual(data, batch_size):
    return run_fasc(
        data,
        np.arange(3),
        np.arange(3, 6),
        0.98,
        0.75,
        3,
        3,
        30,
        "DASS",
        1,
        "dual-cosine",
        lambda_weight=0.125,
        batch_size=batch_size,
        verbose=False,
        return_info=True,
    )


def run_spms(data, batch_size):
    return run_fasc(
        data,
        np.arange(300),
        np.arange(300, 600),
        0.70,
        0.70,
        8,
        50,
        200,
        "DASS",
        2,
        "dual-cosine",
        lambda_weight=1.0,
        batch_size=batch_size,
        verbose=False,
        return_info=True,
    )


class FASCTests(unittest.TestCase):
    def test_permutation_equivariance(self):
        data = dual_fixture()
        permutation = (
            np.array(
                [7, 1, 13, 4, 10, 16, 2, 8, 14, 5, 11, 17, 3, 9, 15, 6, 12, 18]
            )
            - 1
        )
        centers_a, counts_a, labels_a, _ = run_dual(data, 4)
        centers_b, counts_b, permuted_labels, _ = run_dual(
            data[permutation], 4
        )
        labels_b = np.zeros_like(permuted_labels)
        labels_b[permutation] = permuted_labels
        np.testing.assert_array_equal(labels_a == 0, labels_b == 0)
        np.testing.assert_array_equal(
            co_clustering(labels_a), co_clustering(labels_b)
        )
        np.testing.assert_array_equal(counts_a, counts_b)
        np.testing.assert_allclose(centers_a, centers_b, atol=1e-12)

    def test_batch_size_equivalence(self):
        data = dual_fixture()
        centers_a, counts_a, labels_a, info_a = run_dual(data, 2)
        centers_b, counts_b, labels_b, info_b = run_dual(data, 1000)
        np.testing.assert_array_equal(labels_a, labels_b)
        np.testing.assert_array_equal(counts_a, counts_b)
        np.testing.assert_allclose(centers_a, centers_b, atol=1e-12)
        self.assertAlmostEqual(
            info_a["finalObjective"],
            info_b["finalObjective"],
            places=12,
        )

    def test_exact_recovery_under_threshold_separation(self):
        data = np.array(
            [
                [1.00, 0.02, 0.01, 0.01],
                [0.98, 0.03, 0.01, 0.00],
                [1.00, 0.01, 0.02, 0.00],
                [0.99, 0.02, 0.01, 0.01],
                [0.02, 1.00, 0.01, 0.01],
                [0.03, 0.98, 0.01, 0.00],
                [0.01, 1.00, 0.02, 0.00],
                [0.02, 0.99, 0.01, 0.01],
                [0.01, 0.02, 1.00, 0.01],
                [0.01, 0.03, 0.98, 0.00],
                [0.02, 0.01, 1.00, 0.00],
                [0.01, 0.02, 0.99, 0.01],
            ]
        )
        truth = np.repeat(np.arange(3), 4)
        _, counts, labels, info = run_fasc(
            data,
            np.array([0, 1]),
            np.array([2, 3]),
            0.98,
            0.95,
            1,
            6,
            20,
            "DASS",
            2,
            "cosine",
            lambda_weight=1.0,
            batch_size=5,
            verbose=False,
            return_info=True,
        )
        np.testing.assert_array_equal(
            co_clustering(labels),
            truth[:, None] == truth[None, :],
        )
        np.testing.assert_array_equal(np.sort(counts), np.array([4, 4, 4]))
        self.assertEqual(np.sum(labels == 0), 0)
        self.assertEqual(info["terminationReason"], "fixed-point")

    def test_signed_cosine_preprocessing_preserves_direction(self):
        data = np.array(
            [[1.0, -2.0], [2.0, -4.0], [-1.0, 2.0], [-2.0, 4.0]]
        )
        truth = np.array([1, 1, 2, 2])
        _, counts, labels, _ = run_fasc(
            data,
            np.array([0]),
            np.array([1]),
            0.90,
            0.90,
            1,
            2,
            10,
            "SF",
            1,
            "cosine",
            batch_size=2,
            verbose=False,
            return_info=True,
        )
        np.testing.assert_array_equal(
            co_clustering(labels),
            truth[:, None] == truth[None, :],
        )
        np.testing.assert_array_equal(np.sort(counts), np.array([2, 2]))

    def test_final_output_projects_provisional_clusters(self):
        data = np.array(
            [
                [1.00, 0.01],
                [0.99, 0.02],
                [0.01, 1.00],
                [0.02, 0.99],
                [-1.00, 0.00],
            ]
        )
        _, counts, labels, info = run_fasc(
            data,
            np.array([0]),
            np.array([1]),
            0.95,
            0.95,
            1,
            3,
            20,
            "SF",
            2,
            "cosine",
            batch_size=2,
            verbose=False,
            return_info=True,
        )
        np.testing.assert_array_equal(np.sort(counts), np.array([2, 2]))
        self.assertTrue(np.all(counts >= 2))
        self.assertEqual(labels[-1], 0)
        self.assertTrue(info["matureProjectionApplied"])

    def test_cycle_guard_waits_for_fixed_capacity_tail(self):
        data = np.array(
            [
                [-1.00, 0.00],
                [-0.99, 0.01],
                [1.00, 0.00],
                [0.99, 0.01],
            ]
        )
        _, counts, labels, info = run_fasc(
            data,
            np.array([0]),
            np.array([1]),
            0.90,
            0.90,
            1,
            2,
            20,
            "SF",
            1,
            "cosine",
            capacity_schedule=[1, 2, 1, 2],
            batch_size=2,
            verbose=False,
            return_info=True,
        )
        np.testing.assert_array_equal(np.sort(counts), np.array([2, 2]))
        self.assertEqual(np.sum(labels == 0), 0)
        self.assertEqual(info["terminationReason"], "fixed-point")
        self.assertFalse(info["cycleDetected"])
        self.assertTrue(info["capacityScheduleFixed"])
        self.assertGreaterEqual(info["convergeIter"], 5)

    def test_real_spms_permutation_equivariance(self):
        data_path = REPO_DIR / "data" / "dataMatrix.csv"
        data = np.loadtxt(
            data_path, delimiter=",", max_rows=1000, usecols=range(600)
        )
        permutation = np.arange(data.shape[0] - 1, -1, -1)
        centers_a, counts_a, labels_a, info_a = run_spms(data, 137)
        centers_b, counts_b, permuted_labels, info_b = run_spms(
            data[permutation], 1000
        )
        labels_b = np.zeros_like(permuted_labels)
        labels_b[permutation] = permuted_labels
        np.testing.assert_array_equal(labels_a == 0, labels_b == 0)
        np.testing.assert_array_equal(
            co_clustering(labels_a), co_clustering(labels_b)
        )
        np.testing.assert_array_equal(counts_a, counts_b)
        np.testing.assert_allclose(centers_a, centers_b, atol=1e-12)
        self.assertAlmostEqual(
            info_a["finalObjective"],
            info_b["finalObjective"],
            places=8,
        )
        self.assertEqual(len(counts_a), 22)
        self.assertEqual(np.sum(labels_a == 0), 437)
        self.assertTrue(np.all(counts_a >= 2))
        self.assertTrue(info_a["matureProjectionApplied"])
        self.assertEqual(info_a["terminationReason"], "limit-cycle")
        self.assertEqual(info_a["cycleLength"], 2)

    def test_sequential_maximin_promotion(self):
        data = np.array(
            [[0.0, 1.0], [1.0, 0.0], [0.99, 0.10], [0.70, 0.70]]
        )
        centers, counts, labels, _ = run_fasc(
            data,
            np.array([0]),
            np.array([1]),
            1.1,
            0.999,
            1,
            3,
            1,
            "SF",
            1,
            "cosine",
            batch_size=2,
            verbose=False,
            return_info=True,
        )
        diagonal = np.array([1.0, 1.0]) / np.sqrt(2.0)
        near_first = np.array([0.99, 0.10])
        near_first /= np.linalg.norm(near_first)
        self.assertEqual(len(counts), 3)
        self.assertEqual(np.sum(labels == 0), 1)
        self.assertLess(
            np.min(np.linalg.norm(centers - diagonal, axis=1)), 1e-12
        )
        self.assertGreater(
            np.min(np.linalg.norm(centers - near_first, axis=1)), 0.05
        )

    def test_dual_cosine_objective_uses_minimum(self):
        data = dual_fixture()
        lambda_weight = 0.125
        centers, counts, labels, info = run_dual(data, 5)
        normalized = data.copy()
        for columns in (np.arange(3), np.arange(3, 6)):
            norms = np.linalg.norm(
                normalized[:, columns], axis=1, keepdims=True
            )
            norms[norms == 0] = 1
            normalized[:, columns] /= norms
        assigned = labels > 0
        assigned_centers = centers[labels[assigned] - 1]
        positive = np.sum(
            normalized[assigned, :3] * assigned_centers[:, :3], axis=1
        )
        negative = np.sum(
            normalized[assigned, 3:] * assigned_centers[:, 3:], axis=1
        )
        expected = np.sum(np.minimum(positive, negative))
        expected += lambda_weight * np.sum(counts * (counts - 1) / 2)
        self.assertAlmostEqual(
            info["finalObjective"], expected, places=10
        )

    def test_support_callback_names_and_legacy_aliases_match(self):
        data = dual_fixture()
        positional = (
            data,
            np.arange(3),
            np.arange(3, 6),
            0.98,
            0.75,
            3,
            3,
            30,
            "DASS",
            1,
            "dual-cosine",
        )
        support = lambda n: n
        potential = lambda n: n * (n - 1.0) / 2.0
        preferred = run_fasc(
            *positional,
            lambda_weight=0.125,
            support_function=support,
            support_potential_function=potential,
            batch_size=4,
            verbose=False,
            return_info=True,
        )
        legacy = run_fasc(
            *positional,
            lambda_weight=0.125,
            density_function=support,
            density_potential_function=potential,
            batch_size=4,
            verbose=False,
            return_info=True,
        )
        np.testing.assert_allclose(preferred[0], legacy[0], atol=1e-12)
        np.testing.assert_array_equal(preferred[1], legacy[1])
        np.testing.assert_array_equal(preferred[2], legacy[2])
        self.assertAlmostEqual(
            preferred[3]["finalObjective"],
            legacy[3]["finalObjective"],
            places=12,
        )
        with self.assertRaisesRegex(ValueError, "not both"):
            run_fasc(
                *positional,
                support_function=support,
                density_function=support,
                verbose=False,
            )

    def test_all_built_in_metrics(self):
        data = np.array(
            [
                [8, 1, 1, 0, 0],
                [7, 2, 1, 0, 0],
                [1, 8, 1, 0, 0],
                [1, 7, 2, 0, 0],
                [0, 1, 8, 1, 0],
                [0, 1, 7, 2, 0],
                [0, 0, 1, 8, 1],
                [0, 0, 2, 7, 1],
                [1, 0, 0, 1, 8],
                [2, 0, 0, 1, 7],
            ],
            dtype=float,
        )
        metrics = [
            "cosine",
            "euclidean",
            "l1",
            "minimum",
            "maximum",
            "algebraic",
            "logarithmic",
            "geometric",
            "harmonic",
            "enhanced harmonic",
            "entropy",
            "weighted entropy",
            "best average",
            "fitted core",
        ]
        for metric in metrics:
            with self.subTest(metric=metric):
                centers, counts, labels, info = run_fasc(
                    data,
                    np.arange(2),
                    np.arange(2, 5),
                    1.1,
                    -1,
                    2,
                    3,
                    8,
                    "SF",
                    1,
                    metric,
                    batch_size=3,
                    verbose=False,
                    return_info=True,
                )
                self.assertTrue(np.all(np.isfinite(centers)))
                self.assertEqual(np.sum(counts), np.sum(labels > 0))
                self.assertTrue(np.isfinite(info["finalObjective"]))

    def test_custom_kernel_and_representative(self):
        data = np.array([[0, 0], [0.1, 0], [1, 1], [1, 0.9]])

        def similarity(left, right):
            left = np.asarray(left)
            right = np.asarray(right)
            return 1.0 - np.sum(
                np.abs(left[:, None, :] - right[None, :, :]), axis=2
            )

        centers, counts, labels, info = run_fasc(
            data,
            np.array([0]),
            np.array([1]),
            0.95,
            0.5,
            2,
            2,
            10,
            "SF",
            1,
            "custom",
            similarity_function=similarity,
            representative_function=lambda members: np.median(
                np.asarray(members), axis=0
            ),
            verbose=False,
            return_info=True,
        )
        np.testing.assert_array_equal(np.sort(counts), np.array([2, 2]))
        self.assertEqual(np.sum(labels > 0), 4)
        self.assertTrue(np.all(np.isfinite(centers)))
        self.assertIn(
            info["terminationReason"], {"fixed-point", "limit-cycle"}
        )

    def test_legacy_python_label_convention(self):
        data = dual_fixture()
        _, _, matlab_labels, _ = run_dual(data, 4)
        _, _, python_labels, _ = run_fasc(
            data,
            np.arange(3),
            np.arange(3, 6),
            0.98,
            0.75,
            3,
            3,
            30,
            "DASS",
            1,
            "dual-cosine",
            lambda_weight=0.125,
            batch_size=4,
            verbose=False,
            return_info=True,
            label_convention="python",
        )
        np.testing.assert_array_equal(
            python_labels,
            np.where(matlab_labels > 0, matlab_labels - 1, -1),
        )

    def test_sparse_dense_equivalence(self):
        data = dual_fixture()
        dense = run_dual(data, 4)
        sparse = run_dual(sp.csr_matrix(data), 4)
        np.testing.assert_allclose(dense[0], sparse[0], atol=1e-12)
        np.testing.assert_array_equal(dense[1], sparse[1])
        np.testing.assert_array_equal(dense[2], sparse[2])
        self.assertAlmostEqual(
            dense[3]["finalObjective"],
            sparse[3]["finalObjective"],
            places=10,
        )


if __name__ == "__main__":
    unittest.main()
