"""Coverage-aware external-label benchmark for FASC.

This script deliberately contains no SPMS-specific preprocessing. It records
the complete parameter grid so that test labels cannot be used silently to
select a favourable FASC threshold.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.cluster import Birch, HDBSCAN, KMeans, OPTICS
from sklearn.datasets import (
    fetch_openml,
    load_breast_cancer,
    load_digits,
    load_iris,
    load_wine,
)
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    homogeneity_completeness_v_measure,
)
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import (
    LabelEncoder,
    MaxAbsScaler,
    StandardScaler,
    normalize,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = REPO_ROOT / "python"
sys.path.insert(0, str(PYTHON_DIR))

from fasc_core import run_fasc  # noqa: E402


OPENML_DATASETS = {
    "optdigits": 28,
    "pendigits": 32,
    "segment": 36,
    "vehicle": 54,
    "isolet": 300,
    "mfeat-factors": 12,
    "letter": 6,
}

PILOT_DATASETS = (
    "iris",
    "wine",
    "breast-cancer",
    "digits",
)

PUBLIC_DATASETS = (
    "digits",
    "optdigits",
    "pendigits",
    "segment",
    "vehicle",
    "isolet",
    "mfeat-factors",
    "letter",
)

MAXABS_DATASETS = {
    "digits",
    "optdigits",
    "pendigits",
    "mnist-test",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--suite",
        choices=("pilot", "public"),
        default="pilot",
    )
    parser.add_argument(
        "--datasets",
        help="Comma-separated dataset names overriding --suite.",
    )
    parser.add_argument(
        "--mnist-dir",
        type=Path,
        help="Directory containing x.csv and true_labels.csv.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data_cache",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            Path(__file__).resolve().parent
            / "results"
            / "external_label_benchmark.csv"
        ),
    )
    parser.add_argument("--max-samples", type=int, default=2000)
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument(
        "--threshold-quantiles",
        default="0.01",
        help="Quantiles used for the assignment threshold.",
    )
    parser.add_argument(
        "--merge-threshold-quantiles",
        default="0.50,0.75,0.90",
        help=(
            "Merge-threshold quantiles. Assignment and merge grids are "
            "evaluated as a Cartesian product."
        ),
    )
    parser.add_argument(
        "--fasc-capacity-multipliers",
        default="1.0",
        help=(
            "Multipliers applied to the label-blind base capacity "
            "min(100, max(12, ceil(sqrt(N))))."
        ),
    )
    parser.add_argument("--dass-alphas", default="0.05")
    parser.add_argument("--kmeans-repeats", type=int, default=10)
    parser.add_argument(
        "--hdbscan-fractions",
        default="0.005,0.01,0.02",
    )
    parser.add_argument(
        "--hdbscan-min-samples-fraction",
        type=float,
        default=0.005,
        help=(
            "Label-blind min_samples fraction used for every HDBSCAN "
            "minimum-cluster-size setting."
        ),
    )
    parser.add_argument(
        "--birch-radius-quantiles",
        default="0.25,0.50,0.75",
        help=(
            "Quantiles of the unlabeled kth-neighbor Euclidean distance "
            "used as BIRCH radius thresholds."
        ),
    )
    parser.add_argument(
        "--dpmeans-penalty-quantiles",
        default="0.10,0.25,0.50",
        help=(
            "Quantiles of the unlabeled pairwise squared-Euclidean "
            "distribution used as batch DP-means penalties."
        ),
    )
    parser.add_argument("--dpmeans-max-iter", type=int, default=100)
    parser.add_argument(
        "--fasc-only",
        action="store_true",
        help="Skip all reference-method runs.",
    )
    parser.add_argument(
        "--references-only",
        action="store_true",
        help="Skip FASC and run only the reference methods.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output CSV.",
    )
    return parser.parse_args()


def comma_floats(value: str) -> list[float]:
    numbers = [float(item) for item in value.split(",") if item.strip()]
    if not numbers:
        raise ValueError("a parameter grid cannot be empty")
    return numbers


def load_dataset(
    name: str,
    cache_dir: Path,
    mnist_dir: Path | None,
) -> tuple[np.ndarray, np.ndarray]:
    if name == "iris":
        data, target = load_iris(return_X_y=True)
    elif name == "wine":
        data, target = load_wine(return_X_y=True)
    elif name == "breast-cancer":
        data, target = load_breast_cancer(return_X_y=True)
    elif name == "digits":
        data, target = load_digits(return_X_y=True)
    elif name == "mnist-test":
        if mnist_dir is None:
            raise ValueError("mnist-test requires --mnist-dir")
        data_path = mnist_dir / "x.csv"
        label_path = mnist_dir / "true_labels.csv"
        data = np.loadtxt(
            data_path,
            delimiter=",",
            skiprows=60000,
            max_rows=10000,
        )
        target = np.loadtxt(
            label_path,
            dtype=np.int64,
            skiprows=60000,
            max_rows=10000,
        )
    elif name in OPENML_DATASETS:
        fetched = fetch_openml(
            data_id=OPENML_DATASETS[name],
            as_frame=False,
            parser="liac-arff",
            data_home=cache_dir,
        )
        data, target = fetched.data, fetched.target
    else:
        raise ValueError("unknown dataset: {}".format(name))

    data = np.asarray(data, dtype=np.float64)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("{} must be an N-by-D matrix with D >= 2".format(name))
    if np.any(~np.isfinite(data)):
        raise ValueError("{} contains missing or infinite features".format(name))
    encoded = LabelEncoder().fit_transform(np.asarray(target).reshape(-1))
    if data.shape[0] != encoded.size:
        raise ValueError("{} has inconsistent feature and label rows".format(name))
    return data, encoded.astype(np.int64)


def deterministic_subsample(
    data: np.ndarray,
    target: np.ndarray,
    max_samples: int,
) -> tuple[np.ndarray, np.ndarray]:
    if data.shape[0] <= max_samples:
        return data, target
    rng = np.random.default_rng(20260723)
    selected = np.sort(
        rng.choice(data.shape[0], size=max_samples, replace=False)
    )
    return data[selected], target[selected]


def preprocess(
    data: np.ndarray,
    dataset: str,
) -> tuple[np.ndarray, np.ndarray, str]:
    if dataset in MAXABS_DATASETS:
        transformed = MaxAbsScaler().fit_transform(data)
        preprocessing = "maxabs-l2"
    else:
        transformed = StandardScaler().fit_transform(data)
        preprocessing = "standardize-l2"
    unit = normalize(transformed, norm="l2")
    return transformed, np.asarray(unit), preprocessing


def calibrated_thresholds(
    unit_data: np.ndarray,
    neighbor_count: int,
    quantiles: list[float],
) -> list[tuple[float, float]]:
    k = min(max(1, neighbor_count), unit_data.shape[0] - 1)
    model = NearestNeighbors(
        n_neighbors=k + 1,
        metric="cosine",
        algorithm="brute",
        n_jobs=1,
    ).fit(unit_data)
    distances, _ = model.kneighbors(unit_data)
    kth_similarity = 1.0 - distances[:, k]
    return [
        (quantile, float(np.quantile(kth_similarity, quantile)))
        for quantile in quantiles
    ]


def squared_euclidean_matrix(
    data: np.ndarray,
    centers: np.ndarray,
) -> np.ndarray:
    distances = (
        np.sum(data * data, axis=1)[:, None]
        + np.sum(centers * centers, axis=1)[None, :]
        - 2.0 * data @ centers.T
    )
    return np.maximum(distances, 0.0)


def canonical_farthest_index(
    data: np.ndarray,
    minimum_distances: np.ndarray,
) -> int:
    maximum = float(np.max(minimum_distances))
    tolerance = 32.0 * np.finfo(float).eps * max(1.0, abs(maximum))
    candidates = np.flatnonzero(
        np.abs(minimum_distances - maximum) <= tolerance
    )
    if candidates.size == 1:
        return int(candidates[0])
    keys = tuple(
        data[candidates, column]
        for column in range(data.shape[1] - 1, -1, -1)
    )
    return int(candidates[np.lexsort(keys)[0]])


def deterministic_batch_dp_means(
    data: np.ndarray,
    penalty: float,
    max_iter: int,
) -> tuple[np.ndarray, int, str]:
    """Run a deterministic batch DP-means reference configuration.

    New representatives are introduced at the canonically selected farthest
    observation until every observation is within the squared-distance
    penalty. Representatives are then updated from frozen assignments.
    """
    if penalty <= 0 or not np.isfinite(penalty):
        raise ValueError("the DP-means penalty must be positive and finite")
    if max_iter < 1:
        raise ValueError("dpmeans_max_iter must be positive")

    centers = np.mean(data, axis=0, keepdims=True)
    previous_labels = None
    for iteration in range(1, max_iter + 1):
        while True:
            distances = squared_euclidean_matrix(data, centers)
            labels = np.argmin(distances, axis=1)
            minimum_distances = distances[np.arange(data.shape[0]), labels]
            if float(np.max(minimum_distances)) <= penalty:
                break
            farthest = canonical_farthest_index(data, minimum_distances)
            centers = np.vstack((centers, data[farthest]))

        distances = squared_euclidean_matrix(data, centers)
        labels = np.argmin(distances, axis=1)
        active = np.unique(labels)
        remap = np.full(centers.shape[0], -1, dtype=np.int64)
        remap[active] = np.arange(active.size, dtype=np.int64)
        labels = remap[labels]
        centers = np.vstack(
            [np.mean(data[labels == index], axis=0) for index in range(active.size)]
        )
        if previous_labels is not None and np.array_equal(
            labels, previous_labels
        ):
            return labels + 1, iteration, "fixed-point"
        previous_labels = labels.copy()

    return labels + 1, max_iter, "iteration-limit"


def choose_feature_groups(dimension_count: int) -> tuple[np.ndarray, np.ndarray]:
    split = max(1, dimension_count // 2)
    first = np.arange(split, dtype=np.int64)
    second = np.arange(split, dimension_count, dtype=np.int64)
    if second.size == 0:
        second = np.array([dimension_count - 1], dtype=np.int64)
    return first, second


def unique_singleton_noise(labels: np.ndarray) -> np.ndarray:
    output = np.asarray(labels, dtype=np.int64).copy()
    rejected = np.flatnonzero(output == 0)
    first = int(output.max(initial=0)) + 1
    output[rejected] = np.arange(first, first + rejected.size)
    return output


def choose2(values: np.ndarray) -> int:
    integer = np.asarray(values, dtype=np.int64)
    return int(np.sum(integer * (integer - 1) // 2))


def pairwise_scores(
    target: np.ndarray,
    labels: np.ndarray,
) -> tuple[float, float, float]:
    assigned = labels > 0
    class_count = int(np.max(target)) + 1
    active_labels, active_index = np.unique(
        labels[assigned], return_inverse=True
    )
    contingency = np.zeros(
        (class_count, active_labels.size), dtype=np.int64
    )
    np.add.at(
        contingency,
        (target[assigned], active_index),
        1,
    )
    true_positive = choose2(contingency.ravel())
    predicted_positive = choose2(contingency.sum(axis=0))
    actual_positive = choose2(np.bincount(target, minlength=class_count))
    precision = (
        true_positive / predicted_positive if predicted_positive else 0.0
    )
    recall = true_positive / actual_positive if actual_positive else 0.0
    f1 = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return precision, recall, f1


def finite_metric(
    function,
    target: np.ndarray,
    labels: np.ndarray,
) -> float:
    if target.size < 2 or np.unique(labels).size < 1:
        return float("nan")
    return float(function(target, labels))


def evaluate(
    dataset: str,
    method: str,
    target: np.ndarray,
    labels: np.ndarray,
    runtime_seconds: float,
    parameters: dict,
    termination: str = "",
    cycle_length: int = 0,
    projection_applied: bool = False,
) -> dict:
    labels = np.asarray(labels, dtype=np.int64)
    assigned = labels > 0
    all_labels = unique_singleton_noise(labels)
    assigned_target = target[assigned]
    assigned_labels = labels[assigned]
    pair_precision, pair_recall, pair_f1 = pairwise_scores(target, labels)
    homogeneity, completeness, v_measure = (
        homogeneity_completeness_v_measure(target, all_labels)
    )
    class_coverages = [
        float(np.mean(assigned[target == class_index]))
        for class_index in np.unique(target)
    ]
    supports = np.bincount(labels[assigned])
    supports = supports[supports > 0]
    return {
        "dataset": dataset,
        "method": method,
        "sample_count": int(target.size),
        "dimension_count": int(parameters.pop("dimension_count")),
        "reference_class_count": int(np.unique(target).size),
        "cluster_count": int(np.unique(labels[assigned]).size),
        "minimum_support": int(supports.min()) if supports.size else 0,
        "coverage": float(np.mean(assigned)),
        "outlier_fraction": float(np.mean(~assigned)),
        "all_sample_ari": finite_metric(
            adjusted_rand_score, target, all_labels
        ),
        "all_sample_ami": finite_metric(
            adjusted_mutual_info_score, target, all_labels
        ),
        "assigned_only_ari": finite_metric(
            adjusted_rand_score, assigned_target, assigned_labels
        ),
        "assigned_only_ami": finite_metric(
            adjusted_mutual_info_score,
            assigned_target,
            assigned_labels,
        ),
        "pairwise_precision": pair_precision,
        "pairwise_recall": pair_recall,
        "pairwise_f1": pair_f1,
        "homogeneity": float(homogeneity),
        "completeness": float(completeness),
        "v_measure": float(v_measure),
        "macro_class_coverage": float(np.mean(class_coverages)),
        "minimum_class_coverage": float(np.min(class_coverages)),
        "runtime_seconds": runtime_seconds,
        "termination": termination,
        "cycle_length": int(cycle_length),
        "mature_projection_applied": bool(projection_applied),
        "parameters": json.dumps(parameters, sort_keys=True),
    }


def run_fasc_grid(
    dataset: str,
    data: np.ndarray,
    target: np.ndarray,
    threshold_pairs: list[tuple[float, float, float, float]],
    alphas: list[float],
    capacity_multipliers: list[float],
    max_iter: int,
    neighbor_count: int,
    preprocessing: str,
) -> list[dict]:
    sample_count, dimension_count = data.shape
    base_capacity = min(
        100, max(12, int(math.ceil(math.sqrt(sample_count))))
    )
    minimum_support = max(2, int(math.ceil(0.001 * sample_count)))
    idx_pos, idx_neg = choose_feature_groups(dimension_count)
    rows = []
    method_settings = [("fasc-sf", None)]
    method_settings.extend(("fasc-dass", alpha) for alpha in alphas)
    capacities = []
    for multiplier in capacity_multipliers:
        capacity = min(
            sample_count,
            max(2, int(math.ceil(multiplier * base_capacity))),
        )
        if capacity not in [item[1] for item in capacities]:
            capacities.append((multiplier, capacity))
    for capacity_multiplier, maximum_clusters in capacities:
        seed_count = min(8, maximum_clusters)
        for method, alpha in method_settings:
            for (
                assignment_quantile,
                assignment_threshold,
                merge_quantile,
                merge_threshold,
            ) in threshold_pairs:
                started = time.perf_counter()
                _, counts, labels, info = run_fasc(
                    data,
                    idx_pos,
                    idx_neg,
                    merge_threshold,
                    assignment_threshold,
                    seed_count,
                    maximum_clusters,
                    max_iter,
                    "SF" if alpha is None else "DASS",
                    minimum_support,
                    "cosine",
                    lambda_weight=(
                        0.0 if alpha is None else alpha / sample_count
                    ),
                    batch_size=min(4096, sample_count),
                    verbose=False,
                    return_info=True,
                )
                runtime = time.perf_counter() - started
                if counts.size and np.min(counts) < minimum_support:
                    raise AssertionError("FASC returned a provisional cluster")
                rows.append(
                    evaluate(
                        dataset,
                        method,
                        target,
                        labels,
                        runtime,
                        {
                            "dimension_count": dimension_count,
                            "assignment_threshold_quantile": (
                                assignment_quantile
                            ),
                            "merge_threshold_quantile": merge_quantile,
                            "tau_intra": assignment_threshold,
                            "tau_inter": merge_threshold,
                            "knn_k": neighbor_count,
                            "alpha": alpha,
                            "lambda": (
                                0.0
                                if alpha is None
                                else alpha / sample_count
                            ),
                            "base_capacity": base_capacity,
                            "capacity_multiplier": capacity_multiplier,
                            "kmax": maximum_clusters,
                            "seed_count": seed_count,
                            "minimum_support": minimum_support,
                            "preprocessing": preprocessing,
                            "uses_reference_k": False,
                        },
                        termination=info["terminationReason"],
                        cycle_length=info["cycleLength"],
                        projection_applied=info["matureProjectionApplied"],
                    )
                )
    return rows


def run_reference_methods(
    dataset: str,
    unit_data: np.ndarray,
    target: np.ndarray,
    kmeans_repeats: int,
    hdbscan_fractions: list[float],
    hdbscan_min_samples_fraction: float,
    birch_radius_quantiles: list[float],
    dpmeans_penalty_quantiles: list[float],
    dpmeans_max_iter: int,
    neighbor_count: int,
    preprocessing: str,
) -> list[dict]:
    sample_count, dimension_count = unit_data.shape
    reference_k = int(np.unique(target).size)
    rows = []
    for seed in range(kmeans_repeats):
        started = time.perf_counter()
        labels = (
            KMeans(
                n_clusters=reference_k,
                n_init=1,
                random_state=seed,
            )
            .fit_predict(unit_data)
            .astype(np.int64)
            + 1
        )
        rows.append(
            evaluate(
                dataset,
                "kmeans-unit-oracle-k",
                target,
                labels,
                time.perf_counter() - started,
                {
                    "dimension_count": dimension_count,
                    "reference_k": reference_k,
                    "random_seed": seed,
                    "preprocessing": preprocessing,
                    "uses_reference_k": True,
                },
            )
        )

    hdbscan_min_samples = max(
        5, int(math.ceil(hdbscan_min_samples_fraction * sample_count))
    )
    for fraction in hdbscan_fractions:
        minimum_cluster_size = max(
            2, int(math.ceil(fraction * sample_count))
        )
        started = time.perf_counter()
        raw = HDBSCAN(
            min_cluster_size=minimum_cluster_size,
            min_samples=hdbscan_min_samples,
            metric="euclidean",
            n_jobs=1,
        ).fit_predict(unit_data)
        labels = np.where(raw >= 0, raw + 1, 0).astype(np.int64)
        rows.append(
            evaluate(
                dataset,
                "hdbscan",
                target,
                labels,
                time.perf_counter() - started,
                {
                    "dimension_count": dimension_count,
                    "minimum_cluster_fraction": fraction,
                    "minimum_cluster_size": minimum_cluster_size,
                    "minimum_samples_fraction": (
                        hdbscan_min_samples_fraction
                    ),
                    "minimum_samples": hdbscan_min_samples,
                    "preprocessing": preprocessing,
                    "uses_reference_k": False,
                },
            )
        )

        started = time.perf_counter()
        raw = OPTICS(
            min_samples=max(5, int(math.ceil(0.005 * sample_count))),
            min_cluster_size=minimum_cluster_size,
            metric="euclidean",
            cluster_method="xi",
            xi=0.05,
            n_jobs=1,
        ).fit_predict(unit_data)
        labels = np.where(raw >= 0, raw + 1, 0).astype(np.int64)
        rows.append(
            evaluate(
                dataset,
                "optics-xi",
                target,
                labels,
                time.perf_counter() - started,
                {
                    "dimension_count": dimension_count,
                    "minimum_cluster_fraction": fraction,
                    "minimum_cluster_size": minimum_cluster_size,
                    "minimum_samples": max(
                        5, int(math.ceil(0.005 * sample_count))
                    ),
                    "xi": 0.05,
                    "preprocessing": preprocessing,
                    "uses_reference_k": False,
                },
            )
        )

    k = min(max(1, neighbor_count), sample_count - 1)
    neighbor_model = NearestNeighbors(
        n_neighbors=k + 1,
        metric="euclidean",
        algorithm="brute",
        n_jobs=1,
    ).fit(unit_data)
    neighbor_distances, _ = neighbor_model.kneighbors(unit_data)
    kth_distance = neighbor_distances[:, k]
    pairwise_squared = squared_euclidean_matrix(unit_data, unit_data)
    pairwise_squared = pairwise_squared[
        np.triu_indices(sample_count, k=1)
    ]
    for quantile in dpmeans_penalty_quantiles:
        penalty = max(
            float(np.quantile(pairwise_squared, quantile)),
            np.finfo(float).eps,
        )
        started = time.perf_counter()
        labels, iterations, termination = deterministic_batch_dp_means(
            unit_data,
            penalty,
            dpmeans_max_iter,
        )
        rows.append(
            evaluate(
                dataset,
                "dp-means-batch",
                target,
                labels,
                time.perf_counter() - started,
                {
                    "dimension_count": dimension_count,
                    "penalty_quantile": quantile,
                    "squared_distance_penalty": penalty,
                    "penalty_calibration": (
                        "pairwise-squared-euclidean-quantile"
                    ),
                    "maximum_iterations": dpmeans_max_iter,
                    "iterations": iterations,
                    "preprocessing": preprocessing,
                    "uses_reference_k": False,
                },
                termination=termination,
            )
        )
    for quantile in birch_radius_quantiles:
        radius = float(np.quantile(kth_distance, quantile))
        started = time.perf_counter()
        labels = (
            Birch(
                threshold=max(radius, np.finfo(float).eps),
                branching_factor=50,
                n_clusters=None,
            )
            .fit_predict(unit_data)
            .astype(np.int64)
            + 1
        )
        rows.append(
            evaluate(
                dataset,
                "birch-adaptive",
                target,
                labels,
                time.perf_counter() - started,
                {
                    "dimension_count": dimension_count,
                    "radius_quantile": quantile,
                    "radius_threshold": radius,
                    "knn_k": neighbor_count,
                    "preprocessing": preprocessing,
                    "uses_reference_k": False,
                },
            )
        )
    return rows


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.fasc_only and args.references_only:
        raise ValueError(
            "--fasc-only and --references-only cannot be used together"
        )
    assignment_quantiles = comma_floats(args.threshold_quantiles)
    merge_quantiles = (
        comma_floats(args.merge_threshold_quantiles)
        if args.merge_threshold_quantiles
        else None
    )
    alphas = comma_floats(args.dass_alphas)
    capacity_multipliers = comma_floats(args.fasc_capacity_multipliers)
    hdbscan_fractions = comma_floats(args.hdbscan_fractions)
    birch_radius_quantiles = comma_floats(
        args.birch_radius_quantiles
    )
    dpmeans_penalty_quantiles = comma_floats(
        args.dpmeans_penalty_quantiles
    )
    if any(value <= 0 for value in capacity_multipliers):
        raise ValueError("FASC capacity multipliers must be positive")
    if args.hdbscan_min_samples_fraction <= 0:
        raise ValueError("the HDBSCAN min_samples fraction must be positive")
    if args.datasets:
        dataset_names = tuple(
            item.strip() for item in args.datasets.split(",") if item.strip()
        )
    elif args.suite == "pilot":
        dataset_names = PILOT_DATASETS
    else:
        dataset_names = PUBLIC_DATASETS
    if args.mnist_dir is not None and "mnist-test" not in dataset_names:
        dataset_names = tuple(dataset_names) + ("mnist-test",)
    if args.output.exists():
        if args.overwrite:
            args.output.unlink()
        else:
            raise FileExistsError(
                "{} already exists; pass --overwrite or choose another "
                "output path".format(args.output)
            )

    for dataset_name in dataset_names:
        print("Loading {}...".format(dataset_name), flush=True)
        data, target = load_dataset(
            dataset_name, args.cache_dir, args.mnist_dir
        )
        data, target = deterministic_subsample(
            data, target, args.max_samples
        )
        transformed, unit_data, preprocessing = preprocess(
            data, dataset_name
        )
        rows = []
        if not args.references_only:
            assignment_thresholds = calibrated_thresholds(
                unit_data, args.knn_k, assignment_quantiles
            )
            if merge_quantiles is None:
                threshold_pairs = [
                    (
                        quantile,
                        threshold,
                        quantile,
                        threshold,
                    )
                    for quantile, threshold in assignment_thresholds
                ]
            else:
                merge_thresholds = calibrated_thresholds(
                    unit_data, args.knn_k, merge_quantiles
                )
                threshold_pairs = [
                    (
                        assignment_quantile,
                        assignment_threshold,
                        merge_quantile,
                        merge_threshold,
                    )
                    for assignment_quantile, assignment_threshold
                    in assignment_thresholds
                    for merge_quantile, merge_threshold in merge_thresholds
                ]
            rows.extend(
                run_fasc_grid(
                    dataset_name,
                    transformed,
                    target,
                    threshold_pairs,
                    alphas,
                    capacity_multipliers,
                    args.max_iter,
                    args.knn_k,
                    preprocessing,
                )
            )
        if not args.fasc_only:
            rows.extend(
                run_reference_methods(
                    dataset_name,
                    unit_data,
                    target,
                    args.kmeans_repeats,
                    hdbscan_fractions,
                    args.hdbscan_min_samples_fraction,
                    birch_radius_quantiles,
                    dpmeans_penalty_quantiles,
                    args.dpmeans_max_iter,
                    args.knn_k,
                    preprocessing,
                )
            )
        write_rows(args.output, rows)
        print(
            "Completed {} configurations for {}.".format(
                len(rows), dataset_name
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
