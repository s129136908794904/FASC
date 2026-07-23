"""Authoritative Python reference implementation of FASC.

The state transition in this module mirrors ``matlab/src/FASC.m``:
deterministic maximin seeding and promotion, frozen-state blocked
assignment, connected-component merging, canonical cluster identifiers,
and collision-resistant fixed-point/limit-cycle detection after the capacity
schedule enters a certified fixed tail.
"""

import hashlib
import time

import numpy as np
import scipy.sparse as sp


_BUILT_IN_SCORES = {
    "dual-cosine",
    "cosine",
    "euclidean",
    "euclidean-distance",
    "l1",
    "l1-norm",
    "manhattan",
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
}


# ---------------------------------------------------------------------------
# Public compatibility helpers
# ---------------------------------------------------------------------------
def get_row(matrix, index):
    """Return one matrix row as a dense one-dimensional array."""
    row = matrix[index]
    if sp.issparse(row):
        row = row.toarray()
    return np.asarray(row).reshape(-1)


def get_rows(matrix, indices):
    """Return selected rows, densifying a sparse input for compatibility."""
    rows = matrix[indices]
    if sp.issparse(rows):
        return rows.toarray()
    return np.asarray(rows)


def get_mean(matrix, mask, count):
    """Return the dense mean of rows selected by a Boolean mask."""
    members = matrix[mask]
    return _row_sum(members) / count


def normalize_l2(matrix):
    """Return a row-wise L2-normalized matrix."""
    return _normalize_rows_l2(matrix)


def normalize_l1(matrix):
    """Return a row-wise absolute-L1-normalized matrix."""
    return _normalize_rows_l1(matrix)


def normalize_dual_l2(matrix, idx_pos, idx_neg):
    """Normalize positive and negative column groups independently."""
    return _normalize_dual_l2(matrix, idx_pos, idx_neg)


def compute_sim_matrix(
    left,
    right,
    algo_key,
    idx_pos=None,
    idx_neg=None,
    similarity_function=None,
):
    """Evaluate one already-preprocessed similarity block.

    This helper is retained for compatibility. ``run_fasc`` performs the
    required preprocessing before calling the same internal evaluator.
    """
    key, similarity_function = _similarity_specification(
        algo_key, similarity_function
    )
    return _similarity_matrix(
        left,
        right,
        key,
        similarity_function,
        _index_array(idx_pos),
        _index_array(idx_neg),
    )


def calculate_objective(
    data,
    centers,
    counts,
    assignments,
    valid_k,
    strategy,
    algo,
    idx_pos,
    idx_neg,
    lambda_weight=1.0,
):
    """Compatibility wrapper for the former zero-based Python API."""
    legacy_labels = np.asarray(assignments, dtype=np.int64)
    labels = np.where(legacy_labels >= 0, legacy_labels + 1, 0)
    key, similarity_function = _similarity_specification(algo, None)
    options = {
        "batch_size": 100000,
        "lambda_weight": float(lambda_weight),
        "support_potential_function": lambda n: n * (n - 1.0) / 2.0,
        "objective_function": None,
    }
    return _objective(
        data,
        np.asarray(centers)[:valid_k],
        np.asarray(counts)[:valid_k],
        labels,
        str(strategy).upper(),
        key,
        similarity_function,
        _index_array(idx_pos),
        _index_array(idx_neg),
        options,
    )


# ---------------------------------------------------------------------------
# Reference FASC entry point
# ---------------------------------------------------------------------------
def run_fasc(
    data_matrix,
    idx_pos,
    idx_neg,
    sim_inter,
    sim_inner,
    init_limit,
    max_clusters,
    max_iter,
    strategy,
    min_vol,
    algo,
    log_callback=None,
    *,
    lambda_weight=1.0,
    support_function=None,
    support_potential_function=None,
    density_function=None,
    density_potential_function=None,
    batch_size=100000,
    similarity_function=None,
    representative_function=None,
    capacity_schedule=None,
    objective_function=None,
    track_objective=False,
    verbose=True,
    return_info=False,
    label_convention="matlab",
):
    """Run the authoritative Python FASC state transition.

    The first eleven parameters match the historical Python API and the
    positional MATLAB API. Optional keyword arguments map to the MATLAB
    name-value parameters.

    Parameters
    ----------
    label_convention : {"matlab", "python"}
        ``"matlab"`` returns outlier label 0 and cluster labels 1..K.
        ``"python"`` retains the former Python convention: -1 and 0..K-1.
    return_info : bool
        When true, return ``(centers, counts, labels, info)``. Otherwise,
        retain the historical five-value return signature.
    support_function : callable, optional
        Non-decreasing support transformation ``phi(n)`` used by DASS.
        The default is the identity. ``density_function`` is retained as a
        backward-compatible alias; do not supply both names.
    support_potential_function : callable, optional
        Cumulative support reward used to rank states in a detected DASS
        cycle. The default is ``n * (n - 1) / 2``.
        ``density_potential_function`` is a backward-compatible alias.
    """
    options = _parse_options(
        lambda_weight=lambda_weight,
        support_function=support_function,
        support_potential_function=support_potential_function,
        density_function=density_function,
        density_potential_function=density_potential_function,
        batch_size=batch_size,
        similarity_function=similarity_function,
        representative_function=representative_function,
        capacity_schedule=capacity_schedule,
        objective_function=objective_function,
        track_objective=track_objective,
        verbose=verbose,
        log_callback=log_callback,
        label_convention=label_convention,
    )
    idx_pos = _index_array(idx_pos)
    idx_neg = _index_array(idx_neg)
    _validate_inputs(
        data_matrix,
        idx_pos,
        idx_neg,
        sim_inter,
        sim_inner,
        init_limit,
        max_clusters,
        max_iter,
        min_vol,
    )

    strategy_key = str(strategy).upper()
    if strategy_key not in {"SF", "SIMFIRST", "DASS", "DENSITYFIRST"}:
        raise ValueError("strategy must be 'SF' or 'DASS'")

    similarity_key, similarity_function = _similarity_specification(
        algo, options["similarity_function"]
    )
    if similarity_key == "dual-cosine":
        _validate_dual_indices(
            idx_pos, idx_neg, data_matrix.shape[1]
        )

    log = options["log"]
    log("========= FASC Start")
    wall_start = time.perf_counter()
    cpu_start = time.process_time()

    data_work = _prepare_data(
        data_matrix, similarity_key, idx_pos, idx_neg
    )
    sample_count, dimension_count = data_work.shape

    initial_capacity = _capacity_at(
        options["capacity_schedule"], 1, max_clusters
    )
    seed_count = min(init_limit, initial_capacity, sample_count)
    seed_idx = _select_seeds(
        data_work,
        seed_count,
        similarity_key,
        similarity_function,
        idx_pos,
        idx_neg,
        options["batch_size"],
    )

    center_dtype = data_work.dtype
    centers = np.zeros((seed_count, dimension_count), dtype=center_dtype)
    for cluster_index, row_index in enumerate(seed_idx):
        centers[cluster_index] = _representative(
            data_work[row_index : row_index + 1],
            similarity_key,
            similarity_function,
            options["representative_function"],
            idx_pos,
            idx_neg,
            options["batch_size"],
        )
    counts = np.zeros(seed_count, dtype=np.int64)
    labels = np.zeros(sample_count, dtype=np.int64)

    iteration = 0
    previous_labels = None
    first_seen = {}

    iter_data = []
    objective_history = []
    label_agreement_history = []
    state_hash_history = []

    cycle_detected = False
    cycle_length = 0
    cycle_remaining = 0
    cycle_best_objective = -np.inf
    cycle_best_state = None
    termination_reason = "iteration-limit"
    mature_projection_applied = False

    while iteration < max_iter or cycle_remaining > 0:
        iteration_start = time.perf_counter()
        iteration += 1
        capacity = _capacity_at(
            options["capacity_schedule"], iteration, max_clusters
        )
        capacity_schedule_fixed = _capacity_schedule_is_fixed(
            options["capacity_schedule"], iteration
        )
        log("Iteration {} start".format(iteration))

        if centers.shape[0] > capacity:
            labels[labels > capacity] = 0
            centers = centers[:capacity].copy()
            counts = counts[:capacity].copy()

        previous_centers = centers.copy()
        previous_counts = counts.copy()
        labels_before = labels.copy()

        labels, max_affinity = _assign_in_blocks(
            data_work,
            centers,
            counts,
            strategy_key,
            sim_inner,
            similarity_key,
            similarity_function,
            idx_pos,
            idx_neg,
            options,
        )

        holes = max(0, capacity - centers.shape[0])
        promotion_idx = _select_promotions(
            data_work,
            labels,
            max_affinity,
            holes,
            similarity_key,
            similarity_function,
            options["representative_function"],
            idx_pos,
            idx_neg,
            options["batch_size"],
        )
        protected = np.zeros(
            centers.shape[0] + promotion_idx.size, dtype=bool
        )
        if promotion_idx.size:
            first_new = centers.shape[0] + 1
            labels[promotion_idx] = np.arange(
                first_new, first_new + promotion_idx.size
            )
            protected[first_new - 1 :] = True

        centers, counts, labels = _consolidate(
            data_work,
            labels,
            protected,
            min_vol,
            sim_inter,
            similarity_key,
            similarity_function,
            options["representative_function"],
            idx_pos,
            idx_neg,
            options["batch_size"],
        )

        stability, label_agreement = _stability(
            previous_centers,
            previous_counts,
            centers,
            counts,
            labels_before,
            labels,
            similarity_key,
            similarity_function,
            idx_pos,
            idx_neg,
        )

        current_hash = _state_hash(labels, capacity)
        repeated_state = (
            cycle_remaining == 0
            and capacity_schedule_fixed
            and current_hash in first_seen
        )
        cycle_candidate_state = None
        if cycle_remaining > 0 or repeated_state:
            cycle_candidate_state = _project_mature_state(
                centers, counts, labels, min_vol
            )
            mature_projection_applied = (
                mature_projection_applied
                or not np.array_equal(cycle_candidate_state[2], labels)
            )
            current_objective = _objective(
                data_work,
                cycle_candidate_state[0],
                cycle_candidate_state[1],
                cycle_candidate_state[2],
                strategy_key,
                similarity_key,
                similarity_function,
                idx_pos,
                idx_neg,
                options,
            )
        elif options["track_objective"]:
            current_objective = _objective(
                data_work,
                centers,
                counts,
                labels,
                strategy_key,
                similarity_key,
                similarity_function,
                idx_pos,
                idx_neg,
                options,
            )
        else:
            current_objective = np.nan

        outlier_count = int(np.sum(labels == 0))
        iter_data.append([stability, centers.shape[0], outlier_count])
        objective_history.append(current_objective)
        label_agreement_history.append(label_agreement)
        state_hash_history.append(current_hash)

        log("\tCurrent cluster count: {}".format(centers.shape[0]))
        log("\tOutlier count: {}".format(outlier_count))
        if np.isfinite(current_objective):
            log("\tObjective: {:.10g}".format(current_objective))
        log("\tInter-iteration similarity: {:.8f}%".format(stability * 100))
        log(
            "\tCompleted in {:.2f}s".format(
                time.perf_counter() - iteration_start
            )
        )

        if cycle_remaining > 0:
            if current_objective > cycle_best_objective:
                cycle_best_objective = current_objective
                cycle_best_state = cycle_candidate_state
            cycle_remaining -= 1
            if cycle_remaining == 0:
                centers, counts, labels = cycle_best_state
                termination_reason = "limit-cycle"
                break
        else:
            next_capacity = _capacity_at(
                options["capacity_schedule"], iteration + 1, max_clusters
            )
            is_fixed_point = (
                previous_labels is not None
                and np.array_equal(labels, previous_labels)
                and next_capacity == capacity
                and capacity_schedule_fixed
            )
            if is_fixed_point:
                termination_reason = "fixed-point"
                break

            if capacity_schedule_fixed:
                if current_hash in first_seen:
                    cycle_detected = True
                    cycle_length = iteration - first_seen[current_hash]
                    cycle_best_objective = current_objective
                    cycle_best_state = cycle_candidate_state
                    cycle_remaining = cycle_length - 1
                    log(
                        "Repeated label fingerprint detected (period {}); "
                        "replaying cycle.".format(cycle_length)
                    )
                    if cycle_remaining == 0:
                        termination_reason = "limit-cycle"
                        break
                else:
                    first_seen[current_hash] = iteration

        previous_labels = labels.copy()

    labels_before_projection = labels.copy()
    centers, counts, labels = _project_mature_state(
        centers, counts, labels, min_vol
    )
    final_projection_changed = not np.array_equal(
        labels, labels_before_projection
    )
    mature_projection_applied = (
        mature_projection_applied or final_projection_changed
    )
    final_capacity = _capacity_at(
        options["capacity_schedule"], iteration, max_clusters
    )
    final_objective = _objective(
        data_work,
        centers,
        counts,
        labels,
        strategy_key,
        similarity_key,
        similarity_function,
        idx_pos,
        idx_neg,
        options,
    )
    wall_elapsed = time.perf_counter() - wall_start
    cpu_elapsed = time.process_time() - cpu_start

    info = {
        "iterData": np.asarray(iter_data, dtype=float),
        "objective": np.asarray(objective_history, dtype=float),
        "labelAgreement": np.asarray(
            label_agreement_history, dtype=float
        ),
        "stateHash": list(state_hash_history),
        "convergeIter": iteration,
        "outLiersCount": int(np.sum(labels == 0)),
        "terminationReason": termination_reason,
        "cycleDetected": bool(cycle_detected),
        "cycleLength": int(cycle_length),
        "lambda": options["lambda_weight"],
        "batchSize": options["batch_size"],
        "finalObjective": float(final_objective),
        "finalStateHash": _state_hash(labels, final_capacity),
        "finalProjectionChanged": bool(final_projection_changed),
        "matureProjectionApplied": bool(mature_projection_applied),
        "capacityScheduleFixed": bool(
            _capacity_schedule_is_fixed(
                options["capacity_schedule"], iteration
            )
        ),
        "elapsedTime": wall_elapsed,
        "cpuTime": cpu_elapsed,
        "labelConvention": options["label_convention"],
    }

    output_labels = _format_labels(
        labels, options["label_convention"]
    )
    log("")
    log("Termination: {}".format(termination_reason))
    log(
        "Total clock time: {:.2f}s; CPU time: {:.2f}s".format(
            wall_elapsed, cpu_elapsed
        )
    )
    log("Outliers count: {}".format(info["outLiersCount"]))
    log("FASC End =========")

    if return_info:
        return centers, counts, output_labels, info
    return (
        centers,
        counts,
        output_labels,
        info["iterData"],
        info["convergeIter"],
    )


# ---------------------------------------------------------------------------
# Validation and preprocessing
# ---------------------------------------------------------------------------
def _parse_options(
    *,
    lambda_weight,
    support_function,
    support_potential_function,
    density_function,
    density_potential_function,
    batch_size,
    similarity_function,
    representative_function,
    capacity_schedule,
    objective_function,
    track_objective,
    verbose,
    log_callback,
    label_convention,
):
    if not np.isscalar(lambda_weight) or not np.isfinite(lambda_weight):
        raise ValueError("lambda_weight must be one finite scalar")
    if float(lambda_weight) < 0:
        raise ValueError("lambda_weight must be nonnegative")
    if not _is_positive_integer(batch_size):
        raise ValueError("batch_size must be a positive integer")
    if support_function is not None and density_function is not None:
        raise ValueError(
            "specify support_function or density_function, not both"
        )
    if (
        support_potential_function is not None
        and density_potential_function is not None
    ):
        raise ValueError(
            "specify support_potential_function or "
            "density_potential_function, not both"
        )
    if support_function is None:
        support_function = density_function
    if support_potential_function is None:
        support_potential_function = density_potential_function
    if support_function is None:
        support_function = lambda n: n
    if support_potential_function is None:
        support_potential_function = lambda n: n * (n - 1.0) / 2.0
    callbacks = {
        "support_function": support_function,
        "support_potential_function": support_potential_function,
    }
    for name, callback in callbacks.items():
        if not callable(callback):
            raise TypeError("{} must be callable".format(name))
    for name, callback in {
        "similarity_function": similarity_function,
        "representative_function": representative_function,
        "objective_function": objective_function,
    }.items():
        if callback is not None and not callable(callback):
            raise TypeError("{} must be callable or None".format(name))
    if not isinstance(track_objective, (bool, np.bool_)):
        raise TypeError("track_objective must be Boolean")
    if not isinstance(verbose, (bool, np.bool_)):
        raise TypeError("verbose must be Boolean")
    convention = str(label_convention).lower()
    if convention not in {"matlab", "python"}:
        raise ValueError("label_convention must be 'matlab' or 'python'")

    def log(message):
        if not verbose:
            return
        if log_callback is None:
            print(message)
        else:
            log_callback(message)

    return {
        "lambda_weight": float(lambda_weight),
        "support_function": support_function,
        "support_potential_function": support_potential_function,
        "batch_size": int(batch_size),
        "similarity_function": similarity_function,
        "representative_function": representative_function,
        "capacity_schedule": capacity_schedule,
        "objective_function": objective_function,
        "track_objective": bool(track_objective),
        "verbose": bool(verbose),
        "log": log,
        "label_convention": convention,
    }


def _validate_inputs(
    data_matrix,
    idx_pos,
    idx_neg,
    sim_inter,
    sim_inner,
    init_limit,
    max_clusters,
    max_iter,
    min_vol,
):
    if not hasattr(data_matrix, "shape") or len(data_matrix.shape) != 2:
        raise TypeError("data_matrix must be a two-dimensional matrix")
    if data_matrix.shape[0] == 0 or data_matrix.shape[1] == 0:
        raise ValueError("data_matrix must be nonempty")
    values = data_matrix.data if sp.issparse(data_matrix) else np.asarray(data_matrix)
    if not np.issubdtype(values.dtype, np.number):
        raise TypeError("data_matrix must be numeric")
    if np.iscomplexobj(values) or np.any(~np.isfinite(values)):
        raise ValueError("data_matrix must be real and finite")
    for name, value in {
        "sim_inter": sim_inter,
        "sim_inner": sim_inner,
    }.items():
        if not np.isscalar(value) or not np.isfinite(value):
            raise ValueError("{} must be one finite scalar".format(name))
    for name, value in {
        "init_limit": init_limit,
        "max_clusters": max_clusters,
        "max_iter": max_iter,
        "min_vol": min_vol,
    }.items():
        if not _is_positive_integer(value):
            raise ValueError("{} must be a positive integer".format(name))
    if init_limit > max_clusters:
        raise ValueError("init_limit cannot exceed max_clusters")
    if idx_pos.ndim != 1 or idx_neg.ndim != 1:
        raise ValueError("idx_pos and idx_neg must be one-dimensional")


def _validate_dual_indices(idx_pos, idx_neg, dimension_count):
    combined = np.concatenate([idx_pos, idx_neg])
    if idx_pos.size == 0 or idx_neg.size == 0:
        raise ValueError("dual-cosine index groups must be nonempty")
    if np.any(combined < 0) or np.any(combined >= dimension_count):
        raise ValueError("dual-cosine indices are out of range")
    if np.unique(combined).size != combined.size:
        raise ValueError(
            "dual-cosine indices must be unique and disjoint"
        )


def _similarity_specification(algorithm, supplied_function):
    if callable(algorithm):
        key = "custom"
        similarity_function = algorithm
    elif isinstance(algorithm, str):
        key = algorithm.lower()
        similarity_function = supplied_function
    else:
        raise TypeError("algo must be a score name or callable")
    if key == "custom" and similarity_function is None:
        raise ValueError("a custom similarity function is required")
    if (
        key != "custom"
        and key not in _BUILT_IN_SCORES
        and similarity_function is None
    ):
        raise ValueError(
            "unsupported similarity '{}'; supply similarity_function".format(
                key
            )
        )
    if similarity_function is not None and key != "dual-cosine":
        key = "custom"
    return key, similarity_function


def _prepare_data(data_matrix, similarity_key, idx_pos, idx_neg):
    if sp.issparse(data_matrix):
        data_work = data_matrix.astype(
            data_matrix.dtype
            if np.issubdtype(data_matrix.dtype, np.floating)
            else np.float64,
            copy=True,
        ).tocsr()
    else:
        array = np.asarray(data_matrix)
        dtype = (
            array.dtype
            if np.issubdtype(array.dtype, np.floating)
            else np.float64
        )
        data_work = np.array(array, dtype=dtype, copy=True)
    if similarity_key == "dual-cosine":
        return _normalize_dual_l2(data_work, idx_pos, idx_neg)
    if similarity_key == "cosine":
        return _normalize_rows_l2(data_work)
    if similarity_key in _BUILT_IN_SCORES:
        return _normalize_rows_l1(data_work)
    return data_work


def _capacity_at(schedule, iteration, maximum_capacity):
    if schedule is None:
        value = maximum_capacity
    elif callable(schedule):
        value = schedule(iteration)
    else:
        sequence = np.asarray(schedule).reshape(-1)
        if sequence.size == 0:
            raise ValueError("capacity_schedule cannot be empty")
        value = sequence[min(iteration - 1, sequence.size - 1)]
    if not _is_positive_integer(value):
        raise ValueError("capacity schedule values must be positive integers")
    return min(int(value), int(maximum_capacity))


def _capacity_schedule_is_fixed(schedule, iteration):
    """Return whether the implementation can certify a constant future tail."""
    if schedule is None:
        return True
    if callable(schedule):
        return False
    sequence = np.asarray(schedule).reshape(-1)
    if sequence.size == 0:
        raise ValueError("capacity_schedule cannot be empty")
    return iteration >= sequence.size


# ---------------------------------------------------------------------------
# Deterministic state transition
# ---------------------------------------------------------------------------
def _select_seeds(
    data,
    seed_count,
    similarity_key,
    similarity_function,
    idx_pos,
    idx_neg,
    batch_size,
):
    sample_count = data.shape[0]
    seed_idx = np.zeros(seed_count, dtype=np.int64)
    if seed_count == 0:
        return seed_idx
    available = np.ones(sample_count, dtype=bool)
    seed_idx[0] = _lexicographic_min(data, np.flatnonzero(available))
    available[seed_idx[0]] = False
    nearest_similarity = np.full(sample_count, -np.inf)

    for seed_number in range(1, seed_count):
        new_center = _dense_rows(
            data, np.array([seed_idx[seed_number - 1]])
        )
        for first in range(0, sample_count, batch_size):
            last = min(first + batch_size, sample_count)
            similarities = _similarity_matrix(
                data[first:last],
                new_center,
                similarity_key,
                similarity_function,
                idx_pos,
                idx_neg,
            )
            nearest_similarity[first:last] = np.maximum(
                nearest_similarity[first:last], similarities[:, 0]
            )
        candidate_scores = nearest_similarity.copy()
        candidate_scores[~available] = np.inf
        minimum_score = np.min(candidate_scores)
        candidates = np.flatnonzero(
            available & (candidate_scores == minimum_score)
        )
        seed_idx[seed_number] = _lexicographic_min(data, candidates)
        available[seed_idx[seed_number]] = False
    return seed_idx


def _assign_in_blocks(
    data,
    centers,
    counts,
    strategy_key,
    threshold,
    similarity_key,
    similarity_function,
    idx_pos,
    idx_neg,
    options,
):
    sample_count = data.shape[0]
    cluster_count = centers.shape[0]
    labels = np.zeros(sample_count, dtype=np.int64)
    max_affinity = np.full(sample_count, -np.inf)
    if cluster_count == 0:
        return labels, max_affinity

    if strategy_key in {"DASS", "DENSITYFIRST"}:
        support_bonus = np.asarray(
            options["support_function"](counts.astype(float))
        ).reshape(-1)
        if (
            support_bonus.size != cluster_count
            or np.iscomplexobj(support_bonus)
            or np.any(~np.isfinite(support_bonus))
        ):
            raise ValueError(
                "support_function must return one finite value per support"
            )
        support_bonus = options["lambda_weight"] * support_bonus
    else:
        support_bonus = np.zeros(cluster_count)

    batch_size = options["batch_size"]
    for first in range(0, sample_count, batch_size):
        last = min(first + batch_size, sample_count)
        similarities = _similarity_matrix(
            data[first:last],
            centers,
            similarity_key,
            similarity_function,
            idx_pos,
            idx_neg,
        )
        max_affinity[first:last] = np.max(similarities, axis=1)
        scores = similarities + support_bonus[None, :]
        scores[similarities < threshold] = -np.inf
        best_cluster = np.argmax(scores, axis=1)
        best_score = scores[np.arange(last - first), best_cluster]
        accepted = np.isfinite(best_score)
        block_labels = np.zeros(last - first, dtype=np.int64)
        block_labels[accepted] = best_cluster[accepted] + 1
        labels[first:last] = block_labels
    return labels, max_affinity


def _select_promotions(
    data,
    labels,
    max_affinity,
    hole_count,
    similarity_key,
    similarity_function,
    representative_function,
    idx_pos,
    idx_neg,
    batch_size,
):
    outliers = labels == 0
    promotion_count = min(int(hole_count), int(np.sum(outliers)))
    promotion_idx = np.zeros(promotion_count, dtype=np.int64)
    available = outliers.copy()
    affinity = max_affinity.copy()

    for promotion_number in range(promotion_count):
        score = affinity.copy()
        score[~available] = np.inf
        minimum_score = np.min(score)
        candidates = np.flatnonzero(
            available & (score == minimum_score)
        )
        chosen = _lexicographic_min(data, candidates)
        promotion_idx[promotion_number] = chosen
        available[chosen] = False

        representative = _representative(
            data[chosen : chosen + 1],
            similarity_key,
            similarity_function,
            representative_function,
            idx_pos,
            idx_neg,
            batch_size,
        )[None, :]
        available_rows = np.flatnonzero(available)
        for first in range(0, available_rows.size, batch_size):
            rows = available_rows[first : first + batch_size]
            similarities = _similarity_matrix(
                data[rows],
                representative,
                similarity_key,
                similarity_function,
                idx_pos,
                idx_neg,
            )
            affinity[rows] = np.maximum(
                affinity[rows], similarities[:, 0]
            )
    return promotion_idx


def _consolidate(
    data,
    labels,
    protected,
    min_volume,
    merge_threshold,
    similarity_key,
    similarity_function,
    representative_function,
    idx_pos,
    idx_neg,
    batch_size,
):
    cluster_count = int(np.max(labels, initial=0))
    if protected.size < cluster_count:
        protected = np.pad(
            protected, (0, cluster_count - protected.size)
        )
    else:
        protected = protected[:cluster_count].copy()

    centers, counts = _recompute_clusters(
        data,
        labels,
        cluster_count,
        similarity_key,
        similarity_function,
        representative_function,
        idx_pos,
        idx_neg,
        batch_size,
    )

    while cluster_count > 1:
        center_similarity = _similarity_matrix(
            centers,
            centers,
            similarity_key,
            similarity_function,
            idx_pos,
            idx_neg,
        )
        adjacency = center_similarity >= merge_threshold
        adjacency = adjacency | adjacency.T
        np.fill_diagonal(adjacency, True)
        components = _connected_components(adjacency)
        new_count = int(np.max(components))
        if new_count == cluster_count:
            break

        remap = np.concatenate(
            [np.array([0], dtype=np.int64), components]
        )
        labels = remap[labels]
        new_protected = np.zeros(new_count, dtype=bool)
        for old_index, component in enumerate(components):
            new_protected[component - 1] |= protected[old_index]
        protected = new_protected
        cluster_count = new_count
        centers, counts = _recompute_clusters(
            data,
            labels,
            cluster_count,
            similarity_key,
            similarity_function,
            representative_function,
            idx_pos,
            idx_neg,
            batch_size,
        )

    if cluster_count:
        keep = (counts >= min_volume) | protected
        if not np.all(keep):
            remap = np.zeros(cluster_count + 1, dtype=np.int64)
            remap[1:][keep] = np.arange(1, np.sum(keep) + 1)
            labels = remap[labels]
            centers = centers[keep]
            counts = counts[keep]

    return _canonicalize_clusters(centers, counts, labels)


def _connected_components(adjacency):
    cluster_count = adjacency.shape[0]
    components = np.zeros(cluster_count, dtype=np.int64)
    component_count = 0
    for start in range(cluster_count):
        if components[start] != 0:
            continue
        component_count += 1
        queue = [start]
        components[start] = component_count
        head = 0
        while head < len(queue):
            current = queue[head]
            head += 1
            neighbors = np.flatnonzero(
                adjacency[current] & (components == 0)
            )
            if neighbors.size:
                components[neighbors] = component_count
                queue.extend(neighbors.tolist())
    return components


def _recompute_clusters(
    data,
    labels,
    cluster_count,
    similarity_key,
    similarity_function,
    representative_function,
    idx_pos,
    idx_neg,
    batch_size,
):
    centers = np.zeros(
        (cluster_count, data.shape[1]), dtype=data.dtype
    )
    counts = np.zeros(cluster_count, dtype=np.int64)
    for cluster_index in range(1, cluster_count + 1):
        member_idx = np.flatnonzero(labels == cluster_index)
        counts[cluster_index - 1] = member_idx.size
        if member_idx.size:
            centers[cluster_index - 1] = _representative(
                data[member_idx],
                similarity_key,
                similarity_function,
                representative_function,
                idx_pos,
                idx_neg,
                batch_size,
            )
    return centers, counts


def _representative(
    members,
    similarity_key,
    similarity_function,
    representative_function,
    idx_pos,
    idx_neg,
    batch_size,
):
    dimension_count = members.shape[1]
    if representative_function is not None:
        center = np.asarray(representative_function(members))
        center = center.reshape(-1)
        if (
            center.size != dimension_count
            or np.iscomplexobj(center)
            or np.any(~np.isfinite(center))
        ):
            raise ValueError(
                "representative_function must return one finite D-vector"
            )
        return center.astype(members.dtype, copy=False)

    if similarity_key == "dual-cosine":
        center = _row_sum(members) / members.shape[0]
        summed = _row_sum(members)
        center[idx_pos] = _normalize_vector_l2(summed[idx_pos])
        center[idx_neg] = _normalize_vector_l2(summed[idx_neg])
        return center
    if similarity_key == "cosine":
        return _normalize_vector_l2(_row_sum(members))
    if similarity_key in {"euclidean", "euclidean-distance"}:
        return _row_sum(members) / members.shape[0]
    if similarity_key in {"l1", "l1-norm", "manhattan"}:
        median = np.median(_dense_rows(members, slice(None)), axis=0)
        return _normalize_vector_l1(median)
    return _similarity_medoid(
        members,
        similarity_key,
        similarity_function,
        idx_pos,
        idx_neg,
        batch_size,
    )


def _similarity_medoid(
    members,
    similarity_key,
    similarity_function,
    idx_pos,
    idx_neg,
    requested_batch_size,
):
    member_count = members.shape[0]
    if member_count == 1:
        return get_row(members, 0)
    max_entries = 10000000
    batch_size = min(
        requested_batch_size,
        max(1, int(np.floor(max_entries / member_count))),
    )
    scores = np.zeros(member_count)
    for first in range(0, member_count, batch_size):
        last = min(first + batch_size, member_count)
        similarities = _similarity_matrix(
            members[first:last],
            members,
            similarity_key,
            similarity_function,
            idx_pos,
            idx_neg,
        )
        scores[first:last] = np.sum(similarities, axis=1)
    best_score = np.max(scores)
    candidates = np.flatnonzero(scores == best_score)
    return get_row(members, _lexicographic_min(members, candidates))


def _canonicalize_clusters(centers, counts, labels):
    cluster_count = counts.size
    if cluster_count <= 1:
        return centers, counts, labels
    keys = np.column_stack([-counts.astype(float), centers.astype(float)])
    lexicographic_keys = tuple(
        keys[:, index] for index in range(keys.shape[1] - 1, -1, -1)
    )
    order = np.lexsort(lexicographic_keys)
    inverse = np.zeros(cluster_count, dtype=np.int64)
    inverse[order] = np.arange(1, cluster_count + 1)
    remap = np.concatenate(
        [np.array([0], dtype=np.int64), inverse]
    )
    return centers[order], counts[order], remap[labels]


def _project_mature_state(centers, counts, labels, min_volume):
    """Map one-pass protected clusters below min_volume back to outliers."""
    keep = counts >= int(min_volume)
    if np.all(keep):
        return _capture_state(centers, counts, labels)
    remap = np.zeros(counts.size + 1, dtype=np.int64)
    remap[np.flatnonzero(keep) + 1] = np.arange(
        1, int(np.sum(keep)) + 1, dtype=np.int64
    )
    projected_labels = remap[labels]
    projected_centers = centers[keep].copy()
    projected_counts = counts[keep].copy()
    return _canonicalize_clusters(
        projected_centers, projected_counts, projected_labels
    )


# ---------------------------------------------------------------------------
# Objective, diagnostics, and state identity
# ---------------------------------------------------------------------------
def _objective(
    data,
    centers,
    counts,
    labels,
    strategy_key,
    similarity_key,
    similarity_function,
    idx_pos,
    idx_neg,
    options,
):
    objective_function = options["objective_function"]
    if objective_function is not None:
        score = objective_function(data, centers, counts, labels)
        if (
            not np.isscalar(score)
            or np.iscomplexobj(score)
            or not np.isfinite(score)
        ):
            raise ValueError(
                "objective_function must return one finite scalar"
            )
        return float(score)

    cohesion = 0.0
    assigned_rows = np.flatnonzero(labels > 0)
    batch_size = options["batch_size"]
    for first in range(0, assigned_rows.size, batch_size):
        rows = assigned_rows[first : first + batch_size]
        similarities = _similarity_matrix(
            data[rows],
            centers,
            similarity_key,
            similarity_function,
            idx_pos,
            idx_neg,
        )
        cohesion += np.sum(
            similarities[
                np.arange(rows.size), labels[rows].astype(int) - 1
            ]
        )

    if strategy_key in {"DASS", "DENSITYFIRST"}:
        support_potential = np.asarray(
            options["support_potential_function"](
                counts.astype(float)
            )
        ).reshape(-1)
        if (
            support_potential.size != counts.size
            or np.iscomplexobj(support_potential)
            or np.any(~np.isfinite(support_potential))
        ):
            raise ValueError(
                "support_potential_function must return one finite "
                "value per support"
            )
        return float(
            cohesion
            + options["lambda_weight"] * np.sum(support_potential)
        )
    return float(cohesion)


def _stability(
    previous_centers,
    previous_counts,
    centers,
    counts,
    previous_labels,
    labels,
    similarity_key,
    similarity_function,
    idx_pos,
    idx_neg,
):
    if previous_labels is None or previous_labels.size == 0:
        return 0.0, 0.0
    label_agreement = float(np.mean(previous_labels == labels))

    count_length = max(previous_counts.size, counts.size)
    count_a = np.zeros(count_length)
    count_b = np.zeros(count_length)
    count_a[: previous_counts.size] = np.sort(previous_counts)[::-1]
    count_b[: counts.size] = np.sort(counts)[::-1]
    denominator = np.linalg.norm(count_a) * np.linalg.norm(count_b)
    if denominator == 0:
        count_similarity = float(
            np.linalg.norm(count_a - count_b) == 0
        )
    else:
        count_similarity = float(np.dot(count_a, count_b) / denominator)

    if previous_centers.size == 0 or centers.size == 0:
        center_similarity = float(
            previous_centers.size == 0 and centers.size == 0
        )
    else:
        pairwise = _similarity_matrix(
            centers,
            previous_centers,
            similarity_key,
            similarity_function,
            idx_pos,
            idx_neg,
        )
        pair_count = min(pairwise.shape)
        matched = np.zeros(pair_count)
        weights = np.zeros(pair_count)
        used_current = np.zeros(centers.shape[0], dtype=bool)
        used_previous = np.zeros(previous_centers.shape[0], dtype=bool)
        for match_index in range(pair_count):
            available = pairwise.copy()
            available[used_current] = -np.inf
            available[:, used_previous] = -np.inf
            flat_index = np.argmax(available.ravel(order="F"))
            row, column = np.unravel_index(
                flat_index, available.shape, order="F"
            )
            matched[match_index] = available[row, column]
            weights[match_index] = counts[row]
            used_current[row] = True
            used_previous[column] = True
        if np.sum(weights) == 0:
            center_similarity = float(np.mean(matched))
        else:
            center_similarity = float(
                np.sum(matched * weights) / np.sum(weights)
            )
    return (
        (float(count_similarity) + float(center_similarity)) / 2.0,
        label_agreement,
    )


def _capture_state(centers, counts, labels):
    return centers.copy(), counts.copy(), labels.copy()


def _state_hash(labels, capacity):
    digest = hashlib.sha256()
    digest.update(np.asarray([capacity], dtype="<i4").tobytes())
    digest.update(np.asarray(labels, dtype="<i4").tobytes(order="C"))
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# Similarity kernels
# ---------------------------------------------------------------------------
def _similarity_matrix(
    left,
    right,
    similarity_key,
    similarity_function,
    idx_pos,
    idx_neg,
):
    row_count = left.shape[0]
    column_count = right.shape[0]
    if row_count == 0 or column_count == 0:
        return np.zeros((row_count, column_count))

    if similarity_key == "dual-cosine":
        positive = _matrix_product(
            left[:, idx_pos], right[:, idx_pos]
        )
        negative = _matrix_product(
            left[:, idx_neg], right[:, idx_neg]
        )
        similarities = np.minimum(positive, negative)
    elif similarity_key == "custom":
        similarities = _custom_similarity_matrix(
            left, right, similarity_function
        )
    elif similarity_key == "cosine":
        similarities = _matrix_product(
            _normalize_rows_l2(left), _normalize_rows_l2(right)
        )
    elif similarity_key in {"euclidean", "euclidean-distance"}:
        distances = _pairwise_squared_euclidean(left, right)
        similarities = 1.0 - distances / 2.0
    elif similarity_key in {"l1", "l1-norm", "manhattan"}:
        similarities = _pairwise_l1_score(left, right, "l1")
    elif similarity_key == "minimum":
        similarities = _pairwise_l1_score(left, right, "minimum")
    elif similarity_key == "maximum":
        similarities = _pairwise_l1_score(left, right, "maximum")
    else:
        similarities = _exotic_similarity_matrix(
            left, right, similarity_key
        )

    similarities = _dense_output(similarities)
    if (
        similarities.shape != (row_count, column_count)
        or np.iscomplexobj(similarities)
        or np.any(~np.isfinite(similarities))
    ):
        raise ValueError(
            "similarity '{}' returned invalid output".format(
                similarity_key
            )
        )
    return similarities


def _custom_similarity_matrix(left, right, similarity_function):
    try:
        candidate = similarity_function(left, right)
        candidate = _dense_output(candidate)
        if candidate.shape == (left.shape[0], right.shape[0]):
            return candidate
    except Exception:
        pass

    similarities = np.zeros((left.shape[0], right.shape[0]))
    for row_index in range(left.shape[0]):
        left_row = get_row(left, row_index)
        for column_index in range(right.shape[0]):
            right_row = get_row(right, column_index)
            similarities[row_index, column_index] = similarity_function(
                left_row, right_row
            )
    return similarities


def _pairwise_squared_euclidean(left, right):
    left_norm = _row_squared_norm(left)
    right_norm = _row_squared_norm(right)
    distances = (
        left_norm[:, None]
        + right_norm[None, :]
        - 2.0 * _matrix_product(left, right)
    )
    return np.maximum(distances, 0.0)


def _pairwise_l1_score(left, right, mode):
    left_dense = _dense_rows(left, slice(None))
    right_dense = _dense_rows(right, slice(None))
    similarities = np.zeros((left.shape[0], right.shape[0]))
    for column_index, center in enumerate(right_dense):
        if mode == "l1":
            distances = np.sum(np.abs(left_dense - center), axis=1)
            similarities[:, column_index] = 1.0 - distances / 2.0
        elif mode == "minimum":
            similarities[:, column_index] = np.sum(
                np.minimum(left_dense, center), axis=1
            )
        else:
            similarities[:, column_index] = (
                1.0
                - np.sum(np.maximum(left_dense, center), axis=1) / 2.0
            )
    return similarities


def _exotic_similarity_matrix(left, right, algorithm):
    left_dense = _dense_rows(left, slice(None))
    right_dense = _dense_rows(right, slice(None))
    similarities = np.zeros((left.shape[0], right.shape[0]))
    for column_index, center in enumerate(right_dense):
        shared = (left_dense != 0) & (center != 0)
        if algorithm == "algebraic":
            values = (left_dense + center) / 2.0
            similarities[:, column_index] = np.sum(
                np.where(shared, values, 0.0), axis=1
            )
        elif algorithm == "logarithmic":
            values = np.zeros_like(left_dense, dtype=float)
            rows, columns = np.where(shared)
            a = left_dense[rows, columns]
            b = center[columns]
            tolerance = 8.0 * np.spacing(
                np.maximum(np.abs(a), np.abs(b))
            )
            equal = np.abs(a - b) <= tolerance
            terms = np.empty_like(a, dtype=float)
            terms[equal] = (a[equal] + b[equal]) / 2.0
            relative = (b[~equal] - a[~equal]) / a[~equal]
            terms[~equal] = (
                a[~equal] * relative / np.log1p(relative)
            )
            values[rows, columns] = terms
            similarities[:, column_index] = np.sum(values, axis=1)
        elif algorithm == "geometric":
            values = np.sqrt(left_dense * center)
            similarities[:, column_index] = np.sum(
                np.where(shared, values, 0.0), axis=1
            )
        elif algorithm == "harmonic":
            values = np.zeros_like(left_dense, dtype=float)
            numerator = 2.0 * left_dense * center
            denominator = left_dense + center
            np.divide(
                numerator,
                denominator,
                out=values,
                where=shared,
            )
            similarities[:, column_index] = np.sum(values, axis=1)
        elif algorithm == "enhanced harmonic":
            values = np.zeros_like(left_dense, dtype=float)
            left_shared = left_dense[shared]
            center_shared = np.broadcast_to(
                center, left_dense.shape
            )[shared]
            values[shared] = np.sqrt(
                2.0
                / (
                    np.power(left_shared, -2.0)
                    + np.power(center_shared, -2.0)
                )
            )
            similarities[:, column_index] = np.sum(values, axis=1)
        elif algorithm == "entropy":
            entropy_left = _entropy_rows(left_dense)
            entropy_center = _entropy_rows(center[None, :])[0]
            entropy_mean = _entropy_rows(
                (left_dense + center) / 2.0
            )
            similarities[:, column_index] = (
                1.0
                - (
                    2.0 * entropy_mean
                    - entropy_left
                    - entropy_center
                )
                / np.log(4.0)
            )
        elif algorithm == "weighted entropy":
            weighted_left = _entropy_weight_rows(left_dense)
            weighted_center = _entropy_weight_rows(center[None, :])[0]
            entropy_left = _entropy_rows(weighted_left)
            entropy_center = _entropy_rows(
                weighted_center[None, :]
            )[0]
            entropy_mean = _entropy_rows(
                (weighted_left + weighted_center) / 2.0
            )
            similarities[:, column_index] = (
                1.0
                - (
                    2.0 * entropy_mean
                    - entropy_left
                    - entropy_center
                )
                / np.log(4.0)
            )
        elif algorithm in {"best average", "fitted core"}:
            values = np.zeros_like(left_dense, dtype=float)
            left_shared = left_dense[shared]
            center_shared = np.broadcast_to(
                center, left_dense.shape
            )[shared]
            with np.errstate(over="ignore", invalid="ignore"):
                values[shared] = np.power(
                    np.power(left_shared, center_shared)
                    * np.power(center_shared, left_shared),
                    1.0 / (left_shared + center_shared),
                )
            similarities[:, column_index] = np.sum(values, axis=1)
        else:
            raise ValueError(
                "unsupported similarity algorithm: {}".format(algorithm)
            )
    return similarities


def _entropy_rows(matrix):
    positive = matrix > 0
    safe = np.where(positive, matrix, 1.0)
    return -np.sum(
        np.where(positive, matrix * np.log(safe), 0.0), axis=1
    )


def _entropy_weight_rows(matrix):
    entropy = _entropy_rows(matrix)
    weighted = np.array(matrix, copy=True)
    selected = entropy < 3.0
    if np.any(selected):
        exponent = 0.25 + 0.25 * entropy[selected]
        weighted[selected] = np.power(
            weighted[selected], exponent[:, None]
        )
        weighted[selected] = _normalize_rows_l1(
            weighted[selected]
        )
    return weighted


# ---------------------------------------------------------------------------
# Matrix utilities
# ---------------------------------------------------------------------------
def _normalize_rows_l2(matrix):
    if sp.issparse(matrix):
        norms = np.sqrt(
            np.asarray(matrix.multiply(matrix).sum(axis=1)).reshape(-1)
        )
        norms[norms == 0] = 1
        return sp.diags(1.0 / norms).dot(matrix).tocsr()
    array = np.asarray(matrix)
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return array / norms


def _normalize_rows_l1(matrix):
    if sp.issparse(matrix):
        row_sums = np.asarray(abs(matrix).sum(axis=1)).reshape(-1)
        row_sums[row_sums == 0] = 1
        return sp.diags(1.0 / row_sums).dot(matrix).tocsr()
    array = np.asarray(matrix)
    row_sums = np.sum(np.abs(array), axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    return array / row_sums


def _normalize_dual_l2(matrix, idx_pos, idx_neg):
    if sp.issparse(matrix):
        positive_norm = np.sqrt(
            np.asarray(
                matrix[:, idx_pos]
                .multiply(matrix[:, idx_pos])
                .sum(axis=1)
            ).reshape(-1)
        )
        negative_norm = np.sqrt(
            np.asarray(
                matrix[:, idx_neg]
                .multiply(matrix[:, idx_neg])
                .sum(axis=1)
            ).reshape(-1)
        )
        positive_norm[positive_norm == 0] = 1
        negative_norm[negative_norm == 0] = 1
        positive_columns = np.zeros(matrix.shape[1], dtype=bool)
        negative_columns = np.zeros(matrix.shape[1], dtype=bool)
        positive_columns[idx_pos] = True
        negative_columns[idx_neg] = True
        coo = matrix.tocoo(copy=True)
        scale = np.ones(coo.data.size)
        pos_entries = positive_columns[coo.col]
        neg_entries = negative_columns[coo.col]
        scale[pos_entries] = 1.0 / positive_norm[coo.row[pos_entries]]
        scale[neg_entries] = 1.0 / negative_norm[coo.row[neg_entries]]
        coo.data *= scale
        return coo.tocsr()

    normalized = np.array(matrix, copy=True)
    normalized[:, idx_pos] = _normalize_rows_l2(
        normalized[:, idx_pos]
    )
    normalized[:, idx_neg] = _normalize_rows_l2(
        normalized[:, idx_neg]
    )
    return normalized


def _normalize_vector_l2(vector):
    norm = np.linalg.norm(vector)
    if norm == 0:
        norm = 1
    return np.asarray(vector) / norm


def _normalize_vector_l1(vector):
    row_sum = np.sum(np.abs(vector))
    if row_sum == 0:
        row_sum = 1
    return np.asarray(vector) / row_sum


def _row_sum(matrix):
    if sp.issparse(matrix):
        return np.asarray(matrix.sum(axis=0)).reshape(-1)
    return np.sum(np.asarray(matrix), axis=0)


def _row_squared_norm(matrix):
    if sp.issparse(matrix):
        return np.asarray(
            matrix.multiply(matrix).sum(axis=1)
        ).reshape(-1)
    return np.sum(np.asarray(matrix) ** 2, axis=1)


def _matrix_product(left, right):
    product = left @ right.T
    return _dense_output(product)


def _dense_output(value):
    if sp.issparse(value):
        return value.toarray()
    if hasattr(value, "A"):
        return np.asarray(value.A)
    return np.asarray(value)


def _dense_rows(matrix, indices):
    rows = matrix[indices]
    if sp.issparse(rows):
        return rows.toarray()
    return np.asarray(rows)


def _lexicographic_min(data, candidates):
    candidates = np.asarray(candidates, dtype=np.int64).reshape(-1)
    if candidates.size == 0:
        raise RuntimeError("no deterministic candidate is available")
    for dimension_index in range(data.shape[1]):
        if sp.issparse(data):
            values = np.asarray(
                data[candidates, dimension_index].toarray()
            ).reshape(-1)
        else:
            values = np.asarray(
                data[candidates, dimension_index]
            ).reshape(-1)
        minimum_value = np.min(values)
        candidates = candidates[values == minimum_value]
        if candidates.size == 1:
            break
    return int(candidates[0])


def _format_labels(labels, convention):
    if convention == "matlab":
        return labels.copy()
    return np.where(labels > 0, labels - 1, -1)


def _index_array(indices):
    if indices is None:
        return np.array([], dtype=np.int64)
    array = np.asarray(indices)
    if array.ndim != 1:
        raise ValueError("feature indices must be one-dimensional")
    if array.size and (
        not np.issubdtype(array.dtype, np.integer)
        or np.any(~np.isfinite(array))
    ):
        raise ValueError("feature indices must contain integers")
    return array.astype(np.int64, copy=False)


def _is_positive_integer(value):
    return (
        np.isscalar(value)
        and np.isfinite(value)
        and float(value) >= 1
        and float(value).is_integer()
    )
