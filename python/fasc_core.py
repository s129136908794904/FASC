import numpy as np
import scipy.sparse as sp
from scipy.spatial.distance import cdist
import time

# ====================================================================
# Data Extraction & Normalization Helpers
# ====================================================================
def get_row(mat, idx):
    if sp.issparse(mat):
        return mat[idx].toarray().flatten()
    return mat[idx]

def get_rows(mat, idxs):
    if sp.issparse(mat):
        return mat[idxs].toarray()
    return mat[idxs]

def get_mean(mat, mask, count):
    if sp.issparse(mat):
        return np.asarray(mat[mask].sum(axis=0)).flatten() / count
    return np.mean(mat[mask], axis=0)

def normalize_l2(mat):
    """L2 normalization for standard Cosine."""
    if sp.issparse(mat):
        norms = np.array(np.sqrt(mat.multiply(mat).sum(axis=1))).flatten()
        norms[norms == 0] = 1
        return mat.multiply(1.0 / norms[:, None]).tocsr()
    else:
        norms = np.linalg.norm(mat, axis=1, keepdims=True)
        norms[norms == 0] = 1
        return mat / norms

def normalize_dual_l2(mat, idx_pos, idx_neg):
    """Independent L2 normalization for positive and negative modes (Dual-Cosine)."""
    if sp.issparse(mat):
        mat = mat.tolil()
        pos_norms = np.array(np.sqrt(mat[:, idx_pos].multiply(mat[:, idx_pos]).sum(axis=1))).flatten()
        neg_norms = np.array(np.sqrt(mat[:, idx_neg].multiply(mat[:, idx_neg]).sum(axis=1))).flatten()
        pos_norms[pos_norms == 0] = 1
        neg_norms[neg_norms == 0] = 1
        
        # In-place scaling for LIL matrix rows is tricky, easier via CSR format element-wise
        mat = mat.tocsr()
        mat[:, idx_pos] = mat[:, idx_pos].multiply(1.0 / pos_norms[:, None])
        mat[:, idx_neg] = mat[:, idx_neg].multiply(1.0 / neg_norms[:, None])
        return mat
    else:
        mat = mat.copy()
        pos_norms = np.linalg.norm(mat[:, idx_pos], axis=1, keepdims=True)
        neg_norms = np.linalg.norm(mat[:, idx_neg], axis=1, keepdims=True)
        pos_norms[pos_norms == 0] = 1
        neg_norms[neg_norms == 0] = 1
        mat[:, idx_pos] /= pos_norms
        mat[:, idx_neg] /= neg_norms
        return mat

def normalize_l1(mat):
    """L1 normalization for non-cosine metrics (Manhattan, Entropy, Harmonic, etc.)."""
    if sp.issparse(mat):
        norms = np.array(mat.sum(axis=1)).flatten()
        norms[norms == 0] = 1
        return mat.multiply(1.0 / norms[:, None]).tocsr()
    else:
        norms = np.sum(mat, axis=1, keepdims=True)
        norms[norms == 0] = 1
        return mat / norms

# ====================================================================
# Similarity Computation Core
# ====================================================================
def compute_sim_matrix(X, C, algo_key, idx_pos=None, idx_neg=None):
    """
    Computes similarities between Dataset X (N x D) and Centers C (K x D).
    Uses smart chunking for memory safety on massive sparse matrices.
    """
    algo = algo_key.lower()
    
    # --- 1. Dual Cosine ---
    if algo == 'dual-cosine':
        if sp.issparse(X):
            sims_pos = np.asarray(X[:, idx_pos] @ C[:, idx_pos].T)
            sims_neg = np.asarray(X[:, idx_neg] @ C[:, idx_neg].T)
        else:
            sims_pos = X[:, idx_pos] @ C[:, idx_pos].T
            sims_neg = X[:, idx_neg] @ C[:, idx_neg].T
        return np.minimum(sims_pos, sims_neg)
        
    # --- 2. Standard Cosine ---
    elif algo == 'cosine':
        if sp.issparse(X):
            return np.asarray(X @ C.T)
        else:
            return X @ C.T
            
    # --- 3. Distance & Exotic Metrics (Memory-Safe Chunking) ---
    else:
        N = X.shape[0]
        K = C.shape[0]
        sims = np.zeros((N, K))
        
        chunk_size = 10000 # Safely convert to dense in chunks of 10,000 spectra
        eps = 1e-12
        
        for i in range(0, N, chunk_size):
            end = min(i + chunk_size, N)
            A = X[i:end].toarray() if sp.issparse(X) else X[i:end]
            
            if algo in ['euclidean', 'euclidean-distance']:
                D = cdist(A, C, metric='euclidean')
                sims[i:end] = 1 - D / np.sqrt(2)
                
            elif algo in['l1', 'l1-norm', 'manhattan', 'minimum', 'maximum']:
                D = cdist(A, C, metric='cityblock')
                if algo == 'maximum':
                    sims[i:end] = 0.5 - D / 4  # equivalent to 1 - sum(max)/2 for L1 normalized
                else: # minimum or l1
                    sims[i:end] = 1 - D / 2
                    
            else:
                # Dense Broadcasting for Exotic Metrics
                for k in range(K):
                    B = C[k, :]
                    mask = (A > 0) & (B > 0)
                    
                    if algo == 'algebraic':
                        vals = (A + B) / 2
                        sims[i:end, k] = np.sum(np.where(mask, vals, 0), axis=1)
                        
                    elif algo == 'geometric':
                        vals = np.sqrt(A * B)
                        sims[i:end, k] = np.sum(np.where(mask, vals, 0), axis=1)
                        
                    elif algo == 'harmonic':
                        vals = 2 * (A * B) / (A + B + eps)
                        sims[i:end, k] = np.sum(np.where(mask, vals, 0), axis=1)
                        
                    elif algo == 'enhanced harmonic':
                        vals = np.sqrt(2 * (A**2 * B**2) / (A**2 + B**2 + eps))
                        sims[i:end, k] = np.sum(np.where(mask, vals, 0), axis=1)
                        
                    elif algo == 'logarithmic':
                        log_diff = np.log(B + eps) - np.log(A + eps)
                        vals = np.where(np.abs(log_diff) < 1e-6, (A+B)/2, (B - A) / log_diff)
                        sims[i:end, k] = np.sum(np.where(mask, vals, 0), axis=1)
                        
                    elif algo in ['best average', 'fitted core']:
                        vals = np.power(np.power(A, B) * np.power(B, A), 1 / (A + B + eps))
                        sims[i:end, k] = np.sum(np.where(mask, vals, 0), axis=1)
                        
                    elif algo == 'entropy':
                        def H(v):
                            v_nz = np.where(v > 0, v, 1)
                            return -np.sum(np.where(v > 0, v * np.log(v_nz), 0), axis=1)
                        
                        H_A = H(A)
                        B_nz = np.where(B > 0, B, 1)
                        H_B = -np.sum(np.where(B > 0, B * np.log(B_nz), 0))
                        H_mid = H((A + B) / 2)
                        sims[i:end, k] = 1 - (2 * H_mid - H_A - H_B) / np.log(4)
                        
                    elif algo == 'weighted entropy':
                        def H_and_w(v_mat):
                            v_nz = np.where(v_mat > 0, v_mat, 1)
                            ent = -np.sum(np.where(v_mat > 0, v_mat * np.log(v_nz), 0), axis=1, keepdims=True)
                            mask_ent = ent < 3
                            w_mat = np.where(mask_ent, np.power(v_mat, 0.25 + 0.25 * ent), v_mat)
                            sums = np.sum(w_mat, axis=1, keepdims=True)
                            sums[sums == 0] = 1
                            w_mat = w_mat / sums
                            w_nz = np.where(w_mat > 0, w_mat, 1)
                            new_ent = -np.sum(np.where(w_mat > 0, w_mat * np.log(w_nz), 0), axis=1)
                            return w_mat, new_ent
                            
                        A_w, H_A = H_and_w(A)
                        B_w, H_B = H_and_w(B[None, :])
                        B_w, H_B = B_w[0], H_B[0]
                        
                        mid = (A_w + B_w) / 2
                        mid_nz = np.where(mid > 0, mid, 1)
                        H_mid = -np.sum(np.where(mid > 0, mid * np.log(mid_nz), 0), axis=1)
                        sims[i:end, k] = 1 - (2 * H_mid - H_A - H_B) / np.log(4)
                        
                    else:
                        raise ValueError(f"Unsupported similarity algorithm: {algo}")
        return sims

# ====================================================================
# FASC Execution and Objective Function
# ====================================================================
def calculate_objective(data, centers, counts, assignments, valid_k, strategy, algo, idx_pos, idx_neg):
    mask = (assignments >= 0) & (assignments < valid_k)
    if not np.any(mask):
        return -np.inf
    
    # We only care about the similarity to the ASSIGNED center
    assigned_data = data[mask]
    assigned_centers = centers[assignments[mask]]
    
    # Compute similarity purely pair-wise to save memory
    algo_key = algo.lower()
    if algo_key == 'dual-cosine':
        if sp.issparse(assigned_data):
            sims_pos = np.asarray(assigned_data[:, idx_pos].multiply(assigned_centers[:, idx_pos]).sum(axis=1)).flatten()
            sims_neg = np.asarray(assigned_data[:, idx_neg].multiply(assigned_centers[:, idx_neg]).sum(axis=1)).flatten()
        else:
            sims_pos = np.sum(assigned_data[:, idx_pos] * assigned_centers[:, idx_pos], axis=1)
            sims_neg = np.sum(assigned_data[:, idx_neg] * assigned_centers[:, idx_neg], axis=1)
        similarities = np.minimum(sims_pos, sims_neg) 
        
    elif algo_key == 'cosine':
        if sp.issparse(assigned_data):
            similarities = np.asarray(assigned_data.multiply(assigned_centers).sum(axis=1)).flatten()
        else:
            similarities = np.sum(assigned_data * assigned_centers, axis=1)
    else:
        # For complex metrics, extract the diagonal of the matrix computation
        # (Chunked for safety)
        similarities = np.zeros(assigned_data.shape[0])
        chunk = 10000
        for i in range(0, assigned_data.shape[0], chunk):
            end = min(i + chunk, assigned_data.shape[0])
            S_mat = compute_sim_matrix(assigned_data[i:end], assigned_centers[i:end], algo_key)
            similarities[i:end] = np.diag(S_mat)
        
    sum_sim = np.sum(similarities)
    
    if strategy.upper() in["DASS", "DENSITYFIRST"]:
        c = counts[:valid_k]
        density_potential = np.sum(c * (c - 1) / 2) # Matches \binom{n_j}{2}
        return sum_sim + density_potential
    else:
        return sum_sim

def run_fasc(data_matrix, idx_pos, idx_neg, sim_inter, sim_inner, 
             init_limit, max_clusters, max_iter, strategy, min_vol, algo, log_callback=None):
    
    def log(msg):
        if log_callback: log_callback(msg)
        else: print(msg)

    log("========= FASC Start =========")
    start_time = time.time()
    
    N, D = data_matrix.shape
    algo_key = algo.lower()
    is_dual = (algo_key == 'dual-cosine')
    
    # 0. Initial Pre-normalization (Based on selected metric)
    log(f"Pre-normalizing data for '{algo_key}' metric...")
    if algo_key == 'dual-cosine':
        data_matrix = normalize_dual_l2(data_matrix, idx_pos, idx_neg)
    elif algo_key == 'cosine':
        data_matrix = normalize_l2(data_matrix)
    else:
        data_matrix = normalize_l1(data_matrix)
        
    initial_size = 150
    counts = np.zeros(initial_size, dtype=int)
    centers = np.zeros((initial_size, D), dtype=np.float32)
    assignments = np.full(N, -1, dtype=int)
    
    valid_k = 0
    centers[0] = get_row(data_matrix, 0)
    
    iter_sim_stack = np.zeros(1)
    iter_data = np.zeros((max_iter, 3))
    
    best_obj = -np.inf
    best_state = {}
    
    sim_history =[]
    cycle_window = 20
    no_improve_count = 0
    early_stop_limit = 100
    converged = False
    cycle_detected = False
    
    for it in range(max_iter):
        iter_start = time.time()
        counts_last = counts.copy()
        meaningful_k_last = np.sum(counts >= min_vol)
        centers_last = centers.copy()
        valid_k_last = valid_k
        
        log(f"Iteration {it+1} start")
        
        # --- PHASE 1: Assignment ---
        if it == 0:
            np.random.seed(42)
            rand_idx = np.random.permutation(N)
            for i in rand_idx:
                row_i = get_row(data_matrix, i) 
                
                if valid_k == 0:
                    valid_k = 1
                    centers[0] = row_i
                    counts[0] = 1
                    assignments[i] = 0
                    continue
                    
                active_c = centers[:valid_k]
                score_vec = compute_sim_matrix(row_i[None, :], active_c, algo_key, idx_pos, idx_neg).flatten()
                
                mask = (score_vec >= sim_inner)
                if np.any(mask):
                    if strategy.upper() in ['DASS', 'DENSITYFIRST']:
                        score = score_vec + counts[:valid_k]
                    else:
                        score = score_vec.copy()
                        
                    score[~mask] = -np.inf
                    best_c = np.argmax(score)
                    
                    c_old = centers[best_c] * counts[best_c]
                    counts[best_c] += 1
                    c_new = c_old + row_i
                    
                    # Normalize based on algorithm
                    if algo_key == 'dual-cosine':
                        c_new[idx_pos] /= (np.linalg.norm(c_new[idx_pos]) + 1e-12)
                        c_new[idx_neg] /= (np.linalg.norm(c_new[idx_neg]) + 1e-12)
                    elif algo_key == 'cosine':
                        c_new /= (np.linalg.norm(c_new) + 1e-12)
                    else:
                        c_new /= (np.sum(c_new) + 1e-12)
                        
                    centers[best_c] = c_new
                    assignments[i] = best_c
                else:
                    if valid_k < init_limit:
                        if valid_k == len(centers):
                            centers = np.vstack([centers, np.zeros((valid_k, D))])
                            counts = np.append(counts, np.zeros(valid_k, dtype=int))
                        centers[valid_k] = row_i
                        counts[valid_k] = 1
                        assignments[i] = valid_k
                        valid_k += 1
                    else:
                        assignments[i] = -1
            fill_idx =[]
        else:
            active_c = centers[:valid_k]
            
            # Massive matrix computation properly chunked inside helper function
            S_mat = compute_sim_matrix(data_matrix, active_c, algo_key, idx_pos, idx_neg)
            mask = (S_mat >= sim_inner)
            
            if strategy.upper() in['DASS', 'DENSITYFIRST']:
                score = S_mat + counts[:valid_k]
            else:
                score = S_mat.copy()
            
            score[~mask] = -np.inf
            best_scores = np.max(score, axis=1)
            best_c = np.argmax(score, axis=1)
            
            assigned = best_scores != -np.inf
            assignments = np.where(assigned, best_c, -1)
            
            outliers = np.sum(~assigned)
            holes = max_clusters - valid_k
            if holes > 0 and outliers > 0:
                fill_count = min(holes, outliers)
                outlier_idx = np.where(~assigned)[0]
                fill_idx = np.random.choice(outlier_idx, fill_count, replace=False)
            else:
                fill_idx =[]

        # --- PHASE 2: Post-Iteration (Kill & Merge) ---
        # 1. Kill
        old_valid_k = valid_k
        remained = np.zeros(old_valid_k, dtype=bool)
        for i in range(old_valid_k):
            mask = (assignments == i)
            c_count = np.sum(mask)
            counts[i] = c_count
            if c_count >= min_vol:
                new_c = get_mean(data_matrix, mask, c_count)
                
                if algo_key == 'dual-cosine':
                    new_c[idx_pos] /= (np.linalg.norm(new_c[idx_pos]) + 1e-12)
                    new_c[idx_neg] /= (np.linalg.norm(new_c[idx_neg]) + 1e-12)
                elif algo_key == 'cosine':
                    new_c /= (np.linalg.norm(new_c) + 1e-12)
                else:
                    new_c /= (np.sum(new_c) + 1e-12)
                    
                centers[i] = new_c
                remained[i] = True
            else:
                centers[i] = 0
                counts[i] = 0
                assignments[mask] = -1
                
        active_centers = centers[:old_valid_k][remained]
        active_counts = counts[:old_valid_k][remained]
        
        valid_k = np.sum(remained)
        centers[:valid_k] = active_centers
        counts[:valid_k] = active_counts
        
        old_ids = np.where(remained)[0]
        id_map = {old: new for new, old in enumerate(old_ids)}
        id_map[-1] = -1
        assignments = np.vectorize(id_map.get)(assignments)
        
        sort_idx = np.argsort(counts[:valid_k])[::-1]
        centers[:valid_k] = centers[:valid_k][sort_idx]
        counts[:valid_k] = counts[:valid_k][sort_idx]
        inv_sort = {old: new for new, old in enumerate(sort_idx)}
        inv_sort[-1] = -1
        assignments = np.vectorize(inv_sort.get)(assignments)
        
        if len(fill_idx) > 0:
            vk_before = valid_k
            valid_k += len(fill_idx)
            if valid_k > len(centers):
                centers = np.vstack([centers, np.zeros((valid_k, D))])
                counts = np.append(counts, np.zeros(valid_k, dtype=int))
            centers[vk_before:valid_k] = get_rows(data_matrix, fill_idx)
            counts[vk_before:valid_k] = 1
            assignments[fill_idx] = np.arange(vk_before, valid_k)
            
        # 2. Merge
        merged_mask = np.ones(valid_k, dtype=bool)
        curr_idx = 0
        while True:
            rem_idx = np.where(merged_mask)[0]
            if curr_idx >= len(rem_idx): break
            comp_idx = rem_idx[curr_idx:]
            target_c = centers[comp_idx[0]]
            
            s_vec = compute_sim_matrix(centers[comp_idx], target_c[None, :], algo_key, idx_pos, idx_neg).flatten()
            merge_mask = (s_vec >= sim_inter)
            m_idx = comp_idx[merge_mask]
            
            if len(m_idx) > 1:
                m_count = np.sum(counts[m_idx])
                w_centers = centers[m_idx] * counts[m_idx][:, None]
                n_center = np.sum(w_centers, axis=0)
                
                if algo_key == 'dual-cosine':
                    n_center[idx_pos] /= (np.linalg.norm(n_center[idx_pos]) + 1e-12)
                    n_center[idx_neg] /= (np.linalg.norm(n_center[idx_neg]) + 1e-12)
                elif algo_key == 'cosine':
                    n_center /= (np.linalg.norm(n_center) + 1e-12)
                else:
                    n_center /= (np.sum(n_center) + 1e-12)
                    
                centers[m_idx[0]] = n_center
                centers[m_idx[1:]] = 0
                merged_mask[m_idx[1:]] = False
                counts[m_idx[0]] = m_count
                counts[m_idx[1:]] = 0
                
                for idx in m_idx[1:]:
                    assignments[assignments == idx] = m_idx[0]
                curr_idx = 0
            else:
                curr_idx += 1
                
        active_centers = centers[:valid_k][merged_mask]
        active_counts = counts[:valid_k][merged_mask]
        
        rem_count = np.sum(merged_mask)
        valid_k = min(rem_count, max_clusters)
        centers[:valid_k] = active_centers[:valid_k]
        counts[:valid_k] = active_counts[:valid_k]
        
        active_ids = np.unique(assignments[assignments >= 0])
        merge_map = {old: new for new, old in enumerate(active_ids[:valid_k])}
        merge_map[-1] = -1
        assignments = np.vectorize(lambda x: merge_map.get(x, -1))(assignments)

        # --- PHASE 3: Convergence & Objective ---
        iter_sim = 0
        if it > 0:
            c1 = counts_last[:meaningful_k_last]
            c2 = counts[:np.sum(counts >= min_vol)]
            min_len = min(len(c1), len(c2))
            if min_len > 0:
                n1, n2 = np.linalg.norm(c1[:min_len]), np.linalg.norm(c2[:min_len])
                sim_struct = (c1[:min_len] @ c2[:min_len]) / (n1 * n2 + 1e-12)
            else:
                sim_struct = 0
                
            sim_cent = 0
            if valid_k > 0 and valid_k_last > 0:
                S_mat = compute_sim_matrix(centers[:valid_k], centers_last[:valid_k_last], algo_key, idx_pos, idx_neg)
                sim_cent = np.mean(np.max(S_mat, axis=1)) 
            
            iter_sim = (sim_struct + sim_cent) / 2
            
        iter_sim_stack[0] = iter_sim
        outliers = np.sum(assignments == -1)
        iter_data[it] =[iter_sim, valid_k, outliers]
        
        curr_obj = calculate_objective(data_matrix, centers, counts, assignments, valid_k, strategy, algo_key, idx_pos, idx_neg)
        if curr_obj > best_obj:
            best_obj = curr_obj
            best_state = {'c': centers[:valid_k].copy(), 'n': counts[:valid_k].copy(), 'a': assignments.copy(), 'k': valid_k}

        sim_history.append(iter_sim)
        if len(sim_history) > cycle_window:
            recent = sim_history[-cycle_window:-1]
            if any(abs(r - iter_sim) < 1e-7 for r in recent) and iter_sim < 0.999999:
                log("\n>>> Converged! (Limited Cycle Detected)")
                converged = True
                cycle_detected = True
        
        if iter_sim > 0.9999:
            no_improve_count += 1
        else:
            no_improve_count = 0
            
        if no_improve_count >= early_stop_limit:
            log("\n>>> Converged! (Early Stopping)")
            converged = True
            
        log(f"\tOutlier count: {outliers}")
        log(f"\tCompleted in {time.time()-iter_start:.2f}s")
        log(f"\tInter-iteration similarity: {iter_sim*100:.5f}%")
        
        if converged:
            if cycle_detected:
                log(f"Restoring Best State from Cycle History (Objective: {best_obj:.4f})...")
                centers[:best_state['k']] = best_state['c']
                counts[:best_state['k']] = best_state['n']
                assignments = best_state['a']
                valid_k = best_state['k']
            break

    centers = centers[:valid_k]
    counts = counts[:valid_k]
    iter_data = iter_data[:it+1]
    
    log("FASC End =========")
    log(f"Total Time: {time.time()-start_time:.2f}s")
    
    return centers, counts, assignments, iter_data, it+1