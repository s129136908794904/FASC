function [clusterCenterMtx, clusterListCounts, ptcClusterIdx, iterInfo] = FASC(...
    dataMatrix, idx_pos, idx_neg, simThres_inter, simThres_inner, initialClusterCountLimit, ...
    clusterCountLimit, iterLimit, strategy, ...
    minClusterVolume, similarityAlgorithm)

rng(42);
fprintf('=========FASC Start\n\n')

time_clock_start = clock;
time_CPU_start = cputime;
ptcCount = size(dataMatrix, 1);
iterStackVolume = 1;

similarityAlgorithm = string(similarityAlgorithm);
similarityAlgorithmKey = lower(similarityAlgorithm);
isDualCosine = similarityAlgorithmKey == "dual-cosine";
useUnitCenters = isDualCosine;
mergeStrategy = "all iter";

% Pre-normalization for Dual-Cosine to speed up objective calculation
if isDualCosine
    dataMatrix(:,idx_pos) = dataMatrix(:,idx_pos) ./ vecnorm(dataMatrix(:,idx_pos), 2, 2);
    dataMatrix(:,idx_neg) = dataMatrix(:,idx_neg) ./ vecnorm(dataMatrix(:,idx_neg), 2, 2);
end

initialSize = 150;
clusterListCounts = zeros(initialSize, 1);
clusterCenterMtx = zeros(initialSize, size(dataMatrix, 2), 'like', dataMatrix);
ptcClusterIdx = zeros(ptcCount, 1)-1;

validClusterCount = 0;
clusterCenterMtx(1, :) = dataMatrix(1, :);

iterSimStack = zeros(iterStackVolume, 1);
iterData = zeros(iterLimit, 3);
iteration = 0;

% === [Early Stopping Initialization] ===
best_iterSim_for_stop = -1; 
no_improvement_count = 0;   
early_stop_limit = 100;     

% === [Cycle Protection & Optimization Initialization] ===
cycle_check_window = 20;
sim_history_log = []; 

% [MODIFIED] Track global best based on Objective Function (Psi or L), not just stability
global_best_objective = -inf; 
global_best_state.clusterCenterMtx = [];
global_best_state.clusterListCounts = [];
global_best_state.ptcClusterIdx = [];
global_best_state.validClusterCount = 0;
cycle_detected_flag = false;

while iteration < iterLimit
    this_time_clock_start = clock;
    
    clusterListCounts_last = clusterListCounts;
    meaningfulClusterCount_last = sum(clusterListCounts >= minClusterVolume);
    clusterCenterMtx_last = clusterCenterMtx;       
    validClusterCount_last = validClusterCount;     

    iteration = iteration + 1;
    fprintf('Iteration %d start \n', iteration);

    % ================= Clustering Assignment Step =================
    if iteration == 1
        this_idx_random = randperm(ptcCount);
        for i  = this_idx_random
            if isDualCosine
                this_sims_pos = clusterCenterMtx(1:validClusterCount,idx_pos) * dataMatrix(i,idx_pos)';
                this_sims_neg = clusterCenterMtx(1:validClusterCount,idx_neg) * dataMatrix(i,idx_neg)';
            else
                this_sims_full = Matrix_SimilarityAlgorithms(...
                    clusterCenterMtx(1:validClusterCount,:), dataMatrix(i,:), ...
                    similarityAlgorithmKey, 'Parallel', false);
                this_sims_pos = this_sims_full;
                this_sims_neg = this_sims_full;
            end
            idx_found = CandidateFinder_dual( ...
                strategy, this_sims_pos, this_sims_neg, clusterListCounts,validClusterCount,simThres_inner);

            if idx_found ~= 0
                this_cCenter = clusterCenterMtx(idx_found,:);
                existingCount = clusterListCounts(idx_found);
                this_cCenter = this_cCenter * existingCount + dataMatrix(i, :);
                newCount = existingCount + 1;
                if useUnitCenters
                    clusterCenterMtx(idx_found, idx_pos) = this_cCenter(idx_pos) / (norm(this_cCenter(idx_pos)));
                    clusterCenterMtx(idx_found, idx_neg) = this_cCenter(idx_neg) / (norm(this_cCenter(idx_neg)));
                else
                    clusterCenterMtx(idx_found, :) = this_cCenter / newCount;
                end
                clusterListCounts(idx_found) = newCount;
                ptcClusterIdx(i) = idx_found;
            else
                if validClusterCount < initialClusterCountLimit
                    if validClusterCount == size(clusterCenterMtx, 1)
                        newSize = min(ceil(validClusterCount * 2), initialClusterCountLimit);
                        clusterCenterMtx(end+1:newSize, :) = zeros(newSize - size(clusterCenterMtx,1), size(clusterCenterMtx,2), 'like', dataMatrix);
                        clusterListCounts(end+1:newSize) = 0;
                    end
                    validClusterCount = validClusterCount + 1;
                    clusterCenterMtx(validClusterCount, :) = dataMatrix(i, :);
                    clusterListCounts(validClusterCount) = 1;
                    ptcClusterIdx(i) = validClusterCount;
                else
                    ptcClusterIdx(i) = 0;
                end
            end
        end
        fprintf('Initial cluster span in %.2fs\n', etime(clock, this_time_clock_start));
        ptcIdx_fill = [];

    else % iteration > 1
        activeCenters = clusterCenterMtx(1:validClusterCount, :);

        if isDualCosine
            simMatrix_pos = dataMatrix(:,idx_pos) * activeCenters(:,idx_pos)';
            simMatrix_neg = dataMatrix(:,idx_neg) * activeCenters(:,idx_neg)';
        else
            simMatrix_full = Matrix_SimilarityAlgorithms(...
                dataMatrix, activeCenters, similarityAlgorithmKey);
            simMatrix_pos = simMatrix_full;
            simMatrix_neg = simMatrix_full;
        end

        outers = true(ptcCount,1);
        parfor i  = 1 : ptcCount
            this_sims_pos = simMatrix_pos(i,:)';
            this_sims_neg = simMatrix_neg(i,:)';

            idx_found = CandidateFinder_dual( ...
                strategy, this_sims_pos, this_sims_neg, clusterListCounts,validClusterCount,simThres_inner);

            if idx_found ~= 0
                ptcClusterIdx(i) = idx_found;
                outers(i) = false;
            else
                ptcClusterIdx(i) = 0;
            end

        end

        outerCount = sum(outers);
        holeCount = clusterCountLimit - validClusterCount;

        if holeCount > 0 && outerCount > 0
            if  outerCount >= holeCount
                fillCount = holeCount;
            else
                fillCount = outerCount;
            end
            idx_outer = find(ptcClusterIdx == 0);
            this_idx_random = randperm(length(idx_outer));
            ptcIdx_fill = idx_outer(this_idx_random(1:fillCount));
        else
            ptcIdx_fill = [];
        end
    end

    % ================= Post-Iteration Processing (Merge/Kill) =================
    [clusterCenterMtx, clusterListCounts, validClusterCount, ptcClusterIdx] = ...
        FASC_dual_postIterationProcessor_parallel(...
        clusterCenterMtx, clusterListCounts, validClusterCount, dataMatrix, ...
        ptcClusterIdx, iteration, simThres_inter, clusterCountLimit, minClusterVolume, ...
        mergeStrategy, similarityAlgorithm, ptcIdx_fill,idx_pos,idx_neg,useUnitCenters);

    % ================= Convergence Check (Phase 3 in text) =================
    
    % 1. Check Stability (Eq 3 & 4 in text): Monitors inter-iteration similarity
    [iterSimStack, iterData, converged] = FASC_checkConvergence_parallel(...
        iterSimStack, iterData, iteration, ...
        clusterCenterMtx, clusterCenterMtx_last, ...   
        clusterListCounts, clusterListCounts_last, ... 
        validClusterCount, validClusterCount_last, ... 
        meaningfulClusterCount_last, ...
        ptcClusterIdx, ... 
        minClusterVolume, ...
        similarityAlgorithmKey, idx_pos, idx_neg);

    current_iterSim = iterSimStack(end);

    % 2. Calculate Optimization Objective (Psi or Lyapunov)
    % [MODIFIED] Calculate the Energy/Potential of the current state
    current_Objective = FASC_calculate_Objective(...
        dataMatrix, clusterCenterMtx, clusterListCounts, ptcClusterIdx, ...
        validClusterCount, strategy, similarityAlgorithmKey, idx_pos, idx_neg);

    % [MODIFIED] Update Global Best based on Objective Function (Optimization Goal)
    % This ensures that if we hit a limit cycle, we pick the state with best Potential/Energy
    if current_Objective > global_best_objective
        global_best_objective = current_Objective;
        global_best_state.clusterCenterMtx = clusterCenterMtx;
        global_best_state.clusterListCounts = clusterListCounts;
        global_best_state.ptcClusterIdx = ptcClusterIdx;
        global_best_state.validClusterCount = validClusterCount;
        % fprintf('\tNew Best Objective: %.4f\n', current_Objective); % Debug
    end

    % 3. Cycle Detection Logic
    sim_history_log = [sim_history_log; current_iterSim];
    if length(sim_history_log) > cycle_check_window
        recent_history = sim_history_log(end-cycle_check_window : end-1);
        match_idx = find(abs(recent_history - current_iterSim) < 1e-7, 1);
        
        % If we see a repeated stability score, it indicates a cycle
        if ~isempty(match_idx) && current_iterSim < 0.9999999
            dist_back = length(recent_history) - match_idx + 1;
            if dist_back > 1 
                fprintf('\n>>> Converged! (Limited Cycle Detected)\n');
                cycle_detected_flag = true;
                converged = true; 
            end
        end
    end

    % 4. Early Stopping based on Stability
    if current_iterSim > best_iterSim_for_stop
        best_iterSim_for_stop = current_iterSim;
        no_improvement_count = 0;
    else
        no_improvement_count = no_improvement_count + 1;
    end
    
    if no_improvement_count >= early_stop_limit
        fprintf('\n>>> Converged! (Early Stopping)\n');
        converged = true; 
    end

    fprintf('\tOutlier count: %d\n', iterData(iteration,3));
    fprintf('\tCompleted in %.2fs\n', etime(clock, this_time_clock_start));
    fprintf('\tInter iteration similarity: %.5f%%\n', current_iterSim*100); 
    
    if converged
        if cycle_detected_flag
             % [MODIFIED] Restore the state with the Optimal Objective (Energy/Potential)
             fprintf('Restoring Best State from Cycle History (Objective: %.4f)...\n', global_best_objective);
             clusterCenterMtx = global_best_state.clusterCenterMtx;
             clusterListCounts = global_best_state.clusterListCounts;
             ptcClusterIdx = global_best_state.ptcClusterIdx;
             validClusterCount = global_best_state.validClusterCount;
        elseif no_improvement_count < early_stop_limit
            fprintf('Converged!\n');
        end
        break;
    end

end

validClusterCount = validClusterCount - sum(...
    clusterListCounts < minClusterVolume & clusterListCounts > 0);

[clusterCenterMtx, clusterListCounts, iterInfo, ptcClusterIdx] = FASC_finalizeResults(...
    clusterCenterMtx, clusterListCounts, ptcClusterIdx, iterData, validClusterCount, ...
    time_clock_start, time_CPU_start, iteration,iterStackVolume);

fprintf('FASC End=========\n');

end

% ================= New Helper Function for Phase 3 Logic =================

function objScore = FASC_calculate_Objective(dataMatrix, clusterCenterMtx, clusterListCounts, ...
    ptcClusterIdx, validClusterCount, strategy, similarityAlgorithmKey, idx_pos, idx_neg)
% Calculates the optimization objective defined in Eq (1) or (2).
% SF:  Maximizes Sum Similarity (Equivalent to minimizing Eq 1 Lyapunov)
% DASS: Maximizes Eq 2 Potential (Sum Sim + Density Potential)

    % 1. Calculate Sum of Similarities (Internal Cohesion)
    % Identify valid assignments
    msk_assigned = ptcClusterIdx > 0 & ptcClusterIdx <= validClusterCount;
    if ~any(msk_assigned)
        objScore = -inf;
        return;
    end
    
    assigned_ptcs = dataMatrix(msk_assigned, :);
    assigned_centers = clusterCenterMtx(ptcClusterIdx(msk_assigned), :);
    
    % Compute similarity between particles and their assigned centers
    if similarityAlgorithmKey == "dual-cosine"
        % Fast path for Dual Cosine (assuming pre-normalized or unit centers)
        % Note: dataMatrix is pre-normalized in main, centers are maintained normalized
        sims_pos = dot(assigned_ptcs(:, idx_pos), assigned_centers(:, idx_pos), 2);
        sims_neg = dot(assigned_ptcs(:, idx_neg), assigned_centers(:, idx_neg), 2);
        % Use averaged similarity for objective
        similarities = (sims_pos + sims_neg) / 2;
    else
        % Generic path for other algorithms
        % Calculate row-wise similarity efficiently
        nAssigned = size(assigned_ptcs, 1);
        similarities = zeros(nAssigned, 1);
        % Chunking to prevent memory issues if massive
        chunkSize = 10000;
        for i = 1:chunkSize:nAssigned
            idx_end = min(i+chunkSize-1, nAssigned);
            block_ptc = assigned_ptcs(i:idx_end, :);
            block_ctr = assigned_centers(i:idx_end, :);
            
            % Since SimilarityAlgorithms is vector-based or matrix-matrix, 
            % we need a row-wise variant. For now, using loop or simplified metrics.
            % Most metrics in FASC are dot-product based or distance based.
            if contains(similarityAlgorithmKey, ["cosine", "euclidean", "l1"])
                 similarities(i:idx_end) = row_wise_metric(block_ptc, block_ctr, similarityAlgorithmKey);
            else
                 % Fallback to single call (slow but compatible)
                 for k = i:idx_end
                     similarities(k) = SimilarityAlgorithms(assigned_ptcs(k,:), assigned_centers(k,:), similarityAlgorithmKey);
                 end
            end
        end
    end
    
    sum_similarity = sum(similarities);
    
    % 2. Calculate Strategy-Specific Objective
    strategyKey = upper(string(strategy));
    
    if strategyKey == "DASS" || strategyKey == "DENSITYFIRST"
        % Eq (2): Psi = Sum(Sim) + Sum(n*(n-1)/2)
        % The second term represents the density potential
        counts = clusterListCounts(1:validClusterCount);
        density_potential = sum(counts .* (counts - 1) / 2);
        objScore = sum_similarity + density_potential;
    else
        % Eq (1): L_phi ~ Sum(D_phi). Minimizing divergence <=> Maximizing Similarity
        % for standard bounded metrics.
        objScore = sum_similarity;
    end
    
    function s = row_wise_metric(A, B, method)
       switch method
           case 'cosine'
               s = dot(A, B, 2) ./ (vecnorm(A,2,2).*vecnorm(B,2,2));
           case {'euclidean', 'euclidean-distance'}
               d = vecnorm(A - B, 2, 2);
               s = 1 - d/sqrt(2);
           case {'l1', 'l1-norm', 'manhattan'}
               d = sum(abs(A-B), 2);
               s = 1 - d/2;
           otherwise
               s = dot(A, B, 2); % Fallback
       end
    end

end


% ================= Existing Helper functions (Unchanged) =================

function idx_found = CandidateFinder_dual(strategy, sims_pos, sims_neg,clusterListCounts,validClusterCount,simThresh)
strategyKey = upper(string(strategy));
switch strategyKey
    case {"SIMFIRST","SF"}
        msk_found_pos = sims_pos >= simThresh; 
        msk_found_neg = sims_neg >= simThresh; 
        idx_found = find(msk_found_pos & msk_found_neg);
        if ~isempty(idx_found)
            sims_sum = sims_pos + sims_neg;
            sims_sum = sims_sum(idx_found);
            [~, idx_found_sum] = max(sims_sum);
            idx_found = idx_found(idx_found_sum);
        else
            idx_found = 0;
        end
    case {"DENSITYFIRST","DASS"}
        clusterListCounts_valid = clusterListCounts(1:validClusterCount);
        K = find(sims_pos >= simThresh & sims_neg >= simThresh);
        if ~isempty(K)
            sims = sims_pos + sims_neg;
            sims_K = sims(K);
            clusterListCounts_valid_K = clusterListCounts_valid(K);
            [~, idx_found] = max(sims_K + clusterListCounts_valid_K);
            if length(idx_found) > 1
                idx_found = idx_found(1);
            end
            idx_found = K(idx_found);
        else
            idx_found = 0;
        end
end
end

function [clusterCenterMtx, clusterListCounts, validClusterCount, ptcClusterIdx] = ...
    FASC_dual_postIterationProcessor_parallel(clusterCenterMtx, clusterListCounts, validClusterCount, ...
    dataMatrix, ptcClusterIdx, iteration, simThres_inter, clusterCountLimit, ...
    minClusterVolume, mergeStrategy, similarityAlgorithm, ptcIdx_fill,idx_pos,idx_neg,useUnitCenters)

[clusterCenterMtx, clusterListCounts, ptcClusterIdx, validClusterCount] = ...
    FASC_dual_clusterKiller(dataMatrix, clusterCenterMtx, clusterListCounts, ...
    ptcClusterIdx, validClusterCount, minClusterVolume,idx_pos,idx_neg,useUnitCenters);

if validClusterCount == 0
    return
end

[clusterCenterMtx, clusterListCounts, ptcClusterIdx] = ...
    Clustering_resultSorter(clusterCenterMtx, clusterListCounts, ptcClusterIdx, validClusterCount);

fillCount = length(ptcIdx_fill);
validClusterCount_before = validClusterCount;
validClusterCount = validClusterCount + fillCount;
clusterCenterMtx(validClusterCount_before + 1 : validClusterCount, :) = dataMatrix(ptcIdx_fill,:);
clusterListCounts(validClusterCount_before + 1 : validClusterCount) = 1;
ptcClusterIdx(ptcIdx_fill) = validClusterCount_before + 1 : validClusterCount;

switch mergeStrategy
    case "first iter"
        if iteration == 1
            [clusterCenterMtx, clusterListCounts, validClusterCount, ptcClusterIdx] = ...
                Clustering_FASC_dual_clusterMerger(dataMatrix, clusterCenterMtx, clusterListCounts, ...
                simThres_inter, clusterCountLimit, ptcClusterIdx, validClusterCount, similarityAlgorithm,idx_pos,idx_neg);
            [clusterCenterMtx, clusterListCounts, ptcClusterIdx] = ...
                Clustering_resultSorter(clusterCenterMtx, clusterListCounts, ptcClusterIdx, validClusterCount);
        end
    case "all iter"
        [clusterCenterMtx, clusterListCounts, validClusterCount, ptcClusterIdx] = ...
            Clustering_FASC_dual_clusterMerger(dataMatrix, clusterCenterMtx, clusterListCounts, ...
            simThres_inter, clusterCountLimit, ptcClusterIdx, validClusterCount, similarityAlgorithm,idx_pos,idx_neg);
        [clusterCenterMtx, clusterListCounts, ptcClusterIdx] = ...
            Clustering_resultSorter(clusterCenterMtx, clusterListCounts, ptcClusterIdx, validClusterCount);
end

if iteration == 1
    if validClusterCount <= clusterCountLimit
        clusterCenterMtx_next = zeros(clusterCountLimit,size(clusterCenterMtx,2), 'like', dataMatrix);
        clusterCenterMtx_next(1:validClusterCount,:) = clusterCenterMtx(1:validClusterCount,:);
        clusterCenterMtx = clusterCenterMtx_next;
        clusterListCounts_next = zeros(clusterCountLimit,1);
        clusterListCounts_next(1:validClusterCount) = clusterListCounts(1:validClusterCount);
        clusterListCounts = clusterListCounts_next;
    else
        clusterCenterMtx = clusterCenterMtx(1:clusterCountLimit, :);
        clusterListCounts = clusterListCounts(1:clusterCountLimit);
        validClusterCount = clusterCountLimit;
    end
end
end

function [clusterCenterMtx, clusterListCounts, ptcClusterIdx, validClusterCount] = ...
    FASC_dual_clusterKiller(dataMatrix, clusterCenterMtx, clusterListCounts, ptcClusterIdx,...
    validClusterCount, minClusterVolume,idx_pos,idx_neg,useUnitCenters)
msk_remainedID = true(validClusterCount,1);
lastRoundIDs = 1 : validClusterCount;
for i = 1 : validClusterCount
    msk_thisCluster = (ptcClusterIdx == i);
    clusterListCounts(i) = sum(msk_thisCluster);
    if clusterListCounts(i) >= minClusterVolume
        clusterCenterMtx(i, :) = mean(dataMatrix(msk_thisCluster, :), 1);
        if useUnitCenters
            clusterCenterMtx(i, idx_pos) = clusterCenterMtx(i, idx_pos) / (norm(clusterCenterMtx(i, idx_pos)));
            clusterCenterMtx(i, idx_neg) = clusterCenterMtx(i, idx_neg) / (norm(clusterCenterMtx(i, idx_neg)));
        end
    else
        clusterCenterMtx(i,:) = 0;
        msk_remainedID(i) = false;
        clusterListCounts(i) = 0;
        ptcClusterIdx(msk_thisCluster) = 0;
    end
end
validClusterCount = sum(msk_remainedID);
clusterCenterMtx(1:validClusterCount,:) = clusterCenterMtx(msk_remainedID,:);
clusterListCounts(1:validClusterCount) = clusterListCounts(msk_remainedID);
remainedIDs = lastRoundIDs(msk_remainedID);
newIDs = 1 : length(remainedIDs);
for i = 1 : length(remainedIDs)
    if remainedIDs(i) ~= newIDs(i)
        ptcClusterIdx(ptcClusterIdx == remainedIDs(i)) = newIDs(i);
    end
end
end

function [clusterCenterMtx, clusterListCounts, ptcClusterIdx] = ...
    Clustering_resultSorter(clusterCenterMtx, clusterListCounts, ptcClusterIdx, validClusterCount)
[~, Msk_sortID] = sort(clusterListCounts(1:validClusterCount), 'descend');
Msk_sort = [Msk_sortID; (validClusterCount+1:1:length(clusterListCounts))'];
clusterListCounts = clusterListCounts(Msk_sort);
clusterCenterMtx = clusterCenterMtx(Msk_sort, :);
ptcClusterIdx_sort = zeros(size(ptcClusterIdx));
for i = 1 : validClusterCount
    ptcClusterIdx_sort(ptcClusterIdx == Msk_sortID(i)) = i;
end
ptcClusterIdx = ptcClusterIdx_sort;
end

function [clusterCenterMtx, clusterListCounts, validClusterCount, ptcClusterIdx] = ...
    Clustering_FASC_dual_clusterMerger(dataMatrix, clusterCenterMtx, clusterListCounts, sim_inter, ...
    clusterCountLimit, ptcClusterIdx,validClusterCount, similarityAlgorithm,idx_pos,idx_neg)
similarityAlgorithm = string(similarityAlgorithm);
similarityAlgorithmKey = lower(similarityAlgorithm);
isDualCosine = similarityAlgorithmKey == "dual-cosine";
validClusterCount_ini = validClusterCount;
fprintf('\tCurrent cluster count: %d \n', validClusterCount_ini);
msk_remain = true(validClusterCount, 1);
idx_thisC = 1;
while true
    idx_remain = find(msk_remain);
    idx_compute = idx_remain(idx_thisC : end);
    this_Center = clusterCenterMtx(idx_compute(1), :);
    if isDualCosine
        sims_pos = clusterCenterMtx(idx_compute, idx_pos) * this_Center(idx_pos)';
        sims_neg = clusterCenterMtx(idx_compute, idx_neg) * this_Center(idx_neg)';
    else
        sims_full = Matrix_SimilarityAlgorithms(clusterCenterMtx(idx_compute, :), this_Center, similarityAlgorithmKey);
        sims_pos = sims_full;
        sims_neg = sims_full;
    end
    msk_merge = sims_pos >= sim_inter & sims_neg >= sim_inter;
    idx_merge = idx_compute(msk_merge);
    if numel(idx_merge) > 1
        mergedCount = sum(clusterListCounts(idx_merge));
        weightedCenters = clusterCenterMtx(idx_merge, :) .* clusterListCounts(idx_merge);
        newCenter = sum(weightedCenters, 1);
        if isDualCosine
            newCenter(idx_pos) = newCenter(idx_pos) / (norm(newCenter(idx_pos)));
            newCenter(idx_neg) = newCenter(idx_neg) / (norm(newCenter(idx_neg)));
        else
            newCenter = newCenter / mergedCount;
        end
        clusterCenterMtx(idx_merge(1), :) = newCenter;
        clusterCenterMtx(idx_merge(2:end), :) = 0;
        msk_remain(idx_merge(2:end)) = false;
        clusterListCounts(idx_merge(1)) = mergedCount;
        clusterListCounts(idx_merge(2:end)) = 0;
        ptcClusterIdx(logical(sum(ptcClusterIdx' == idx_merge))) = idx_merge(1);
        idx_thisC = 1;
    else
        idx_thisC = idx_thisC + 1;
    end
    if ~((idx_thisC <= clusterCountLimit) && (idx_thisC < length(idx_remain)))
        break
    end
    if idx_thisC > clusterCountLimit || idx_thisC > length(idx_remain)
        break
    end
end
msk_remain = clusterListCounts(1:validClusterCount_ini) > 0;
remainedClusterCount = sum(msk_remain);
if remainedClusterCount == validClusterCount_ini
    fprintf('\tNo cluster merged \n');
else
    if clusterCountLimit > remainedClusterCount
        validClusterCount = remainedClusterCount;
        remainedIDs = unique(ptcClusterIdx);
        remainedIDs = sort(remainedIDs(remainedIDs > 0));
        for i = 1 : length(remainedIDs)
            if remainedIDs(i) ~= i
                ptcClusterIdx(ptcClusterIdx == remainedIDs(i)) = i;
            end
        end
        clusterCenterMtx(1:remainedClusterCount,:) = clusterCenterMtx(msk_remain,:);
        clusterListCounts(1:remainedClusterCount) = clusterListCounts(msk_remain);
        clusterCenterMtx(remainedClusterCount+1:end,:) = 0;
        clusterListCounts(remainedClusterCount+1:end) = 0;
    else
        validClusterCount = clusterCountLimit;
        remainedIDs = unique(ptcClusterIdx);
        remainedIDs = sort(remainedIDs(remainedIDs > 0));
        for i = 1 : validClusterCount
            if remainedIDs(i) ~= i
                ptcClusterIdx(ptcClusterIdx == remainedIDs(i)) = i;
            end
        end
        ptcClusterIdx(ptcClusterIdx > validClusterCount) = 0;
        clusterCenterMtx(1:remainedClusterCount,:) = clusterCenterMtx(msk_remain,:);
        clusterListCounts(1:remainedClusterCount) = clusterListCounts(msk_remain);
        clusterCenterMtx(validClusterCount+1:end,:) = 0;
        clusterListCounts(validClusterCount+1:end) = 0;
    end
    for i = 1 : validClusterCount
        msk_thisCluster = (ptcClusterIdx == i);
        clusterListCounts(i) = sum(msk_thisCluster);
        clusterCenterMtx(i, :) = mean(dataMatrix(msk_thisCluster, :), 1);
        if isDualCosine
            clusterCenterMtx(i, idx_pos) = clusterCenterMtx(i, idx_pos) / (norm(clusterCenterMtx(i, idx_pos)));
            clusterCenterMtx(i, idx_neg) = clusterCenterMtx(i, idx_neg) / (norm(clusterCenterMtx(i, idx_neg)));
        end
    end
    fprintf('\tMerged cluster count: %d \n', remainedClusterCount);
end
end

function [iterSimStack, iterData, converged] = FASC_checkConvergence_parallel(...
    iterSimStack, iterData, iteration, ...
    clusterCenterMtx, clusterCenterMtx_last, ...
    clusterListCounts, clusterListCounts_last, ...
    validClusterCount, validClusterCount_last, ...
    meaningfulClusterCount_last, ...
    ptcClusterIdx, ...
    minClusterVolume, ... 
    similarityAlgorithmKey, idx_pos, idx_neg)

iterSim = 0;
sim_structure = 0;
sim_centroid = 0;

if iteration > 1
    meaningfulClusterCount = sum(clusterListCounts >= minClusterVolume);
    iterCountsMtx = zeros(max([meaningfulClusterCount_last meaningfulClusterCount]), 2);
    iterCountsMtx(1:meaningfulClusterCount_last, 1) = clusterListCounts_last(1:meaningfulClusterCount_last);
    iterCountsMtx(1:meaningfulClusterCount, 2) = clusterListCounts(1:meaningfulClusterCount);
    
    norm1 = norm(iterCountsMtx(:, 1));
    norm2 = norm(iterCountsMtx(:, 2));
    if norm1 > 0 && norm2 > 0
        sim_structure = (iterCountsMtx(:, 1)' * iterCountsMtx(:, 2)) / (norm1 * norm2);
    else
        sim_structure = 0;
    end

    if validClusterCount > 0 && validClusterCount_last > 0
        centers_curr = clusterCenterMtx(1:validClusterCount, :);
        centers_last = clusterCenterMtx_last(1:validClusterCount_last, :);
        
        if similarityAlgorithmKey == "dual-cosine"
            sims_pos = centers_curr(:, idx_pos) * centers_last(:, idx_pos)';
            sims_neg = centers_curr(:, idx_neg) * centers_last(:, idx_neg)';
            simMatrix = (sims_pos + sims_neg) / 2;
        else
            simMatrix = Matrix_SimilarityAlgorithms(centers_curr, centers_last, similarityAlgorithmKey, 'Parallel', false);
        end
        
        num_matches = min(validClusterCount, validClusterCount_last);
        matched_sims = zeros(num_matches, 1);
        matched_counts = zeros(num_matches, 1);
        
        [sortedSims, sortIdx] = sort(simMatrix(:), 'descend');
        [rows, cols] = ind2sub(size(simMatrix), sortIdx);
        
        used_curr = false(validClusterCount, 1);
        used_last = false(validClusterCount_last, 1);
        match_count = 0;
        
        for k = 1:length(sortedSims)
            r = rows(k);
            c = cols(k);
            if ~used_curr(r) && ~used_last(c)
                match_count = match_count + 1;
                matched_sims(match_count) = sortedSims(k);
                matched_counts(match_count) = clusterListCounts(r); 
                used_curr(r) = true;
                used_last(c) = true;
                if match_count >= num_matches
                    break;
                end
            end
        end
        total_weight = sum(matched_counts);
        if total_weight > 0
            sim_centroid = sum(matched_sims .* matched_counts) / total_weight;
        else
            sim_centroid = 0;
        end
    else
        sim_centroid = 0;
    end
    iterSim = (sim_structure + sim_centroid) / 2;
end

iterSimStack = circshift(iterSimStack, -1);
iterSimStack(end) = iterSim;
outlierCount = sum(ptcClusterIdx==0);
iterData(iteration, :) = [iterSim, validClusterCount,outlierCount];

thres_iterSim = 1.0;
tol_iterSim = 1e-6;
effectiveThres_iterSim = thres_iterSim - tol_iterSim;
converged = all(iterSimStack >= effectiveThres_iterSim);
end

function [clusterCenterMtx, clusterListCounts, iterInfo, ptcClusterIdx] = FASC_finalizeResults(...
    clusterCenterMtx, clusterListCounts, ptcClusterIdx, iterData, validClusterCount, ...
    time_clock_start, time_CPU_start, iteration, iterStackVolume)
clusterCenterMtx = clusterCenterMtx(1:validClusterCount, :);
clusterListCounts = clusterListCounts(1:validClusterCount, :);
ptcClusterIdx(ptcClusterIdx > validClusterCount) = 0;
iterInfo.iterData = iterData(1:iteration, :);
iterInfo.convergeIter = iteration - (iterStackVolume-1) * (iteration >= iterStackVolume);
iterInfo.outLiersCount = sum(ptcClusterIdx == 0);
time_clock_end = clock;
time_CPU_end = cputime;
fprintf('\nTotal clock time: %.2fs\nCPU time: %.2fs\n', ...
    etime(time_clock_end, time_clock_start), time_CPU_end - time_CPU_start);
fprintf('\nOutliers count: %d\n', ...
    iterInfo.outLiersCount);
end

% Matrix_SimilarityAlgorithms
function sims = Matrix_SimilarityAlgorithms(mtx1, mtx2, algorithmName, varargin)
p = inputParser;
addParameter(p, 'Order', 20, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Parallel', true, @islogical);
addParameter(p, 'BatchSize', 100000, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Normalize', true, @islogical);
parse(p, varargin{:});
order = p.Results.Order;
useParallel = p.Results.Parallel;
batchSize = p.Results.BatchSize;
doNormalize = p.Results.Normalize;
rowCount = size(mtx1, 1);
colCount = size(mtx2, 1);
algorithmKey = lower(string(algorithmName));
if isa(mtx1,'single') || isa(mtx2,'single')
    mtx1 = single(mtx1);
    mtx2 = single(mtx2);
end

switch algorithmKey
    case 'cosine'
        sims = fast_cosine_similarity(mtx1, mtx2, doNormalize);
        return;
    case {'euclidean-distance','euclidean'}
        sims = fast_euclidean_similarity(mtx1, mtx2, doNormalize);
        return;
    case {'l1-norm','l1','manhattan'}
        sims = fast_l1_similarity(mtx1, mtx2, doNormalize);
        return;
    case 'minimum'
        sims = fast_minimum_similarity(mtx1, mtx2, doNormalize);
        return;
    case 'maximum'
        sims = fast_maximum_similarity(mtx1, mtx2, doNormalize);
        return;
end
if doNormalize
    mtx1_norm = normalize_matrix(mtx1);
    mtx2_norm = normalize_matrix(mtx2);
else
    mtx1_norm = mtx1;
    mtx2_norm = mtx2;
end
sims = zeros(rowCount, colCount);

if useParallel && (rowCount * colCount > 10000)
    sims = parallel_batch_computation(mtx1_norm, mtx2_norm, algorithmKey, order, batchSize);
elseif rowCount * colCount > 1000
    sims = batch_computation(mtx1_norm, mtx2_norm, algorithmKey, order, batchSize);
else
    sims = direct_computation(mtx1_norm, mtx2_norm, algorithmKey, order);
end
end

function sims = fast_cosine_similarity(mtx1, mtx2, doNormalize)

    if doNormalize
        denom1 = vecnorm(mtx1, 2, 2);
        denom2 = vecnorm(mtx2, 2, 2);
        mtx1_norm = mtx1 ./ denom1;
        mtx2_norm = mtx2 ./ denom2;
    else
        norms1 = vecnorm(mtx1, 2, 2);
        norms2 = vecnorm(mtx2, 2, 2);
        mtx1_norm = mtx1 ./ (norms1);
        mtx2_norm = mtx2 ./ (norms2);
    end
    mtx1_norm = cast(mtx1_norm, 'like', mtx1);
    mtx2_norm = cast(mtx2_norm, 'like', mtx2);
    sims = mtx1_norm * mtx2_norm';
end

function sims = fast_euclidean_similarity(mtx1, mtx2, doNormalize)
    if doNormalize
        mtx1 = normalize_matrix(mtx1);
        mtx2 = normalize_matrix(mtx2);
    end
    if isa(mtx1,'single') || isa(mtx2,'single')
        mtx1 = single(mtx1);
        mtx2 = single(mtx2);
    end
    distances = pdist2(mtx1, mtx2, 'euclidean');
    sims = 1 - distances / sqrt(2);
end

function sims = fast_l1_similarity(mtx1, mtx2, doNormalize)
    if doNormalize
        mtx1 = normalize_matrix(mtx1);
        mtx2 = normalize_matrix(mtx2);
    end
    if isa(mtx1,'single') || isa(mtx2,'single')
        mtx1 = single(mtx1);
        mtx2 = single(mtx2);
    end
    distances = pdist2(mtx1, mtx2, 'cityblock');
    sims = 1 - distances / 2;
end

function sims = fast_minimum_similarity(mtx1, mtx2, doNormalize)
    if doNormalize
        mtx1 = normalize_matrix(mtx1);
        mtx2 = normalize_matrix(mtx2);
    end
    rowCount = size(mtx1, 1);
    colCount = size(mtx2, 1);
    sims = zeros(rowCount, colCount);
    for i = 1:rowCount
        combined = cat(3, repmat(mtx1(i, :), colCount, 1), mtx2);
        sims(i, :) = sum(min(combined, [], 3), 2);
    end
end

function sims = fast_maximum_similarity(mtx1, mtx2, doNormalize)
    if doNormalize
        mtx1 = normalize_matrix(mtx1);
        mtx2 = normalize_matrix(mtx2);
    end
    rowCount = size(mtx1, 1);
    colCount = size(mtx2, 1);
    sims = zeros(rowCount, colCount);
    for i = 1:rowCount
        combined = cat(3, repmat(mtx1(i, :), colCount, 1), mtx2);
        sims(i, :) = 1 - sum(max(combined, [], 3), 2) / 2;
    end
end

function mtx_norm = normalize_matrix(mtx)
    origType = class(mtx);
    row_sums = sum(mtx, 2, 'native');
    row_sums(row_sums == 0) = 1;
    mtx_norm = mtx ./ row_sums;
    if ~strcmp(class(mtx_norm), origType)
        mtx_norm = cast(mtx_norm, origType);
    end
end

function sims = parallel_batch_computation(mtx1, mtx2, algorithmName, order, batchSize)
    rowCount = size(mtx1, 1);
    colCount = size(mtx2, 1);
    numBatches = ceil(rowCount / batchSize);
    batchResults = cell(numBatches, 1);
    batchRowIndices = cell(numBatches, 1);
    for batchIdx = 1:numBatches
        startRow = (batchIdx - 1) * batchSize + 1;
        endRow = min(batchIdx * batchSize, rowCount);
        batchRowIndices{batchIdx} = startRow:endRow;
    end
    parfor batchIdx = 1:numBatches
        batchRows = batchRowIndices{batchIdx};
        batchSims = zeros(length(batchRows), colCount);
        for i = 1:length(batchRows)
            for j = 1:colCount
                batchSims(i, j) = SimilarityAlgorithms(...
                    mtx1(batchRows(i), :), mtx2(j, :), algorithmName, order);
            end
        end
        batchResults{batchIdx} = batchSims;
    end
    sims = zeros(rowCount, colCount);
    for batchIdx = 1:numBatches
        batchRows = batchRowIndices{batchIdx};
        sims(batchRows, :) = batchResults{batchIdx};
    end
end

function sims = batch_computation(mtx1, mtx2, algorithmName, order, batchSize)
    rowCount = size(mtx1, 1);
    colCount = size(mtx2, 1);
    sims = zeros(rowCount, colCount);
    for startRow = 1:batchSize:rowCount
        endRow = min(startRow + batchSize - 1, rowCount);
        batchRows = startRow:endRow;
        for i = 1:length(batchRows)
            for j = 1:colCount
                sims(batchRows(i), j) = SimilarityAlgorithms(...
                    mtx1(batchRows(i), :), mtx2(j, :), algorithmName, order);
            end
        end
    end
end

function sims = direct_computation(mtx1, mtx2, algorithmName, order)
    rowCount = size(mtx1, 1);
    colCount = size(mtx2, 1);
    sims = zeros(rowCount, colCount);
    for i = 1:rowCount
        for j = 1:colCount
            sims(i, j) = SimilarityAlgorithms(...
                mtx1(i, :), mtx2(j, :), algorithmName, order);
        end
    end
end

function sim = SimilarityAlgorithms(vector1,vector2,algorithmName)
vector1 = vector1 / sum(vector1);
vector2 = vector2 / sum(vector2);
algorithmKey = lower(string(algorithmName));
switch algorithmKey
    case "algebraic"
        vectLen = length(vector1);
        sim = 0;
        for i = 1 : vectLen
            if vector1(i)&&vector2(i)
                a = vector1(i);
                b = vector2(i);
                sim = sim + (a+b)/2;
            end
        end
    case "cosine"
        sim = vector1 * vector2'/norm(vector1)/norm(vector2);
    case {"euclidean-distance","euclidean"}
        sim = 1-norm(vector1-vector2)/sqrt(2);
    case {"l1-norm","l1","manhattan"}
        sim = 1 - sum(abs(vector1 - vector2)) / 2;
    case "minimum"
        combination = [vector1;vector2];
        sim = sum(min(combination));
    case "maximum"
        combination = [vector1;vector2];
        sim = 1 - sum(max(combination))/2;
    case "logarithmic"
        vectLen = length(vector1);
        sim = 0;
        for i = 1 : vectLen
            if vector1(i)&&vector2(i)
                a = vector1(i);
                b = vector2(i);
                sim = sim + (b-a)/(log(b)-log(a));
            end
        end
    case "geometric"
        vectLen = length(vector1);
        sim = 0;
        for i = 1 : vectLen
            if vector1(i)&&vector2(i)
                a = vector1(i);
                b = vector2(i);
                sim = sim + sqrt(a*b);
            end
        end
    case "harmonic"
        vectLen = length(vector1);
        sim = 0;
        for i = 1 : vectLen
            if vector1(i)&&vector2(i)
                sim = sim + 2*(vector1(i)*vector2(i))/(vector1(i)+vector2(i));
            end
        end
    case "enhanced harmonic"
        vectLen = length(vector1);
        sim = 0;
        for i = 1 : vectLen
            if vector1(i)&&vector2(i)
                sim = sim + sqrt(2 / ((vector1(i)^(-2)) + (vector2(i)^(-2))));
            end
        end
    case "entropy"
        entropy1 = entropyCalculater(vector1);
        entropy2 = entropyCalculater(vector2);
        entropy1P2 = entropyCalculater((vector1+vector2)/2);
        sim = 1 - (2 * entropy1P2 - entropy1 - entropy2) / 1.3863; %ln(4)
    case "weighted entropy"
        entropy1 = entropyCalculater(vector1);
        if entropy1 < 3
            vector1 = vector1.^(0.25 + 0.25 * entropy1);
            vector1 = vector1 / sum(vector1);
            entropy1 = entropyCalculater(vector1);
        end
        entropy2 = entropyCalculater(vector2);
        if entropy2 < 3
            vector2 = vector2.^(0.25 + 0.25 * entropy2);
            vector2 = vector2 / sum(vector2);
            entropy2 = entropyCalculater(vector2);
        end
        entropy1P2 = entropyCalculater((vector1+vector2)/2);
        sim = 1 - (2 * entropy1P2 - entropy1 - entropy2) / 1.3863; %ln(4)
    case "best average"
        vectLen = length(vector1);
        sim = 0;
        for i = 1 : vectLen
            if vector1(i)&&vector2(i)
                a = vector1(i);
                b = vector2(i);
                sim = sim + (a^b*b^a)^(1/(a+b));
            end
        end
    case "fitted core"
        vectLen = length(vector1);
        sim = 0;
        for i = 1 : vectLen
            if vector1(i)&&vector2(i)
                a = vector1(i);
                b = vector2(i);
                sim = sim + (a^b*b^a)^(1/(a+b));
            end
        end
    otherwise
        error('Unsupported similarity algorithm: %s', char(algorithmName));
end
    function entropy = entropyCalculater(vector)
        entropy = 0;
        for j = 1:length(vector)
            if vector(j)
                entropy = entropy - (vector(j) * log(vector(j)));
            end
        end
    end
end