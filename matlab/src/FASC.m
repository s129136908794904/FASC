function [clusterCenterMtx, clusterListCounts, ptcClusterIdx, iterInfo] = FASC(...
    dataMatrix, idx_pos, idx_neg, simThres_inter, simThres_inner, initialClusterCountLimit, ...
    clusterCountLimit, iterLimit, strategy, minClusterVolume, similarityAlgorithm, varargin)
%FASC Flexible Adaptive Stable Clustering.
%
% The first eleven inputs retain the original public API. Optional name-value
% arguments expose the quantities used in the manuscript:
%   Lambda                    Support weight (default: 1)
%   SupportFunction           phi(n) in the DASS rule (default: @(n)n)
%   SupportPotentialFunction  cumulative support reward
%                             (default: n(n-1)/2)
%   DensityFunction           Backward-compatible alias for SupportFunction
%   DensityPotentialFunction  Backward-compatible alias for
%                             SupportPotentialFunction
%   BatchSize                 Maximum rows in a similarity block
%   SimilarityFunction        Custom pairwise similarity function
%   RepresentativeFunction    Custom deterministic representative function
%   CapacitySchedule          Numeric K_t vector or function handle
%                             (cycle detection starts only in a certified
%                             fixed tail; numeric schedules hold their last
%                             value, while function handles are not assumed
%                             to become constant)
%   ObjectiveFunction         Optional cycle-scoring function
%   TrackObjective            Score every iteration instead of cycles only
%   Verbose                   Print iteration diagnostics

opts = FASC_parseOptions(varargin{:});
FASC_validateInputs(dataMatrix, idx_pos, idx_neg, simThres_inter, simThres_inner, ...
    initialClusterCountLimit, clusterCountLimit, iterLimit, minClusterVolume);

strategyKey = upper(string(strategy));
if ~ismember(strategyKey, ["SF", "SIMFIRST", "DASS", "DENSITYFIRST"])
    error('FASC:InvalidStrategy', 'Strategy must be ''SF'' or ''DASS''.');
end

[similarityKey, similarityFcn] = FASC_similaritySpecification(...
    similarityAlgorithm, opts.SimilarityFunction);
isDualCosine = similarityKey == "dual-cosine";
if isDualCosine
    FASC_validateDualIndices(idx_pos, idx_neg, size(dataMatrix, 2));
end

if opts.Verbose
    fprintf('=========FASC Start\n\n');
end
timeClockStart = tic;
timeCPUStart = cputime;

% All assignment passes use the same preprocessed representation. Dual
% cosine is normalized independently by polarity, ordinary cosine uses
% row-wise L2 normalization, and composition scores use absolute-L1
% preprocessing so signed rows are not reversed by an algebraic row sum.
dataWork = FASC_prepareData(dataMatrix, similarityKey, idx_pos, idx_neg);
sampleCount = size(dataWork, 1);

initialCapacity = FASC_capacityAt(opts.CapacitySchedule, 1, clusterCountLimit);
seedCount = min([initialClusterCountLimit, initialCapacity, sampleCount]);
seedIdx = FASC_selectSeeds(dataWork, seedCount, similarityKey, similarityFcn, ...
    idx_pos, idx_neg, opts.BatchSize);

clusterCenterMtx = zeros(seedCount, size(dataWork, 2), 'like', dataWork);
for j = 1:seedCount
    clusterCenterMtx(j, :) = FASC_representative(dataWork(seedIdx(j), :), ...
        similarityKey, similarityFcn, opts.RepresentativeFunction, ...
        idx_pos, idx_neg, opts.BatchSize);
end
clusterListCounts = zeros(seedCount, 1);
ptcClusterIdx = zeros(sampleCount, 1);

iteration = 0;
previousLabels = [];
firstSeen = containers.Map('KeyType', 'char', 'ValueType', 'double');

iterData = zeros(iterLimit, 3);
objectiveHistory = nan(iterLimit, 1);
labelAgreementHistory = nan(iterLimit, 1);
stateHashHistory = cell(iterLimit, 1);

cycleDetected = false;
cycleLength = 0;
cycleRemaining = 0;
cycleBestObjective = -inf;
cycleBestState = struct('centers', [], 'counts', [], 'labels', []);
terminationReason = "iteration-limit";
matureProjectionApplied = false;

while iteration < iterLimit || cycleRemaining > 0
    iterationStart = tic;
    iteration = iteration + 1;
    capacity = FASC_capacityAt(opts.CapacitySchedule, iteration, clusterCountLimit);
    capacityScheduleFixed = FASC_capacityScheduleIsFixed( ...
        opts.CapacitySchedule, iteration);

    if opts.Verbose
        fprintf('Iteration %d start\n', iteration);
    end

    % A decreasing schedule removes the lowest-priority canonical clusters.
    if size(clusterCenterMtx, 1) > capacity
        removed = ptcClusterIdx > capacity;
        ptcClusterIdx(removed) = 0;
        clusterCenterMtx = clusterCenterMtx(1:capacity, :);
        clusterListCounts = clusterListCounts(1:capacity);
    end

    previousCenters = clusterCenterMtx;
    previousCounts = clusterListCounts;
    labelsBefore = ptcClusterIdx;

    % Phase 1: representatives and supports are frozen for the whole pass.
    [ptcClusterIdx, maxAffinity] = FASC_assignInBlocks(dataWork, ...
        clusterCenterMtx, clusterListCounts, strategyKey, simThres_inner, ...
        similarityKey, similarityFcn, idx_pos, idx_neg, opts);

    % Phase 2: promote a deterministic maximin subset of outliers.
    holes = max(0, capacity - size(clusterCenterMtx, 1));
    promotionIdx = FASC_selectPromotions(dataWork, ptcClusterIdx, ...
        maxAffinity, holes, similarityKey, similarityFcn, ...
        opts.RepresentativeFunction, idx_pos, idx_neg, opts.BatchSize);
    protected = false(size(clusterCenterMtx, 1) + numel(promotionIdx), 1);
    if ~isempty(promotionIdx)
        firstNew = size(clusterCenterMtx, 1) + 1;
        ptcClusterIdx(promotionIdx) = (firstNew:firstNew + numel(promotionIdx) - 1)';
        protected(firstNew:end) = true;
    end

    [clusterCenterMtx, clusterListCounts, ptcClusterIdx] = FASC_consolidate(...
        dataWork, ptcClusterIdx, protected, minClusterVolume, simThres_inter, ...
        similarityKey, similarityFcn, opts.RepresentativeFunction, ...
        idx_pos, idx_neg, opts.BatchSize);

    [stability, labelAgreement] = FASC_stability(previousCenters, previousCounts, ...
        clusterCenterMtx, clusterListCounts, labelsBefore, ptcClusterIdx, ...
        similarityKey, similarityFcn, idx_pos, idx_neg, opts.BatchSize);

    currentHash = FASC_stateHash(ptcClusterIdx, capacity);
    repeatedState = cycleRemaining == 0 && capacityScheduleFixed && ...
        isKey(firstSeen, currentHash);
    cycleCandidateState = struct('centers', [], 'counts', [], 'labels', []);
    if cycleRemaining > 0 || repeatedState
        [candidateCenters, candidateCounts, candidateLabels] = ...
            FASC_projectMatureState(clusterCenterMtx, clusterListCounts, ...
            ptcClusterIdx, minClusterVolume);
        matureProjectionApplied = matureProjectionApplied || ...
            ~isequal(candidateLabels, ptcClusterIdx);
        cycleCandidateState = FASC_captureState(candidateCenters, ...
            candidateCounts, candidateLabels);
        currentObjective = FASC_objective(dataWork, candidateCenters, ...
            candidateCounts, candidateLabels, strategyKey, similarityKey, ...
            similarityFcn, idx_pos, idx_neg, opts);
    elseif opts.TrackObjective
        currentObjective = FASC_objective(dataWork, clusterCenterMtx, ...
            clusterListCounts, ptcClusterIdx, strategyKey, similarityKey, ...
            similarityFcn, idx_pos, idx_neg, opts);
    else
        currentObjective = nan;
    end

    iterData(iteration, :) = [stability, size(clusterCenterMtx, 1), ...
        sum(ptcClusterIdx == 0)];
    objectiveHistory(iteration, 1) = currentObjective;
    labelAgreementHistory(iteration, 1) = labelAgreement;
    stateHashHistory{iteration, 1} = currentHash;

    if opts.Verbose
        fprintf('\tCurrent cluster count: %d\n', size(clusterCenterMtx, 1));
        fprintf('\tOutlier count: %d\n', sum(ptcClusterIdx == 0));
        if isfinite(currentObjective)
            fprintf('\tObjective: %.10g\n', currentObjective);
        end
        fprintf('\tInter-iteration similarity: %.8f%%\n', stability * 100);
        fprintf('\tCompleted in %.2fs\n', toc(iterationStart));
    end

    if cycleRemaining > 0
        if currentObjective > cycleBestObjective
            cycleBestObjective = currentObjective;
            cycleBestState = cycleCandidateState;
        end
        cycleRemaining = cycleRemaining - 1;
        if cycleRemaining == 0
            clusterCenterMtx = cycleBestState.centers;
            clusterListCounts = cycleBestState.counts;
            ptcClusterIdx = cycleBestState.labels;
            terminationReason = "limit-cycle";
            break;
        end
    else
        nextCapacity = FASC_capacityAt(opts.CapacitySchedule, iteration + 1, clusterCountLimit);
        isFixedPoint = ~isempty(previousLabels) && ...
            isequal(ptcClusterIdx, previousLabels) && nextCapacity == capacity && ...
            capacityScheduleFixed;

        if isFixedPoint
            terminationReason = "fixed-point";
            break;
        end

        if capacityScheduleFixed
            if isKey(firstSeen, currentHash)
                cycleDetected = true;
                cycleLength = iteration - firstSeen(currentHash);
                cycleBestObjective = currentObjective;
                cycleBestState = cycleCandidateState;
                cycleRemaining = cycleLength - 1;
                if opts.Verbose
                    fprintf(['Repeated label fingerprint detected (period %d); ' ...
                        'replaying cycle.\n'], cycleLength);
                end
                if cycleRemaining == 0
                    terminationReason = "limit-cycle";
                    break;
                end
            else
                firstSeen(currentHash) = iteration;
            end
        end
    end

    previousLabels = ptcClusterIdx;
end

labelsBeforeProjection = ptcClusterIdx;
[clusterCenterMtx, clusterListCounts, ptcClusterIdx] = ...
    FASC_projectMatureState(clusterCenterMtx, clusterListCounts, ...
    ptcClusterIdx, minClusterVolume);
finalProjectionChanged = ~isequal(ptcClusterIdx, labelsBeforeProjection);
matureProjectionApplied = matureProjectionApplied || finalProjectionChanged;

iterInfo.iterData = iterData(1:iteration, :);
iterInfo.objective = objectiveHistory(1:iteration, :);
iterInfo.labelAgreement = labelAgreementHistory(1:iteration, :);
iterInfo.stateHash = stateHashHistory(1:iteration, :);
iterInfo.convergeIter = iteration;
iterInfo.outLiersCount = sum(ptcClusterIdx == 0);
iterInfo.terminationReason = char(terminationReason);
iterInfo.cycleDetected = cycleDetected;
iterInfo.cycleLength = cycleLength;
iterInfo.lambda = opts.Lambda;
iterInfo.batchSize = opts.BatchSize;
iterInfo.finalObjective = FASC_objective(dataWork, clusterCenterMtx, ...
    clusterListCounts, ptcClusterIdx, strategyKey, similarityKey, ...
    similarityFcn, idx_pos, idx_neg, opts);
finalCapacity = FASC_capacityAt(opts.CapacitySchedule, iteration, clusterCountLimit);
iterInfo.finalStateHash = FASC_stateHash(ptcClusterIdx, finalCapacity);
iterInfo.finalProjectionChanged = finalProjectionChanged;
iterInfo.matureProjectionApplied = matureProjectionApplied;
iterInfo.capacityScheduleFixed = FASC_capacityScheduleIsFixed( ...
    opts.CapacitySchedule, iteration);
iterInfo.elapsedTime = toc(timeClockStart);
iterInfo.cpuTime = cputime - timeCPUStart;

if opts.Verbose
    fprintf('\nTermination: %s\n', iterInfo.terminationReason);
    fprintf('Total clock time: %.2fs\nCPU time: %.2fs\n', ...
        iterInfo.elapsedTime, iterInfo.cpuTime);
    fprintf('Outliers count: %d\n', iterInfo.outLiersCount);
    fprintf('FASC End=========\n');
end
end

function opts = FASC_parseOptions(varargin)
p = inputParser;
p.FunctionName = 'FASC';
addParameter(p, 'Lambda', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'SupportFunction', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'SupportPotentialFunction', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'DensityFunction', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'DensityPotentialFunction', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'BatchSize', 100000, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && fix(x) == x);
addParameter(p, 'SimilarityFunction', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'RepresentativeFunction', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'CapacitySchedule', [], ...
    @(x) isempty(x) || isa(x, 'function_handle') || ...
    (isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x >= 1) && all(fix(x) == x)));
addParameter(p, 'ObjectiveFunction', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'TrackObjective', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;
supportSupplied = ~any(strcmpi('SupportFunction', p.UsingDefaults));
legacySupportSupplied = ~any(strcmpi('DensityFunction', p.UsingDefaults));
potentialSupplied = ~any(strcmpi('SupportPotentialFunction', p.UsingDefaults));
legacyPotentialSupplied = ~any(strcmpi( ...
    'DensityPotentialFunction', p.UsingDefaults));
if supportSupplied && legacySupportSupplied
    error('FASC:ConflictingSupportFunction', ...
        ['Specify SupportFunction or the legacy DensityFunction alias, ' ...
         'not both.']);
end
if potentialSupplied && legacyPotentialSupplied
    error('FASC:ConflictingSupportPotentialFunction', ...
        ['Specify SupportPotentialFunction or the legacy ' ...
         'DensityPotentialFunction alias, not both.']);
end
if supportSupplied
    opts.SupportFunction = opts.SupportFunction;
elseif legacySupportSupplied
    opts.SupportFunction = opts.DensityFunction;
else
    opts.SupportFunction = @(n) n;
end
if potentialSupplied
    opts.SupportPotentialFunction = opts.SupportPotentialFunction;
elseif legacyPotentialSupplied
    opts.SupportPotentialFunction = opts.DensityPotentialFunction;
else
    opts.SupportPotentialFunction = @(n) n .* (n - 1) ./ 2;
end
if ~isa(opts.SupportFunction, 'function_handle')
    error('FASC:InvalidSupportFunction', ...
        'SupportFunction must be a function handle.');
end
if ~isa(opts.SupportPotentialFunction, 'function_handle')
    error('FASC:InvalidSupportPotentialFunction', ...
        'SupportPotentialFunction must be a function handle.');
end
opts.BatchSize = double(opts.BatchSize);
end

function FASC_validateInputs(dataMatrix, idx_pos, idx_neg, simInter, simInner, ...
    initialLimit, clusterLimit, iterLimit, minVolume)
validateattributes(dataMatrix, {'numeric'}, {'2d', 'nonempty', 'real', 'finite'}, ...
    mfilename, 'dataMatrix');
validateattributes(simInter, {'numeric'}, {'scalar', 'real', 'finite'}, ...
    mfilename, 'simThres_inter');
validateattributes(simInner, {'numeric'}, {'scalar', 'real', 'finite'}, ...
    mfilename, 'simThres_inner');
validateattributes(initialLimit, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'initialClusterCountLimit');
validateattributes(clusterLimit, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'clusterCountLimit');
validateattributes(iterLimit, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'iterLimit');
validateattributes(minVolume, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'minClusterVolume');
if initialLimit > clusterLimit
    error('FASC:InvalidCapacity', ...
        'initialClusterCountLimit cannot exceed clusterCountLimit.');
end
if ~(isnumeric(idx_pos) && isvector(idx_pos) && isnumeric(idx_neg) && isvector(idx_neg))
    error('FASC:InvalidPolarityIndices', 'idx_pos and idx_neg must be numeric vectors.');
end
end

function FASC_validateDualIndices(idxPos, idxNeg, dimensionCount)
idx = [idxPos(:); idxNeg(:)];
if isempty(idxPos) || isempty(idxNeg) || any(~isfinite(idx)) || ...
        any(fix(idx) ~= idx) || any(idx < 1) || any(idx > dimensionCount) || ...
        numel(unique(idx)) ~= numel(idx)
    error('FASC:InvalidPolarityIndices', ...
        'Dual-cosine indices must be nonempty, unique, disjoint, and in range.');
end
end

function [key, similarityFcn] = FASC_similaritySpecification(algorithm, suppliedFcn)
if isa(algorithm, 'function_handle')
    key = "custom";
    similarityFcn = algorithm;
elseif ischar(algorithm) || (isstring(algorithm) && isscalar(algorithm))
    key = lower(string(algorithm));
    similarityFcn = suppliedFcn;
else
    error('FASC:InvalidSimilarity', ...
        'similarityAlgorithm must be a scalar name or function handle.');
end

if key == "custom" && isempty(similarityFcn)
    error('FASC:InvalidSimilarity', 'A custom similarity function is required.');
end
if key ~= "custom" && ~FASC_isBuiltInSimilarity(key) && isempty(similarityFcn)
    error('FASC:UnsupportedSimilarity', ...
        'Unsupported similarity ''%s''. Supply SimilarityFunction for a custom kernel.', key);
end
if ~isempty(similarityFcn) && key ~= "dual-cosine"
    key = "custom";
end
end

function tf = FASC_isBuiltInSimilarity(key)
tf = ismember(key, ["dual-cosine", "cosine", "euclidean", ...
    "euclidean-distance", "l1", "l1-norm", "manhattan", "minimum", ...
    "maximum", "algebraic", "logarithmic", "geometric", "harmonic", ...
    "enhanced harmonic", "entropy", "weighted entropy", "best average", ...
    "fitted core"]);
end

function dataWork = FASC_prepareData(dataMatrix, similarityKey, idxPos, idxNeg)
if ~isfloat(dataMatrix)
    dataWork = double(dataMatrix);
else
    dataWork = dataMatrix;
end

if similarityKey == "dual-cosine"
    dataWork(:, idxPos) = FASC_normalizeRowsL2(dataWork(:, idxPos));
    dataWork(:, idxNeg) = FASC_normalizeRowsL2(dataWork(:, idxNeg));
elseif similarityKey == "cosine"
    dataWork = FASC_normalizeRowsL2(dataWork);
elseif FASC_isBuiltInSimilarity(similarityKey)
    dataWork = FASC_normalizeRowsL1(dataWork);
end
end

function capacity = FASC_capacityAt(schedule, iteration, maximumCapacity)
if isempty(schedule)
    capacity = maximumCapacity;
elseif isa(schedule, 'function_handle')
    capacity = schedule(iteration);
else
    capacity = schedule(min(iteration, numel(schedule)));
end
validateattributes(capacity, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
capacity = min(double(capacity), double(maximumCapacity));
end

function tf = FASC_capacityScheduleIsFixed(schedule, iteration)
% A finite numeric schedule holds its last value. The future of a function
% handle cannot be certified from observed values, so recurrence detection
% remains disabled for that interface.
if isempty(schedule)
    tf = true;
elseif isa(schedule, 'function_handle')
    tf = false;
else
    tf = iteration >= numel(schedule);
end
end

function seedIdx = FASC_selectSeeds(dataMatrix, seedCount, similarityKey, ...
    similarityFcn, idxPos, idxNeg, batchSize)
sampleCount = size(dataMatrix, 1);
seedIdx = zeros(seedCount, 1);
if seedCount == 0
    return;
end

available = true(sampleCount, 1);
seedIdx(1) = FASC_lexicographicMin(dataMatrix, find(available));
available(seedIdx(1)) = false;
nearestSeedSimilarity = -inf(sampleCount, 1);

for s = 2:seedCount
    newCenter = dataMatrix(seedIdx(s - 1), :);
    for first = 1:batchSize:sampleCount
        last = min(first + batchSize - 1, sampleCount);
        rows = first:last;
        sims = FASC_similarityMatrix(dataMatrix(rows, :), newCenter, ...
            similarityKey, similarityFcn, idxPos, idxNeg, batchSize);
        nearestSeedSimilarity(rows) = max(nearestSeedSimilarity(rows), sims(:, 1));
    end
    candidateScores = nearestSeedSimilarity;
    candidateScores(~available) = inf;
    minimumScore = min(candidateScores);
    candidates = find(available & candidateScores == minimumScore);
    seedIdx(s) = FASC_lexicographicMin(dataMatrix, candidates);
    available(seedIdx(s)) = false;
end
end

function idx = FASC_lexicographicMin(dataMatrix, candidates)
if isempty(candidates)
    error('FASC:InternalSelectionError', 'No deterministic candidate is available.');
end
candidates = candidates(:);
for d = 1:size(dataMatrix, 2)
    values = dataMatrix(candidates, d);
    minimumValue = min(values);
    candidates = candidates(values == minimumValue);
    if isscalar(candidates)
        break;
    end
end
idx = candidates(1);
end

function [labels, maxAffinity] = FASC_assignInBlocks(dataMatrix, centers, counts, ...
    strategyKey, threshold, similarityKey, similarityFcn, idxPos, idxNeg, opts)
sampleCount = size(dataMatrix, 1);
clusterCount = size(centers, 1);
labels = zeros(sampleCount, 1);
maxAffinity = -inf(sampleCount, 1);
if clusterCount == 0
    return;
end

if ismember(strategyKey, ["DASS", "DENSITYFIRST"])
    supportBonus = opts.SupportFunction(double(counts(:)));
    if ~isnumeric(supportBonus) || numel(supportBonus) ~= clusterCount || ...
            any(~isfinite(supportBonus(:)))
        error('FASC:InvalidSupportFunction', ...
            'SupportFunction must return one finite value per cluster support.');
    end
    supportBonus = opts.Lambda .* reshape(double(supportBonus), 1, []);
else
    supportBonus = zeros(1, clusterCount);
end

for first = 1:opts.BatchSize:sampleCount
    last = min(first + opts.BatchSize - 1, sampleCount);
    rows = first:last;
    sims = FASC_similarityMatrix(dataMatrix(rows, :), centers, ...
        similarityKey, similarityFcn, idxPos, idxNeg, opts.BatchSize);
    maxAffinity(rows) = max(sims, [], 2);
    scores = sims + supportBonus;
    scores(sims < threshold) = -inf;
    [bestScore, bestCluster] = max(scores, [], 2);
    accepted = isfinite(bestScore);
    blockLabels = zeros(numel(rows), 1);
    blockLabels(accepted) = bestCluster(accepted);
    labels(rows) = blockLabels;
end
end

function promotionIdx = FASC_selectPromotions(dataMatrix, labels, maxAffinity, ...
    holeCount, similarityKey, similarityFcn, representativeFcn, idxPos, idxNeg, ...
    batchSize)
outliers = labels == 0;
promotionCount = min(holeCount, sum(outliers));
promotionIdx = zeros(promotionCount, 1);
available = outliers;
for p = 1:promotionCount
    score = maxAffinity;
    score(~available) = inf;
    minimumScore = min(score);
    candidates = find(available & score == minimumScore);
    promotionIdx(p) = FASC_lexicographicMin(dataMatrix, candidates);
    available(promotionIdx(p)) = false;

    representative = FASC_representative(dataMatrix(promotionIdx(p), :), ...
        similarityKey, similarityFcn, representativeFcn, idxPos, idxNeg, ...
        batchSize);
    availableRows = find(available);
    for first = 1:batchSize:numel(availableRows)
        last = min(first + batchSize - 1, numel(availableRows));
        rows = availableRows(first:last);
        similarity = FASC_similarityMatrix(dataMatrix(rows, :), representative, ...
            similarityKey, similarityFcn, idxPos, idxNeg, batchSize);
        maxAffinity(rows) = max(maxAffinity(rows), similarity(:, 1));
    end
end
end

function [centers, counts, labels] = FASC_consolidate(dataMatrix, labels, protected, ...
    minVolume, mergeThreshold, similarityKey, similarityFcn, representativeFcn, ...
    idxPos, idxNeg, batchSize)
clusterCount = max([0; labels(:)]);
if numel(protected) < clusterCount
    protected(end + 1:clusterCount, 1) = false;
else
    protected = protected(1:clusterCount);
end

[centers, counts] = FASC_recomputeClusters(dataMatrix, labels, clusterCount, ...
    similarityKey, similarityFcn, representativeFcn, idxPos, idxNeg, batchSize);

% Merge connected components simultaneously, then repeat because a newly
% computed representative can create additional merge edges.
while clusterCount > 1
    centerSimilarity = FASC_similarityMatrix(centers, centers, similarityKey, ...
        similarityFcn, idxPos, idxNeg, batchSize);
    adjacency = centerSimilarity >= mergeThreshold;
    adjacency = adjacency | adjacency';
    adjacency(1:clusterCount + 1:end) = true;
    components = FASC_connectedComponents(adjacency);
    newCount = max(components);
    if newCount == clusterCount
        break;
    end

    labels = FASC_remapLabels(labels, components);
    protected = accumarray(components, double(protected), [newCount, 1], @max, 0) > 0;
    clusterCount = newCount;
    [centers, counts] = FASC_recomputeClusters(dataMatrix, labels, clusterCount, ...
        similarityKey, similarityFcn, representativeFcn, idxPos, idxNeg, batchSize);
end

if clusterCount > 0
    keep = counts >= minVolume | protected;
    if ~all(keep)
        remap = zeros(clusterCount, 1);
        remap(keep) = 1:sum(keep);
        labels = FASC_remapLabels(labels, remap);
        centers = centers(keep, :);
        counts = counts(keep);
    end
end

[centers, counts, labels] = FASC_canonicalizeClusters(centers, counts, labels);
end

function labels = FASC_remapLabels(labels, remap)
assigned = labels > 0;
oldLabels = labels(assigned);
labels(assigned) = remap(oldLabels);
end

function components = FASC_connectedComponents(adjacency)
clusterCount = size(adjacency, 1);
components = zeros(clusterCount, 1);
componentCount = 0;
for start = 1:clusterCount
    if components(start) ~= 0
        continue;
    end
    componentCount = componentCount + 1;
    queue = start;
    components(start) = componentCount;
    head = 1;
    while head <= numel(queue)
        current = queue(head);
        head = head + 1;
        neighbours = find(adjacency(current, :) & components' == 0);
        if ~isempty(neighbours)
            components(neighbours) = componentCount;
            queue = [queue, neighbours]; %#ok<AGROW>
        end
    end
end
end

function [centers, counts] = FASC_recomputeClusters(dataMatrix, labels, clusterCount, ...
    similarityKey, similarityFcn, representativeFcn, idxPos, idxNeg, batchSize)
centers = zeros(clusterCount, size(dataMatrix, 2), 'like', dataMatrix);
counts = zeros(clusterCount, 1);
for j = 1:clusterCount
    members = dataMatrix(labels == j, :);
    counts(j) = size(members, 1);
    if counts(j) > 0
        centers(j, :) = FASC_representative(members, similarityKey, ...
            similarityFcn, representativeFcn, idxPos, idxNeg, batchSize);
    end
end
end

function center = FASC_representative(members, similarityKey, similarityFcn, ...
    representativeFcn, idxPos, idxNeg, batchSize)
if ~isempty(representativeFcn)
    center = representativeFcn(members);
    if ~isnumeric(center) || ~isreal(center) || any(~isfinite(center(:))) || ...
            ~isequal(size(center), [1, size(members, 2)])
        error('FASC:InvalidRepresentativeFunction', ...
            'RepresentativeFunction must return one finite 1-by-D vector.');
    end
    center = cast(center, 'like', members);
    return;
end

switch similarityKey
    case "dual-cosine"
        center = mean(members, 1);
        center(idxPos) = FASC_normalizeRowsL2(sum(members(:, idxPos), 1));
        center(idxNeg) = FASC_normalizeRowsL2(sum(members(:, idxNeg), 1));
    case "cosine"
        center = FASC_normalizeRowsL2(sum(members, 1));
    case {"euclidean", "euclidean-distance"}
        center = mean(members, 1);
    case {"l1", "l1-norm", "manhattan"}
        center = FASC_normalizeRowsL1(median(members, 1));
    otherwise
        % A deterministic empirical Frechet medoid is valid for any bounded
        % symmetric similarity. Large custom problems should provide a
        % domain-specific RepresentativeFunction to avoid quadratic work.
        center = FASC_similarityMedoid(members, similarityKey, similarityFcn, ...
            idxPos, idxNeg, batchSize);
end
end

function center = FASC_similarityMedoid(members, similarityKey, similarityFcn, ...
    idxPos, idxNeg, requestedBatchSize)
memberCount = size(members, 1);
if memberCount == 1
    center = members(1, :);
    return;
end
maxEntries = 1e7;
batchSize = min(requestedBatchSize, max(1, floor(maxEntries / memberCount)));
scores = zeros(memberCount, 1);
for first = 1:batchSize:memberCount
    last = min(first + batchSize - 1, memberCount);
    rows = first:last;
    sims = FASC_similarityMatrix(members(rows, :), members, similarityKey, ...
        similarityFcn, idxPos, idxNeg, batchSize);
    scores(rows) = sum(sims, 2);
end
bestScore = max(scores);
candidates = find(scores == bestScore);
best = FASC_lexicographicMin(members, candidates);
center = members(best, :);
end

function [centers, counts, labels] = FASC_canonicalizeClusters(centers, counts, labels)
clusterCount = numel(counts);
if clusterCount <= 1
    return;
end
keys = [-double(counts(:)), double(full(centers))];
[~, order] = sortrows(keys, 1:size(keys, 2));
inverse = zeros(clusterCount, 1);
inverse(order) = (1:clusterCount)';
labels = FASC_remapLabels(labels, inverse);
centers = centers(order, :);
counts = counts(order);
end

function [centers, counts, labels] = FASC_projectMatureState( ...
    centers, counts, labels, minVolume)
% Remove one-pass protected clusters that have not reached minimum support.
keep = counts >= minVolume;
if all(keep)
    return;
end
remap = zeros(numel(counts), 1);
remap(keep) = (1:sum(keep))';
labels = FASC_remapLabels(labels, remap);
centers = centers(keep, :);
counts = counts(keep);
[centers, counts, labels] = FASC_canonicalizeClusters(centers, counts, labels);
end

function score = FASC_objective(dataMatrix, centers, counts, labels, strategyKey, ...
    similarityKey, similarityFcn, idxPos, idxNeg, opts)
if ~isempty(opts.ObjectiveFunction)
    score = opts.ObjectiveFunction(dataMatrix, centers, counts, labels);
    validateattributes(score, {'numeric'}, {'scalar', 'real', 'finite'});
    score = double(score);
    return;
end

cohesion = 0;
assignedRows = find(labels > 0);
for first = 1:opts.BatchSize:numel(assignedRows)
    last = min(first + opts.BatchSize - 1, numel(assignedRows));
    rows = assignedRows(first:last);
    sims = FASC_similarityMatrix(dataMatrix(rows, :), centers, similarityKey, ...
        similarityFcn, idxPos, idxNeg, opts.BatchSize);
    linearIdx = sub2ind(size(sims), (1:numel(rows))', labels(rows));
    cohesion = cohesion + sum(sims(linearIdx));
end

if ismember(strategyKey, ["DASS", "DENSITYFIRST"])
    supportPotential = opts.SupportPotentialFunction(double(counts(:)));
    if ~isnumeric(supportPotential) || numel(supportPotential) ~= numel(counts) || ...
            any(~isfinite(supportPotential(:)))
        error('FASC:InvalidSupportPotentialFunction', ...
            'SupportPotentialFunction must return one finite value per support.');
    end
    score = double(cohesion) + opts.Lambda * sum(double(supportPotential(:)));
else
    score = double(cohesion);
end
end

function [stability, labelAgreement] = FASC_stability(previousCenters, previousCounts, ...
    centers, counts, previousLabels, labels, similarityKey, similarityFcn, ...
    idxPos, idxNeg, batchSize)
if isempty(previousLabels)
    stability = 0;
    labelAgreement = 0;
    return;
end
labelAgreement = mean(previousLabels == labels);

countLength = max(numel(previousCounts), numel(counts));
countA = zeros(countLength, 1);
countB = zeros(countLength, 1);
countA(1:numel(previousCounts)) = sort(previousCounts, 'descend');
countB(1:numel(counts)) = sort(counts, 'descend');
denominator = norm(countA) * norm(countB);
if denominator == 0
    countSimilarity = double(norm(countA - countB) == 0);
else
    countSimilarity = (countA' * countB) / denominator;
end

if isempty(previousCenters) || isempty(centers)
    centerSimilarity = double(isempty(previousCenters) && isempty(centers));
else
    pairwise = FASC_similarityMatrix(centers, previousCenters, similarityKey, ...
        similarityFcn, idxPos, idxNeg, batchSize);
    pairCount = min(size(pairwise));
    matched = zeros(pairCount, 1);
    weights = zeros(pairCount, 1);
    usedCurrent = false(size(centers, 1), 1);
    usedPrevious = false(size(previousCenters, 1), 1);
    for match = 1:pairCount
        available = pairwise;
        available(usedCurrent, :) = -inf;
        available(:, usedPrevious) = -inf;
        [value, linearIdx] = max(available(:));
        [row, col] = ind2sub(size(available), linearIdx);
        matched(match) = value;
        weights(match) = counts(row);
        usedCurrent(row) = true;
        usedPrevious(col) = true;
    end
    if sum(weights) == 0
        centerSimilarity = mean(matched);
    else
        centerSimilarity = sum(matched .* weights) / sum(weights);
    end
end
stability = (double(countSimilarity) + double(centerSimilarity)) / 2;
end

function state = FASC_captureState(centers, counts, labels)
state.centers = centers;
state.counts = counts;
state.labels = labels;
end

function hash = FASC_stateHash(labels, capacity)
payload = [int32(capacity); int32(labels(:))];
bytes = typecast(payload, 'uint8');
try
    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(typecast(bytes, 'int8'));
    digestBytes = typecast(digest.digest(), 'uint8');
    hash = lower(reshape(dec2hex(digestBytes, 2).', 1, []));
catch exception
    error('FASC:StateHashUnavailable', ...
        'SHA-256 state hashing is unavailable: %s', exception.message);
end
end

function sims = FASC_similarityMatrix(mtx1, mtx2, similarityKey, similarityFcn, ...
    idxPos, idxNeg, batchSize)
if isempty(mtx1) || isempty(mtx2)
    sims = zeros(size(mtx1, 1), size(mtx2, 1));
    return;
end

if similarityKey == "dual-cosine"
    simsPos = mtx1(:, idxPos) * mtx2(:, idxPos)';
    simsNeg = mtx1(:, idxNeg) * mtx2(:, idxNeg)';
    sims = min(simsPos, simsNeg);
elseif similarityKey == "custom"
    sims = FASC_customSimilarityMatrix(mtx1, mtx2, similarityFcn);
else
    sims = Matrix_SimilarityAlgorithms(mtx1, mtx2, similarityKey, ...
        'Parallel', false, 'BatchSize', batchSize, 'Normalize', false);
end

if ~isnumeric(sims) || ~isreal(sims) || ...
        ~isequal(size(sims), [size(mtx1, 1), size(mtx2, 1)]) || ...
        any(~isfinite(sims(:)))
    error('FASC:InvalidSimilarityOutput', ...
        ['Similarity ''%s'' returned size %s for expected size %s ', ...
        '(%d nonfinite values).'], char(similarityKey), mat2str(size(sims)), ...
        mat2str([size(mtx1, 1), size(mtx2, 1)]), sum(~isfinite(sims(:))));
end
end

function sims = FASC_customSimilarityMatrix(mtx1, mtx2, similarityFcn)
matrixCallSucceeded = false;
try
    candidate = similarityFcn(mtx1, mtx2);
    matrixCallSucceeded = isnumeric(candidate) && ...
        isequal(size(candidate), [size(mtx1, 1), size(mtx2, 1)]);
catch
    candidate = [];
end
if matrixCallSucceeded
    sims = candidate;
    return;
end

sims = zeros(size(mtx1, 1), size(mtx2, 1));
for i = 1:size(mtx1, 1)
    for j = 1:size(mtx2, 1)
        sims(i, j) = similarityFcn(mtx1(i, :), mtx2(j, :));
    end
end
end

function normalized = FASC_normalizeRowsL2(matrix)
norms = vecnorm(matrix, 2, 2);
norms(norms == 0) = 1;
normalized = matrix ./ norms;
end

function normalized = FASC_normalizeRowsL1(matrix)
rowSums = sum(abs(matrix), 2);
rowSums(rowSums == 0) = 1;
normalized = matrix ./ rowSums;
end

function sims = Matrix_SimilarityAlgorithms(mtx1, mtx2, algorithmName, varargin)
p = inputParser;
addParameter(p, 'Order', 20, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Parallel', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'BatchSize', 100000, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'Normalize', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
order = p.Results.Order;
useParallel = p.Results.Parallel;
batchSize = p.Results.BatchSize;
doNormalize = p.Results.Normalize;
rowCount = size(mtx1, 1);
colCount = size(mtx2, 1);
algorithmKey = lower(string(algorithmName));

switch algorithmKey
    case "cosine"
        sims = fast_cosine_similarity(mtx1, mtx2);
        return;
    case {"euclidean-distance", "euclidean"}
        sims = fast_euclidean_similarity(mtx1, mtx2, doNormalize);
        return;
    case {"l1-norm", "l1", "manhattan"}
        sims = fast_l1_similarity(mtx1, mtx2, doNormalize);
        return;
    case "minimum"
        sims = fast_minimum_similarity(mtx1, mtx2, doNormalize);
        return;
    case "maximum"
        sims = fast_maximum_similarity(mtx1, mtx2, doNormalize);
        return;
end

if doNormalize
    mtx1 = FASC_normalizeRowsL1(mtx1);
    mtx2 = FASC_normalizeRowsL1(mtx2);
end
if useParallel && rowCount * colCount > 10000
    sims = parallel_batch_computation(mtx1, mtx2, algorithmKey, order, batchSize);
elseif rowCount * colCount > 1000
    sims = batch_computation(mtx1, mtx2, algorithmKey, order, batchSize);
else
    sims = direct_computation(mtx1, mtx2, algorithmKey, order);
end
end

function sims = fast_cosine_similarity(mtx1, mtx2)
mtx1 = FASC_normalizeRowsL2(mtx1);
mtx2 = FASC_normalizeRowsL2(mtx2);
sims = mtx1 * mtx2';
end

function sims = fast_euclidean_similarity(mtx1, mtx2, doNormalize)
if doNormalize
    mtx1 = FASC_normalizeRowsL1(mtx1);
    mtx2 = FASC_normalizeRowsL1(mtx2);
end
distances = pdist2(mtx1, mtx2, 'squaredeuclidean');
sims = 1 - distances / 2;
end

function sims = fast_l1_similarity(mtx1, mtx2, doNormalize)
if doNormalize
    mtx1 = FASC_normalizeRowsL1(mtx1);
    mtx2 = FASC_normalizeRowsL1(mtx2);
end
distances = pdist2(mtx1, mtx2, 'cityblock');
sims = 1 - distances / 2;
end

function sims = fast_minimum_similarity(mtx1, mtx2, doNormalize)
if doNormalize
    mtx1 = FASC_normalizeRowsL1(mtx1);
    mtx2 = FASC_normalizeRowsL1(mtx2);
end
sims = zeros(size(mtx1, 1), size(mtx2, 1));
for i = 1:size(mtx1, 1)
    sims(i, :) = sum(min(mtx1(i, :), mtx2), 2)';
end
end

function sims = fast_maximum_similarity(mtx1, mtx2, doNormalize)
if doNormalize
    mtx1 = FASC_normalizeRowsL1(mtx1);
    mtx2 = FASC_normalizeRowsL1(mtx2);
end
sims = zeros(size(mtx1, 1), size(mtx2, 1));
for i = 1:size(mtx1, 1)
    sims(i, :) = (1 - sum(max(mtx1(i, :), mtx2), 2) / 2)';
end
end

function sims = parallel_batch_computation(mtx1, mtx2, algorithmName, order, batchSize)
rowCount = size(mtx1, 1);
colCount = size(mtx2, 1);
batchStarts = 1:batchSize:rowCount;
batchResults = cell(numel(batchStarts), 1);
parfor batchIdx = 1:numel(batchStarts)
    first = batchStarts(batchIdx);
    last = min(first + batchSize - 1, rowCount);
    block = zeros(last - first + 1, colCount);
    for i = 1:size(block, 1)
        for j = 1:colCount
            block(i, j) = SimilarityAlgorithms(...
                mtx1(first + i - 1, :), mtx2(j, :), algorithmName, order); %#ok<PFBNS>
        end
    end
    batchResults{batchIdx} = block;
end
sims = zeros(rowCount, colCount);
for batchIdx = 1:numel(batchStarts)
    first = batchStarts(batchIdx);
    last = min(first + batchSize - 1, rowCount);
    sims(first:last, :) = batchResults{batchIdx};
end
end

function sims = batch_computation(mtx1, mtx2, algorithmName, order, batchSize)
rowCount = size(mtx1, 1);
colCount = size(mtx2, 1);
sims = zeros(rowCount, colCount);
for first = 1:batchSize:rowCount
    last = min(first + batchSize - 1, rowCount);
    for i = first:last
        for j = 1:colCount
            sims(i, j) = SimilarityAlgorithms(mtx1(i, :), mtx2(j, :), ...
                algorithmName, order);
        end
    end
end
end

function sims = direct_computation(mtx1, mtx2, algorithmName, order)
sims = zeros(size(mtx1, 1), size(mtx2, 1));
for i = 1:size(mtx1, 1)
    for j = 1:size(mtx2, 1)
        sims(i, j) = SimilarityAlgorithms(mtx1(i, :), mtx2(j, :), ...
            algorithmName, order);
    end
end
end

function sim = SimilarityAlgorithms(vector1, vector2, algorithmName, ~)
vector1 = FASC_normalizeRowsL1(vector1);
vector2 = FASC_normalizeRowsL1(vector2);
algorithmKey = lower(string(algorithmName));
shared = vector1 ~= 0 & vector2 ~= 0;

switch algorithmKey
    case "algebraic"
        sim = sum((vector1(shared) + vector2(shared)) / 2);
    case "cosine"
        denominator = norm(vector1) * norm(vector2);
        if denominator == 0
            sim = 0;
        else
            sim = vector1 * vector2' / denominator;
        end
    case {"euclidean-distance", "euclidean"}
        sim = 1 - sum((vector1 - vector2).^2) / 2;
    case {"l1-norm", "l1", "manhattan"}
        sim = 1 - sum(abs(vector1 - vector2)) / 2;
    case "minimum"
        sim = sum(min(vector1, vector2));
    case "maximum"
        sim = 1 - sum(max(vector1, vector2)) / 2;
    case "logarithmic"
        a = vector1(shared);
        b = vector2(shared);
        equalValues = abs(a - b) <= 8 * eps(max(abs(a), abs(b)));
        terms = zeros(size(a));
        terms(equalValues) = (a(equalValues) + b(equalValues)) / 2;
        relativeDifference = (b(~equalValues) - a(~equalValues)) ./ a(~equalValues);
        terms(~equalValues) = a(~equalValues) .* relativeDifference ./ ...
            log1p(relativeDifference);
        sim = sum(terms);
    case "geometric"
        sim = sum(sqrt(vector1(shared) .* vector2(shared)));
    case "harmonic"
        sim = sum(2 .* vector1(shared) .* vector2(shared) ./ ...
            (vector1(shared) + vector2(shared)));
    case "enhanced harmonic"
        sim = sum(sqrt(2 ./ (vector1(shared).^(-2) + vector2(shared).^(-2))));
    case "entropy"
        entropy1 = entropyCalculater(vector1);
        entropy2 = entropyCalculater(vector2);
        entropyMean = entropyCalculater((vector1 + vector2) / 2);
        sim = 1 - (2 * entropyMean - entropy1 - entropy2) / log(4);
    case "weighted entropy"
        vector1 = entropyWeight(vector1);
        vector2 = entropyWeight(vector2);
        entropy1 = entropyCalculater(vector1);
        entropy2 = entropyCalculater(vector2);
        entropyMean = entropyCalculater((vector1 + vector2) / 2);
        sim = 1 - (2 * entropyMean - entropy1 - entropy2) / log(4);
    case {"best average", "fitted core"}
        a = vector1(shared);
        b = vector2(shared);
        sim = sum((a.^b .* b.^a).^(1 ./ (a + b)));
    otherwise
        error('FASC:UnsupportedSimilarity', ...
            'Unsupported similarity algorithm: %s', char(algorithmName));
end

    function entropy = entropyCalculater(vector)
        positive = vector > 0;
        entropy = -sum(vector(positive) .* log(vector(positive)));
    end

    function weighted = entropyWeight(vector)
        entropy = entropyCalculater(vector);
        weighted = vector;
        if entropy < 3
            weighted = vector .^ (0.25 + 0.25 * entropy);
            weighted = FASC_normalizeRowsL1(weighted);
        end
    end
end
