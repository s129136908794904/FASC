function [mixingPairs, simMtx] = Clustering_clusterListCountsHister(...
    clusterListCounts, clusterCenterMtx, sim_inter, similarityAlgorithm, idx_pos, idx_neg, ifPlot, ifSave, saveFileName, saveDir)
% Clustering_clusterListCountsHister
% Updated: Fixed colorbar tick offset by increasing colormap resolution.

colorMap = sanitizeBR1ColorMap(myColorMap("BR1"));
clusterCount = length(clusterListCounts);

similarityAlgorithm = string(similarityAlgorithm);
similarityAlgorithmKey = lower(similarityAlgorithm);
isDualCosine = similarityAlgorithmKey == "dual-cosine";

% Initialize the matrix with NaNs (Upper triangle logic for visualization)
simMtx = NaN(clusterCount, clusterCount);

if isDualCosine
    % === Dual-Cosine Logic ===
    simMtx_pos = clusterCenterMtx(:,idx_pos) * clusterCenterMtx(:,idx_pos)';
    simMtx_neg = clusterCenterMtx(:,idx_neg) * clusterCenterMtx(:,idx_neg)';
    
    for i = 1 : clusterCount-1
        for j = i+1 : clusterCount
            val_pos = simMtx_pos(i,j);
            val_neg = simMtx_neg(i,j);
            simMtx(i,j) = min(val_pos, val_neg);
        end
    end
else
    % === General Logic ===
    fullSimMtx = Matrix_SimilarityAlgorithms(clusterCenterMtx, clusterCenterMtx, similarityAlgorithmKey, 'Parallel', false);
    for i = 1 : clusterCount-1
        for j = i+1 : clusterCount
            simMtx(i,j) = fullSimMtx(i,j);
        end
    end
end

if ifPlot || ifSave
    stepwiseColors = buildStepwiseColormap(colorMap, sim_inter);
    renderMap = expandDiscreteColormap(stepwiseColors);
    
    validMask = ~isnan(simMtx);
    alphaData = double(validMask);
    colorData = simMtx;
    colorData(validMask) = min(max(colorData(validMask),0),1);
    colorData(~validMask) = 0;

    % Create Figure (Use Normalized Units for GUI compatibility)
    f = figure('Visible', 'off'); 
    
    % Create Axes with Square Aspect Ratio
    ax = axes(Parent=f);
    set(ax,'FontName','arial','FontSize',8,'LineWidth',0.5);
    
    % Plot Image
    image(ax, colorData, "AlphaData", alphaData, "CDataMapping","scaled");
    
    % --- KEY LAYOUT SETTINGS ---
    axis(ax, 'square');       % Force 1:1 Aspect Ratio
    axis(ax, 'tight');        
    set(ax,...
        'YDir','reverse', ...
        'XAxisLocation','top', ...
        'YAxisLocation','right', ...
        'TickDir','out', ...
        'Layer','top', ...
        'Color','none');
    box(ax,'off');
    
    colormap(ax,renderMap);
    
    % --- Automatic Colorbar (Right Side) ---
    c = colorbar(ax);
    c.Location = 'eastoutside'; 
    c.FontName = 'arial';
    c.FontSize = 8;
    c.LineWidth = 0.5;
    
    % Force limits to match data range [0 1] exactly
    caxis(ax,[0 1]);
    c.Limits = [0 1];
    
    tickValues = 0:0.1:1;
    tickLabels = arrayfun(@(x) sprintf('%.1f', x), tickValues, 'UniformOutput', false);
    c.Ticks = tickValues;
    c.TickLabels = tickLabels;

    % Ticks Logic for Matrix Axis
    maxTickCount = 6;
    if clusterCount <= maxTickCount
        heatmapTicks = 1:clusterCount;
    else
        targetIntervals = maxTickCount - 1;
        approxStep = clusterCount / targetIntervals;
        magnitude = 10^floor(log10(approxStep));
        fraction = approxStep / magnitude;
        if fraction <= 1
            niceFraction = 1;
        elseif fraction <= 2
            niceFraction = 2;
        elseif fraction <= 5
            niceFraction = 5;
        else
            niceFraction = 10;
        end
        niceStep = niceFraction * magnitude;
        startTick = ceil(1 / niceStep) * niceStep;
        if startTick < 1
            startTick = 1;
        end
        heatmapTicks = startTick:niceStep:clusterCount;
        if isempty(heatmapTicks)
            heatmapTicks = startTick;
        end
        if heatmapTicks(end) ~= clusterCount
            heatmapTicks = [heatmapTicks clusterCount];
        end
    end
    heatmapTicks = unique(round(heatmapTicks));
    heatmapTicks = heatmapTicks(heatmapTicks >= 1 & heatmapTicks <= clusterCount);
    if isempty(heatmapTicks)
        heatmapTicks = 1:clusterCount;
    end
    heatmapTickLabels = arrayfun(@(x) num2str(x), heatmapTicks, 'UniformOutput', false);
    ax.XTick = heatmapTicks;
    ax.XTickLabel = heatmapTickLabels;
    ax.YTick = heatmapTicks;
    ax.YTickLabel = heatmapTickLabels;

    if ifSave
        heatmapBaseName = saveDir + "heatMap" + saveFileName;
        exportgraphics(f, heatmapBaseName + ".svg", "ContentType","vector");
    end
end


idx_sel = simMtx >= sim_inter;
idx_sel_scalar = find(idx_sel==1);
[idx_sel_is, idx_sel_js] = ind2sub(size(simMtx), idx_sel_scalar);

mixingPairs = [idx_sel_is idx_sel_js];
clusterListFractions = clusterListCounts/sum(clusterListCounts);


if ifPlot || ifSave
    x_plot = 1 : clusterCount;
    f = figure('Visible', 'off'); 
    axBar = axes(Parent=f);
    
    set(axBar,'FontName','arial','FontSize',8,'LineWidth',0.5);
    hold(axBar,"on")

    bar(axBar,x_plot,clusterListFractions,FaceColor=[0.82 0.82 0.82],FaceAlpha=0.9,EdgeColor="none")

    axBar.XLim = [0.5 clusterCount+0.5];
    set(axBar,'FontSize',8,'FontName','arial');
    set(axBar,'ycolor','k');
    ylabel(axBar,"Fraction",VerticalAlignment="baseline",FontSize=8,FontName='arial',Color='k')
    xlabel(axBar,"Cluster",VerticalAlignment="top",FontSize=8,FontName='arial')
    set(axBar,'tickdir','out')
    grid(axBar,'off')
    box(axBar,'off')

    axBar.XTick = heatmapTicks;
    axBar.XTickLabel = heatmapTickLabels;
    ax.FontName = 'arial';
    ax.FontSize = 8;

    for i = 1:length(idx_sel_scalar)
        this_x_hi = idx_sel_is(i);
        this_x_lo = idx_sel_js(i);
        plot(axBar,[this_x_hi this_x_lo this_x_lo],...
            [clusterListFractions(this_x_hi) clusterListFractions(this_x_hi) ...
            clusterListFractions(this_x_lo)],"k-",LineWidth=0.5)
    end

    for i = 1:length(idx_sel_scalar)
        this_x_hi = idx_sel_is(i);
        this_x_lo = idx_sel_js(i);
        scatter(axBar,[this_x_hi this_x_lo this_x_lo],...
            [clusterListFractions(this_x_hi) clusterListFractions(this_x_hi) clusterListFractions(this_x_lo)],5,...
            Marker="o",MarkerEdgeColor="k",MarkerFaceColor="w")
    end

    hold(axBar,"off")

    if ifSave
        histBaseName = saveDir + "hist_clusterListFraction_" + saveFileName;
        exportgraphics(f, histBaseName + ".svg", "ContentType","vector");
    end

end
end

% =========================================================================
% Helper Functions (Visualization)
% =========================================================================
function stepwiseColors = buildStepwiseColormap(baseMap, simThres)
stepEdges = 0:0.1:1;
stepwiseColors = zeros(numel(stepEdges)-1,3);
[bluePalette, redPalette] = splitPalette(baseMap);

tol = 1e-12;
blueMask = stepEdges(2:end) <= (simThres + tol);
redMask = ~blueMask;

blueCount = nnz(blueMask);
redCount = nnz(redMask);

if blueCount > 0
    stepwiseColors(blueMask,:) = samplePaletteSegment(bluePalette, blueCount);
end
if redCount > 0
    stepwiseColors(redMask,:) = samplePaletteSegment(redPalette, redCount);
end
end

function discreteMap = expandDiscreteColormap(stepwiseColors)
% INCREASED RESOLUTION: Changed from 10 to 100 to fix tick alignment offset.
repeatFactor = 100; 

numBins = size(stepwiseColors,1);
% Create a high-resolution map by repeating each color 'repeatFactor' times.
% This ensures that if repeatFactor is 100, the map has 1000 entries (for 10 bins).
% Value 0.1 maps to index 101 (start of bin 2), creating a sharp transition exactly at 0.1.
discreteMap = zeros(numBins * repeatFactor, 3);
for i = 1:numBins
    rows = (i-1)*repeatFactor + (1:repeatFactor);
    discreteMap(rows,:) = repmat(stepwiseColors(i,:), repeatFactor, 1);
end
end

function colors = samplePaletteSegment(segment, count)
if count <= 0
    colors = zeros(0,3);
    return
end

if isempty(segment)
    colors = repmat([0 0 0], count,1);
    return
end

if size(segment,1) == 1
    colors = repmat(segment, count,1);
    return
end

queryPts = linspace(1, size(segment,1), count);
colors = interp1(1:size(segment,1), segment, queryPts, "linear");
end

function [bluePalette, redPalette] = splitPalette(baseMap)
if isempty(baseMap)
    bluePalette = [0.4 0.6 0.9];
    redPalette = [0.9 0.4 0.4];
    return
end

midIdx = floor(size(baseMap,1)/2);
if midIdx < 1
    midIdx = 1;
end
bluePalette = baseMap(1:midIdx,:);
if midIdx < size(baseMap,1)
    redPalette = baseMap(midIdx+1:end,:);
else
    redPalette = baseMap(midIdx:end,:);
end
end

% =========================================================================
% Helper Functions (Similarity Calculation)
% =========================================================================
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
if useParallel && isempty(gcp('nocreate'))
    try
        parpool('local');
    catch
        warning('');
        useParallel = false;
    end
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
    epsVal = eps(cast(1, 'like', mtx1));
    if doNormalize
        denom1 = vecnorm(mtx1, 2, 2) + epsVal;
        denom2 = vecnorm(mtx2, 2, 2) + epsVal;
        mtx1_norm = mtx1 ./ denom1;
        mtx2_norm = mtx2 ./ denom2;
    else
        norms1 = vecnorm(mtx1, 2, 2);
        norms2 = vecnorm(mtx2, 2, 2);
        mtx1_norm = mtx1 ./ (norms1 + epsVal);
        mtx2_norm = mtx2 ./ (norms2 + epsVal);
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