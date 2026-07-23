function tests = test_FASC
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
addpath(fullfile(testDir, '..', 'src'));
testCase.TestData.dualData = dualFixture();
repoDir = fileparts(fileparts(testDir));
testCase.TestData.spmsDataPath = fullfile(repoDir, 'data', 'dataMatrix.csv');
end

function testPermutationEquivariance(testCase)
data = testCase.TestData.dualData;
permutation = [7 1 13 4 10 16 2 8 14 5 11 17 3 9 15 6 12 18];

[centersA, countsA, labelsA] = runDual(data, 4);
[centersB, countsB, labelsPermuted] = runDual(data(permutation, :), 4);
labelsB = zeros(size(labelsPermuted));
labelsB(permutation) = labelsPermuted;

verifyEqual(testCase, labelsA == 0, labelsB == 0);
verifyEqual(testCase, coClustering(labelsA), coClustering(labelsB));
verifyEqual(testCase, countsA, countsB);
verifyEqual(testCase, centersA, centersB, 'AbsTol', 1e-12);
end

function testBatchSizeEquivalence(testCase)
data = testCase.TestData.dualData;
[centersA, countsA, labelsA, infoA] = runDual(data, 2);
[centersB, countsB, labelsB, infoB] = runDual(data, 1000);

verifyEqual(testCase, labelsA, labelsB);
verifyEqual(testCase, countsA, countsB);
verifyEqual(testCase, centersA, centersB, 'AbsTol', 1e-12);
verifyEqual(testCase, infoA.finalObjective, infoB.finalObjective, 'AbsTol', 1e-12);
end

function testExactRecoveryUnderThresholdSeparation(testCase)
data = [
    1.00 0.02 0.01 0.01;
    0.98 0.03 0.01 0.00;
    1.00 0.01 0.02 0.00;
    0.99 0.02 0.01 0.01;
    0.02 1.00 0.01 0.01;
    0.03 0.98 0.01 0.00;
    0.01 1.00 0.02 0.00;
    0.02 0.99 0.01 0.01;
    0.01 0.02 1.00 0.01;
    0.01 0.03 0.98 0.00;
    0.02 0.01 1.00 0.00;
    0.01 0.02 0.99 0.01];
truth = repelem((1:3)', 4);

[~, counts, labels, info] = FASC(data, 1:2, 3:4, ...
    0.98, 0.95, 1, 6, 20, 'DASS', 2, 'cosine', ...
    'Lambda', 1, 'BatchSize', 5, 'Verbose', false);

verifyEqual(testCase, coClustering(labels), truth == truth');
verifyEqual(testCase, sort(counts), [4; 4; 4]);
verifyEqual(testCase, sum(labels == 0), 0);
verifyEqual(testCase, string(info.terminationReason), "fixed-point");
end

function testSignedCosinePreprocessingPreservesDirection(testCase)
data = [
    1 -2;
    2 -4;
    -1 2;
    -2 4];
truth = [1; 1; 2; 2];

[~, counts, labels] = FASC(data, 1, 2, ...
    0.90, 0.90, 1, 2, 10, 'SF', 1, 'cosine', ...
    'BatchSize', 2, 'Verbose', false);

verifyEqual(testCase, coClustering(labels), truth == truth');
verifyEqual(testCase, sort(counts), [2; 2]);
end

function testFinalOutputProjectsProvisionalClusters(testCase)
data = [
    1.00 0.01;
    0.99 0.02;
    0.01 1.00;
    0.02 0.99;
    -1.00 0.00];

[~, counts, labels, info] = FASC(data, 1, 2, ...
    0.95, 0.95, 1, 3, 20, 'SF', 2, 'cosine', ...
    'BatchSize', 2, 'Verbose', false);

verifyEqual(testCase, sort(counts), [2; 2]);
verifyTrue(testCase, all(counts >= 2));
verifyEqual(testCase, labels(end), 0);
verifyTrue(testCase, info.matureProjectionApplied);
end

function testCycleGuardWaitsForFixedCapacityTail(testCase)
data = [
    -1.00 0.00;
    -0.99 0.01;
     1.00 0.00;
     0.99 0.01];

[~, counts, labels, info] = FASC(data, 1, 2, ...
    0.90, 0.90, 1, 2, 20, 'SF', 1, 'cosine', ...
    'CapacitySchedule', [1 2 1 2], ...
    'BatchSize', 2, 'Verbose', false);

verifyEqual(testCase, sort(counts), [2; 2]);
verifyEqual(testCase, sum(labels == 0), 0);
verifyEqual(testCase, string(info.terminationReason), "fixed-point");
verifyFalse(testCase, info.cycleDetected);
verifyTrue(testCase, info.capacityScheduleFixed);
verifyGreaterThanOrEqual(testCase, info.convergeIter, 5);
end

function testRealSpmsPermutationEquivariance(testCase)
dataPath = testCase.TestData.spmsDataPath;
verifyEqual(testCase, exist(dataPath, 'file'), 2, ...
    'The repository SPMS fixture is missing.');
data = readmatrix(dataPath, 'Range', [1 1 1000 600]);
permutation = size(data, 1):-1:1;

[centersA, countsA, labelsA, infoA] = runSpms(data, 137);
[centersB, countsB, labelsPermuted, infoB] = ...
    runSpms(data(permutation, :), 1000);
labelsB = zeros(size(labelsPermuted));
labelsB(permutation) = labelsPermuted;

verifyEqual(testCase, labelsA == 0, labelsB == 0);
verifyEqual(testCase, coClustering(labelsA), coClustering(labelsB));
verifyEqual(testCase, countsA, countsB);
verifyEqual(testCase, centersA, centersB, 'AbsTol', 1e-12);
verifyEqual(testCase, infoA.finalObjective, infoB.finalObjective, ...
    'AbsTol', 1e-10);
verifyEqual(testCase, numel(countsA), 22);
verifyEqual(testCase, sum(labelsA == 0), 437);
verifyTrue(testCase, all(countsA >= 2));
verifyTrue(testCase, infoA.matureProjectionApplied);
verifyEqual(testCase, string(infoA.terminationReason), "limit-cycle");
verifyEqual(testCase, infoA.cycleLength, 2);
end

function testPromotionUpdatesMaximinAffinity(testCase)
data = [0 1; 1 0; 0.99 0.10; 0.70 0.70];
[centers, counts, labels] = FASC(data, 1, 2, ...
    1.1, 0.999, 1, 3, 1, 'SF', 1, 'cosine', ...
    'BatchSize', 2, 'Verbose', false);

diagonalRepresentative = [1 1] / sqrt(2);
nearFirstPromotion = [0.99 0.10] / norm([0.99 0.10]);
verifyEqual(testCase, numel(counts), 3);
verifyEqual(testCase, sum(labels == 0), 1);
verifyLessThan(testCase, min(vecnorm(centers - diagonalRepresentative, 2, 2)), ...
    1e-12);
verifyGreaterThan(testCase, min(vecnorm(centers - nearFirstPromotion, 2, 2)), ...
    0.05);
end

function testDualCosineObjectiveUsesMinimum(testCase)
data = testCase.TestData.dualData;
lambda = 0.125;
[centers, counts, labels, info] = FASC(data, 1:3, 4:6, ...
    0.98, 0.75, 3, 3, 30, 'DASS', 1, 'dual-cosine', ...
    'Lambda', lambda, 'BatchSize', 5, 'Verbose', false);

normalized = data;
normalized(:, 1:3) = normalizeRows(normalized(:, 1:3));
normalized(:, 4:6) = normalizeRows(normalized(:, 4:6));
assigned = labels > 0;
assignedCenters = centers(labels(assigned), :);
positive = sum(normalized(assigned, 1:3) .* assignedCenters(:, 1:3), 2);
negative = sum(normalized(assigned, 4:6) .* assignedCenters(:, 4:6), 2);
expected = sum(min(positive, negative)) + ...
    lambda * sum(counts .* (counts - 1) / 2);

verifyEqual(testCase, info.finalObjective, expected, 'AbsTol', 1e-10);
end

function testSupportCallbackNamesAndLegacyAliasesMatch(testCase)
data = testCase.TestData.dualData;
support = @(n) n;
potential = @(n) n .* (n - 1) ./ 2;

[centersA, countsA, labelsA, infoA] = FASC(data, 1:3, 4:6, ...
    0.98, 0.75, 3, 3, 30, 'DASS', 1, 'dual-cosine', ...
    'Lambda', 0.125, ...
    'SupportFunction', support, ...
    'SupportPotentialFunction', potential, ...
    'BatchSize', 4, 'Verbose', false);
[centersB, countsB, labelsB, infoB] = FASC(data, 1:3, 4:6, ...
    0.98, 0.75, 3, 3, 30, 'DASS', 1, 'dual-cosine', ...
    'Lambda', 0.125, ...
    'DensityFunction', support, ...
    'DensityPotentialFunction', potential, ...
    'BatchSize', 4, 'Verbose', false);

verifyEqual(testCase, centersA, centersB, 'AbsTol', 1e-12);
verifyEqual(testCase, countsA, countsB);
verifyEqual(testCase, labelsA, labelsB);
verifyEqual(testCase, infoA.finalObjective, infoB.finalObjective, ...
    'AbsTol', 1e-12);
end

function testAllBuiltInMetricsRun(testCase)
data = [
    8 1 1 0 0;
    7 2 1 0 0;
    1 8 1 0 0;
    1 7 2 0 0;
    0 1 8 1 0;
    0 1 7 2 0;
    0 0 1 8 1;
    0 0 2 7 1;
    1 0 0 1 8;
    2 0 0 1 7];
metrics = ["cosine", "euclidean", "l1", "minimum", "maximum", ...
    "algebraic", "logarithmic", "geometric", "harmonic", ...
    "enhanced harmonic", "entropy", "weighted entropy", ...
    "best average", "fitted core"];

for metric = metrics
    fprintf('Built-in metric smoke test: %s\n', metric);
    [centers, counts, labels, info] = FASC(data, 1:2, 3:5, ...
        1.1, -1, 2, 3, 8, 'SF', 1, metric, ...
        'BatchSize', 3, 'Verbose', false);
    verifyTrue(testCase, all(isfinite(centers(:))), char(metric));
    verifyEqual(testCase, sum(counts), sum(labels > 0), char(metric));
    verifyTrue(testCase, isfinite(info.finalObjective), char(metric));
end
end

function testCustomKernelAndRepresentative(testCase)
data = [0 0; 0.1 0; 1 1; 1 0.9];
similarity = @(left, right) 1 - pdist2(left, right, 'cityblock');
representative = @(members) median(members, 1);

[centers, counts, labels, info] = FASC(data, 1, 2, ...
    0.95, 0.5, 2, 2, 10, 'SF', 1, 'custom', ...
    'SimilarityFunction', similarity, ...
    'RepresentativeFunction', representative, ...
    'Verbose', false);

verifyEqual(testCase, sort(counts), [2; 2]);
verifyEqual(testCase, sum(labels > 0), 4);
verifyTrue(testCase, all(isfinite(centers(:))));
verifyTrue(testCase, ismember(string(info.terminationReason), ...
    ["fixed-point", "limit-cycle"]));
end

function [centers, counts, labels, info] = runDual(data, batchSize)
[centers, counts, labels, info] = FASC(data, 1:3, 4:6, ...
    0.98, 0.75, 3, 3, 30, 'DASS', 1, 'dual-cosine', ...
    'Lambda', 0.125, 'BatchSize', batchSize, 'Verbose', false);
end

function [centers, counts, labels, info] = runSpms(data, batchSize)
[centers, counts, labels, info] = FASC(data, 1:300, 301:600, ...
    0.70, 0.70, 8, 50, 200, 'DASS', 2, 'dual-cosine', ...
    'Lambda', 1, 'BatchSize', batchSize, 'Verbose', false);
end

function relation = coClustering(labels)
relation = labels == labels';
assigned = labels > 0;
relation = relation & assigned & assigned';
end

function normalized = normalizeRows(matrix)
norms = vecnorm(matrix, 2, 2);
norms(norms == 0) = 1;
normalized = matrix ./ norms;
end

function data = dualFixture()
groupA = [
    1.00 0.02 0.01  0.98 0.03 0.01;
    0.97 0.05 0.02  1.00 0.02 0.01;
    0.99 0.03 0.02  0.97 0.05 0.02;
    1.00 0.01 0.03  0.99 0.03 0.02;
    0.98 0.04 0.01  0.98 0.04 0.01;
    0.96 0.06 0.02  0.97 0.04 0.03];
groupB = groupA(:, [2 1 3 5 4 6]);
groupC = groupA(:, [3 2 1 6 5 4]);
data = [groupA; groupB; groupC];
end
