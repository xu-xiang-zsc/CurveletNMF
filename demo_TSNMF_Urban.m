%% Demo for HYDICE Urban data set.

close all;
clear;
clc;


% load image data set
load Urban_R162;
load Urban_groundTruth4;

norm_y = sqrt(mean(mean(Y.^2)));
% rescale M and Y and lambda
M = M/norm_y;
Y = Y/norm_y;

numSources = size(A,1);
nb = size(M,1);

% reference abundances and endmembers
reference.Y = Y;
reference.S = A;
reference.A = M;
reference.sourcesShape = [nRow, nCol];
sourcesShape = [nRow, nCol]; % 2D shape of images


%% parameter setting

clear parameters;
parameters.reference = reference;
parameters.rank = numSources; % number of sources to recover
parameters.Y = Y; % noisy data
parameters.MaximumIteration = 600;
parameters.phaseRatio = 0.15; % ratio of iterations before refinement phase

% parameters for the update of A
parameters.A.MaximumIteration = 80; % subiterations for A
parameters.A.RelativeDifferenceTolerance = 0.000001;
parameters.A.normConstrained = 0;

% parameters for the update of S
parameters.S.tau_MAD = 15;
parameters.S.MaximumIteration = 30; % subiterations for S
parameters.S.sourcesShape = sourcesShape; % 2D shape of the images
parameters.S.directSparsity = 1; % whether the images are sparse in the direct domain or not
parameters.S.reweightedL1 = 0; % use reweigthed L1
parameters.S.scales = 3; % nscale of Curvelet
parameters.S.angles = 16; % directions of Curvelet
parameters.S.uniformFirstThreshold = 0; % use the same threshold for each source at the
% first iteration (default: 1). If the sources have very different dynamics, better set to 0.

parameters.coeff_len = 49;  % S3D16-49,s3d32-97, s5d16-145


%%

rng(12, 'twister');

result = TSNMF_core(parameters);


%% save results

tau_str = num2str(parameters.S.tau_MAD);
eval(sprintf('save TSNMF_Urban%d_tau%s_S3D16_pca.mat result;', numSources, tau_str));  
eval(sprintf('save TSNMF_Urban%d_tau%s_S3D16_pca.mat reference -append;', numSources, tau_str)); 