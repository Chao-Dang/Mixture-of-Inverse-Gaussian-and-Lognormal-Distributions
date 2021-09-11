function [x_lpss] = LPSS(lpss_design,lpss_strata)

%% Copyright (C) Shields Uncertainty Research Group (SURG)
% All Rights Reserved
% Johns Hopkins University
% Department of Civil Engineering
% Updated: 7 January 2016
% Michael D. Shields

%% Latinized Partially Stratified Sampling
% Generates a latinized partially stratified sample set on U(0,1) as described in:
% Shields, M.D. and Zhang, J. "The generalization of Latin hypercube 
% sampling" Reliability Engineering and System Safety. 148: 96-108

% Input variables
% lpss_design: Vector defining the subdomains to be used                   
%       Example: 5D problem with 2x2D + 1x1D subdomains                         
%                pss_design = [2,2,1]                                     
%       Note: The sum of the values in the pss_design vector equals the   
%             dimension of the problem.   
% lpss_strata: Vector defining how each dimension should be stratified     
%       Example: 5D problem with 2x2D + 1x1D subdomains with 625 samples  
%                pss_strata = [25,25,625]                                 
%       Note: pss_strata(i)^pss_design(i) = number of samples (for all i) 
% Output variables: 
% x_lpss: Generated samples
%       Array (nSamples x nDim) - nDim = vector dimension

%% Check that the PSS design is valid
if length(lpss_design) ~= length(lpss_strata)
    error('Input vectors "lpss_design" and "lpss_strata" must be the same length')
end

sample_check = lpss_strata.^lpss_design;
if range(sample_check) ~=0
    error('All dimensions must have the same number of samples/strata. Check to ensure that all values of lpss_strata.^lpss_design are equal.')
end

%% Latinized Partially Stratified Sampling

nDim = sum(lpss_design); %总的维数 
nSamples = lpss_strata(1).^lpss_design(1);%总点个数  

col = 0;
for i = 1:length(lpss_design) %从1到分层个数 
% for i = 1:1; %从1到分层个数
    nStrata = lpss_strata(i)*ones(1,lpss_design(i));
    x_lpss(:,col+1:col+lpss_design(i)) = LSS(nStrata); %第一列 第二列 =nstrate
    x_lpss(:,col+1:col+lpss_design(i)) = x_lpss(randperm(nSamples),col+1:col+lpss_design(i));
    col = col+lpss_design(i);
end

    
