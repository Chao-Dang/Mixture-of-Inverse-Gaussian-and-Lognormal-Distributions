function x = LSS(nStrata)

%% Copyright (C) Shields Uncertainty Research Group (SURG)
% All Rights Reserved
% Johns Hopkins University
% Department of Civil Engineering
% Updated: 16 July 2015
% Michael D. Shields

%% "Latinized" Stratified Sampling
% Generates a sample set on U(0,1) that is simultaneously a Latin hypercube 
% sample and a stratified sample as described in:
% Shields, M.D. and Zhang, J. "The generalization of Latin hypercube 
% sampling" Reliability Engineering and System Safety. (in review)

% Input:
% nStrata ----- Defines the design stratification  
%               Vector (1 x nRVs)
%               Each column defines the number of stratifications to make
%                   in that direction.
% Example: nStrata = [25,25] produces a 2D design with 25x25 
% stratification and 625 total samples
% 
% Output: 
% x ----------- Generated samples
%               Array (nSamples x nRVs)
%
%% Initialize variables
nSamples = prod(nStrata);
nRVs = length(nStrata);

%% Draw Latinized Stratified Samples

% Draw a Latin hypercube sample
x_temp = lhsdesign(nSamples,nRVs);
x = zeros(size(x_temp));

% Array to identify the candidate strata in each dimension for the LHS
% samples
strata = zeros(size(x_temp));

% Identify the candidate strata in each dimension for the LHS samples
for i = 1:nRVs
    y = 1/nStrata(i):1/nStrata(i):1;
    for j = 1:nSamples
        for k = 1:length(y)
            if x_temp(j,i) < y(k) && strata(j,i) == 0
                strata(j,i) = k;
                continue
            end
        end    
    end
end

% Place the components of the LHS into the strata one-by-one
for i = 1:length(nStrata)
    l = 1;
    for j = 1:prod(nStrata(1:i))
        if l > nStrata(i)
            l = 1;
        end
        for k = 1:prod(nStrata)/prod(nStrata(1:i))
            for m = 1:length(strata)
                if strata(m,i) == l
                    x((j-1)*(prod(nStrata)/prod(nStrata(1:i)))+k,i) = x_temp(m,i);
                    strata(m,i) = 0;
                    break
                end
            end
        end
        l = l+1;
    end
end