%% QM_COMPUTE_SCAT_AVG_WEIGHTS
%
% Computes the average scattering weights for the orthonormal scattering
% coefficients, taken over numerous draws of the training set. These can
% then be used in QM_DISPLAY_SCAT_AVG_WEIGHTS to generate several figures.
% Scattering coefficients must first be computed with QM_COMPUTE_SCAT_2D.
%
% This file is part of ScatNet_QM_2D.
%
% Author: Matthew Hirn
% Email: mhirn@msu.edu
%
% Copyright 2016 Matthew J. Hirn
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%% Settings

% Folder for saving the outputs
foldername = [];

% Data
data_opt.data_set = 'qm_2d';                        % qm7_2d or qm_2d
data_opt.N = [512 398];                             % Number of pixels along x-axis and y-axis

% Densities
densities = {'dirac', 'core', 'valence'};           % Set which densities to use:
                                                    %   (a) dirac (b) atomic (c) core/valence
                                                    %   (d) dirac/atomic (e) dirac/core/valence

% Scattering
scat_opt.M = 2;                                     % Number of scattering layers (max is 2)
scat_opt.O2_invariant = true;                       % If true, scattering coefficients invariant to reflections

% Filter
filt_opt.n_wavelet_per_octave = 1;                  % Number of wavelets per octave/scale                         
filt_opt.J = log2(data_opt.N(1));                   % Max scale 
filt_opt.L = 16;                                    % Number of angles in [0,pi)
filt_opt.xi_psi = 3*pi/4;                           % Central frequency of the mother wavelet
filt_opt.slant_psi = 2/3;                           % Slant of ellipsoidal envelope

% Regression
max_p_init = 2004;                                  % Max number of regression terms
num_bags = 1000;                                    % Number of bagging iterations
bag_per = 0.80;                                     % Size of training bag relative to data set.

%% Load and unwind scattering

% Load scattering
if scat_opt.M == 1
    M1 = true;
    scat_opt.M = 2;
else
    M1 = false;
end
filename = qm_scat_save_name(data_opt, scat_opt, filt_opt);
load(filename);
if M1
    scat_opt.M = 1;
end

% Unwind scattering
num_densities = length(densities);
EU_cell = cell(1, 2 * num_densities);
scat_type = cell(1, 2 * num_densities);
field_val = cell(1, 2 * num_densities);
for i=1:num_densities
    switch densities{i}
        case 'dirac'
            EU_cell{2*i-1} = EU.dirac.L1(1:(scat_opt.M+1));
            EU_cell{2*i} = EU.dirac.L2(1:(scat_opt.M+1));
            field_val{2*i-1} = 'dirac';
            field_val{2*i} = 'dirac';
        case 'atomic'
            EU_cell{2*i-1} = EU.atomic.L1(1:(scat_opt.M+1));
            EU_cell{2*i} = EU.atomic.L2(1:(scat_opt.M+1));
            field_val{2*i-1} = 'atomic';
            field_val{2*i} = 'atomic';
        case 'core'
            EU_cell{2*i-1} = EU.core.L1(1:(scat_opt.M+1));
            EU_cell{2*i} = EU.core.L2(1:(scat_opt.M+1));
            field_val{2*i-1} = 'core';
            field_val{2*i} = 'core';
        case 'valence'
            EU_cell{2*i-1} = EU.valence.L1(1:(scat_opt.M+1));
            EU_cell{2*i} = EU.valence.L2(1:(scat_opt.M+1));
            field_val{2*i-1} = 'valence';
            field_val{2*i} = 'valence';
    end
    scat_type{2*i-1} = 't';
    scat_type{2*i} = 't';
end

[X, coeff_names, coeff_pars] = ...
    unwind_exp_scat_2d(EU_cell, scat_type, field_val, 'density');

clear i filename EU EU_cell scat_type field_val M1 num_densities

%% Load the data (energies)

load(data_opt.data_set, 'T');
T = T';

%% OLS regression (looped over multiple training sets)

% Initialize
Aq_mean_bag = cell(num_bags, 1);
pterm_ind = cell(num_bags, 1);
max_p = max_p_init * ones(num_bags, 1);
ind_bag = cell(num_bags, 1);
j = 1;
pathname = [];

% OLS
while j <= num_bags
    
    % Display counter
    display(['Bag number = ', num2str(j)]);
    
    % Get a subset of the data
    ind_bag{j} = randperm(size(X,1));
    ind_bag{j} = ind_bag{j}(1:round(bag_per * size(X,1)));
    pterm_ind{j} = ols(T(ind_bag{j}), X(ind_bag{j},:), max_p(j));
    max_p(j) = length(pterm_ind{j});
    Xp = X(ind_bag{j}, pterm_ind{j});
    
    % QR
    [Q, ~] = qr(Xp);
    Q = Q(:, 1:max_p(j));
    Aq = Q' * T(ind_bag{j});
    
    % Regression
    T_reg_q = Q * Aq;
    
    % Regression error
    TQ_err = mean(abs(T(ind_bag{j}) - T_reg_q));
    display(['Q avg regression deviation: ', num2str(TQ_err)]);
    
    % Squared amplitude average weights in Q
    % A priori, we compute:
    %   Aq_mean_bag{j}(p) = sqrt( mean((Aq(p) * Q(:,p)).^2) );
    % But since Q is orthonormal, this reduces to:
    %   Aq_mean_bag{j}(p) = (1/sqrt(length(ind_bag{j}))) * abs(Aq(p))
    if TQ_err < 3 % If max_p_init is small, this may need to be increased
        Aq_mean_bag{j} = (1/sqrt(length(ind_bag{j}))) * abs(Aq);
        j = j + 1;
    else
        max_p(j) = max_p_init;
    end
    
    % Save every 10
    if mod(j-1,10) == 0
        
        % Add up over bags
        Aq_mean = cell(min(max_p), 1);
        for p=1:min(max_p)
            Aq_mean{p} = zeros(size(X,2), 1);
            for k=1:(j-1)
                Aq_mean{p}(pterm_ind{k}(1:p)) = Aq_mean{p}(pterm_ind{k}(1:p)) + ...
                    Aq_mean_bag{k}(1:p);
            end
            Aq_mean{p} = Aq_mean{p} / (j-1);
        end
        clear Aq Q R T_reg_q Xp TQ_err
        
        % Path
        old_pathname = pathname;
        filename = qm_regression_scat_save_name(data_opt, scat_opt, filt_opt, densities);
        filename(end-3:end) = [];
        filename = strcat(filename, '_num_bags=', num2str(j-1), '_max_p=', num2str(max_p_init), '.mat');
        foldername = '/Users/mhirn/Dropbox/Mathematics/Matlab_files/work/ENS/DFT/qm/scatnet_qm_2d_internal/results/regression/analysis_display/';
        pathname = strcat(foldername, filename);
        clear filename foldername
        
        % Save workspace
        save(pathname);
        
        % Delete old file
        if ~isempty(old_pathname)
            delete(old_pathname);
        end
    end
    
end

% Add up over bags
Aq_mean = cell(min(max_p), 1);
for p=1:min(max_p)
    Aq_mean{p} = zeros(size(X,2), 1);
    for j=1:num_bags
        Aq_mean{p}(pterm_ind{j}(1:p)) = Aq_mean{p}(pterm_ind{j}(1:p)) + ...
            Aq_mean_bag{j}(1:p);
    end
    Aq_mean{p} = Aq_mean{p} / num_bags;
end

clear Aq i j k p Q T_reg_q Xp TQ_err old_pathname

%% Save

% Path
filename = qm_regression_scat_save_name(data_opt, scat_opt, filt_opt, densities);
filename(end-3:end) = [];
filename = strcat(filename, '_num_bags=', num2str(num_bags), '_max_p=', num2str(max_p_init), '.mat');
pathname = strcat(foldername, filename);
clear filename foldername

% Save workspace
save(pathname);