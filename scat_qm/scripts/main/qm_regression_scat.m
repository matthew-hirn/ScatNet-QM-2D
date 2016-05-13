%% QM_REGRESSION_SCAT
%
% Computes scattering orthogonal least squares regression for the QM data
% sets. Cross validates the optimal number of regression terms and uses
% bagging to reduce the variance.
%
% Important variables created by this script, which are saved:
%   1.) MAE_fold (vector): Vector of mean absolute errors (MAE) indexed by
%       testing fold. For each fold, the MAE value recorded is the average
%       over numerous "bags" (models), each of which uses cross validation
%       to determine the number of regression terms.
%   2.) RMSE_fold (vector): Vector of root mean square errors (RMSE) 
%       indexed by testing fold. For each fold, the RMSE value recorded is 
%       the average over numerous "bags" (models), each of which uses cross
%       validation to determine the number of regression terms.
%   3.) p_mae_fold (matrix): Matrix of number of regression terms, indexed
%       by testing fold and bag. The (k,l) entry is the number of 
%       regression terms that yields the MAE on the kth fold and lth bag.
%   4.) p_rmse_fold (matrix): Matrix of number of regression terms, indexed
%       by testing fold and bag. The (k,l) entry is the number of 
%       regression terms that yields the RMSE on the kth fold and lth bag.
%   5.) res_err_fold_mae (cell):  Residual errors on the database organized
%       by fold. The kth cell contains the residual errors for the
%       molecules in fold k, optimized against the MAE. 
%   6.) res_err_fold_rmse (cell): Residual errors on the database organized
%       by fold. The kth cell contains the residual errors for the
%       molecules in fold k, optimized against the RMSE. 
%   7.) T_reg_fold_mae (cell): Regressed energies on the database 
%       organized by fold. The kth cell contains the regressed energies for
%       the molecules in fold k, optimized against the MAE. 
%   8.) T_reg_fold_rmse (cell): Regressed energies on the database 
%       organized by fold. The kth cell contains the regressed energies for
%       the molecules in fold k, optimized against the RMSE. 
%   9.) coeff_names (cell): List of the scattering coefficient names.
%  10.) coeff_pars (struct): Struct in which each field corresponds to one
%       of the parameters of the scattering invariant coefficients. Each
%       field contains a 1xD vector or 1xD cell, where D is the number of
%       scattering features, and lists the value of the corresponding
%       parameter for each coefficient.
%
% The script prints the key error statistics in addition to saving
% variables.
%
% The user must enter an appropriate folder name in order to save the
% outputs. 
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
scat_opt.M = 2;                                     % Number of scattering layers (M=1 is wavelet; M=2 is scattering)
scat_opt.O2_invariant = true;                       % If true, scattering coefficients invariant to reflections

% Filter
filt_opt.n_wavelet_per_octave = 1;                  % Number of wavelets per octave/scale                         
filt_opt.J = log2(data_opt.N(1));                   % Max scale 
filt_opt.L = 16;                                    % Number of angles in [0,pi)
filt_opt.xi_psi = 3*pi/4;                           % Central frequency of the mother wavelet
filt_opt.slant_psi = 2/3;                           % Slant of ellipsoidal envelope

% Regression
max_p = 2^9 + 2^7;                                  % Max number of regression terms
num_bags = 10;                                      % Number of bags
bag_per = 0.80;                                     % Size of training bag relative to full training set

%% Load scattering

% Load the computed scattering coefficients
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

clear M1 filename

%% Unwind scattering

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

%% Load the data (energies and folds)

load(data_opt.data_set, 'T', 'P');
T = T';

%% Fully learned regression

[p_mae_fold, p_rmse_fold, T_reg_fold_mae, T_reg_fold_rmse, res_err_fold_mae, ...
    res_err_fold_rmse] = qm_kfold_regression(T, X, P, max_p, num_bags, bag_per);

clear X

%% Errors

% MAE and RMSE by fold
MAE_fold = zeros(1, 5);
RMSE_fold = zeros(1, 5);
for i=1:5
    MAE_fold(i) = mean(abs(res_err_fold_mae{i}));
    RMSE_fold(i) = sqrt(mean(res_err_fold_rmse{i}.^2));
end

clear i

%% Print the error statistics

fprintf('\n');
display(['Mean M: ', num2str(mean([p_mae_fold(:); p_rmse_fold(:)]))]);
display(['Std Dev M: ', num2str(std([p_mae_fold(:); p_rmse_fold(:)]))]);
display(['Mean MAE: ', num2str(mean(MAE_fold))]);
display(['Std Dev MAE: ', num2str(std(MAE_fold))]);
display(['Mean RMSE: ', num2str(mean(RMSE_fold))]);
display(['Std Dev RMSE: ', num2str(std(RMSE_fold))]);

%% Save

% Path
filename = qm_regression_scat_save_name(data_opt, scat_opt, filt_opt, densities);
pathname = strcat(foldername, filename);
clear filename foldername

% Save workspace
save(pathname);