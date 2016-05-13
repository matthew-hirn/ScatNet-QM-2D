%% QM_REGRESSION_FOURIER
%
% Computes Fourier orthogonal least squares regression for the QM data
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
%   9.) coeff_names (cell): List of the Fourier coefficient names.
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
data_opt.N = [2^9 2^9];                             % Number of pixels along x-axis and y-axis

% Densities
densities = {'dirac', 'core', 'valence'};           % Set which densities to use:
                                                    %   (a) dirac (b) atomic (c) core/valence
                                                    %   (d) dirac/atomic (e) dirac/core/valence

% Regression
max_p = 2^10;                                       % Max number of regression terms
num_bags = 10;                                      % Number of bags
bag_per = 0.80;                                     % Size of training bag relative to full training set

%% Load and unwind Fourier features

% Load Fourier
filename = qm_fourier_save_name(data_opt);
load(filename);

% Unwind Fourier
num_densities = length(densities);
X = [];
for i=1:num_densities
    switch densities{i}
        case 'dirac'
            X = cat(2, X, Fk.dirac.L1, Fk.dirac.L2);
        case 'atomic'
            X = cat(2, X, Fk.atomic.L1, Fk.atomic.L2);
        case 'core'
            X = cat(2, X, Fk.core.L1, Fk.core.L2);
        case 'valence'
            X = cat(2, X, Fk.valence.L1, Fk.valence.L2);
    end
end

% Coefficient names
N = min(data_opt.N);
coeff_names = cell(N/2, 2 * num_densities);
for i=1:(2 * num_densities)
    
    % L1 or L2
    switch mod(i,2)
        case 1
            c_name_beg = 'L1^1_r=';
        case 0
            c_name_beg = 'L2^2_r=';
    end

    % Density
    c_name_end = densities{ceil(i/2)};
    
    % Loop through radii
    for r=0:(N/2 - 1)
        coeff_names{r+1,i} = strcat(c_name_beg, num2str(r), '_', c_name_end);
    end
  
end
coeff_names = coeff_names(:);

clear filename Fk N i c_name_beg c_name_end r num_densities

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
filename = qm_regression_fourier_save_name(data_opt, densities);
pathname = strcat(foldername, filename);
clear filename foldername

% Save workspace
save(pathname);