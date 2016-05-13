%% QM_REGRESSION_FOURIER_OLS_ONLY
%
% Computes Fourier orthogonal least squares regression for the QM data
% sets. Does not cross validate the optimal number of regression terms or
% use bagging to reduce the variance.
%
% Important variables created by this script, which are saved:
%   1.) MAE (vector): Vector of mean absolute errors (MAE) indexed by
%       the number of regression terms. The mean is taken over the results
%       of all cross validation folds.
%   2.) RMSE (vector): Vector of root mean square errors (RMSE) indexed by
%       the number of regression terms. The mean is taken over the results
%       of all cross validation folds.
%   3.) p_mae (integer): The number of regression terms that yields the
%       minimum of MAE.
%   4.) p_rmse (integer): The number of regression terms that yields the
%       minimum of RMSE.
%   5.) MAE_fold (vector): Vector of mean absolute errors (MAE) indexed by
%       testing fold. For each fold, the optimal number of regression terms
%       that yields the smallest MAE is selected. 
%   6.) RMSE_fold (vector): Vector of root mean square errors (RMSE) 
%       indexed by testing fold. For each fold, the optimal number of 
%       regression terms that yields the smallest RMSE is selected. 
%   7.) p_mae_fold (vector): Vector of number of regression terms, indexed
%       by testing fold. The kth entry is the number of regression terms
%       that yields the minimum MAE on the kth fold (which is the value 
%       MAE_fold (k)). 
%   8.) p_rmse_fold (vector): Vector of number of regression terms, indexed
%       by testing fold. The kth entry is the number of regression terms
%       that yields the minimum RMSE on the kth fold (which is the value
%       RMSE_fold(k)).
%   9.) res_err_pterm (cell): Residual errors on the database organized by
%       the number of regression terms. The pth cell is a vector with the
%       residual errors between the p-term Fourier regression and the true
%       energies for each molecule.
%  10.) res_err_fold (cell): Residual errors on the database organized by
%       fold. The kth cell contains max_p cells; the pth subcell is a
%       vector with the residual errors between the p-term Fourier
%       regression and the true energies for each molecule in the kth fold.
%  11.) coeff_names (cell): List of the Fourier coefficient names.
%  12.) pterm_ind (cell): Cell of coefficient indices, organized by testing
%       fold. Each cell contains the indices, in order of selection, of the
%       Fourier coefficients selected by the orthogonal least squares
%       algorithm. 
%
% This script will also display the MAE and RMSE as a function of the
% number of regression terms, on a log2-log2 plot, in addition to a stem
% plot of the residual errors for the number of regression terms that
% minimizes the RMSE error. 
%
% The script saves the above variables in addition to some others, as well
% as a PNG of the displayed figure. The user must enter an appropriate
% folder name in order to save the results. 
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

%% OLS regression

% OLS 5-fold regression
[pterm_ind, max_p_fold, T_reg_fold, res_err_fold] = ...
    qm_ols_kfold_regression(T, X, P, max_p);

% Update max_p and get min(max_p_fold)
max_p = max(max_p_fold);
min_max_p = min(max_p_fold);

clear X

%% Errors

% Residual errors organized by number of regression terms
res_err_pterm = cell(1,min_max_p);
for p=1:min_max_p
    res_err_pterm{p} = zeros(1,length(T));
    for i=1:5
        res_err_pterm{p}(P{i}) = res_err_fold{i}{p};
    end
end

% Mean absolute error (MAE) and root mean square error (RMSE) by number of regression terms
MAE = zeros(1,min_max_p);
RMSE = zeros(1,min_max_p);
for p=1:min_max_p
    MAE(p) = mean(abs(res_err_pterm{p}));
    RMSE(p) = sqrt(mean(res_err_pterm{p}.^2));
end

% Minimum mean absolute error and root mean square error by fold
MAE_fold = zeros(1, 5);
RMSE_fold = zeros(1, 5);
p_mae_fold = zeros(1, 5);
p_rmse_fold = zeros(1, 5);
for i=1:5
    mae_fold_i = zeros(1, min_max_p);
    rmse_fold_i = zeros(1, min_max_p);
    for p=1:min_max_p
        mae_fold_i(p) = mean(abs(res_err_fold{i}{p}));
        rmse_fold_i(p) = sqrt(mean(res_err_fold{i}{p}.^2));
    end
    [MAE_fold(i), p_mae_fold(i)] = min(mae_fold_i);
    [RMSE_fold(i), p_rmse_fold(i)] = min(rmse_fold_i);
end

clear p i mae_fold_i rmse_fold_i

%% Plots

% New figure
fig = figure('units','normalized','outerposition',[0 0 1 1]);

% Locations of minima
[~,p_mae] = min(MAE);
[~,p_rmse] = min(RMSE);

% Plot MAE
subplot(2,2,1);
hold on;
plot(log2(1:length(MAE)), log2(MAE), 'LineWidth', 2);
plot(log2(p_mae),log2(MAE(p_mae)), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'w', 'LineWidth', 2);
grid on;
xlabel('Log2 of number of regression terms (M)');
ylabel('Log2 of MAE in kcal/mol');
title('Mean absolute error as a function of number of regression terms');
legend('MAE', ['M=',num2str(p_mae), ', Min(MAE)=',num2str(MAE(p_mae)), ' kcal/mol']);
hold off;

% Plot RMSE
subplot(2,2,2);
hold on;
plot(log2(1:length(RMSE)), log2(RMSE), 'LineWidth', 2);
plot(log2(p_rmse), log2(RMSE(p_rmse)), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'w', 'LineWidth', 2);
grid on;
xlabel('Log2 of number of regression terms (M)');
ylabel('Log2 of RMSE in kcal/mol');
title('Root mean square error as a function of number of regression terms');
legend('RMSE', ['M=',num2str(p_rmse), ', Min(RMSE)=', num2str(RMSE(p_rmse)), ' kcal/mol']);
hold off;

% Plot residual errors
subplot(2,2,[3,4]);
hold on;
stem(1:length(T), res_err_pterm{p_rmse});
v = axis;
v(2) = length(T) + 1;
xlabel('Configurations');
ylabel('Residual error (kcal/mol)');
title(['Residual errors for M=',num2str(p_rmse)]);
axis(v);
hold off;

clear v

%% Save

% Path
filename = qm_regression_fourier_ols_only_save_name(data_opt, densities);
pathname = strcat(foldername, filename);
clear filename foldername

% Save figure
set(gcf,'PaperPositionMode','auto')
print(fig, strcat(pathname(1:end-4),'.png'), '-dpng', '-r0');
clear fig

% Save workspace
save(pathname);