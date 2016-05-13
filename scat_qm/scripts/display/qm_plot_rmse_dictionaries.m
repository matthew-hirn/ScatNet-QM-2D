%% QM_PLOT_RMSE_DICTIONARY
%
% Plots the RMSE errors of the Fourier, wavelet, and scattering
% dictionaries on a log2-log2 plot as a function of M, the number of
% regression terms.
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

% Data
data_opt_scat.data_set = 'qm_2d';                   % qm7_2d or qm_2d
data_opt_scat.N = [512 398];                        % Number of pixels along x-axis and y-axis
data_opt_fourier.data_set = data_opt_scat.data_set; % Same data set for both wavelet/scattering and Fourier
data_opt_fourier.N = [2^9 2^9];                     % Number of pixels along x-axis and y-axis

% Densities
densities = {'dirac', 'core', 'valence'};           % Set which densities to use:
                                                    %   (a) dirac (b) atomic (c) core/valence
                                                    %   (d) dirac/atomic (e) dirac/core/valence

% Scattering
scat_opt.O2_invariant = true;                       % If true, makes scattering coefficients invariant to reflections               
                                                    
% Filter
filt_opt.n_wavelet_per_octave = 1;                  % Number of wavelets per octave/scale                         
filt_opt.J = log2(data_opt_scat.N(1));              % Max scale 
filt_opt.L = 16;                                    % Number of angles in [0,pi)
filt_opt.xi_psi = 3*pi/4;                           % Central frequency of the mother wavelet
filt_opt.slant_psi = 2/3;                           % Slant of ellipsoidal envelope

%% Set coulomb level of performance

if strcmpi(data_opt_scat.data_set, 'qm7_2d')
    coulomb_err = 14.8;
elseif strcmpi(data_opt_scat.data_set, 'qm_2d')
    coulomb_err = 5.4;
end

%% Load results

% Fourier
filename = qm_regression_fourier_ols_only_save_name(data_opt_fourier, densities);
S = load(filename, 'RMSE', 'p_rmse');
RMSE_fourier = S.RMSE;
p_rmse_fourier = S.p_rmse;

% Wavelet
scat_opt.M = 1;
filename = qm_regression_scat_ols_only_save_name(data_opt_scat, scat_opt, filt_opt, densities);
S = load(filename, 'RMSE', 'p_rmse');
RMSE_wavelet = S.RMSE;
p_rmse_wavelet = S.p_rmse;

% Scattering
scat_opt.M = 2;
filename = qm_regression_scat_ols_only_save_name(data_opt_scat, scat_opt, filt_opt, densities);
S = load(filename, 'RMSE', 'p_rmse');
RMSE_scat = S.RMSE;
p_rmse_scat = S.p_rmse;

clear filename S

%% Plot

hold on;

% Fourier
p1 = plot(log2(1:length(RMSE_fourier)), log2(RMSE_fourier), ...
    'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2);
plot(log2(p_rmse_fourier), log2(RMSE_fourier(p_rmse_fourier)), ...
    'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [0, 0.4470, 0.7410], ...
    'MarkerFaceColor', 'w', ...
    'LineWidth', 2);

% Wavelet
p2 = plot(log2(1:length(RMSE_wavelet)), log2(RMSE_wavelet), ...
    'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2);
plot(log2(p_rmse_wavelet), log2(RMSE_wavelet(p_rmse_wavelet)), ...
    'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', 'w', ...
    'LineWidth', 2);

% Scattering
p3 = plot(log2(1:length(RMSE_scat)), log2(RMSE_scat), ...
    'Color', [0.4660, 0.6740, 0.1880], ...
    'LineWidth', 2);
plot(log2(p_rmse_scat), log2(RMSE_scat(p_rmse_scat)), ...
    'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], ...
    'MarkerFaceColor', 'w', ...
    'LineWidth', 2);

% Coulomb
RMSE_coul = coulomb_err * ones(1, length(RMSE_scat));
p4 = plot(log2(1:length(RMSE_scat)), log2(RMSE_coul), 'k--', 'LineWidth', 2);

% Tighten axis
grid on;
v = axis;
v(2) = log2(length(RMSE_scat));
axis(v);

% Axis ticks
set(gca,'Xtick',0:log2(length(RMSE_scat)),'XTickLabel',{0:log2(length(RMSE_scat))});

% Labels
set(gca,'fontsize',16)
xlabel('log_2 M (model complexity)');
ylabel('log_2 RMSE');
h_leg = legend([p1,p2,p3,p4], 'Fourier', 'Wavelet', 'Scattering', 'Coulomb Matrix');
set(h_leg, 'FontSize', 16);

clear v p1 p2 p3 p4