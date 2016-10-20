%% QM_SCAT_TIME_BENCHMARK_2D
%
% Benchmarks the time it takes to compute the scattering coefficients of 2D
% planar molecules from QM data bases.
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
data_opt.data_set = 'qm_2d';                        % qm7_2d or qm_2d
data_opt.x_box = [-9,9];                            % Bounding box along x-axis in angstroms (qm7_2d: [-10,11]; qm_2d: [-9,9])
data_opt.y_box = [-7,7];                            % Bounding box along y-axis in angstroms (qm7_2d: [-7,7]; qm_2d: [-7,7])
Nx = 2^9;                                           % Number of pixels along x-axis (y-axis is set to have same sampling rate)

% Scattering
scat_opt.M = 2;                                     % Number of scattering layers (max is 2)
scat_opt.O2_invariant = true;                       % If true, coefficients are O(2) invariant, if false, they are SO(2) invariant
scat_opt.translation.oversampling = 1;              % Oversampling for wavelet transform (can take integer values in [0,Inf]; 0: critical sampling rate; Inf: no downsampling)
scat_opt.translation.compute_low_pass = false;      % Do not compute the low pass on top of the scattering features, because we will do a global average anyway

% Filter
filt_opt.translation.n_wavelet_per_octave = 1;      % Number of wavelets per octave/scale                         
filt_opt.translation.J = log2(Nx);                  % Max scale 
filt_opt.translation.L = 16;                        % Number of angles in [0,pi)
filt_opt.translation.xi_psi = 3*pi/4;               % Central frequency of the mother wavelet
filt_opt.translation.slant_psi = 2/3;               % Slant of ellipsoidal envelope
filt_opt.translation.filter_type = 'morlet';        % morlet (for Morlet wavelet; leave as is)

%% Initialize time struct

time_benchmark = struct;
time_benchmark.filters = zeros(1,1);
time_benchmark.rho = zeros(size(R, 1), 1);
time_benchmark.U = zeros(size(R, 1), 1);
time_benchmark.EU = zeros(size(R, 1), 1);
time_benchmark.total = zeros(size(R, 1), 1);

%% Filters

% Determine Ny
Ny = round((diff(data_opt.y_box) / diff(data_opt.x_box)) * Nx);
data_opt.N = [Nx Ny];
clear Nx Ny

% Construct filterbank
tstart = tic;
filters = filters_factory_2d(data_opt.N, filt_opt);
filt_opt = filters.translation.meta;
time_benchmark.filters = toc(tstart);

%% Load data

% Load
load(data_opt.data_set, 'R', 'Z');

% Convert coordinates from Bohr to Angstrom
R = 0.529177 * R;

%% Compute scattering and time it

% Indices
batch_ind = create_batch_ind(1, size(R, 1));
num_batch = length(batch_ind);

% Loop through batches (can make this loop parfor)
EU_L1 = cell(1, num_batch);
EU_L2 = cell(1, num_batch);
for i=1:size(Z,1)
    
    % Display batch number
    display(strjoin({'Molecule number:', num2str(i)}));
        
    % Approximate densities
    tstart = tic;
    rho = qm_approximate_density_2d_atomic_only(data_opt, R(i,:,:), Z(i,:));
    time_benchmark.rho(i) = toc(tstart);
    
    % Wavelet operators and scattering U transform
    tstart = tic;
    Wop = wavelet_operator_2d(filters, scat_opt);
    [~,U] = scat(rho.atomic, Wop);
    time_benchmark.U(i) = toc(tstart);
    
    % Expected scattering coefficients
    tstart = tic;
    EU_L1{i} = expected_scat_light_2d(U, 't', 1);
    EU_L2{i} = expected_scat_light_2d(U, 't', 2);
    time_benchmark.EU(i) = toc(tstart);
        
end
time_benchmark.total = time_benchmark.rho + time_benchmark.U + time_benchmark.EU;

clear R Z i tstart rho Wop U EU

%% Display average times

display(['Average time to compute rho: ', num2str(mean(time_benchmark.rho))]);
display(['Average time to compute scattering propogator: ', num2str(mean(time_benchmark.U))]);
display(['Average time to compute invariant scattering coefficients: ', num2str(mean(time_benchmark.EU))]);
display(['Average total time: ', num2str(mean(time_benchmark.total))]);