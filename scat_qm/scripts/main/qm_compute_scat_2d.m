%% QM_COMPUTE_SCAT_2D
%
% Computes the 2D scattering transform for the planar QM data sets. There
% are two such data sets:
%   1.) qm7_2d: The 454 nearly planar molecules from the QM7 data set.
%   2.) qm_2d:  A new data set of 4357 organic, planar molecules.
%
% The outputs of this script which are saved are:
%   1.) data_opt (struct): The inputted data options.
%   2.) scat_opt (struct): The inputted scattering options.
%   3.) filt_opt (struct): The filter options (inputted and default)
%   4.) EU (struct): The invariant scattering coefficients. EU has four
%       fields: dirac, atomic, core, valence, corresponding to the type of
%       approximate non-interacting density used. Each density field is a
%       struct, with two fields: L1, L2, for the type of coefficients.
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

% Parallelization
num_conf_batch = 1;                                 % Number of simultaneous scat computations per batch (batches can be run in parallel)

%% Filters

% Determine Ny
Ny = round((diff(data_opt.y_box) / diff(data_opt.x_box)) * Nx);
data_opt.N = [Nx Ny];
clear Nx Ny

% Construct filterbank
filters = filters_factory_2d(data_opt.N, filt_opt);
filt_opt = filters.translation.meta;

%% Load data

% Load
load(data_opt.data_set, 'R', 'Z');

% Convert coordinates from Bohr to Angstrom
R = 0.529177 * R;

%% Partition the data for parallelization 

% Batch indices and number of batches
batch_ind = create_batch_ind(num_conf_batch, size(R, 1));
num_batch = length(batch_ind);

% Batch positions and charges
R_batch = cell(1,num_batch);
Z_batch = cell(1,num_batch);
for i=1:num_batch
    R_batch{i} = R(batch_ind{i},:,:);
    Z_batch{i} = Z(batch_ind{i},:,:);
end

clear i

%% Compute scattering

% Scattering batches initialization
EU_batch_dirac_L1 = cell(1, num_batch);
EU_batch_dirac_L2 = cell(1, num_batch);
EU_batch_atomic_L1 = cell(1, num_batch);
EU_batch_atomic_L2 = cell(1, num_batch);
EU_batch_core_L1 = cell(1, num_batch);
EU_batch_core_L2 = cell(1, num_batch);
EU_batch_valence_L1 = cell(1, num_batch);
EU_batch_valence_L2 = cell(1, num_batch);

% Loop through batches (can make this loop parfor)
for i=1:num_batch
    
    % Display batch number
    display(strjoin({'Batch number:', num2str(i)}));
        
    % Approximate densities
    rho = qm_approximate_density_2d(data_opt, R_batch{i}, Z_batch{i});
    
    % Wavelet operators
    Wop = wavelet_operator_2d(filters, scat_opt);
    
    % Scattering dirac
    [~,U] = scat(rho.dirac, Wop);
    EU_batch_dirac_L1{i} = expected_scat_light_2d(U, 't', 1);
    EU_batch_dirac_L2{i} = expected_scat_light_2d(U, 't', 2);
    
    % Scattering atomic
    [~,U] = scat(rho.atomic, Wop);
    EU_batch_atomic_L1{i} = expected_scat_light_2d(U, 't', 1);
    EU_batch_atomic_L2{i} = expected_scat_light_2d(U, 't', 2);
    
    % Scattering core
    [~,U] = scat(rho.core, Wop);
    EU_batch_core_L1{i} = expected_scat_light_2d(U, 't', 1);
    EU_batch_core_L2{i} = expected_scat_light_2d(U, 't', 2);
    
    % Scattering valence
    [~,U] = scat(rho.valence, Wop);
    EU_batch_valence_L1{i} = expected_scat_light_2d(U, 't', 1);
    EU_batch_valence_L2{i} = expected_scat_light_2d(U, 't', 2);
        
end

% Get the outputted scattering options
[~, scat_opt] = wavelet_operator_2d(filters, scat_opt);

% Put together scattering from individual batches
EU = struct;
EU.dirac = struct;
EU.atomic = struct;
EU.core = struct;
EU.valence = struct;
EU.dirac.L1 = combine_batch_exp_scat(EU_batch_dirac_L1, batch_ind);
EU.dirac.L2 = combine_batch_exp_scat(EU_batch_dirac_L2, batch_ind);
EU.atomic.L1 = combine_batch_exp_scat(EU_batch_atomic_L1, batch_ind);
EU.atomic.L2 = combine_batch_exp_scat(EU_batch_atomic_L2, batch_ind);
EU.core.L1 = combine_batch_exp_scat(EU_batch_core_L1, batch_ind);
EU.core.L2 = combine_batch_exp_scat(EU_batch_core_L2, batch_ind);
EU.valence.L1 = combine_batch_exp_scat(EU_batch_valence_L1, batch_ind);
EU.valence.L2 = combine_batch_exp_scat(EU_batch_valence_L2, batch_ind);

clear filters i Ri Zi rho Wop U EU_batch_dirac_L1 EU_batch_dirac_L2 ...
    EU_batch_atomic_L1 EU_batch_atomic_L2 EU_batch_core_L1 ...
    EU_batch_core_L2 EU_batch_valence_L1 R_batch EU_batch_valence_L2 ...
    R Z batch_ind num_batch num_conf_batch Z_batch

%% Make 2nd order coefficients invariant to reflections

if scat_opt.O2_invariant
    EU.dirac.L1 = qm_scat_invariant_reflect(EU.dirac.L1);
    EU.dirac.L2 = qm_scat_invariant_reflect(EU.dirac.L2);
    EU.atomic.L1 = qm_scat_invariant_reflect(EU.atomic.L1);
    EU.atomic.L2 = qm_scat_invariant_reflect(EU.atomic.L2);
    EU.core.L1 = qm_scat_invariant_reflect(EU.core.L1);
    EU.core.L2 = qm_scat_invariant_reflect(EU.core.L2);
    EU.valence.L1 = qm_scat_invariant_reflect(EU.valence.L1);
    EU.valence.L2 = qm_scat_invariant_reflect(EU.valence.L2);
end

%% Save

filename = qm_scat_save_name(data_opt, scat_opt, filt_opt);
pathname = strcat(foldername, filename);
save(pathname, 'data_opt', 'scat_opt', 'filt_opt', 'EU');

clear filename foldername pathname