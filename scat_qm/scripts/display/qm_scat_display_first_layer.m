%% QM_SCAT_DISPLAY_FIRST_LAYER
%
% Displays the covariant part of the first layer of the scattering 
% transform, which are the wavelet modulus coefficients. 
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

% Data
data_opt.data_set = 'qm_2d';                        % qm7_2d or qm_2d
data_opt.x_box = [-9,9];                            % Bounding box along x-axis in angstroms (qm7_2d: [-10,11]; qm_2d: [-9,9])
data_opt.y_box = [-9,9];                            % Bounding box along y-axis in angstroms (qm7_2d: [-7,7]; qm_2d: [-7,7])
Nx = 2^9;                                           % Number of pixels along x-axis (y-axis is set to have same sampling rate)
mol_ind = 1297;                                     % qm_2d, mol_ind = 1297 is in the paper

% Scattering
scat_opt.M = 1;                                     % Number of scattering layers (max is 2)
scat_opt.translation.oversampling = Inf;            % Oversampling for wavelet transform (can take integer values in [0,Inf]; 0: critical sampling rate; Inf: no downsampling)
scat_opt.translation.compute_low_pass = false;      % Do not compute the low pass on top of the scattering features, because we will do a global average anyway

% Filter
filt_opt.translation.n_wavelet_per_octave = 1;      % Number of wavelets per octave/scale                         
filt_opt.translation.J = log2(Nx);                  % Max scale 
filt_opt.translation.L = 8;                         % Number of angles in [0,pi)
filt_opt.translation.xi_psi = 3*pi/4;               % Central frequency of the mother wavelet
filt_opt.translation.slant_psi = 2/3;               % Slant of ellipsoidal envelope
filt_opt.translation.filter_type = 'morlet';        % morlet (for Morlet wavelet; leave as is)

%% Filters

% Determine Ny
Ny = round((diff(data_opt.y_box) / diff(data_opt.x_box)) * Nx);
data_opt.N = [Nx Ny];
clear Nx Ny

% Construct filterbank
filters = filters_factory_2d(data_opt.N, filt_opt);
filt_opt = filters.translation.meta;

% Wavelet operators
Wop = wavelet_operator_2d(filters, scat_opt);

%% Load data

% Load
load(data_opt.data_set, 'R', 'Z');

% Convert coordinates from Bohr to Angstrom
R = 0.529177 * R;

%% Compute densities

R = R(mol_ind,:,:);
Z = Z(mol_ind,:);
rho = qm_approximate_density_2d(data_opt, R, Z);

%% Compute scattering

% Initialize
U = struct;

% Dirac
[~, U.dirac] = scat(rho.dirac, Wop);

% Atonic
[~, U.atomic] = scat(rho.atomic, Wop);

% Valence
[~, U.valence] = scat(rho.valence, Wop);

% Core
[~, U.core] = scat(rho.core, Wop);

%% Create matrix for display

% Pick which density to use to display
U_disp = U.dirac{2}.signal;

% Initialize image
img = zeros(data_opt.N(1) * filt_opt.J, data_opt.N(2) * filt_opt.L);

% Fill in image
for j=1:filt_opt.J
    for l=1:filt_opt.L
        
        img((j-1) * data_opt.N(1) + 1 : j * data_opt.N(1), (l-1) * data_opt.N(2) + 1 : l * data_opt.N(2)) = ...
            U_disp{j}(end:-1:1,end:-1:1,l);
        
    end
end

clear l j

%% Display

imagesc(log2(img+1).^0.5);
axis('image');
xlabel('Angle indices');
ylabel('log_2 Scales');
set(gca,'Xtick',(data_opt.N(2)/2):data_opt.N(2):(data_opt.N(2) * filt_opt.L),'XTickLabel',{1:filt_opt.L});
set(gca,'Ytick',(data_opt.N(1)/2):data_opt.N(1):(data_opt.N(1) * filt_opt.J),'YTickLabel',{0:filt_opt.J-1});