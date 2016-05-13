%% QM_COMPUTE_FOURIER_2D
%
% Computes the 2D Fourier modulus circular average features for planar QM 
% data sets. There are two such data sets:
%   1.) qm7_2d: The 454 nearly planar molecules from the QM7 data set.
%   2.) qm_2d:  A new data set of 4357 organic, planar molecules.
%
% The outputs of this script which are saved are:
%   1.) data_opt (struct): The inputted data options.
%   2.) Fk (struct): The invariant Fourier coefficients. Fk has four
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
data_opt.x_box = [-9,9];                            % Bounding box along x-axis in angstroms (qm7_2d: [-11,11]; qm_2d: [-9,9])
data_opt.y_box = [-9,9];                            % Bounding box along y-axis in angstroms (qm7_2d: [-11,11]; qm_2d: [-9,9])
Nx = 2^9;                                           % Number of pixels along x-axis (y-axis is set to have same sampling rate)

%% Fourier setup

% Determine Ny
Ny = round((diff(data_opt.y_box) / diff(data_opt.x_box)) * Nx);
data_opt.N = [Nx Ny];
clear Nx Ny

% Grid spacing
x_delta = (data_opt.x_box(2) - data_opt.x_box(1)) / (data_opt.N(1) - 1);
y_delta = (data_opt.y_box(2) - data_opt.y_box(1)) / (data_opt.N(2) - 1);

% Center
xc = floor(data_opt.N(1)/2) + 1;
yc = floor(data_opt.N(1)/2) + 1;

% Engine
[X,Y] = meshgrid(1:data_opt.N(1), 1:data_opt.N(2));
R2 = (X-xc).^2 + (Y-yc).^2;

clear X Y

%% Load data

% Load
load(data_opt.data_set, 'R', 'Z');

% Convert coordinates from Bohr to Angstrom
R = 0.529177 * R;

%% Compute Fourier

% Initialize
num_conf = size(R,1);
N = min(data_opt.N);
Fk_dirac_L1 = zeros(num_conf, N/2);
Fk_dirac_L2 = zeros(num_conf, N/2);
Fk_atomic_L1 = zeros(num_conf, N/2);
Fk_atomic_L2 = zeros(num_conf, N/2);
Fk_core_L1 = zeros(num_conf, N/2);
Fk_core_L2 = zeros(num_conf, N/2);
Fk_valence_L1 = zeros(num_conf, N/2);
Fk_valence_L2 = zeros(num_conf, N/2);

% Loop through configurations
for i=1:num_conf
    
    % Display configuration number
    display(strjoin({'Configuration:', num2str(i)}));
    
    % Approximate densities
    rho = qm_approximate_density_2d(data_opt, R(i,:,:), Z(i,:));
    
    % --- Dirac ---
    
    % Fourier modulus transform
    Fki = fft2(rho.dirac);
    Fki = abs(fftshift(Fki));
    Fki = Fki * x_delta * y_delta;
    
    % Radius zero (average)
    Fkic_L1 = zeros(1,N/2);
    Fkic_L2 = zeros(1,N/2);
    Fkic_L1(1) = mean(rho.dirac(:));
    Fkic_L2(1) = mean(rho.dirac(:).^2);
    
    % Loop over radii
    for r=1:(N/2 - 1)
        c = contourc(1:N, 1:N, R2, [0 0]+r^2);
        c = round(c(:,2:end));                          % pixels located ~ on circle
        Fkir = Fki(sub2ind(size(Fki),c(2,:),c(1,:)));   % extract value
        Fkic_L1(r+1) = mean(Fkir);                      % L1 mean        
        Fkic_L2(r+1) = mean(Fkir.^2);                   % L2 mean
    end
    Fk_dirac_L1(i,:) = Fkic_L1;
    Fk_dirac_L2(i,:) = Fkic_L2;
    
    % --- Atomic ---
    
    % Fourier modulus transform
    Fki = fft2(rho.atomic);
    Fki = abs(fftshift(Fki));
    Fki = Fki * x_delta * y_delta;
    
    % Radius zero (average)
    Fkic_L1 = zeros(1,N/2);
    Fkic_L2 = zeros(1,N/2);
    Fkic_L1(1) = mean(rho.atomic(:));
    Fkic_L2(1) = mean(rho.atomic(:).^2);
    
    % Loop over radii
    for r=1:(N/2 - 1)
        c = contourc(1:N, 1:N, R2, [0 0]+r^2);
        c = round(c(:,2:end));                          % pixels located ~ on circle
        Fkir = Fki(sub2ind(size(Fki),c(2,:),c(1,:)));   % extract value
        Fkic_L1(r+1) = mean(Fkir);                      % L1 mean        
        Fkic_L2(r+1) = mean(Fkir.^2);                   % L2 mean
    end
    Fk_atomic_L1(i,:) = Fkic_L1;
    Fk_atomic_L2(i,:) = Fkic_L2;
    
    % --- Core ---
    
    % Fourier modulus transform
    Fki = fft2(rho.core);
    Fki = abs(fftshift(Fki));
    Fki = Fki * x_delta * y_delta;
    
    % Radius zero (average)
    Fkic_L1 = zeros(1,N/2);
    Fkic_L2 = zeros(1,N/2);
    Fkic_L1(1) = mean(rho.core(:));
    Fkic_L2(1) = mean(rho.core(:).^2);
    
    % Loop over radii
    for r=1:(N/2 - 1)
        c = contourc(1:N, 1:N, R2, [0 0]+r^2);
        c = round(c(:,2:end));                          % pixels located ~ on circle
        Fkir = Fki(sub2ind(size(Fki),c(2,:),c(1,:)));   % extract value
        Fkic_L1(r+1) = mean(Fkir);                      % L1 mean        
        Fkic_L2(r+1) = mean(Fkir.^2);                   % L2 mean
    end
    Fk_core_L1(i,:) = Fkic_L1;
    Fk_core_L2(i,:) = Fkic_L2;
    
    % --- Valence ---
    
    % Fourier modulus transform
    Fki = fft2(rho.valence);
    Fki = abs(fftshift(Fki));
    Fki = Fki * x_delta * y_delta;
    
    % Radius zero (average)
    Fkic_L1 = zeros(1,N/2);
    Fkic_L2 = zeros(1,N/2);
    Fkic_L1(1) = mean(rho.valence(:));
    Fkic_L2(1) = mean(rho.valence(:).^2);
    
    % Loop over radii
    for r=1:(N/2 - 1)
        c = contourc(1:N, 1:N, R2, [0 0]+r^2);
        c = round(c(:,2:end));                          % pixels located ~ on circle
        Fkir = Fki(sub2ind(size(Fki),c(2,:),c(1,:)));   % extract value
        Fkic_L1(r+1) = mean(Fkir);                      % L1 mean        
        Fkic_L2(r+1) = mean(Fkir.^2);                   % L2 mean
    end
    Fk_valence_L1(i,:) = Fkic_L1;
    Fk_valence_L2(i,:) = Fkic_L2;
    
end

% Put together Fourier
Fk = struct;
Fk.dirac = struct;
Fk.atomic = struct;
Fk.core = struct;
Fk.valence = struct;
Fk.dirac.L1 = Fk_dirac_L1;
Fk.dirac.L2 = Fk_dirac_L2;
Fk.atomic.L1 = Fk_atomic_L1;
Fk.atomic.L2 = Fk_atomic_L2;
Fk.core.L1 = Fk_core_L1;
Fk.core.L2 = Fk_core_L2;
Fk.valence.L1 = Fk_valence_L1;
Fk.valence.L2 = Fk_valence_L2;

clear num_conf N i rho Fki Fkic_L1 Fkic_L2 r c Fkir Fk_dirac_L1 ...
    Fk_dirac_L2 Fk_atomic_L1 Fk_atomic_L2 Fk_core_L1 Fk_core_L2 ...
    Fk_valence_L1 Fk_valence_L2 x_delta y_delta xc yc R2 R Z

%% Save

filename = qm_fourier_save_name(data_opt);
pathname = strcat(foldername, filename);
save(pathname, 'data_opt', 'Fk');

clear filename foldername pathname