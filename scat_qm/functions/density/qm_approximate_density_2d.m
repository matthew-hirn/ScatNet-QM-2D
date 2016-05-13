% QM_APPROXIMATE_DENSITY_2D
% Computes the approximate densities for the 2D QM molecules.
%
% Usage:
%   rho = QM_APPROXIMATE_DENSITY_2D(data_opt, R, Z)
%
% Inputs:
%   1.) data_opt (struct): Data options.
%   2.) R (matrix): Positions of the atoms, can include multiple molecules.
%   3.) Z (matrix): Charges of the atoms, can include multiple molecules.
%
% Outputs:
%   1.) rho (struct): Struct with four fields: dirac, atomic, core, 
%       valence, for each of the approximate density types.
%
% See also:
%   QM_COMPUTE_FOURIER_2D, QM_COMPUTE_SCAT_2D
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

function rho = qm_approximate_density_2d(data_opt, R, Z)

%% Grid for the molecule

x = linspace(data_opt.x_box(1), data_opt.x_box(2), data_opt.N(1));
y = linspace(data_opt.y_box(1), data_opt.y_box(2), data_opt.N(2));

%% List of atoms, charges, valence charges

% List of atom names
atom_list = cell(1,6);
atom_list{1} = 'Hydrogen';
atom_list{2} = 'Carbon';
atom_list{3} = 'Nitrogen';
atom_list{4} = 'Oxygen';
atom_list{5} = 'Sulfur';
atom_list{6} = 'Chlorine';

% List of atom total charges
Z_atom = zeros(1,6);
Z_atom(1) = 1;
Z_atom(2) = 6;
Z_atom(3) = 7;
Z_atom(4) = 8;
Z_atom(5) = 16;
Z_atom(6) = 17;

% List of atom core charges
Z_atom_cor = zeros(1,6);
Z_atom_cor(1) = 0;
Z_atom_cor(2) = 2;
Z_atom_cor(3) = 2;
Z_atom_cor(4) = 2;
Z_atom_cor(5) = 10;
Z_atom_cor(6) = 10;

% List of atom valence charges
Z_atom_val = zeros(1,6);
Z_atom_val(1) = 1;
Z_atom_val(2) = 4;
Z_atom_val(3) = 5;
Z_atom_val(4) = 6;
Z_atom_val(5) = 6;
Z_atom_val(6) = 7;

%% Load atomic densities

rho_atom = cell(1,6);
for i=1:6
    
    % Initialize
    rho_atom{i} = struct;
    
    % Dirac
    filename = strcat(atom_list{i}, '_dirac_density.mat');
    S = load(filename, 'dirac_density');
    rho_atom{i}.dirac = S.dirac_density;
    
    % Atomic
    filename = strcat(atom_list{i}, '_cylindrical_density.mat');
    S = load(filename, 'cylindrical_density');
    rho_atom{i}.atomic = S.cylindrical_density;
    
    % Core
    filename = strcat(atom_list{i}, '_cylindrical_density_core.mat');
    S = load(filename, 'cylindrical_density');
    rho_atom{i}.core = S.cylindrical_density;
    
    % Valence
    filename = strcat(atom_list{i}, '_cylindrical_density_valence.mat');
    S = load(filename, 'cylindrical_density');
    rho_atom{i}.valence = S.cylindrical_density;
    
end

clear S

%% Reference world coordinates

% Reference world coordinates for atomic densities
xWorldLimits = [-2.5 2.5];
yWorldLimits = [-2.5 2.5];
ref_rho_atom = imref2d([501 501], xWorldLimits, yWorldLimits);

% Reference world coordinates for molecules
ref_rho_molecule = imref2d(data_opt.N, data_opt.y_box, data_opt.x_box);

%% Compute densities for each molecule

% Pre-allocate rho
rho = struct;
rho.dirac = zeros(data_opt.N(1), data_opt.N(2), size(R,1));
rho.atomic = rho.dirac;
rho.core = rho.dirac;
rho.valence = rho.dirac;

% Loop over molecules
for i=1:size(R,1)
    
    % Positions and charges of atoms in molecule
    Ri = squeeze(R(i,:,:));
    Zi = Z(i,:);
    ind = Zi > 0;
    Zi = Zi(ind);
    Ri = Ri(ind,:);
    num_atoms = size(Ri,1);
    
    % Loop over atoms
    for j=1:num_atoms
        
        % Z index
        Z_ind = find(Z_atom == Zi(j));
        
        % Translation (encoded as affine transformation)
        pos = Ri(j,:);
        Ab = zeros(3,3);
        Ab(1:2,1:2) = eye(2);
        Ab(3,1:2) = [pos(2) pos(1)];
        Ab(3,3) = 1;
        tform = affine2d(Ab);
        
        %--- Dirac ---
        
        % Index position of dirac
        [~, x_ind] = min(abs(x - pos(1)));
        [~, y_ind] = min(abs(y - pos(2)));
        
        % Translated dirac
        rho_atom_translate = zeros(data_opt.N);
        rho_atom_translate(x_ind,y_ind) = Z_atom(Z_ind);
        
        % Normalize to make sure integral is unchanged
        I = trapz(y, rho_atom_translate, 2);
        I = trapz(x, I);
        rho_atom_translate = (Z_atom(Z_ind)/I) * rho_atom_translate;
        
        % Add into dirac density
        rho.dirac(:,:,i) = rho.dirac(:,:,i) + rho_atom_translate;
        
        % --- Atomic ---
        
        % Translate atomic density to position in molecule
        rho_atom_translate = imwarp(rho_atom{Z_ind}.atomic, ...
            ref_rho_atom, tform, 'OutputView', ref_rho_molecule);
        
        % Normalize to make sure integral is unchanged
        I = trapz(y, rho_atom_translate, 2);
        I = trapz(x, I);
        rho_atom_translate = (Z_atom(Z_ind)/I) * rho_atom_translate;
        
        % Add into atomic density
        rho.atomic(:,:,i) = rho.atomic(:,:,i) + rho_atom_translate;
        
        % --- Valence ---
        
        % Translate valence atomic density to position in molecule
        rho_atom_translate = imwarp(rho_atom{Z_ind}.valence, ...
            ref_rho_atom, tform, 'OutputView', ref_rho_molecule);
        
        % Normalize to make sure integral is unchanged
        I = trapz(y, rho_atom_translate, 2);
        I = trapz(x, I);
        rho_atom_translate = (Z_atom_val(Z_ind)/I) * rho_atom_translate;
        
        % Add into valence density
        rho.valence(:,:,i) = rho.valence(:,:,i) + rho_atom_translate;
        
        % --- Core ---
        
        if Z_atom_cor(Z_ind) > 0
            
            % Translate core atomic density to position in molecule
            rho_atom_translate = imwarp(rho_atom{Z_ind}.core, ...
                ref_rho_atom, tform, 'OutputView', ref_rho_molecule);
            
            % Normalize to make sure integral is unchanged
            I = trapz(y, rho_atom_translate, 2);
            I = trapz(x, I);
            rho_atom_translate = (Z_atom_cor(Z_ind)/I) * rho_atom_translate;
            
            % Add into valence density
            rho.core(:,:,i) = rho.core(:,:,i) + rho_atom_translate;
            
        end
        
    end
    
end

end