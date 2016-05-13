% QM_SCAT_INVARIANT_REFLECT
% Takes in 2D scattering coefficients invariant to SO(2) actions, and
% returns scattering coefficients invariant to O(2) actions.
%
% Usage:
%   EU_inv_O2 = QM_SCAT_INVARIANT_REFLECT(EU_inv_SO2)
%
% Inputs:
%   1.) EU_inv_SO2 (cell): 1x3 cell of SO(2) invariant scattering 
%       coefficients, as outputted by EXPECTED_SCAT_LIGHT_2D. 
%
% Outputs:
%   1.) EU_inv_O2 (cell): 1x3 cell of O(2) invariant scatteirng
%       coefficients, having the same structure as EU_inv_SO2.
%
% See also:
%   QM_COMPUTE_SCAT, EXPECTED_SCAT_LIGHT_2D
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

function EU_inv_O2 = qm_scat_invariant_reflect(EU_inv_SO2)

% Initialize
EU_inv_O2 = cell(size(EU_inv_SO2));
EU_inv_O2{1} = EU_inv_SO2{1};
EU_inv_O2{2} = EU_inv_SO2{2};
EU_inv_O2{3} = struct;
EU_inv_O2{3}.signal = [];
EU_inv_O2{3}.meta = struct;
EU_inv_O2{3}.meta.j = [];
EU_inv_O2{3}.meta.theta2_minus_theta1 = [];
EU_inv_O2{3}.meta.resolution_t = [];
EU_inv_O2{3}.meta.resolution_r = [];

% Grab settings
L = max(EU_inv_SO2{3}.meta.theta2_minus_theta1) + 1;

% Compute invariant coefficients
ind = 1:size(EU_inv_SO2{3}.signal,2);
while ~isempty(ind)
    
    % Get next index and then remove it
    i1 = ind(1);
    ind(1) = [];
    
    % Get new coefficients
    if EU_inv_SO2{3}.meta.theta2_minus_theta1(i1) == 0 || EU_inv_SO2{3}.meta.theta2_minus_theta1(i1) == L/2
        
        % Only one index
        coeffs = EU_inv_SO2{3}.signal(:,i1);
        theta = EU_inv_SO2{3}.meta.theta2_minus_theta1(i1);
    else
        
        % Get second index and remove it
        j = EU_inv_SO2{3}.meta.j(:,i1);
        theta21 = EU_inv_SO2{3}.meta.theta2_minus_theta1(i1);
        i2 = find(EU_inv_SO2{3}.meta.j(1,:) == j(1) & EU_inv_SO2{3}.meta.j(2,:) == j(2) & ...
            EU_inv_SO2{3}.meta.theta2_minus_theta1 == L - theta21);
        ind(ind == i2) = [];
        
        % Compute coefficients
        coeffs = 0.5 * (EU_inv_SO2{3}.signal(:,i1) + EU_inv_SO2{3}.signal(:,i2));
        theta = min(EU_inv_SO2{3}.meta.theta2_minus_theta1(i1), EU_inv_SO2{3}.meta.theta2_minus_theta1(i2));
        
    end
    
    % Add coefficients to existing ones
    EU_inv_O2{3}.signal = cat(2, EU_inv_O2{3}.signal, coeffs);
    
    % Add meta information
    EU_inv_O2{3}.meta.j = cat(2, EU_inv_O2{3}.meta.j, EU_inv_SO2{3}.meta.j(:,i1));
    EU_inv_O2{3}.meta.theta2_minus_theta1 = cat(2, EU_inv_O2{3}.meta.theta2_minus_theta1, theta);
    EU_inv_O2{3}.meta.resolution_t = cat(2, EU_inv_O2{3}.meta.resolution_t, EU_inv_SO2{3}.meta.resolution_t(i1));
    EU_inv_O2{3}.meta.resolution_r = cat(2, EU_inv_O2{3}.meta.resolution_r, EU_inv_SO2{3}.meta.resolution_r(i1));
end

% Last meta information
EU_inv_O2{3}.meta.layer = EU_inv_SO2{3}.meta.layer;
EU_inv_O2{3}.meta.moment = EU_inv_SO2{3}.meta.moment;
EU_inv_O2{3}.meta.reflect_invariant = true;

end