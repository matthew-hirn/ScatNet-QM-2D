% QM_REGRESSION_FOURIER_OLS_ONLY_SAVE_NAME
% Create a file name for saving Fourier regression computations that use
% only orthogonal least squares (OLS).
%
% Usage:
%   filename = QM_REGRESSION_FOURIER_OLS_ONLY_SAVE_NAME(data_opt, densities)
%
% Inputs:
%   1.) data_opt (struct): Data options
%   2.) densities (cell): Approximate densities used
%
% Outputs:
%   1.) filename (string): The file name for saving.
%
% See also:
%   QM_REGRESSION_FOURIER_OLS_ONLY
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

function filename = qm_regression_fourier_ols_only_save_name(data_opt, densities)

filename = strcat(data_opt.data_set, '_regression_ols_only_fourier_N=', ...
    num2str(data_opt.N(1)), 'x', num2str(data_opt.N(2)));
for i=1:length(densities)
    filename = strcat(filename, '_', densities{i});
end
filename = strcat(filename, '.mat');

end