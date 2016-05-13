% QM_FOURIER_SAVE_NAME
% Create a file name for saving fourier features.
%
% Usage:
%   filename = QM_FOURIER_SAVE_NAME(data_opt)
%
% Inputs:
%   1.) data_opt (struct): Data options
%
% Outputs:
%   1.) filename (string): The file name for saving.
%
% See also:
%   QM_COMPUTE_FOURIER_2D
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

function filename = qm_fourier_save_name(data_opt)

filename = strcat(data_opt.data_set, '_fourier_N=', ...
    num2str(data_opt.N(1)), 'x', num2str(data_opt.N(2)), '.mat');

end