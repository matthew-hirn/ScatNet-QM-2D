% QM_OLS_KFOLD_REGRESSION
% k-fold regression using orthogonal least squares only. 
%
% Usage:
%   [pterm_ind,max_p_fold,T_reg_fold,res_err_fold] = QM_OLS_KFOLD_REGRESSION(T,X,P,max_p)
%
% Inputs:
%   1.) T (column vector): Function values to regress.
%   2.) X (matrix): Features of data to use for regression.
%   3.) P (cell): Folds, where each cell specifies the indices of the data
%       points in that fold. k = length(P) is the number of folds.
%   4.) max_p (integer): Maximum number of regression terms.
%
% Outputs:
%   1.) pterm_ind (cell): Indices of features for each fold as greedily
%       selected by OLS.
%   2.) max_p_fold (vector): The actual max_p for each fold.
%   3.) T_reg_fold (cell): The regressed T over each test fold.
%   4.) res_err_fold (cell): The residual errors over each test fold.
%
% See also:
%   QM_REGRESSION_FOURIER_OLS_ONLY, QM_REGRESSION_SCAT_OLS_ONLY, OLS
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

function [pterm_ind, max_p_fold, T_reg_fold, res_err_fold] = qm_ols_kfold_regression(T, X, P, max_p)

% Initialize
k = length(P);                  % Number of folds
pterm_ind = cell(1,k);          % Indices of regression terms in order
max_p_fold = zeros(1,k);        % Max number of regression terms by fold
T_reg_fold = cell(1,k);         % Regressed energy by fold
res_err_fold = cell(1,k);       % Residual errors by fold

% Loop through the folds
for i=1:k
    
    % Display the fold
    display(strjoin({'Fold number', num2str(i)}));
    
    % i is the testing fold
    noti = 1:k;
    noti(i) = [];
    
    % Training and testing data
    ind = P(noti);
    ind = cat(2,ind{:});
    trainX = X(ind,:);
    testX = X(P{i},:);
    trainT = T(ind);
    testT = T(P{i});
    
    % OLS training
    pterm_ind{i} = ols(trainT, trainX, max_p);
    max_p_fold(i) = length(pterm_ind{i});
    trainXp = trainX(:,pterm_ind{i});
    testXp = testX(:,pterm_ind{i});
    
    % QR decomposition
    [Qtrain, Rtrain] = qr(trainXp);
    Qtrain = Qtrain(:,1:max_p_fold(i));
    Rtrain = Rtrain(1:max_p_fold(i),1:max_p_fold(i));
    A = Qtrain' * trainT;
    Qtest = testXp / Rtrain;
    
    % OLS regression
    T_reg_fold{i} = cell(1,max_p_fold(i));
    res_err_fold{i} = cell(1,max_p_fold(i));
    for p=1:max_p_fold(i)
        T_reg_fold{i}{p} = Qtest(:,1:p) * A(1:p);
        res_err_fold{i}{p} = testT - T_reg_fold{i}{p};
    end
        
end

end