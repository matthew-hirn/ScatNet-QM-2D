% QM_KFOLD_REGRESSION
% Fully learned k-fold regression on the QM data bases, using orthogonal 
% least squares, cross validation, and bagging.
%
% Usage: [p_mae_fold, p_rmse_fold, T_reg_fold_mae, T_reg_fold_rmse, ...
%           res_err_fold_mae, res_err_fold_rmse, T_reg_fold_bag_mae, ...
%           T_reg_fold_bag_rmse] = QM_KFOLD_REGRESSION(T, X, P, max_p, num_bags, bag_per)
%
% Inputs:
%   1.) T (column vector): Function values to regress.
%   2.) X (matrix): Features of data to use for regression.
%   3.) P (cell): Folds, where each cell specifies the indices of the data
%       points in that fold. k = length(P) is the number of folds.
%   4.) max_p (integer): Maximum number of regression terms [optional; 
%       default is max_p = Inf].
%   5.) num_bags (integer): Number of regression models to learn and
%       average the regressed energy over [optional; default is num_bags =
%       10].
%   6.) bag_per (number in [0,1]): Percentage of the training set to use
%       when learning a regression model. The remainder is used for cross
%       validation [optional; default is bag_per = 0.75].
%
% Outputs:
%   1.) p_mae_fold (matrix): Matrix of the maximum number of regression
%       terms when optimized against mean absolute error (MAE), organized
%       by fold (rows) and bag (columns).
%   2.) p_rmse_fold (matrix): Matrix of the maximum number of regression
%       terms when optimized against root mean square error (RMSE), 
%       organized by fold (rows) and bag (columns).
%   3.) T_reg_fold_mae (cell): The regressed energy values by testing fold,
%       averaged over all of the bags, and optimized for mean absolute
%       error (MAE).
%   4.) T_reg_fold_rmse (cell): The regressed energy values by fold,
%       averaged over all of the bags, and optimized for root mean square
%       error (RMSE).
%   5.) res_err_fold_mae (cell): The residual errors between the regressed
%       energies and the true energies, organized by testing fold, and
%       optimzed for mean absolute error (MAE).
%   6.) res_err_fold_rmse (cell): The residual errors between the regressed
%       energies and the true energies, organized by testing fold, and
%       optimzed for root mean square error (MAE).
%   7.) T_reg_fold_bag_mae (cell): The regressed energy values by testing 
%       fold and bag, and optimized for mean absolute error (MAE).
%   8.) T_reg_fold_bag_rmse (cell): The regressed energy values by testing 
%       fold and bag, and optimized for root mean square error (RMSE).
%
% See also:
%   QM_REGRESSION_FOURIER, QM_REGRESSION_SCAT, OLS
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

function [p_mae_fold, p_rmse_fold, T_reg_fold_mae, T_reg_fold_rmse, res_err_fold_mae, ...
    res_err_fold_rmse, T_reg_fold_bag_mae, T_reg_fold_bag_rmse] = ...
    qm_kfold_regression(T, X, P, max_p, num_bags, bag_per)

% Optional inputs
if nargin < 6
    bag_per = 0.75;
end
if nargin < 5
    num_bags = 10;
end
if nargin < 4
    max_p = Inf;
end

% Initialize
k = length(P);                              % Number of folds
pterm_ind = cell(k, num_bags);              % Indices of regression terms in order by fold/bag
max_p_fold = zeros(k, num_bags);            % Max number of regression terms by fold/bag
p_mae_fold = zeros(k, num_bags);            % Best p per fold/bag for MAE
p_rmse_fold = zeros(k, num_bags);           % Best p per fold/bag for RMSE
T_reg_fold_bag_mae = cell(k, num_bags);     % Regressed energy by fold/training bag for MAE
T_reg_fold_bag_rmse = cell(k, num_bags);    % Regressed energy by fold/training bag for RMSE
T_reg_fold_mae = cell(1,k);                 % Regressed energy by fold for MAE
T_reg_fold_rmse = cell(1,k);                % Regressed energy by fold for RMSE
res_err_fold_mae = cell(1, k);              % Residual errors by fold for MAE
res_err_fold_rmse = cell(1, k);             % Residual errors by fold for RMSE

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
    num_train = size(trainX,1);
    
    % Bagging / learning M
    num_bag = round(bag_per * num_train);
    for b=1:num_bags
        
        % Get bag_per points uniformly randomly
        bag_ind = randperm(num_train);
        not_bag_ind = bag_ind((num_bag+1):end);
        bag_ind = bag_ind(1:num_bag);
        bagTrainX = trainX(bag_ind,:);
        bagTestX = trainX(not_bag_ind,:);
        bagTrainT = trainT(bag_ind);
        bagTestT = trainT(not_bag_ind);
        
        % OLS training
        pterm_ind{i,b} = ols(bagTrainT, bagTrainX, max_p);
        max_p_fold(i,b) = length(pterm_ind{i,b});
        bagTrainXp = bagTrainX(:,pterm_ind{i,b});
        bagTestXp = bagTestX(:,pterm_ind{i,b});
        
        % QR decomposition for the bag train/test
        [bagQtrain, bagRtrain] = qr(bagTrainXp);
        bagQtrain = bagQtrain(:,1:max_p_fold(i,b));
        bagRtrain = bagRtrain(1:max_p_fold(i,b),1:max_p_fold(i,b));
        A = bagQtrain' * bagTrainT;
        bagQtest = bagTestXp / bagRtrain;
        
        % OLS regression on bag test data
        T_reg_bag = cell(1,max_p_fold(i,b));
        res_err_bag = cell(1,max_p_fold(i,b));
        for p=1:max_p_fold(i,b)
            T_reg_bag{p} = bagQtest(:,1:p) * A(1:p);
            res_err_bag{p} = bagTestT - T_reg_bag{p};
        end
        
        % MAE and RMSE on bag test as a function of number of regression terms
        bagMAE = zeros(1, max_p_fold(i,b));
        bagRMSE = zeros(1, max_p_fold(i,b));
        for p=1:max_p_fold(i,b)
            bagMAE(p) = mean(abs(res_err_bag{p}));
            bagRMSE(p) = sqrt(mean(res_err_bag{p}.^2));
        end 
        [~, p_mae_fold(i,b)] = min(bagMAE);
        [~, p_rmse_fold(i,b)] = min(bagRMSE);
        
        % QR decomposition for the full train/test
        trainXp = trainX(:,pterm_ind{i,b});
        testXp = testX(:,pterm_ind{i,b});
        [Qtrain, Rtrain] = qr(trainXp);
        Qtrain = Qtrain(:,1:max_p_fold(i,b));
        Rtrain = Rtrain(1:max_p_fold(i,b),1:max_p_fold(i,b));
        A = Qtrain' * trainT;
        Qtest = testXp / Rtrain;
        
        % OLS regression on the test data
        T_reg_fold_bag_mae{i,b} = Qtest(:,1:p_mae_fold(i,b)) * A(1:p_mae_fold(i,b));
        T_reg_fold_bag_rmse{i,b} = Qtest(:,1:p_rmse_fold(i,b)) * A(1:p_rmse_fold(i,b));               
        
    end
    
    % Average regressions over bags to get regression over fold
    T_reg_fold_mae{i} = zeros(length(testT), 1);
    T_reg_fold_rmse{i} = zeros(length(testT), 1);
    for b=1:num_bags
        T_reg_fold_mae{i} = T_reg_fold_mae{i} + T_reg_fold_bag_mae{i,b};
        T_reg_fold_rmse{i} = T_reg_fold_rmse{i} + T_reg_fold_bag_rmse{i,b};
    end
    T_reg_fold_mae{i} = T_reg_fold_mae{i} / num_bags;
    T_reg_fold_rmse{i} = T_reg_fold_rmse{i} / num_bags;
    
    % Compute residual errors
    res_err_fold_mae{i} = testT - T_reg_fold_mae{i};
    res_err_fold_rmse{i} = testT - T_reg_fold_rmse{i};
    
end

end