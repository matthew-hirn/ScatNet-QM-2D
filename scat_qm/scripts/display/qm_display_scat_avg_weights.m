%% QM_DISPLAY_SCAT_AVG_WEIGHTS
%
% Displays the average scattering weights for the orthonormal scattering
% coefficients, taken over numerous draws of the training set. The weights
% must first be computed with QM_COMPUTE_SCAT_AVG_WEIGHTS.
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
data_opt.N = [512 398];                             % Number of pixels along x-axis and y-axis

% Densities
densities = {'dirac', 'core', 'valence'};           % Set which densities to use:
                                                    %   (a) dirac (b) atomic (c) core/valence
                                                    %   (d) dirac/atomic (e) dirac/core/valence

% Scattering
scat_opt.M = 2;                                     % Number of scattering layers (max is 2)
scat_opt.O2_invariant = true;                       % If true, scattering coefficients invariant to reflections

% Filter
filt_opt.n_wavelet_per_octave = 1;                  % Number of wavelets per octave/scale                         
filt_opt.J = log2(data_opt.N(1));                   % Max scale 
filt_opt.L = 16;                                    % Number of angles in [0,pi)
filt_opt.slant_psi = 2/3;                           % Slant of ellipsoidal envelope

% Regression
max_p_init = 2004;                                  % Max number of regression terms
num_bags = 1000;                                    % Number of bagging iterations
bag_per = 0.80;                                     % Size of training bag relative to data set.

%% Load weights and other information

% Path
filename = qm_regression_scat_save_name(data_opt, scat_opt, filt_opt, densities);
filename(end-3:end) = [];
filename = strcat(filename, '_num_bags=', num2str(num_bags), '_max_p=', num2str(max_p_init), '.mat');
load(filename);

clear filename

%% Display magnitude of weights by scale


% Set which coefficients to use
A = Aq_mean{512};

% Initialize histograms for each layer
layer_0 = zeros(1,1);
layer_1 = zeros(1, filt_opt.J);
layer_2 = zeros(filt_opt.J-1, filt_opt.J-1);

% Compute L1 coefficient histograms
for p=1:size(X,2)
    
    switch coeff_pars.m(p)
        
        case 0
            layer_0 = layer_0 + abs(A(p));
            
        case 1
            j1 = coeff_pars.j1(p);
            layer_1(j1+1) = layer_1(j1+1) + abs(A(p));
            
        case 2
            j1 = coeff_pars.j1(p);
            j2 = coeff_pars.j2(p);
            layer_2(j2, j1+1) = layer_2(j2, j1+1) + abs(A(p));
            
    end
    
end

% Compression exponent
alpha = 0.80;
layer_1 = layer_1 .^ alpha;
layer_2 = layer_2 .^ alpha;

% Display histograms on same scale
m = min([layer_1(:); layer_2(:)]);
M = max([layer_1(:); layer_2(:)]);

% New figure
figure;

% Wavelet 1st layer
subplot(filt_opt.J+1, filt_opt.J, 1:filt_opt.J);
imagesc(layer_1, [m M]);
set(gca,'fontsize',20)
colormap('pink');
set(gca,'Xtick',1:filt_opt.J,'XTickLabel',{0:filt_opt.J-1});
set(gca,'Ytick',[]);
xlabel('Wavelet scale j');

% Scattering 2nd layer
subplot(filt_opt.J+1, filt_opt.J, 2*filt_opt.J+1:(filt_opt.J+1)*filt_opt.J);
imagesc(layer_2, [m M]);
set(gca,'fontsize',20)
colormap('pink');
set(gca,'Xtick',1:filt_opt.J-1,'XTickLabel',{0:filt_opt.J-2});
xlabel('Scattering first scale j');
ylabel('Scattering second scale j''');
c = colorbar;
tcklbls = linspace(m, M, 6);            % The number 6 needs to be adjusted depending on alpha
tcklbls = round(tcklbls.^(1/alpha));
set(c, 'TickLabels', {tcklbls});
c.Label.String = 'kcal/mol';

clear j1 j2 p A c m M alpha layer_0 layer_1 layer_2 tcklbls
    
%% Plot magnitude of L1 and L2 weights as function of number of regression terms
    
% Set which coefficients to use
A = Aq_mean;
max_p_plot = 512;

% Initialize
l1_mag = zeros(max_p_plot, 1);
l2_mag = zeros(max_p_plot, 1);
l10_mag = zeros(max_p_plot, 1);

% Get percentage L1/L2 for each number of regression terms
for p=1:max_p_plot
    ind = find(A{p} > 0);
    ind1 = coeff_pars.moment(ind) == 1;
    ind2 = coeff_pars.moment(ind) == 2;
    ind0 = ind;
    ind0(ind == 1337) = [];
    ind1_not0 = coeff_pars.moment(ind0) == 1;
    
    % Log2
    l1_mag(p) = log2(sum(A{p}(ind(ind1))));
    l2_mag(p) = log2(sum(A{p}(ind(ind2))));
    l10_mag(p) = log2(sum(A{p}(ind0(ind1_not0))));
end

% Display
figure;
hold on;
p1 = plot(log2(1:max_p_plot), l1_mag, 'LineWidth', 2);
p2 = plot(log2(1:max_p_plot), l10_mag, 'LineWidth', 2);
p3 = plot(log2(1:max_p_plot), l2_mag, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]);
v = axis;
v(2) = log2(max_p_plot);
axis(v);
h_leg = legend([p1, p2, p3], 'L^1', 'L^1 except first weight', 'L^2', 'Location', 'best');
set(h_leg, 'FontSize', 16);
set(gca,'fontsize',16)
xlabel('log_2 M (model complexity)');
ylabel('Aggregate magnitude weights (log_2 scale)')
grid on;

clear A p ind ind0 ind1 ind2 ind1_not0 total p1 p2 p3 h_leg v l1_mag l2_mag l10_mag max_p_plot

%% Plot magnitude of layer 0, layer 1, and layer 2 weights as function of number of regression terms
    
% Set which coefficients to use
A = Aq_mean;
max_p_plot = 512;

% Initialize
layer_0_mag = zeros(max_p_plot, 1);
layer_1_mag = zeros(max_p_plot, 1);
layer_2_mag = zeros(max_p_plot, 1);

% Get avg/wavelet/scattering percentage for each number of regression terms
for p=1:max_p_plot
    ind = find(A{p} > 0);
    ind0 = coeff_pars.m(ind) == 0;
    ind1 = coeff_pars.m(ind) == 1;
    ind2 = coeff_pars.m(ind) == 2;
    
    % Log2
    layer_0_mag(p) = log2(sum(A{p}(ind(ind0))));
    layer_1_mag(p) = log2(sum(A{p}(ind(ind1))));
    layer_2_mag(p) = log2(sum(A{p}(ind(ind2))));
end

% Display
figure;
hold on;
p1 = plot(log2(1:max_p_plot), layer_0_mag, 'LineWidth', 2);
p2 = plot(log2(1:max_p_plot), layer_1_mag, 'LineWidth', 2);
p3 = plot(log2(1:max_p_plot), layer_2_mag, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]);
v = axis;
v(2) = log2(max_p_plot);
axis(v);
h_leg = legend([p1, p2, p3], 'Order 0 (Average)', 'Order 1 (Wavelet)', 'Order 2 (Scattering)', 'Location', 'best');
set(h_leg, 'FontSize', 16);
set(gca,'fontsize', 16)
xlabel('log_2 M (model complexity)');
ylabel('Aggregate magnitude weights (log_2 scale)')
grid on;

clear p ind ind0 ind1 ind2 total p1 p2 p3 h_leg v A layer_0_mag layer_1_mag layer_2_mag max_p_plot

%% Plot magnitude of dirac, core, valence weights as function of number of regression terms
    
% Set which coefficients to use
A = Aq_mean;
max_p_plot = 512;

% Initialize
dirac_mag = zeros(max_p_plot, 1);
core_mag = zeros(max_p_plot, 1);
valence_mag = zeros(max_p_plot, 1);
valence0_mag = zeros(max_p_plot, 1);

% Get avg/wavelet/scattering percentage for each number of regression terms
for p=1:max_p_plot
    ind = find(A{p} > 0);
    ind0 = strcmpi(coeff_pars.density(ind), 'dirac');
    ind1 = strcmpi(coeff_pars.density(ind), 'core');
    ind2 = strcmpi(coeff_pars.density(ind), 'valence');
    indb = ind;
    indb(ind == 1337) = [];
    ind2_not0 = coeff_pars.moment(indb) == 1;
    
    % Log2
    dirac_mag(p) = log2(sum(A{p}(ind(ind0))));
    core_mag(p) = log2(sum(A{p}(ind(ind1))));
    valence_mag(p) = log2(sum(A{p}(ind(ind2))));
    valence0_mag(p) = log2(sum(A{p}(indb(ind2_not0))));
end

% Display
figure;
hold on;
p1 = plot(log2(1:max_p_plot), dirac_mag, 'LineWidth', 2);
p2 = plot(log2(1:max_p_plot), core_mag, 'LineWidth', 2);
p3 = plot(log2(1:max_p_plot), valence_mag, 'LineWidth', 2, 'Color', [0.4660, 0.6740, 0.1880]);
p4 = plot(log2(1:max_p_plot), valence0_mag, 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560]);
v = axis;
v(2) = log2(max_p_plot);
axis(v);
h_leg = legend([p1, p2, p3, p4], 'Dirac', 'Core', 'Valence', 'Valence except first weight', 'Location', 'best');
set(h_leg, 'FontSize', 16);
set(gca,'fontsize',16)
xlabel('log_2 M (model complexity)');
ylabel('Aggregate magnitude weights (log_2 scale)')
grid on;

clear p ind ind0 ind1 ind2 indb ind2_not0 total p1 p2 p3 p4 h_leg v A dirac_mag core_mag valence_mag valence0_mag max_p_plot

%% Plot decay of Aq

% Set how many regression coefficients to plot
max_p_plot = 2^10 + 2^9;

% --- Figure 1 ---
figure;
vals = zeros(num_bags, max_p_plot);
for j=1:num_bags
    vals(j,:) = Aq_mean_bag{j}(1:max_p_plot);
end
x = log2(1:max_p_plot);
y = zeros(3, max_p_plot);
y(1,:) = log2(min(vals));
y(2,:) = log2(mean(vals));
y(3,:) = log2(max(vals));
P2 = polyfit(x, y(2,:), 2);
P2x = P2(1) * x.^2 + P2(2) * x + P2(3);

plot_shaded(x, y, 'b');
hold on;
plot(x, P2x, 'r--', 'LineWidth', 2);
set(gca,'fontsize',16)
xlabel('log_2 m');
ylabel('Magnitude of weight (log_2 scale)')
grid on;
v = axis;
v(2) = log2(max_p_plot);
axis(v);
[~,b] = legend('Range of values', 'Mean value', 'Best fit quadratic polynomial', 'Location', 'best');
PatchInLegend = findobj(b, 'type', 'patch');
set(PatchInLegend, 'facea', 0.2);
set(gca,'Xtick',0:log2(max_p_plot),'XTickLabel',{0:log2(max_p_plot)});

clear vals j x y P2x v a b PatchInLegend