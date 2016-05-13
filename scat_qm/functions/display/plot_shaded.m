% PLOT_SHADED
% Creates a plot with a shaded region.
%
% Usage:
%   PLOT_SHADED(x, y, fstr)
%
% Inputs:
%   1.) x: 1xN vector of x coordinates
%   2.) y: Either 1xN, 2xN, or 3xN matrix of y data to plot. The three
%       cases yield the following results:
%       a.) 1xN: Plot of y against x.
%       b.) 2xN: Plot of shaded region between y(1,:) and y(2,:) against x.
%       c.) 3xN: Plot of y(2,:) against x, and plot of shaded region
%           between y(1,:) and y(3,:) against x.
%   3.) fstr: String with plotting format using in PLOT, for example 'r' or
%       'b--', etc.
%
% See also:
%   QM_DISPLAY_SCAT_AVG_WEIGHTS
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

function plot_shaded(x,y,fstr)
 
% Transpose y if necessary
if size(y,1)>size(y,2)
    y=y';
end;

% Just plot one line
if size(y,1)==1 
    plot(x,y,fstr);
end;

% Plot shaded area
if size(y,1)==2 
    px=[x,fliplr(x)]; % make closed patch
    py=[y(1,:), fliplr(y(2,:))];
    patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
end;

% Plot line and shaded area
if size(y,1)==3 
    px=[x,fliplr(x)];
    py=[y(1,:), fliplr(y(3,:))];
    hold on;
    patch(px,py,1,'FaceColor',fstr,'EdgeColor','none','FaceAlpha',0.2);
    plot(x,y(2,:),fstr,'LineWidth',2);
    hold off;
end;