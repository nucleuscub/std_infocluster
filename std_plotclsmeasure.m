% std_plotclsmeasure() - Plot cluster measures
%
% Usage:
%   >>  ;
%
% Inputs:
%      datmeasures      - Structure with statas for each measure. See
%                         output from std_clustmat.m
%      iplot            - Index od the measure to plot
%    
% Outputs:
%    
%
% See also: 
%   std_infocluster, std_clustmat
%
% Author: Ramon Martinez-Cancino, SCCN, 2014
%
% Copyright (C) 2014  Ramon Martinez-Cancino, 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function std_plotclsmeasure(datmeasures,iplot)

% Settings
try
    icadefs;
catch
    colormap('jet');                                                           % Def colormap
    LABEL_FONTSIZE = 10;
    BACKEEGLABCOLOR =[.66 .76 1];
    THRESHOLD_BIN = 16;
end
SLABEL_FONTSIZE = 12;

fignames = {'Mean Distance from Centroid', 'Summatory of Distances to Centroid'};

figure('name',fignames{iplot}, ...
       'color', [.66 .76 1],...
       'Tag','clusterinfo_plot2',...
       'numbertitle', 'off');
hold on;

namemeasure = fieldnames(datmeasures);
nplots = length(namemeasure);

switch iplot
    case 1 % centroid_distmean and stdv
        for i = 1:nplots
            subplot(nplots,1,i);
            errorbar([datmeasures.(namemeasure{i}).centroid_distmean],[datmeasures.(namemeasure{i}).centroid_diststd],'rx');
            xlim([1,length([datmeasures.(namemeasure{i}).centroid_distmean])]);
            ylabel({namemeasure{i};'a/u'},'fontsize', SLABEL_FONTSIZE,'fontweight','bold');
            if i == nplots
                xlabel('Cluster','fontsize', SLABEL_FONTSIZE,'fontweight','bold');
            end
           grid on; 
        end
    case 2
        % under development
end