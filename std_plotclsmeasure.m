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

function std_plotclsmeasure(STUDY,parentcluster,datmeasures,iplot)

% Settings
try
    icadefs;
catch
    colormap('jet');                                                           % Def colormap
    LABEL_FONTSIZE = 10;
    BACKEEGLABCOLOR =[.66 .76 1];
    THRESHOLD_BIN = 16;
end
LABEL_FONTSIZE  = 14;
SLABEL_FONTSIZE = 12;
zflag = 1;
   parentcluster = deblank(parentcluster);

% Getting clusters names
hits_temp = cellfun(@(x)strcmp(x,parentcluster),{STUDY.cluster.name});
parent_indx = find(hits_temp);

% Getting cls
for i = 1:length({STUDY.cluster.name})
    tmpval = STUDY.cluster(i).parent;
    if isempty(tmpval)
        hits_tmp(i) = 0;
    else
        hits_tmp(i) = strcmp(tmpval{1},deblank(parentcluster));
    end
end
cls = find(hits_tmp); clear hits_tmp tmpval


% cls = (parent_indx+1):(parent_indx + length(STUDY.cluster(parent_indx).child));
clsname = {STUDY.cluster(cls).name}';

% for i=1:length(cls)
%     clsname{i} = num2str(cls(i));
% end


fignames = {'Mean Distance from Centroid (z-score)', 'Summatory of Distances to Centroid'};

figure('name',fignames{iplot}, ...
       'color', [.66 .76 1],...s
       'Tag','clusterinfo_plot2',...
       'numbertitle', 'off',...
       'Units', 'Normalized',...
       'Position', [0.0925 0.1308 0.6769 0.7342]);
hold on;


namemeasure = fieldnames(datmeasures);
nplots = length(namemeasure);
plotorder = [find(~strcmp(namemeasure,'jointspace'))' find(strcmp(namemeasure,'jointspace'))];

switch iplot
    case 1 % centroid_distmean and stdv
        for i = 1:nplots
            handles.(['axes_' num2str(i)]) = subplot(nplots,1,i);
            m1  = [datmeasures.(namemeasure{plotorder(i)}).centroid_distmean];
            m2  = [datmeasures.(namemeasure{plotorder(i)}).centroid_diststd];
            h   = errorbar(m1/std(m1),m2/std(m2),'rx');

            if i == 1, set(get(gca,'Title'),'String','Mean Distance from Centroid \pm Std', 'Interpreter', 'tex','Fontsize',LABEL_FONTSIZE); end
            set(gca, 'XTickLabel', []);
            set(gca, 'XTick', [1:length([datmeasures.(namemeasure{plotorder(i)}).centroid_distmean])]);
            hold on;
            plot([0.5:length([datmeasures.(namemeasure{plotorder(i)}).centroid_distmean])+0.5],zeros(1,length([datmeasures.(namemeasure{plotorder(i)}).centroid_distmean])+1),'LineStyle',':','LineWidth',0.1,'Color','b');
            set(h,'LineWidth',1,'MarkerSize',6, 'Marker', 's', 'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);
            xlim([0.5,length([datmeasures.(namemeasure{plotorder(i)}).centroid_distmean])+0.5]);
            if i ~=nplots
            ylabel({['Measure: ' upper(namemeasure{plotorder(i)})];'Mean Distance'},'fontsize', LABEL_FONTSIZE,'fontweight','bold');
            end
            set(gca,'fontsize',SLABEL_FONTSIZE);
            if i == nplots
                ylabel('Joint measures space','fontsize', LABEL_FONTSIZE-2,'fontweight','bold');
                set(gca, 'XTickLabel', clsname,'fontsize',SLABEL_FONTSIZE);
                
                %--
                % Rotating Xticklabels
                lasthandle = handles.(['axes_' num2str(nplots)]);
                Xt = get(lasthandle,'XTick');
                Xl = get(lasthandle,'XLim');
                ax = axis;                 % Current axis limits
                axis(axis);                % Set the axis limit modes (e.g. XLimMode) to manual
                Yl = ax(3:4);              % Y-axis limits
                
                % Place the text labels
                t = text(Xt,Yl(1)*ones(1,length(Xt)),get(lasthandle,'xticklabel'));
                set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
                
                set(lasthandle,'XTickLabel','');  % Remove the default labels
                
                % Get the Extent of each text object.This loop is unavoidable.
                for i = 1:length(t)
                    ext(i,:) = get(t(i),'Extent');
                end
                % Determine the lowest point.  The X-label will be placed so that the top is aligned with this point.
                LowYPoint = min(ext(:,2));
                
                % Place the axis label at this point
                XMidPoint = Xl(1)+abs(diff(Xl))/2;
                text (XMidPoint,LowYPoint,'Cluster','VerticalAlignment','top',...
                    'HorizontalAlignment','center',...
                    'FontWeight','bold',...
                    'Tag','Xtext',...
                    'FontSize',LABEL_FONTSIZE);
                handles.Xtext = findobj(0,'Tag','Xtext');
                %--
                
            end
            grid on;
        end
    case 2
        % under development
end