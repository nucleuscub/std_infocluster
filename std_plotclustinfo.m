% std_plotclustinfo() - Plot relation of Clusters  Vs Subject/Dataset Vs IC
%
% Usage:
%   >>  [mat2plot, UniqSubjInd] = std_plotclustinfo(STUDY,SubjClusIC_Matrix,parentcluster)
%
% Inputs:
%      STUDY               - studyset structure containing some or all files in ALLEEG
%      SubjClusIC_Matrix   - 3D binary matrix of incidence of Clusters Vs Datasets/Subjects Vs IC's
%                            see parse_clustinfo.m
%
% Optional inputs:
%
%
% plot         - [0/1] 1 to plot matrix , otherwise just provide the matrix
%                 to be plotted for later analysis. Default 1.
% figlabel     - [0/1] 1 to chow "ICs" in the label's cell, otherwise does
%                not show the label. Default 1
%
%
% Outputs:
%    mat2plot  - 2D matrix with the info in SubjClusIC_Matrix but with the
%                dimension of IC collapsed. Also, it only account for Subjects/Datasets
%                from different session. Means, that if we have one subject with 2
%                datasets from the same session, it will only show the Subject/Dataset only 1
%                time.
%  UniqSubjInd - Index of the Subjects/Datasets used for the analysis,
%                after detect sessions
%
% See also:
%   std_plotclustinfo

%
% Author: Ramon Martinez-Cancino, SCCN, 2014
%
% Copyright (C) 2014  Ramon Martinez-Cancino,INC, SCCN
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
function [mat2plot, UniqSubjInd] = std_plotclustinfo(STUDY,SubjClusIC_Matrix,parentcluster,varargin)

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('std_infocluster() error: calling convention {''key'', value, ... } error'); return;
end;

try g.plot;              catch, g.plot            = 1;        end; % By default plot
try g.figlabel;          catch, g.figlabel        = 1;        end; % Show labels in fig
try g.colormapstyle;     catch, g.colormapstyle   = 'Reds';   end; %Selct colormap style from brewermap
try g.sortdim1flag ;     catch, g.sortdim1flag    = 0;       end
try g.sortdim2flag ;     catch, g.sortdim2flag    = 0;       end

% Settings
icadefs;
LABEL_FONTSIZE = 15;
SLABEL_FONTSIZE = 12;
color = BACKEEGLABCOLOR;
THRESHOLD_BIN = 30;
parentcluster = deblank(parentcluster);

% Getting clusters names
hits_tmp = cellfun(@(x)strcmp(x,parentcluster),{STUDY.cluster.name});
parent_indx = find(hits_tmp); clear hits_tmp;

% Getting cls
for i = 1:length({STUDY.cluster.name})
    tmpval = STUDY.cluster(i).parent;
    if isempty(tmpval)
        hits_tmp(i) = 0;
    else
        hits_tmp(i) = strcmp(tmpval{1},parentcluster);
    end
end
cls = find(hits_tmp); clear hits_tmp tmpval;
cls_names = {STUDY.cluster(cls).name}';

% Check session in data or take from input (now we assume data is comming with session info)
sessions  = cellfun(@(x) x,{STUDY.datasetinfo.session});
if iscell(sessions)
    sessions = cellfun(@(x) x,sessions);
end

% ------------- Check for data in different sessions -------------

% Getting Sessions index
for i = 1: max(sessions)
    SessionIndex{i}   = find(i == sessions);
end

% Getting index for each session
for session_i = 1:size(SessionIndex,2)
    
    SessionIndx_tmp =  SessionIndex{session_i};
    
    uniqsubj = unique({STUDY.datasetinfo(SessionIndx_tmp).subject});
    
    % Only 1 subject set / session
    if isequal(length(uniqsubj),length({STUDY.datasetinfo((SessionIndx_tmp)).subject}))
        UniqSubjInd_tmp{session_i} = cell2mat({STUDY.datasetinfo(SessionIndx_tmp).index});
        mat2plot_tmp{session_i}    = squeeze(sum(SubjClusIC_Matrix(:,SessionIndx_tmp,:),3));
        
        % More than 1 subject set / session
    else
        for i = 1: length(uniqsubj)
            booleanIndex               = strcmp(uniqsubj(i), {STUDY.datasetinfo(SessionIndx_tmp).subject});
            integerIndex               = find(booleanIndex);
            Ind_tmp(i)                 = integerIndex(1);
            UniqSubjInd_tmp{session_i} = Ind_tmp;
        end
        
        mat2plot_tmp{session_i} = squeeze(sum(SubjClusIC_Matrix(:,UniqSubjInd_tmp{session_i},:),3));
    end
end

mat2plot    = cell2mat(mat2plot_tmp);
UniqSubjInd = cell2mat(UniqSubjInd_tmp);
 if g.sortdim1flag == 1
     [trash,sortdim1_indx] = sort(sum(mat2plot,2));
     mat2plot              = mat2plot(sortdim1_indx,:);
     cls_names             = cls_names(sortdim1_indx);
 end
 
  if g.sortdim2flag == 1
      [trash,sortdim2_indx] = sort(sum(mat2plot,1));
      mat2plot              = mat2plot(:,sortdim2_indx);
      UniqSubjInd           = UniqSubjInd(sortdim2_indx);
      sessions              = sessions(sortdim2_indx);
 end

%**************************************************************************
%**************************************************************************
%**************************************************************************
if g.plot
    ncmap = max(max(mat2plot));
    [mycolormap,scheme] = brewermap(ncmap,g.colormapstyle);
    
    % Looking for the figure if open
    htmp = findall(0,'Type','Figure', 'Tag','clusterinfo_plot1');
    
    % Making FIG
     if isempty(htmp)
    handles.main = figure('Name'       ,['Cluster Info: ' parentcluster],...
                          'numbertitle','off',...
                          'Color'      , color,...
                          'Tag'        ,'clusterinfo_plot1',...
                          'Units'      ,'Normalized',...
                          'Position'   , [0.1125 0.0444 0.6778 0.8300],...
                          'ToolBar'    ,'none');
     else
         handles.main  = htmp;
         axeslist = findobj(handles.main,'Type','axes');
          for iaxes = 1: length(axeslist)
              cla(axeslist(iaxes));
          end
     end
    %
    
    left = .02; bottom = 0.20; width = .80; height = 0.8;
    pos  = get(gca,'position');
    q    = [pos(1) pos(2) 0 0];
    s    = [pos(3) pos(4) pos(3) pos(4)];
    axis off;
    
    % Saving data
    setappdata(handles.main,'STUDY',STUDY);
    setappdata(handles.main,'SubjClusIC_Matrix',SubjClusIC_Matrix);
    setappdata(handles.main,'g',g);
    setappdata(handles.main,'parentcluster',parentcluster);
    
    % FIGURE 1 ----------------------------------------------------------------
    htmp = findobj(handles.main,'Tag','axes1');
    if isempty(htmp)
        handles.fig_1 = axes('Position',[left bottom width height].*s+q);
    else
        handles.fig_1 = htmp;
        axes(handles.fig_1 );
    end
        
    imagesc(mat2plot); % plot fig1
    set( handles.fig_1,'tag','axes1');
    set(get(handles.fig_1,'Title'),'String', {'\fontsize{15}#ICs/Subject/Cluster'; ...
                                             ['\fontsize{12}#ICs: ' num2str(sum(sum(sum(SubjClusIC_Matrix)))) ' components'];...
                                             ['\fontsize{12}#Clusters: ' num2str(size(SubjClusIC_Matrix,1)) ' clusters (cls)']},'FontSize',LABEL_FONTSIZE); % Title
   
    if max(size(mat2plot))< THRESHOLD_BIN && g.figlabel
        % Creating labels for grid
        for i= 1:size(mat2plot,1)
            for j = 1:size(mat2plot,2)
                if (mat2plot(i,j) ~= 0)
                    gridlabel{i,j} = [num2str(mat2plot(i,j),'%d')];
                else
                    gridlabel{i,j} = num2str(' ','%d ');
                end
            end
        end
        
        gridlabel = strtrim(cellstr(gridlabel));                                                % Remove any space padding
        [x,y]     = meshgrid(1:size(mat2plot,2),1:size(mat2plot,1));                            % Create x and y coordinates for the strings
        
        text(x(:),y(:),gridlabel(:),'HorizontalAlignment','center','FontSize',SLABEL_FONTSIZE); % Plot the strings (,'FontWeight','bold')
    end
    
    handles.Ytext = ylabel('Cluster','FontWeight','bold',...
                                    'FontSize',LABEL_FONTSIZE,...
                                    'Units','Normalized');                       % Y Labels
    
    set(handles.fig_1,'XLim',[0.5 size(mat2plot,2) + 0.5]);                                           % Setting the X limits
    set(handles.fig_1,'YLim',[0.5 size(mat2plot,1) + 0.5]);                                           % Setting the Y limits
    
    set(handles.fig_1, 'Xtick',1:size(mat2plot,2));                                                   % Setting the X ticks
    set(handles.fig_1, 'Ytick',1:size(mat2plot,1)) ;                                                  % Setting the Y ticks
    
    set(handles.fig_1,'YTickLabel',cls_names,'FontSize',SLABEL_FONTSIZE);                              % Note: need to pick the real plotted subj
    set(handles.fig_1,'XTickLabel',''); 
    
    % Axis
    tmp_caxis = caxis;                                                     % Color plot
    caxis([0 max(tmp_caxis)]);                                             % ...
    colormap(mycolormap);                                                  % ...
    if max(size(mat2plot))< THRESHOLD_BIN
        grid on;
    else
        grid off;
    end
    
    % Colorbar
    pos=get(handles.fig_1,'pos');
    handles.colorbar = colorbar('location','westoutside','position',[0.04 pos(2) 0.02 pos(4)]);
    cbartick =get(handles.colorbar,'Ticks');
    count = 1;
    for i= 1: length(cbartick)
        if floor(cbartick(i)) == cbartick(i)
            cbarticknew(count) = cbartick(i);
             count = count+1;
        end
    end
    set(handles.colorbar,'Ticks',cbarticknew);
    title(handles.colorbar, '#ICs','FontSize',LABEL_FONTSIZE,'FontWeight','bold');
    axcopy;
    % FIGURE 2 ----------------------------------------------------------------
    htmp = findobj(handles.main,'Tag','axes2');
    if isempty(htmp)
        handles.fig_2 = axes('Position',[left bottom-0.25 width .24].*s+q);
    else
        handles.fig_2 = htmp;
        axes(handles.fig_2 );
    end
    bar(sum(mat2plot,1));
    set(handles.fig_2,'tag','axes2');
    hold on;
    ylabel('#ICs/Subject','FontWeight','bold','FontSize',LABEL_FONTSIZE);             % Y Labels
    
    set(handles.fig_2,'XLim',[0.5 size(mat2plot,2) + 0.5]);                           % Setting the X limits

    ymaxval1 = max(sum(mat2plot,1))*(1 + 0.1);
    ymaxval2 = median(sum(mat2plot,1),2)*1.50;
    ymaxval = ymaxval2;
    if ymaxval1 > ymaxval2, ymaxval= ymaxval1; end;
    set(handles.fig_2,'YLim',[0 ymaxval]); % Setting the Y limits
    
    % Setting the X ticks with sessions
    for i=1: length(STUDY.datasetinfo(UniqSubjInd))
        if iscell(STUDY.datasetinfo(UniqSubjInd(1)).session)
            tmpvalnum(i) = cell2mat(STUDY.datasetinfo(UniqSubjInd(i)).session);
            tmpvalP{i}   = num2str(cell2mat(STUDY.datasetinfo(UniqSubjInd(i)).session));
        else
            tmpvalnum(i) = STUDY.datasetinfo(UniqSubjInd(i)).session;
            tmpval{i} = num2str(STUDY.datasetinfo(UniqSubjInd(i)).session);
        end
%         xticklabel_value{i} = [STUDY.datasetinfo(UniqSubjInd(i)).subject '_' tmpval];
    end
    
    if unique(tmpvalnum) ~= 1
        for i =1: length(STUDY.datasetinfo(UniqSubjInd))
            xticklabel_value{i} = [STUDY.datasetinfo(UniqSubjInd(i)).subject '_' tmpval{i}];
        end
        xlabeltext = ('Subject_{(Session)}');  % X Labels
    else
        for i =1: length(STUDY.datasetinfo(UniqSubjInd))
            xticklabel_value{i} = STUDY.datasetinfo(UniqSubjInd(i)).subject;
        end
        xlabeltext = ('Subject');              % X Labels
    end
    
    % Rotating Xticklabels
    Xt = get(handles.fig_2,'XTick');
    Xl = get(handles.fig_2,'XLim');
    ax = axis;                 % Current axis limits
    axis(axis);                % Set the axis limit modes (e.g. XLimMode) to manual
    Yl = ax(3:4);              % Y-axis limits
    
    % Place the text labels
    t = text(Xt,Yl(1)*ones(1,length(Xt)),xticklabel_value);
    set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
    
    set(handles.fig_2,'XTickLabel','');  % Remove the default labels
    
    % Get the Extent of each text object.This loop is unavoidable.
    for i = 1:length(t)
        ext(i,:) = get(t(i),'Extent');
    end
    % Determine the lowest point.  The X-label will be placed so that the top is aligned with this point.
    LowYPoint = min(ext(:,2));
    
    % Place the axis label at this point
    XMidPoint = Xl(1)+abs(diff(Xl))/2;
    text(XMidPoint,LowYPoint,xlabeltext,'VerticalAlignment','top',...
                                        'HorizontalAlignment','center',...
                                        'FontWeight','bold',...
                                        'Tag','Xtext',...
                                        'FontSize',LABEL_FONTSIZE);
    handles.Xtext = findobj(0,'Tag','Xtext');                                   
    lineplots(1) = ceil(prctile(sum(mat2plot,1),25));
    lineplots(2) = ceil(prctile(sum(mat2plot,1),75));
    lineplots(3) = max(sum(mat2plot,1));
    lineplots    = unique(lineplots);
    xval = 0:size(mat2plot,2)+1;
    for i = 1:length(lineplots)
        yval = repmat(lineplots(i),[1 size(mat2plot,2)+2]);
        plot(xval,yval,'LineStyle',':','LineWidth',1,'Color','b');
    end
    
    set(handles.fig_2,'YTick',lineplots,'FontSize',SLABEL_FONTSIZE-1);
    % box off;
    axcopy;
    % FIGURE 3 ----------------------------------------------------------------
    htmp = findobj(handles.main,'Tag','axes3');
    if isempty(htmp)
        handles.fig_3 = axes('Position',[left+width+0.01 bottom .25 height].*s+q);
    else
        handles.fig_3 = htmp;
        axes(handles.fig_3 );
    end
    bar(sum(mat2plot,2));
    set(gca,'yaxislocation','right','tag','axes3');
    hold on;
    
    %     nlinestmp = floor(max(sum(mat2plot,2))/size(mat2plot,2));
    lineplots(1) = ceil(prctile(sum(mat2plot,2),25));
    lineplots(2) = ceil(prctile(sum(mat2plot,2),75));
    lineplots(3) = max(sum(mat2plot,2));
    lineplots    = unique(lineplots);
    xval         = 0.5:size(mat2plot,1)+0.5;
    for i = 1:length(lineplots)
        yval = repmat(lineplots(i),[1 size(mat2plot,1)+1]);
        plot(xval,yval,'LineStyle',':','LineWidth',0.1,'Color','b');
    end
    
    
    set(handles.fig_3,'View',[90 90]);                                              % Rotate Axis
    set(handles.fig_3,'XLim',[0.5 size(mat2plot,1) + 0.5]);                         % Setting the X limits
    set(handles.fig_3,'YLim',[0 max(sum(mat2plot,2))*(1 + 0.1)]);                   % Setting the Y limits
    set(handles.fig_3,'XTickLabel','');                                             % Setting the X ticks
    set(handles.fig_3,'YTick',lineplots,'FontSize',SLABEL_FONTSIZE-1);
    colormap(mycolormap(2,:))                                                   % Color plot
    ylabel('#ICs/Cluster','FontWeight','bold','FontSize',LABEL_FONTSIZE);         % Y Labels
    
    % FIGURE 4 ----------------------------------------------------------------
    htmp = findobj(handles.main,'Tag','axes4');
    if isempty(htmp)
        handles.fig_4 = axes('Position',[left+width+0.01 bottom-0.25 .25 .24].*s+q);
    else
        handles.fig_4 = htmp;
        axes(handles.fig_4 );
    end
    [nbins edges] = histcounts(sum(mat2plot,1));
    histogram(sum(mat2plot,1),edges,'FaceColor',[0.0080    0.0080    1.0000]);
    set(handles.fig_4,'XLim',[0 ymaxval],'tag','axes4');
    
    set(handles.fig_4,'View',[90 90]);
    set(handles.fig_4,'XDir','reverse');
    set(handles.fig_4,'XTickLabel',''); 
    
    % color of hist
    colormap(mycolormap);
    histhandle = get(handles.fig_4,'Children');
    set(histhandle,'FaceColor',mycolormap(1,:));
    get(handles.fig_4,'Children');
    ylabel('ICs/Subject PDF','FontWeight','bold','FontSize',LABEL_FONTSIZE);   
    tmpylim = minmax(nbins);
    set(handles.fig_4,'Ylim',[tmpylim(1) tmpylim(2)+0.1*tmpylim(2)]);
    axcopy;
    % Setting callback ----------------------------------------------------
    % Adding sort functionality
    htmp1 = findobj(handles.main,'tag','checkbox_sort1');
    htmp2 = findobj(handles.main,'tag','checkbox_sort2');
    if any([isempty(htmp1),isempty(htmp2)])
        handles.checkbox_sort1 = uicontrol('parent',handles.main,...
            'style','checkbox',...
            'units','normalized',...
            'position',[left*s(1)+q(1)- 0.06 bottom.*s(2)+q(2)- 0.02 0.05    0.02],...
            'string','Sort',...
            'backgroundcolor',color,...
            'FontWeight','bold',...
            'Tag','checkbox_sort1',...
            'Callback',{@sortgraphs,1,handles.main});
        
        handles.checkbox_sort2 = uicontrol('parent',handles.main,...
            'style','checkbox',...
            'units','normalized',...
            'position',[width*s(3)+q(3)+0.1  height*s(4)+q(4)-0.635 0.06 0.02],...
            'string','Sort',...
            'backgroundcolor',color,...
            'FontWeight','bold',...
            'Tag','checkbox_sort2',...
            'Callback',{@sortgraphs,2,handles.main});
    end 
end
end
%--------------------------------------------------------------------------
function sortgraphs(src,evt,dim2sort,handle)
vals = getappdata(handle);

htmp1 = findobj(handle,'tag','checkbox_sort1');
htmp2 = findobj(handle,'tag','checkbox_sort2');

if dim2sort == 1
    vals.g.sortdim1flag = get(htmp1,'value');
elseif dim2sort == 2
    vals.g.sortdim2flag = get(htmp2,'value');
end
  
std_plotclustinfo(vals.STUDY,vals.SubjClusIC_Matrix,vals.parentcluster,...
    'plot'          ,vals.g.plot,...
    'figlabel'      ,vals.g.figlabel,...
    'colormapstyle' ,vals.g.colormapstyle,...
    'sortdim1flag'  ,vals.g.sortdim1flag,...
    'sortdim2flag'  ,vals.g.sortdim2flag)
end