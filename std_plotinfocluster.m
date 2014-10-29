% std_plotinfocluster() - Plot relation of Clusters  Vs Subject/Dataset Vs IC
%
% Usage:
%   >>  [mat2plot, UniqSubjInd] = std_plotinfocluster(STUDY,SubjClusIC_Matrix,parentcluster)
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
%   std_plotinfocluster
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
function [mat2plot, UniqSubjInd] = std_plotinfocluster(STUDY,SubjClusIC_Matrix,parentcluster,varargin)

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

% Checking dependencies for std_plotinfocluster
B2rExist =  exist('b2r', 'file');

if B2rExist ~= 2
    WhereIsInfocluster = which('std_infocluster');
    [path,name, ext]   = fileparts(WhereIsInfocluster);
    addpath([path filesep 'dependencies/b2r']);
end


% Settings
% icadefs;
colormap('jet');                                                           % Def colormap
% color = BACKEEGLABCOLOR;

% Getting clusters names
hits_temp = cellfun(@(x)strcmp(x,parentcluster),{STUDY.cluster.name});
parent_indx = find(hits_temp);
cls = (parent_indx+1):(parent_indx + length(STUDY.cluster(parent_indx).child));
cls_names = {STUDY.cluster(cls).name}';

% Check session in data or take from input (now we assume data is comming with session info)

sessions  = cellfun(@(x) x,{STUDY.datasetinfo.session});

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

% Getting names



%**************************************************************************
%**************************************************************************
%**************************************************************************
if g.plot
    % Making FIG
    h = figure('Name', ['Cluster Info: ' parentcluster],...
               'numbertitle', 'off');                                      
    
    set(h, 'ToolBar', 'none');
    left = 0.05; bottom = 0.3; width = 0.8; height = 0.7;
    pos  = get(gca,'position');
    q    = [pos(1) pos(2) 0 0];
    s    = [pos(3) pos(4) pos(3) pos(4)];
    axis off;
    
    % FIGURE 1 ----------------------------------------------------------------
    h(1) = axes('Position',[left bottom width height].*s+q);
    
    imagesc(mat2plot); % plot fig1
    
    
    % Creating labels for grid
    for i= 1:size(mat2plot,1)
        for j = 1:size(mat2plot,2)
            if (mat2plot(i,j) ~= 0)
                if g.figlabel
                    gridlabel{i,j} = [num2str(mat2plot(i,j),'%d') ' ICs'];
                else
                    gridlabel{i,j} = [num2str(mat2plot(i,j),'%d')];
                end
                
            else
                gridlabel{i,j} = num2str(' ','%d ');
            end
        end
    end
    
    gridlabel = strtrim(cellstr(gridlabel));                               % Remove any space padding
    [x,y]     = meshgrid(1:size(mat2plot,2),1:size(mat2plot,1));           % Create x and y coordinates for the strings
    
    text(x(:),y(:),gridlabel(:),'HorizontalAlignment','center');           % Plot the strings (,'FontWeight','bold')
    ylabel('Clusters','FontWeight','bold');                                % Y Labels
    
    set(gca,'XLim',[0.5 size(mat2plot,2) + 0.5]);                          % Setting the X limits
    set(gca,'YLim',[0.5 size(mat2plot,1) + 0.5]);                          % Setting the Y limits
    
    set(gca, 'Xtick',1:size(mat2plot,2));                                  % Setting the X ticks
    set(gca, 'Ytick',1:size(mat2plot,1)) ;                                 % Setting the Y ticks
    
    set(gca,'YTickLabel',cls_names);                                       % Note: need to pick the real plotted subj
                                
    tmp_caxis = caxis;                                                     % Color plot
    caxis([-max(tmp_caxis) max(tmp_caxis)]);                               % ...
    colormap(b2r(-max(tmp_caxis), max(tmp_caxis)));                        % ...                     
    
    grid on;
    % box off;
    
    
    % FIGURE 2 ----------------------------------------------------------------
    h(2) = axes('Position',[left bottom-0.25 width .24].*s+q);
    
    bar(sum(mat2plot,1));
    
    ylabel('No ICs','FontWeight','bold');                                  % Y Labels
    xlabel('Subjects/Datasets','FontWeight','bold');                       % X Labels
    set(gca,'XLim',[0.5 size(mat2plot,2) + 0.5]);                          % Setting the X limits
    set(gca,'YLim',[0 max(sum(mat2plot,1))*(1 + 0.1)]);                    % Setting the Y limits
    
    colormap(b2r(-max(tmp_caxis), max(tmp_caxis)));                        % Color plot
    
    set(gca,'XTickLabel',{STUDY.datasetinfo(UniqSubjInd).subject});        % Setting the X ticks
    % box off;
    
    % FIGURE 3 ----------------------------------------------------------------
    h(3) = axes('Position',[left+width+0.01 bottom .1 height].*s+q);
    
    bar(sum(mat2plot,2));
    
    set(h(3),'View',[90 90]);                                              % Rotate Axis
    set(gca,'XLim',[0.5 size(mat2plot,1) + 0.5]);                          % Setting the X limits
    set(gca,'YLim',[0 max(sum(mat2plot,2))*(1 + 0.1)]);                    % Setting the Y limits
    set(gca,'XTickLabel','');                                              % Setting the X ticks
    
    colormap(b2r(-max(tmp_caxis), max(tmp_caxis)));                        % Color plot
    
    ylabel('No ICs','FontWeight','bold');                                  % Y Labels
    
    % box off;
end
