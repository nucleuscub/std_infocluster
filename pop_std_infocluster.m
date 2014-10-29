% pop_std_infocluster() - sumarize relation of Cluster Vs Subject/Dataset Vs IC
%
% Usage:
%   >>  clustinfo = pop_std_infocluster(STUDY, ALLEEG);
%
% Inputs:
%      STUDY    - studyset structure containing some or all files in ALLEEG
%      ALLEEG   - vector of loaded EEG datasets
%
% Optional inputs:
%
% parentcluster    - Parentcluster to analyze. By default the 1st cluster
% keepsession      - Keep original session info [0,1]. Default 1
% plotinfo         - Plot Cluster Vs Subject/Dataset Vs No of IC [0,1]. Default 1
% csvsave          - Save as CSV [0,1]. Default 0
% verbose          - [0,1] Default 1
% calc             - NOT used right now, but for future purposes
% figlabel         - [0/1] 1 to chow "ICs" in the label's cell, otherwise does
%                    not show the label. Default 1
%
%
% Outputs:
%      STUDY    - studyset structure containing some or all files in ALLEEG
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

function [STUDY, com] = pop_std_infocluster(STUDY, ALLEEG, varargin)

com = ''; % this initialization ensure that the function will return something
% if the user press the cancel button

% addpath([pwd filesep 'dependencies']);  % Comment this 

% display help if not enough arguments

if nargin < 2
    help pop_std_infocluster;
    return;
end

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

pop_flag = 0;
try g.parentcluster;     catch, g.parentcluster   = STUDY.cluster(1).name; pop_flag = 1;         end; % By default the 1st cluster
try g.keepsession;       catch, g.keepsession     = 1;                     pop_flag = 1;         end; % Keep original session info
try g.plotinfo;          catch, g.plotinfo        = 1;                     pop_flag = 1;         end; % Plot by default
try g.csvsave;           catch, g.csvsave         = 0;                     pop_flag = 1;         end; % Save as CSV
try g.savepath;          catch, g.savepath        = STUDY.filepath;        pop_flag = 1;         end; % Path to save
try g.filename;          catch, g.filename        = STUDY.name;            pop_flag = 1;         end; % Path to save
try g.verbose;           catch, g.verbose         = 1;                     pop_flag = 1;         end; % verbose
try g.calc;              catch, g.calc            = 1;                     pop_flag = 1;         end; % Choose what to calc
try g.figlabel;          catch, g.figlabel        = 1;                     pop_flag = 1;         end; % Show labels in fig

% Pop up window
if nargin < 3 || pop_flag ==1
    
    % Checking dependencies for GUI
    % NOTE: Move this to the eepluging_std_infocluster once done ----- CHECK THIS
     PropeditExist =  exist('PropertyEditor', 'file');
     JidePropExist =  exist('JidePropertyGridField', 'file');
     
     if PropeditExist ~= 2 || JidePropExist ~= 2
         WhereIsInfocluster = which('std_infocluster');
         [path,name, ext]   = fileparts(WhereIsInfocluster);
         addpath([path filesep 'dependencies/PropertyGrid-2010-09-16-mod']);
     end
     
    
    icadefs;
    color = BACKEEGLABCOLOR;
    
    setappdata(0,'STUDY',STUDY);   %handles.STUDY  = STUDY;
    setappdata(0,'opts',g);        %handles.opts   = g;
    setappdata(0,'ALLEEG',ALLEEG); %handles.ALLEEG = ALLEEG;
    
    handles.mainfig = figure('MenuBar','none',...
                             'Name','pop_std_infocluster',...
                             'NumberTitle','off',...
                             'Units', 'Normalized',...
                             'Color', color,...
                             'Position',[0.43,0.54,0.24,0.27]);
    
    % Panel 1
    % .........................................................................
    handles.Panel1 = uipanel('Parent', handles.mainfig,'Units', 'Normalized','BackgroundColor',color,...
        'Position',[0.032,0.575,0.936,0.376]);
    
    % Texts
    handles.Text2 = uicontrol('Parent', handles.Panel1,'Style','Text');
    set(handles.Text2,'String','Parent Cluster','Units', 'Normalized','BackgroundColor',color,'Position',[.035 .626 .239 .143]);
    
    handles.Text1 = uicontrol('Parent', handles.Panel1,'Style','Text');
    set(handles.Text1,'String','Properties','Units', 'Normalized','BackgroundColor',color,'Position',[.035 .231 .175 .165]);
    
    % Popupmenus
    poplist1 = unique([STUDY.cluster.parent]);
    handles.popup_parent = uicontrol('Parent', handles.Panel1,'Style','Popupmenu');
    set(handles.popup_parent,'String',poplist1,'Units', 'Normalized','Position',[.345 .582 .624 .204]); % insert list of clusters
    
    poplist2 = {'No. ICs';'No. ERP'};
    handles.popup_prop = uicontrol('Parent', handles.Panel1,'Style','Popupmenu');
    set(handles.popup_prop,'String',poplist2,'Units', 'Normalized','Position',[.345 .206 .624 .204]);   % insert new metrics
    
    
    % Panel 2
    % .........................................................................
    handles.Panel2 = uipanel('Parent', handles.mainfig,'Units', 'Normalized','BackgroundColor',color,...
        'Position',[0.032,0.286,0.936,0.270]);
    
    % Edit
    handles.edit_save = uicontrol('Parent', handles.Panel2,'Style','Edit');
    set(handles.edit_save,'Units', 'Normalized','BackgroundColor',[1 1 1],'enable', 'off','Position',[.05 .258 .765 .303]);
    
    % Button
    handles.button_path = uicontrol('Parent', handles.Panel2,'Style','PushButton');
    set(handles.button_path,'String','...','Units', 'Normalized','enable', 'off','Position',[.85 .242 .091 .318],'CallBack', {@callback_button_path,handles});
    
    %Checkox
    handles.checkbox_save = uicontrol('Parent', handles.Panel2,'Style','checkbox');
    set(handles.checkbox_save,'String','Save Analysis','Units', 'Normalized','BackgroundColor',color,'Position',[0.05 .621 .345 .273],...
        'CallBack', {@callback_chckbutton_save,handles});
    
    % Buttoms
    % .........................................................................
    button_help = uicontrol('Parent', handles.mainfig,'Style','PushButton');
    set(button_help,'String','Help','Units', 'Normalized','Position',[0.035 0.044 0.2 0.083],'CallBack', @callback_button_help);
    
    handles.button_cancel = uicontrol('Parent', handles.mainfig,'Style','PushButton');
    set(handles.button_cancel,'String','Cancel','Units', 'Normalized','Position',[0.49 0.044 0.2 0.083],'CallBack', {@callback_buton_cancel,handles.mainfig});
    
    handles.button_ok = uicontrol('Parent', handles.mainfig,'Style','PushButton');
    set(handles.button_ok,'String','OK','Units', 'Normalized','Position',[0.765 0.044 0.2 0.083],'CallBack', {@callback_button_ok,handles});
    
    handles.button_opt = uicontrol('Parent', handles.mainfig,'Style','PushButton');
    set(handles.button_opt,'String','Opt','Units', 'Normalized','Position',[0.765 0.151 0.2 0.083],'CallBack',@callback_button_opts);
    
    % .........................................................................
    uiwait(handles.mainfig);
    
end
% OUTPUTS
STUDY          = getappdata(0,'STUDY');
clust_statout  = getappdata(0,'clust_statout');

args = getappdata(0,'args');
% com = ['[STUDY , clust_statout] = std_infocluster(STUDY,ALLEEG,' args{:}  ');']
com            = 'Need to be updated';

%**************************************************************************
%**************************************************************************
% AUXILIAR FUNCTIONS

% _________________________________________________________________________
function callback_buton_cancel(src,eventdata,h)
close(h);

% _________________________________________________________________________
function callback_button_path(src,eventdata,handles)

[FolderName] = uigetdir('','Select path to save results');

opts = getappdata(0,'opts');
opts.savepath = FolderName; % setappdata(0,'handles.opts.savepath',FolderName); % Sending output to main
setappdata(0,'opts',opts);

% _________________________________________________________________________
function callback_chckbutton_save(src,eventdata,handles)

chck_state = get(src,'Value');
opts = getappdata(0,'opts');   % retrievineg opts

if chck_state == 1
    set(handles.edit_save,   'enable', 'on');
    set(handles.button_path, 'enable', 'on');
    
    DefaultFilename = [opts.filename '_IC_Subj'];
    
    if size(DefaultFilename) <= 50
        set(handles.edit_save, 'String',DefaultFilename);
    else
        set(handles.edit_save, 'String',DefaultFilename(end-51:end));
    end
    
    opts.csvsave = 1;
    
else
    set(handles.edit_save,  'enable', 'off');
    set(handles.button_path,'enable', 'off');
    set(handles.edit_save,  'String', '');
end

setappdata(0,'opts',opts); % Storing opts

% _________________________________________________________________________
function handles = callback_button_ok(src,eventdata,handles)

options = '';

%Check  and retrieving GUI inputs
ALLEEG = getappdata(0,'ALLEEG');                                      % retrievineg ALLEEG
STUDY  = getappdata(0,'STUDY');                                       % retrievineg STUDY
opts   = getappdata(0,'opts');                                        % retrievineg opts

opts.parentcluster   = char(get(handles.popup_parent, 'String'));
opts.calc            = get(handles.popup_prop, 'Value');
opts.csvsave         = get(handles.checkbox_save,'Value');

if ~isempty(get(handles.edit_save,'String')), opts.filename = get(handles.edit_save,'String'); end;% csvsave


optionsname = fieldnames(opts);
args = '';
c = 1;

for i = 1:length(optionsname)
    args{c}   = char(optionsname(i));
    %     args{c+1} = char(eval(['handles.opts.' char(optionsname(i))]));
    args{c+1} = eval(['opts.' char(optionsname(i))]);
    c = c + 2;
end

% Evaluating function
[STUDY , clust_statout] = std_infocluster(STUDY,ALLEEG,args{:});

% Saving data

guidata(src, handles);          % Sending output to main
setappdata(0,'ALLEEG',ALLEEG);               % Storing ALLEEG
setappdata(0,'STUDY',STUDY);                 % Storing STUDY
setappdata(0,'opts',opts);                   % Storing opts
setappdata(0,'clust_statout',clust_statout); % Storing output
setappdata(0,'args',args); % Storing output

close(handles.mainfig); % Closing figure

% _________________________________________________________________________
function callback_button_help(src,eventdata)

pophelp('pop_std_infocluster')

function callback_button_opts(src,eventdata)

opts   = getappdata(0,'opts');                                        % retrievineg opts

% create figure
f = figure( ...
    'MenuBar', 'none', ...
    'Name', 'Options: std_infocluster', ...
    'NumberTitle', 'off', ...
    'Toolbar', 'none', ...
    'Units', 'normalized', ...
    'Position',[0.43,0.54,0.20,0.20]);

items  = {opts} ;
editor = PropertyEditor(f, 'Items', items);