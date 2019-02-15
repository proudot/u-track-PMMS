function varargout = motionAnalysisVisualGUI(varargin)
% MOTIONANALYSISVISUALGUI M-file for motionAnalysisVisualGUI.fig
%      MOTIONANALYSISVISUALGUI, by itself, creates a new MOTIONANALYSISVISUALGUI or raises the existing
%      singleton*.
%
%      H = MOTIONANALYSISVISUALGUI returns the handle to a new MOTIONANALYSISVISUALGUI or the handle to
%      the existing singleton*.
%
%      MOTIONANALYSISVISUALGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOTIONANALYSISVISUALGUI.M with the given input arguments.
%
%      MOTIONANALYSISVISUALGUI('Property','Value',...) creates a new MOTIONANALYSISVISUALGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before motionAnalysisVisualGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to motionAnalysisVisualGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help motionAnalysisVisualGUI

% Last Modified by GUIDE v2.5 02-Apr-2012 11:44:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @motionAnalysisVisualGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @motionAnalysisVisualGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before motionAnalysisVisualGUI is made visible.
function motionAnalysisVisualGUI_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.text_copyright, 'String', getLCCBCopyright());

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

% Choose default command line output for detectionVisualGUI
handles.output = hObject;

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);
userData.userData_main = get(userData.handles_main.figure1, 'UserData');

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.MD = userData_main.MD(userData_main.id);  % Get the current Movie Data
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get icon infomation
userData.questIconData = userData.userData_main.questIconData;
userData.colormap = userData.userData_main.colormap;


% Get first channel index
chan = find(userData.crtProc.checkChannelOutput,1);
assert(~isempty(chan), 'User-defined: the process does not have output yet.')

% Make sure detection output is valid
tracksFinal = userData.crtPackage.processes_{2}.loadChannelOutput(chan);
if isempty(tracksFinal)
   error('User-defined: there is no detection information in the output variable.') 
end

allEvents = vertcat(tracksFinal.seqOfEvents);
userData.firstframe = min(allEvents(:,1));
userData.lastframe = max(allEvents(:,1));


% Set frames
set(handles.text_framenum, 'String', ['( Tracks availabe from frame ',num2str(userData.firstframe),' to ',num2str(userData.lastframe),' )'])

set(handles.edit_min, 'String', userData.firstframe)
set(handles.edit_max, 'String', userData.lastframe)
set(handles.checkbox_showConf, 'Value', 0)
set(handles.checkbox_simplifyLin1, 'Value', 0)
set(handles.checkbox_simplifyLin2, 'Value', 0,'Enable','off')

set([handles.edit_offset1 handles.edit_offset2], 'String', '0')

% ----------------------Set up help icon------------------------    
% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = motionAnalysisVisualGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1)


% --- Executes on button press in pushbutton_display.
function pushbutton_display_Callback(hObject, eventdata, handles)
% Tool 1: plotTrakcs2D

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;


% Frame range
min = str2double(get(handles.edit_min, 'String'));
max = str2double(get(handles.edit_max, 'String'));
if isnan(min) || isnan(max) || min<userData.firstframe || ...
        max >userData.lastframe || min>max
    errordlg('Please provide a valid value to parameter "Frames to Include in Plot:".','Error','modal')
    return
end    
    
% Offset
dx = str2double(get(handles.edit_offset1, 'String'));
dy = str2double(get(handles.edit_offset2, 'String'));
if isnan(dx) || isnan(dy) || dx <0 || dy <0
    errordlg('Please provide a valid value to parameter "Offset for Plotting:".','Error','modal')
    return
end

showConf = get(handles.checkbox_showConf,'Value');
simplifyLin = get(handles.checkbox_simplifyLin1,'Value')+...
    get(handles.checkbox_simplifyLin2,'Value');

% Load image
imagePath = get(handles.edit_image,'String');
if ~isempty(imagePath)
    I=imread(imagePath);
else 
    I=[];
end

% Load output
chan = find(userData.crtProc.checkChannelOutput,1);
tracksFinal = userData.crtPackage.processes_{2}.loadChannelOutput(chan);
diffAnalysisRes = userData.crtPackage.processes_{3}.loadChannelOutput(chan);

% Call function
plotTracksDiffAnalysis2D(tracksFinal,diffAnalysisRes,[min max],1,I,...
    showConf,simplifyLin,[dx dy]);

% --- Executes on button press in pushbutton_image.
function pushbutton_image_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');
if(isempty(userData)), userData = struct(); end;

[file,path] = uigetfile('*.*','Select an Image',...
             get(handles.edit_image, 'String'));
        
if ~any([file,path]), return; end

set(handles.edit_image, 'String', [path file])


% --- Executes on button press in checkbox_simplifyLin1.
function checkbox_simplifyLin1_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.checkbox_simplifyLin2,'Enable','on');
else
    set(handles.checkbox_simplifyLin2,'Value',0,'Enable','off');
end
