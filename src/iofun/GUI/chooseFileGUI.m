function varargout = chooseFileGUI(varargin)
%CHOOSEFILEGUI allows the user to select from a list among the choices specified in the input cell array
% 
%SYNOPSIS  selectedItem = chooseFileGUI(specificationCell)
%
%INPUT     1d-cell array with all the choices
%
%OUTPUT    index of selected list entry or empty matrix if cancelled
%
%REMARKS   can be used for more than choosing files!
%
%c: 02/03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last Modified by GUIDE v2.5 06-Feb-2003 11:02:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chooseFileGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @chooseFileGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before chooseFileGUI is made visible.
function chooseFileGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chooseFileGUI (see VARARGIN)

% Choose default command line output for chooseFileGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%write the inputstring into the listbox
try
    listCell = varargin{1};
catch
    if findstr(lasterr,'Index exceeds matrix dimensions.')
        error('input to chooseFileGUI must not be empty!')
    else
        rethrow(lasterr);
    end
end

%test input
if isempty(listCell)|~iscell(listCell)
    error('the input has to be a non-empty cell array!')
end

set(handles.cfg_fileListBox,'String',listCell);

%---set the right size of listbox and figure

%get handles
listBoxH = handles.cfg_fileListBox;
cfgH = handles.chooseFileGUI;

%remember old sizes
listBoxPos = get(listBoxH,'Position');
cfgPos = get(cfgH,'Position');

%find longest string
for i = 1:length(listCell)
    strLength(i) = length(listCell{i});
end
[maxLength] = max(strLength);

%make Listbox width 1.3x the lenght of the longest string (but no lower than oldWidth)
if 1.3*maxLength > listBoxPos(3)
    listBoxPos(3) = 1.3*maxLength;
    cfgPos(3) = listBoxPos(3)+4;
    set(cfgH,'Position',cfgPos);
    set(listBoxH,'Position',listBoxPos);    
end

%----end set right size


% UIWAIT makes chooseFileGUI wait for user response (see UIRESUME)
uiwait(handles.chooseFileGUI);

% --- Outputs from this function are returned to the command line.
function varargout = chooseFileGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles) %if closed other than by button
    handles.chooseFileGUI=[];
    handles.output=[];
end
% Get default command line output from handles structure
varargout{1} = handles.output;
%delete listData.mat;
close(handles.chooseFileGUI);

% --- Executes during object creation, after setting all properties.
function cfg_fileListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cfg_fileListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% load listData;
% set(hObject,'String',listData);

% --- Executes on selection change in cfg_fileListBox.
function cfg_fileListBox_Callback(hObject, eventdata, handles)
% hObject    handle to cfg_fileListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cfg_fileListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cfg_fileListBox


% --- Executes on button press in cfg_OK_PB.
function cfg_OK_PB_Callback(hObject, eventdata, handles)
% hObject    handle to cfg_OK_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output=get(handles.cfg_fileListBox,'Value');
% Update handles structure
guidata(hObject, handles);
% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.chooseFileGUI);

% --- Executes on button press in cfg_cancel_PB.
function cfg_cancel_PB_Callback(hObject, eventdata, handles)
% hObject    handle to cfg_cancel_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output=[];
% Update handles structure
guidata(hObject, handles);
% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.chooseFileGUI);
