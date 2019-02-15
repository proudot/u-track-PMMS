function varargout = CropPanel(varargin)
% CropPanel M-file for CropPanel.fig
%      CropPanel, by itself, creates a new CropPanel or raises the existing
%      singleton*.
%
%      H = CropPanel returns the handle to a new CropPanel or the handle to
%      the existing singleton*.
%
%      CropPanel('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CropPanel.M with the given input arguments.
%
%      CropPanel('Property','Value',...) creates a new CropPanel or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before The_Crapper_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CropPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CropPanel

% Last Modified by GUIDE v2.5 10-Apr-2009 16:32:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CropPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @CropPanel_OutputFcn, ...
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

% --- Executes just before CropPanel is made visible.
function CropPanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CropPanel (see VARARGIN)

val = get(handles.menu_color,'Value');
switch val
case 1
    handles.colormap_type = 'gray';
case 2
    handles.colormap_type = 'Jet';
case 3
    handles.colormap_type = 'HSV';
case 4
    handles.colormap_type = 'HOT';    
end

% Choose default command line output for CropPanel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projSetupGUI wait for user response (see UIRESUME)
%uiwait(handles.output);

initialize_gui(hObject, handles, false);

% UIWAIT makes CropPanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CropPanel_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function density_Callback(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.metricdata.density = density;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function volume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function volume_Callback(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volume as text
%        str2double(get(hObject,'String')) returns contents of volume as a double
volume = str2double(get(hObject, 'String'));
if isnan(volume)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new volume value
handles.metricdata.volume = volume;
guidata(hObject,handles)

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mass = handles.metricdata.density * handles.metricdata.volume;
set(handles.mass, 'String', mass);

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --------------------------------------------------------------------
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to unitgroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, croppanel means
% we are we are just re-initializing a GUI by calling croppanel from the cmd line
% while croppanel is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

handles.metricdata.density = 0;
handles.metricdata.volume  = 0;

%set(handles.density, 'String', handles.metricdata.density);
% set(handles.volume,  'String', handles.metricdata.volume);
% set(handles.mass, 'String', 0);
% 
% set(handles.unitgroup, 'SelectedObject', handles.english);
% 
% set(handles.text4, 'String', 'lb/cu.in');
% set(handles.text5, 'String', 'cu.in');
% set(handles.text6, 'String', 'lb');

% Update handles structure
guidata(handles.figure1, handles);





% --- Executes on button press in pushLoadStacks.
function pushLoadStacks_Callback(hObject, eventdata, handles)
% hObject    handle to pushLoadStacks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%allow user to choose directory

try
    basePath = uigetdir('Please select directory of movies to be collected');
    if basePath == 0
        disp('Interrupted by user.');
        return;
    end
end

%find all makiAnalysis files in chosen directory
fileList = searchFiles('tif',[],basePath,1);
channelDir = unique(fileList(:,2));
% if strcmp(channelDir{1},pwd) == 1
%    channelDir(1) = [];
% end
    
%allow user to choose files
selectIdx = listSelectGUI(channelDir,[],'move');

%get number of chosen files
channelDirs = channelDir(selectIdx,:);
numChannels = length(selectIdx);

if numChannels == 0 %no directory has been specified
    f = warndlg('Please specify at least one channel', 'Warning', 'modal');
    return;
end

clear fileList
for n_idx = 1:numChannels
    fileList{n_idx} = searchFiles('tif',[],channelDirs{n_idx},1);
    numImages(n_idx) = length(fileList{n_idx});
end

if numChannels > 1
    if  length(unique(numImages)) > 1
        f = warndlg(['The selected directores have different number of images.  Please check again before proceed.'], 'Warning', 'modal');
        return;
    end
    for c_idx = 1:numChannels
        [this_row this_col] = size(imread([fileList{c_idx}{1,2} filesep fileList{c_idx}{1,1}]));
        r(c_idx) = this_row;
        c(c_idx) = this_col;
    end
    if sum(abs(diff(r))) > 0 | sum(abs(diff(c))) > 0
        f = warndlg(['The selected directores have images with different dimensions.  Please check again before proceed.'], 'Warning', 'modal');
        return;        
    end
    
    slider_channel_step(1) = 1/(numChannels - 1);
    slider_channel_step(2) = 1/(numChannels - 1);
    set(handles.slider_channel, 'Visible', 'on','Max',numChannels,...
        'Value',1,'Min',1,'sliderstep',slider_channel_step);    
else
    set(handles.slider_channel, 'Visible', 'off');
end
    
mask_idx = 1;
slider_time_step(1) = 1/numImages(1);
slider_time_step(2) = 1/numImages(1);

channel_idx = 1;
thismask = imread([fileList{channel_idx}{mask_idx,2} filesep fileList{channel_idx}{mask_idx,1}]);
set(handles.slider_time, 'Visible', 'on','Max',numImages(1),...
    'Value',1,'Min',1,'sliderstep',slider_time_step);
 
axes(handles.axes_image);    imshow(thismask, []);   colormap(handles.colormap_type);
handles.thismask = thismask;
handles.mask_idx = mask_idx;

handles.rowSize = size(thismask,1);
handles.colSize = size(thismask,2);
handles.numChannels = numChannels;
handles.numImages = numImages;
handles.channel_idx = channel_idx;
handles.fileList = fileList;
handles.channelDirs = channelDirs;

%text(handles.colSize+5,20,['Image Index:' num2str(mask_idx)],'BackgroundColor',[.7 .9 .7]);
set(handles.text_raw, 'String', ['Raw Data: Image ' num2str(mask_idx) ' in Channel ' num2str(channel_idx)]);    

set(handles.text_raw , 'Visible', 'on');
set(handles.push_selectroi , 'Visible', 'on');
set(handles.push_selectroifreehand , 'Visible', 'on');

handles.iscrop = 0;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in backbutton.
function backbutton_Callback(hObject, eventdata, handles)
% hObject    handle to backbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%delete(hObject);
delete(handles.figure1);
%uiresume(handles.output);




function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_selectroi.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to push_selectroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_play.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to push_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushExtractMeta.
function pushExtractMeta_Callback(hObject, eventdata, handles)
% hObject    handle to pushExtractMeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stk2tif('default');

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_time_crop_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.crop_mask_idx = get(handles.slider_time_crop,'Value');
handles.crop_mask_idx = round(handles.crop_mask_idx);
set(handles.slider_time_crop,'Value',handles.crop_mask_idx);

handles.crop_channel_idx = get(handles.slider_channel_crop,'Value');
if handles.crop_channel_idx == 0
    handles.crop_channel_idx = 1;
end
thiscropmask = imread([handles.cropfileList{handles.crop_channel_idx}{handles.crop_mask_idx}]);
axes(handles.axes_cropimage);
imshow(thiscropmask, []);   colormap(handles.colormap_type);
set(handles.text_cropped, 'String', ['Raw Data: Image ' num2str(handles.crop_mask_idx) ' in Channel ' num2str(handles.crop_channel_idx)]);    



handles.mask_idx = handles.crop_mask_idx;
thismask = imread([handles.fileList{handles.channel_idx}{handles.mask_idx,2} filesep handles.fileList{handles.channel_idx}{handles.mask_idx,1}]);
axes(handles.axes_image);    imshow(thismask, []);   colormap(handles.colormap_type);
if handles.iscrop
    hold on;
    plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
    hold off;
end

handles.thismask = thismask;
set(handles.text_raw, 'String', ['Raw Data: Image ' num2str(handles.mask_idx) ' in Channel ' num2str(handles.channel_idx)]);   
set(handles.slider_time,'Value',handles.mask_idx);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_time_crop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_time_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_crop2.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in push_crop3.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in push_crop4.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_channel_crop_Callback(hObject, eventdata, handles)
% hObject    handle to slider_channel_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.crop_channel_idx = get(handles.slider_channel_crop,'Value');
handles.crop_channel_idx = round(handles.crop_channel_idx);
set(handles.slider_channel_crop, 'String', num2str(handles.crop_channel_idx));

thiscropmask = imread([handles.cropfileList{handles.crop_channel_idx}{handles.crop_mask_idx}]);
axes(handles.axes_cropimage);
imshow(thiscropmask, []);   colormap(handles.colormap_type);
set(handles.text_cropped, 'String', ['Cropped Data: Image ' num2str(handles.crop_mask_idx) ' in Channel ' num2str(handles.crop_channel_idx)]);    

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_channel_crop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_channel_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
% Hint: place code in OpeningFcn to populate axes4




% --- Executes during object creation, after setting all properties.
function axes_cropimage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_cropimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
% Hint: place code in OpeningFcn to populate axes_cropimage




% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'Visible','off');
% Hint: place code in OpeningFcn to populate axes5


% --- Executes on selection change in menu_color.
function menu_color_Callback(hObject, eventdata, handles)
% hObject    handle to menu_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menu_color contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_color

val = get(handles.menu_color,'Value');
switch val
case 1
    handles.colormap_type = 'gray';
case 2
    handles.colormap_type = 'Jet';
case 3
    handles.colormap_type = 'HSV';
case 4
    handles.colormap_type = 'HOT';    
end

guidata(hObject, handles);



% --- Executes on slider movement.
function slider_time_Callback(hObject, eventdata, handles)
% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% hObject    handle to slider_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.mask_idx = get(handles.slider_time,'Value');
handles.mask_idx = round(handles.mask_idx);
set(handles.slider_time,'Value',handles.mask_idx);
thismask = imread([handles.fileList{handles.channel_idx}{handles.mask_idx,2} filesep handles.fileList{handles.channel_idx}{handles.mask_idx,1}]);
axes(handles.axes_image);
imshow(thismask, []);   colormap(handles.colormap_type);
set(handles.text_raw, 'String', ['Raw Data: Image ' num2str(handles.mask_idx) ' in Channel ' num2str(handles.channel_idx)]);    

if handles.iscrop
    hold on;
    plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
    hold off;
end

handles.thismask = thismask;

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider_channel_Callback(hObject, eventdata, handles)
% hObject    handle to slider_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.channel_idx = get(handles.slider_channel,'Value');
handles.channel_idx = round(handles.channel_idx);
set(handles.slider_channel, 'String', num2str(handles.channel_idx));

thismask = imread([handles.fileList{handles.channel_idx}{handles.mask_idx,2} filesep handles.fileList{handles.channel_idx}{handles.mask_idx,1}]);
axes(handles.axes_image);
imshow(thismask, []);   colormap(handles.colormap_type);
set(handles.text_raw, 'String', ['Raw Data: Image ' num2str(handles.mask_idx) ' in Channel ' num2str(handles.channel_idx)]);    

if handles.iscrop
    hold on;
    plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
    hold off;
end

handles.thismask = thismask;

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in push_selectroi.
function push_selectroi_Callback(hObject, eventdata, handles)
% hObject    handle to push_selectroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_image);
img = handles.thismask;
imshow(img, []);    colormap(handles.colormap_type);
%text(handles.colSize+5,20,['Image Index:' num2str(handles.mask_idx)],'BackgroundColor',[.7 .9 .7]);

% % Update handles structure
% guidata(hObject, handles);
% img_show(:,:,1) = zeros(size(img));
% img_show(:,:,2) = mat2gray(img);
% img_show(:,:,3) = mat2gray(img);

[segimg, rect] = imcrop;
if isempty(rect)
   handles.iscrop = 0;
   guidata(hObject, handles);
   return;
end

x0 = round(rect(1));
x1 = round(rect(1) + rect(3));
y0 = round(rect(2));
y1 = round(rect(2) + rect(4));

% Check boundaries
if x0<=0, x0=1; end
if y0<=0, y0=1; end
if y1>handles.rowSize, y1=handles.rowSize; end
if x1>handles.colSize, x1=handles.colSize; end

valid_position_x = [x0, x0, x1, x1, x0];
valid_position_y = [y0, y1, y1, y0, y0];
segimg = roipoly(img,valid_position_x,valid_position_y);
c = contourc(double(segimg), [0 0]);
handles.pixel_list = [c(1,2:end)' c(2,2:end)'];

hold on;
plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
hold off;
handles.iscrop = 1;


% Select output directory
thispath=uigetdir('','Select an output directory');
if thispath==0
    guidata(hObject, handles);    
    thispath = warndlg('Please specify an output directory then press any crop button again', 'Warning', 'modal');
    return;
end

for c_idx=1:handles.numChannels
    [fpath,dirName,fno,fext]=getFilenameBody(handles.channelDirs{c_idx});
    if ~exist([thispath filesep dirName num2str(fno)])    
        mkdir(thispath,[dirName num2str(fno)]);
    end
    cropDir{c_idx} = [dirName num2str(fno)];
end

% Initializing waitbar
h=waitbar(0,'Processing...');

% Processing files
for c_idx=1:handles.numChannels
    for n_idx=1:handles.numImages(1)

        thismask = imread([handles.fileList{c_idx}{n_idx,2} filesep handles.fileList{c_idx}{n_idx,1}]);

        % Crop image
        imgC=thismask(y0:y1,x0:x1);

        % Prepare filename with path
        filename=[thispath,filesep,cropDir{c_idx},filesep,'crop_',handles.fileList{c_idx}{n_idx,1}];
        handles.cropfileList{c_idx}{n_idx} = filename;
        
        % Write file to disk
        imwrite(imgC,filename, 'Compression', 'none');

        % Update waitbar
        waitbar(((c_idx-1)*handles.numImages(1) + n_idx)/(handles.numChannels*handles.numImages(1)),h);
    end
end

% Close waitbar
close(h);

% Display the crop images
if handles.numChannels > 1
    slider_channel_step(1) = 1/(handles.numChannels - 1);
    slider_channel_step(2) = 1/(handles.numChannels - 1);
    set(handles.slider_channel_crop, 'Visible', 'on','Max',handles.numChannels,...
        'Value',1,'Min',1,'sliderstep',slider_channel_step);    
else
    set(handles.slider_channel_crop, 'Visible', 'off');
end

crop_mask_idx = 1;
slider_time_step(1) = 1/handles.numImages(1);
slider_time_step(2) = 1/handles.numImages(1);

crop_channel_idx = 1;
thiscropmask = imread([handles.cropfileList{crop_channel_idx}{crop_mask_idx}]);
set(handles.slider_time_crop, 'Visible', 'on','Max',handles.numImages(1),...
    'Value',1,'Min',1,'sliderstep',slider_time_step);
 
axes(handles.axes_cropimage);    imshow(thiscropmask, []);   colormap(handles.colormap_type);
handles.thiscropmask = thiscropmask;
handles.crop_mask_idx = crop_mask_idx;
handles.crop_channel_idx = crop_channel_idx;

set(handles.text_cropped, 'String', ['Cropped Data: Image ' num2str(crop_mask_idx) ' in Channel ' num2str(crop_channel_idx)]);    
set(handles.text_cropped, 'Visible', 'on');

mask_idx = 1;
channel_idx = handles.channel_idx;
thismask = imread([handles.fileList{channel_idx}{mask_idx,2} filesep handles.fileList{channel_idx}{mask_idx,1}]);
axes(handles.axes_image);    imshow(thismask, []);   colormap(handles.colormap_type);
if handles.iscrop
    hold on;
    plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
    hold off;
end

handles.thismask = thismask;
handles.mask_idx = mask_idx;
set(handles.text_raw, 'String', ['Raw Data: Image ' num2str(mask_idx) ' in Channel ' num2str(channel_idx)]);   
set(handles.slider_time,'Value',handles.mask_idx);

set(handles.push_play, 'Visible', 'on');
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in push_play.
function push_play_Callback(hObject, eventdata, handles)
% hObject    handle to push_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

channel_idx = handles.channel_idx;
crop_channel_idx = handles.crop_channel_idx;

for mask_idx = 1:handles.numImages

    thismask = imread([handles.fileList{channel_idx}{mask_idx,2} filesep handles.fileList{channel_idx}{mask_idx,1}]);
    axes(handles.axes_image);    imshow(thismask, []);   colormap(handles.colormap_type);

    if handles.iscrop
        hold on;
        plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
        hold off;
    end

    thiscropmask = imread([handles.cropfileList{crop_channel_idx}{mask_idx}]);
    axes(handles.axes_cropimage);    imshow(thiscropmask, []);   colormap(handles.colormap_type);

    set(handles.slider_time,'Value',mask_idx);
    set(handles.slider_time_crop,'Value',mask_idx);

    set(handles.text_cropped, 'String', ['Cropped Data: Image ' num2str(mask_idx) ' in Channel ' num2str(crop_channel_idx)]);
    set(handles.text_raw, 'String', ['Raw Data: Image ' num2str(mask_idx) ' in Channel ' num2str(channel_idx)]);
    pause(0.1);
end
handles.crop_mask_idx = mask_idx;
handles.mask_idx = mask_idx;

guidata(hObject, handles);

% --- Executes on button press in push_crop2.
function push_crop2_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in push_crop3.
function push_crop3_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in push_crop4.
function push_crop4_Callback(hObject, eventdata, handles)
% hObject    handle to push_crop4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes_logo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_logo

imagePath=['hello_kitty.tif'];
axes(hObject);
imshow(imread(imagePath), []);



% --- Executes on button press in push_selectroifreehand.
function push_selectroifreehand_Callback(hObject, eventdata, handles)
% hObject    handle to push_selectroifreehand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_image);
img = handles.thismask;
imshow(img, []);    colormap(handles.colormap_type);
%text(handles.colSize+5,20,['Image Index:' num2str(handles.mask_idx)],'BackgroundColor',[.7 .9 .7]);

h = imfreehand;
position = wait(h);     

% set(handles.startbutton, 'Visible', 'on');
% set(handles.savebutton, 'Visible', 'on');
if isempty(position)
   handles.iscrop = 0;
   guidata(hObject, handles);
   return;
end

valid_position = position;
valid_position(find(position(:,1)<0.5),1) = 1;
valid_position(find(position(:,2)<0.5),2) = 1;
valid_position(find(position(:,1)>handles.colSize),1) = handles.colSize;
valid_position(find(position(:,2)>handles.rowSize),2) = handles.rowSize;
segimg = roipoly(img,valid_position(:,1),valid_position(:,2));
c = contourc(double(segimg), [0 0]);
pixel_list = [c(1,2:end)' c(2,2:end)'];
       
axes(handles.axes_image);
img = handles.thismask;
imshow(img, []);    colormap(handles.colormap_type);
handles.pixel_list = pixel_list;
handles.segimg = segimg;
hold on;
plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
hold off;
handles.iscrop = 1;

% Select output directory
thispath=uigetdir('','Select an output directory');
if thispath==0
    guidata(hObject, handles);    
    thispath = warndlg('Please specify an output directory then press any crop button again', 'Warning', 'modal');
    return;
end


for c_idx=1:handles.numChannels
    [fpath,dirName,fno,fext]=getFilenameBody(handles.channelDirs{c_idx});
    if ~exist([thispath filesep dirName num2str(fno)])    
        mkdir(thispath,[dirName num2str(fno)]);
    end
    cropDir{c_idx} = [dirName num2str(fno)];
end

% Initializing waitbar
h=waitbar(0,'Processing...');

% Processing files
for c_idx=1:handles.numChannels
    for n_idx=1:handles.numImages(1)

        thisImg = imread([handles.fileList{c_idx}{n_idx,2} filesep handles.fileList{c_idx}{n_idx,1}]);
        if size(thisImg,3) == 3
            continue;
        end
        % Crop image        
        cropimg = uint16(thisImg).*uint16(segimg);
        imgC = cropimg(min(pixel_list(:,2)):max(pixel_list(:,2)), min(pixel_list(:,1)):max(pixel_list(:,1)));
%         common = ml_imgcommonpixel(thisImg);
%         imgC(find(imgC == 0)) = common;

        % Prepare filename with path
        filename=[thispath,filesep,cropDir{c_idx},filesep,'crop_',handles.fileList{c_idx}{n_idx,1}];
        handles.cropfileList{c_idx}{n_idx} = filename;
        
        % Write file to disk
        imwrite(imgC,filename, 'Compression', 'none');

        % Update waitbar
        waitbar(((c_idx-1)*handles.numImages(1) + n_idx)/(handles.numChannels*handles.numImages(1)),h);
    end
end

% Close waitbar
close(h);

% Display the crop images
if handles.numChannels > 1
    slider_channel_step(1) = 1/(handles.numChannels - 1);
    slider_channel_step(2) = 1/(handles.numChannels - 1);
    set(handles.slider_channel_crop, 'Visible', 'on','Max',handles.numChannels,...
        'Value',1,'Min',1,'sliderstep',slider_channel_step);    
else
    set(handles.slider_channel_crop, 'Visible', 'off');
end

crop_mask_idx = 1;
slider_time_step(1) = 1/handles.numImages(1);
slider_time_step(2) = 1/handles.numImages(1);

crop_channel_idx = 1;
thiscropmask = imread([handles.cropfileList{crop_channel_idx}{crop_mask_idx}]);
set(handles.slider_time_crop, 'Visible', 'on','Max',handles.numImages(1),...
    'Value',1,'Min',1,'sliderstep',slider_time_step);
 
axes(handles.axes_cropimage);    imshow(thiscropmask, []);   colormap(handles.colormap_type);
handles.thiscropmask = thiscropmask;
handles.crop_mask_idx = crop_mask_idx;
handles.crop_channel_idx = crop_channel_idx;

set(handles.text_cropped, 'String', ['Cropped Data: Image ' num2str(crop_mask_idx) ' in Channel ' num2str(crop_channel_idx)]);    
set(handles.text_cropped, 'Visible', 'on');

mask_idx = 1;
channel_idx = handles.channel_idx;
thismask = imread([handles.fileList{channel_idx}{mask_idx,2} filesep handles.fileList{channel_idx}{mask_idx,1}]);
axes(handles.axes_image);    imshow(thismask, []);   colormap(handles.colormap_type);
if handles.iscrop
    hold on;
    plot(handles.pixel_list(:,1),handles.pixel_list(:,2),'r--');
    hold off;
end

handles.thismask = thismask;
handles.mask_idx = mask_idx;
set(handles.text_raw, 'String', ['Raw Data: Image ' num2str(mask_idx) ' in Channel ' num2str(channel_idx)]);   
set(handles.slider_time,'Value',handles.mask_idx);

set(handles.push_play, 'Visible', 'on');
% Update handles structure
guidata(hObject, handles);


