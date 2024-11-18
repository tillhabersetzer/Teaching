function varargout = infra_meg_GUI(varargin)
% INFRA_MEG_GUI MATLAB code for infra_meg_GUI.fig
%      INFRA_MEG_GUI, by itself, creates a new INFRA_MEG_GUI or raises the existing
%      singleton*.
%
%      H = INFRA_MEG_GUI returns the handle to a new INFRA_MEG_GUI or the handle to
%      the existing singleton*.
%
%      INFRA_MEG_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INFRA_MEG_GUI.M with the given input arguments.
%
%      INFRA_MEG_GUI('Property','Value',...) creates a new INFRA_MEG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before infra_meg_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to infra_meg_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help infra_meg_GUI

% Last Modified by GUIDE v2.5 12-Dec-2018 09:24:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @infra_meg_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @infra_meg_GUI_OutputFcn, ...
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


% --- Executes just before infra_meg_GUI is made visible.
function infra_meg_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to infra_meg_GUI (see VARARGIN)

% Choose default command line output for infra_meg_GUI
handles.output = handles;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes infra_meg_GUI wait for user response (see UIRESUME)
% uiwait(handles.Infra_MEG);


% --- Outputs from this function are returned to the command line.
function varargout = infra_meg_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;


% --- Executes on button press in pushbutton_play_stop.
function pushbutton_play_stop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global InfraMEGButtonStatus
InfraMEGButtonStatus = get(hObject,'String');
switch InfraMEGButtonStatus
    case {'play'}
        set(hObject,'String','stop');
    case {'stop'}
        set(hObject,'String','play');
end
        


% --- Executes on button press in pushbutton_shift_count.
function pushbutton_shift_count_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_shift_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global InfraMEGNumShiftCount
c = round( str2double(get(handles.pushbutton_shift_count,'String') ) );
if c > 10
    c = 10;
elseif c < -10
    c = -10;
end
InfraMEGNumShiftCount = c;

function edit_shift_count_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shift_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shift_count as text
%        str2double(get(hObject,'String')) returns contents of edit_shift_count as a double


% --- Executes during object creation, after setting all properties.
function edit_shift_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shift_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
