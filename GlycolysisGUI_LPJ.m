function varargout = GlycolysisGUI_LPJ(varargin)
% GLYCOLYSISGUI_LPJ MATLAB code for GlycolysisGUI_LPJ.fig
%      GLYCOLYSISGUI_LPJ, by itself, creates a new GLYCOLYSISGUI_LPJ or raises the existing
%      singleton*.
%
%      H = GLYCOLYSISGUI_LPJ returns the handle to a new GLYCOLYSISGUI_LPJ or the handle to
%      the existing singleton*.
%
%      GLYCOLYSISGUI_LPJ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GLYCOLYSISGUI_LPJ.M with the given input arguments.
%
%      GLYCOLYSISGUI_LPJ('Property','Value',...) creates a new GLYCOLYSISGUI_LPJ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GlycolysisGUI_LPJ_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GlycolysisGUI_LPJ_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GlycolysisGUI_LPJ

% Last Modified by GUIDE v2.5 13-Oct-2014 11:07:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GlycolysisGUI_LPJ_OpeningFcn, ...
                   'gui_OutputFcn',  @GlycolysisGUI_LPJ_OutputFcn, ...
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


% --- Executes just before GlycolysisGUI_LPJ is made visible.
function GlycolysisGUI_LPJ_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GlycolysisGUI_LPJ (see VARARGIN)

% Choose default command line output for GlycolysisGUI_LPJ
handles.output = hObject;

% Update handles structure
handles.a = 1;
handles.b = 1;
handles.initialConditions = [2 2];
handles.results = {};
handles.timeSpan = 50;
guidata(hObject, handles);

% UIWAIT makes GlycolysisGUI_LPJ wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Plot stability boundary
a = linspace(0, 0.2, 10000);
b1=(2^(1/2)*((1 - 8*a).^(1/2) - 2*a + 1).^(1/2))/2;
b2=(2^(1/2)*(1 - (1 - 8*a).^(1/2) - 2*a).^(1/2))/2;
a1 = a(b1 == real(b1));
b1 = b1(b1 == real(b1));
a2 = a(b2 == real(b2));
b2 = b2(b2 == real(b2));
axes(handles.stabilityBoundary);
plot(a1, b1)
hold on
plot(a2, b2)
hold off
title('Stability Boundary')
xlabel(handles.stabilityBoundary, 'a');
ylabel(handles.stabilityBoundary, 'b');




% --- Outputs from this function are returned to the command line.
function varargout = GlycolysisGUI_LPJ_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% Updates phase plot when anything changes
function handles = updatePhasePlot(handles)
syms x1 x2;
a = handles.a;
b = handles.b;
x1dot = -x1+a * x2 + x1^2 * x2;
x2dot = b - a*x2 - x1^2 * x2;
[xstar,ystar]=solve(x1dot,x2dot);

% Define axes
axes(handles.phasePlot)

% Plot fixed points
plot(xstar,ystar,'ro')
hold on
title('Phase Portrait');
xlabel ('x_1 (Concentration of ADP)')
ylabel ('x_2 (Concentration of F6P)')

% Solve with initial conditions
tspan=handles.timeSpan;
x0=handles.initialConditions';
fun = sprintf('[-x(1)+ %g*x(2) + x(1)^2*x(2); %g - %g*x(2) - x(1)^2*x(2)]', a, b, a);
dx=inline(fun,'t','x');
[~,x]=ode45(dx,tspan,x0);
plot(x(:,1),x(:,2))
hold off

% Update results table
% trace was found using jacobian at fixed point [b, b^2/(a+b^2)]
% det always > 0
trace = -1 - a - b^2 + 2*b^2/(a+b^2);

if trace > 0 % unstable, limit cycle predicted
    det = b^2 + a;
    Tcycle = 2 * pi / sqrt(det);
    str = sprintf('Limit Cycle with period %gs', Tcycle);
    handles.results = [{handles.a, handles.b, str}; handles.results];
else % stable, check whether eigenvalues are real or complex
    jac = [ (2*b^3)/(b^2 + a) - 1,   b^2 + a; ...
           -(2*b^3)/(b^2 + a), - b^2 - a];
    eigen = eig(jac);
    if eigen == real(eigen)
        str = sprintf('Stable Node');
    else
        str = sprintf('Stable Focus');
    end
    handles.results = [{handles.a, handles.b, str}; handles.results];
end
set(handles.resultsTable,'Data',handles.results)


% --- Executes on mouse press over axes background.
function stabilityBoundary_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to stabilityBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Updates values of a and b (reaction constants)
coords = get(hObject, 'CurrentPoint');
handles.a = coords(1,1);
handles.b = coords(1,2);
astr = sprintf('Current Value of a: %g', handles.a);
bstr = sprintf('Current Value of b: %g', handles.b);
set(handles.adisp, 'String', astr);
set(handles.bdisp, 'String', bstr);
handles = updatePhasePlot(handles);
guidata(hObject, handles)





function adisp_Callback(hObject, eventdata, handles)
% hObject    handle to adisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of adisp as text
%        str2double(get(hObject,'String')) returns contents of adisp as a double


% --- Executes during object creation, after setting all properties.
function adisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to adisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bdisp_Callback(hObject, eventdata, handles)
% hObject    handle to bdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bdisp as text
%        str2double(get(hObject,'String')) returns contents of bdisp as a double


% --- Executes during object creation, after setting all properties.
function bdisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function initCondDisp_Callback(hObject, eventdata, handles)
% hObject    handle to initCondDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initCondDisp as text
%        str2double(get(hObject,'String')) returns contents of initCondDisp as a double


% --- Executes during object creation, after setting all properties.
function initCondDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initCondDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in initCondEdit.
function initCondEdit_Callback(hObject, eventdata, handles)
% hObject    handle to initCondEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Intial ADP Concentration','Initial F6P Concentration'};
dlg_title = 'Initial Conditions';
num_lines = 1;
ic = handles.initialConditions;
def = {num2str(ic(1)), num2str(ic(2))};
answer = inputdlg(prompt,dlg_title,num_lines,def);
assignin('base', 'answer', answer)
if isempty(answer)
    %do nothing
else
    % Update initial conditions
    handles.initialConditions = [eval(answer{1}), eval(answer{2})];
    ic = handles.initialConditions;
    icstr = sprintf('Initial Conditions: [%g %g]', ic(1), ic(2));
    set(handles.initCondDisp, 'String', icstr);
end
handles = updatePhasePlot(handles);
guidata(hObject, handles)



function timeSpanDisp_Callback(hObject, eventdata, handles)
% hObject    handle to timeSpanDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeSpanDisp as text
%        str2double(get(hObject,'String')) returns contents of timeSpanDisp as a double


% --- Executes during object creation, after setting all properties.
function timeSpanDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeSpanDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in timeSpanButton.
function timeSpanButton_Callback(hObject, eventdata, handles)
% hObject    handle to timeSpanButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Time Span (s)'};
dlg_title = 'Time Span';
num_lines = 1;
tspan = handles.timeSpan;
def = {num2str(tspan)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
assignin('base', 'answer', answer)
if isempty(answer)
    %do nothing
else
    % Update time span
    handles.timeSpan = eval(answer{1});
    tspan = handles.timeSpan;
    tstr = sprintf('Time Span: %gs', tspan);
    set(handles.timeSpanDisp, 'String', tstr);
end
handles = updatePhasePlot(handles);
guidata(hObject, handles)
