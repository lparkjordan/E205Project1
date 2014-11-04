function varargout = E205PredatorPrey(varargin)
% E205PREDATORPREY MATLAB code for E205PredatorPrey.fig
%      E205PREDATORPREY, by itself, creates a new E205PREDATORPREY or raises the existing
%      singleton*.
%
%      H = E205PREDATORPREY returns the handle to a new E205PREDATORPREY or the handle to
%      the existing singleton*.
%
%      E205PREDATORPREY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in E205PREDATORPREY.M with the given input arguments.
%
%      E205PREDATORPREY('Property','Value',...) creates a new E205PREDATORPREY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before E205PredatorPrey_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to E205PredatorPrey_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help E205PredatorPrey

% Last Modified by GUIDE v2.5 03-Nov-2014 17:11:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @E205PredatorPrey_OpeningFcn, ...
                   'gui_OutputFcn',  @E205PredatorPrey_OutputFcn, ...
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


% --- Executes just before E205PredatorPrey is made visible.
function E205PredatorPrey_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to E205PredatorPrey (see VARARGIN)

% Choose default command line output for E205PredatorPrey
handles.output = hObject;

% Update handles structure
handles.mu = 1;
handles.sigma = 1;
handles.initialConditions = [0.5 0.5];
handles.results = {};
handles.timeSpan = 50;

guidata(hObject, handles);

% UIWAIT makes E205PredatorPrey wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = E205PredatorPrey_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Updates phase plot when anything changes
function handles = updatePhasePlot(handles)
syms x1 x2;
mu = handles.mu;
sigma = handles.sigma;
x1dot = x1 - x1 * x2 - mu * x1^2;
x2dot = x1 * x2 - x2 - sigma * x2 * x1dot;
[xstar,ystar]=solve(x1dot,x2dot);

% Define axes
axes(handles.phase_plot_axes)

% Plot fixed points
plot(xstar,ystar,'ro')
hold on
title('Phase Portrait');
xlabel ('x_1 (Prey)')
ylabel ('x_2 (Predator)')

% Solve with initial conditions
tspan=handles.timeSpan;
x0=handles.initialConditions';
fun = sprintf('[x(1) - x(1)*x(2) -  %g*x(1)^2; x(1)*x(2) - x(2) - %g*x(2) * (x(1) - x(1)*x(2) -  %g*x(1)^2) ]', mu, sigma, mu);
dx=inline(fun,'t','x');
[~,x]=ode45(dx,tspan,x0);
plot(x(:,1),x(:,2))
hold off

% % Update results table
% % trace was found using jacobian at fixed point [b, b^2/(a+b^2)]
% % det always > 0
% trace = -1 - a - b^2 + 2*b^2/(a+b^2);
% 
% if trace > 0 % unstable, limit cycle predicted
%     det = b^2 + a;
%     Tcycle = 2 * pi / sqrt(det);
%     str = sprintf('Limit Cycle with period %gs', Tcycle);
%     handles.results = [{handles.a, handles.b, str}; handles.results];
% else % stable, check whether eigenvalues are real or complex
%     jac = [ (2*b^3)/(b^2 + a) - 1,   b^2 + a; ...
%            -(2*b^3)/(b^2 + a), - b^2 - a];
%     eigen = eig(jac);
%     if eigen == real(eigen)
%         str = sprintf('Stable Node');
%     else
%         str = sprintf('Stable Focus');
%     end
%     handles.results = [{handles.a, handles.b, str}; handles.results];
% end
% set(handles.resultsTable,'Data',handles.results)


% --- Executes on mouse press over axes background.
function parameter_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to stabilityBoundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Updates values of a and b (reaction constants)
coords = get(hObject, 'CurrentPoint');
handles.mu = coords(1,1);
handles.sigma = coords(1,2);
mu_str = sprintf('Current Value of mu: %g', handles.mu);
sigma_str = sprintf('Current Value of sigma: %g', handles.sigma);
set(handles.mu_disp, 'String', mu_str);
set(handles.sigma_disp, 'String', sigma_str);
handles = updatePhasePlot(handles);
guidata(hObject, handles)

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function mu_Callback(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu_disp as text
%        str2double(get(hObject,'String')) returns contents of mu_disp as a double


% --- Executes during object creation, after setting all properties.
function mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_disp as text
%        str2double(get(hObject,'String')) returns contents of sigma_disp as a double


% --- Executes during object creation, after setting all properties.
function sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mu_disp_Callback(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu_disp as text
%        str2double(get(hObject,'String')) returns contents of mu_disp as a double


% --- Executes during object creation, after setting all properties.
function mu_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_disp_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_disp as text
%        str2double(get(hObject,'String')) returns contents of sigma_disp as a double


% --- Executes during object creation, after setting all properties.
function sigma_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_disp (see GCBO)
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
