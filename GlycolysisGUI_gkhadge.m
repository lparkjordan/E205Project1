function varargout = GlycolysisGUI_gkhadge(varargin)
% GLYCOLYSISGUI_GKHADGE MATLAB code for GlycolysisGUI_gkhadge.fig
%      GLYCOLYSISGUI_GKHADGE, by itself, creates a new GLYCOLYSISGUI_GKHADGE or raises the existing
%      singleton*.
%
%      H = GLYCOLYSISGUI_GKHADGE returns the handle to a new GLYCOLYSISGUI_GKHADGE or the handle to
%      the existing singleton*.
%
%      GLYCOLYSISGUI_GKHADGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GLYCOLYSISGUI_GKHADGE.M with the given input arguments.
%
%      GLYCOLYSISGUI_GKHADGE('Property','Value',...) creates a new GLYCOLYSISGUI_GKHADGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GlycolysisGUI_gkhadge_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GlycolysisGUI_gkhadge_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GlycolysisGUI_gkhadge

% Last Modified by GUIDE v2.5 13-Oct-2014 10:33:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GlycolysisGUI_gkhadge_OpeningFcn, ...
                   'gui_OutputFcn',  @GlycolysisGUI_gkhadge_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT


% --- Executes just before GlycolysisGUI_gkhadge is made visible.
function GlycolysisGUI_gkhadge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GlycolysisGUI_gkhadge (see VARARGIN)

% Choose default command line output for GlycolysisGUI_gkhadge
handles.output = hObject;

% handles.limitCyclePeriod=0.2;
% set(handles.limitCyclePeriod,'String','cycle');
axes(handles.axes2)
ezplot('a+b^2-2*b^2 + (a+b^2)^2');
xlim([0 2]);
ylim([0 2]);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GlycolysisGUI_gkhadge wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% set(handles.axes2,'ButtonDownFcn',@axes2_ButtonDownFcn);
end

% --- Outputs from this function are returned to the command line.
function varargout = GlycolysisGUI_gkhadge_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% function axes2_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

dx = @(x1,x2,a,b) deal(-x1+a*x2+x1.^2.*x2, b-a*x2-x1.^2.*x2);

x1max = 1.5;
x2max = 1.5;
x1v = 0:x1max/10:x1max;
x2v = 0:x2max/10:x2max;
[x1n, x2n] = meshgrid(x1v,x2v);
x1v = x1n(:);
x2v = x2n(:);


axes(handles.axes1) % select figure handle
hold off

mouse=1;
tspan=[0 50];
x0 = [1,0]';
while(mouse==1)    
%     [t,x]=ode45(dx,tspan,x0);                 
%     plot(x(:,1),x(:,2))               
%%select new values of a and b and press the left mousebutton
%press the right mouse button if you want to stop.
        [av,bv,mouse]=ginput(1);
        [v1, v2] = dx(x1v,x2v, av,bv);
        axes(handles.axes1)
        hold off
%         plot(t,t);
        quiver(x1v,x2v,v1,v2);
        hold all
        plot(bv, bv/(av+bv^2),'o'); %stable point
        xlim([0 x1max]);
        ylim([0 x2max]);
        
        dxplot = @(t,x) [-x(1)+av*x(2)+x(1).^2.*x(2); bv-av*x(2)-x(1).^2.*x(2)];
        [t,x]=ode45(dxplot,tspan,x0);
        plot(x(:,1),x(:,2));
        
        
        trace = 1-2*bv^2/(av+bv^2)+av+bv^2;
        if (trace > 0)
            set(handles.limitCyclePeriod,'String','N/A');
            
        else
            det = av+bv^2; % product of eigenvalues (magnitude squared)
            %Trace is 2*Re(eigenvalues)
            w = sqrt(det^2 - (trace/2)^2); %Imaginary part of eigenvalues
            T = 2*pi/w;
            set(handles.limitCyclePeriod,'String',T);            
        end
        
        set(handles.avalue,'String',av);
        set(handles.bvalue,'String',bv);
        set(handles.stablepointx1,'String', bv);
        set(handles.stablepointx2,'String', bv/(av+bv^2));
end    
guidata(hObject,handles);

end
