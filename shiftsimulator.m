function varargout = shiftsimulator(varargin)
% SHIFTSIMULATOR MATLAB code for shiftsimulator.fig
%      SHIFTSIMULATOR, by itself, creates a new SHIFTSIMULATOR or raises the existing
%      singleton*.
%
%      H = SHIFTSIMULATOR returns the handle to a new SHIFTSIMULATOR or the handle to
%      the existing singleton*.
%
%      SHIFTSIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHIFTSIMULATOR.M with the given input arguments.
%
%      SHIFTSIMULATOR('Property','Value',...) creates a new SHIFTSIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before shiftsimulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to shiftsimulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help shiftsimulator

% Last Modified by GUIDE v2.5 10-Aug-2016 21:47:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @shiftsimulator_OpeningFcn, ...
                   'gui_OutputFcn',  @shiftsimulator_OutputFcn, ...
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


% --- Executes just before shiftsimulator is made visible.
function shiftsimulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to shiftsimulator (see VARARGIN)

% Choose default command line output for shiftsimulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes shiftsimulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = shiftsimulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function rhod_Callback(hObject, eventdata, handles)
% hObject    handle to rhod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhod as text
%        str2double(get(hObject,'String')) returns contents of rhod as a double


% --- Executes during object creation, after setting all properties.
function rhod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rhog_Callback(hObject, eventdata, handles)
% hObject    handle to rhog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhog as text
%        str2double(get(hObject,'String')) returns contents of rhog as a double


% --- Executes during object creation, after setting all properties.
function rhog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rhog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phi_Callback(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phi as text
%        str2double(get(hObject,'String')) returns contents of phi as a double


% --- Executes during object creation, after setting all properties.
function phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delf1_Callback(hObject, eventdata, handles)
% hObject    handle to delf1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delf1 as text
%        str2double(get(hObject,'String')) returns contents of delf1 as a double


% --- Executes during object creation, after setting all properties.
function delf1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delf1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delf3_Callback(hObject, eventdata, handles)
% hObject    handle to delf3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delf3 as text
%        str2double(get(hObject,'String')) returns contents of delf3 as a double


% --- Executes during object creation, after setting all properties.
function delf3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delf3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delf5_Callback(hObject, eventdata, handles)
% hObject    handle to delf5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delf5 as text
%        str2double(get(hObject,'String')) returns contents of delf5 as a double


% --- Executes during object creation, after setting all properties.
function delf5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delf5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delg1_Callback(hObject, eventdata, handles)
% hObject    handle to delg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delg1 as text
%        str2double(get(hObject,'String')) returns contents of delg1 as a double


% --- Executes during object creation, after setting all properties.
function delg1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delg1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delg3_Callback(hObject, eventdata, handles)
% hObject    handle to delg3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delg3 as text
%        str2double(get(hObject,'String')) returns contents of delg3 as a double


% --- Executes during object creation, after setting all properties.
function delg3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delg3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delg5_Callback(hObject, eventdata, handles)
% hObject    handle to delg5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delg5 as text
%        str2double(get(hObject,'String')) returns contents of delg5 as a double


% --- Executes during object creation, after setting all properties.
function delg5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delg5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f1=5e6;        % Hz, Fundamental frequency
zq=8.84e6;     % Modulus of quartz, kg/m^2-s

sauerbrey=@(n,drho) 2*n*f1^2*drho/zq;
delfstarliq=@(n,etarho) (f1^1.5/zq)*(n*etarho/pi)^0.5*(1i-1);

dgliq =1300; %Shift in g3 only due to the liquid (assuming water)
 
if get(handles.filminliquid,'value')==1
    etarhoexpt=pi*dgliq^2*zq^2/(3*f1^3);
    Rliq=@(n,drho) delfstarliq(n,etarhoexpt)./sauerbrey(n,drho);  % defined in Elizabth's 2-layer paper
elseif get(handles.filminliquid,'value')==0
    etarhoexpt = 0;
    Rliq=@(n,drho) 0+0i;
end

d=@(n,d1,phi) d1.*n.^(1-phi./180);    % d/lambda
Dn=@(n,d1,phi)   2.*pi.*d(n,d1,phi).*(1-1i.*tand(phi./2));    % defined in Elizabth's 2015 Viscoelastic paper
delfstardn=@(Dn,Rliq)  -(Dn.^-2+Rliq.^2)./((cot(Dn)./Dn)+Rliq);
delfstar2layer=@(n,d1,phi,drho) delfstardn(Dn(n,d1,phi),Rliq(n,drho));

pd = 1e-6*str2double(get(handles.rhod,'string')); %mg/m^2
rhog3 = 1000*str2double(get(handles.rhog,'string')); %Pa-g/cm^3
phi = str2double(get(handles.phi,'string')); %deg.

% d1 = f1*pd./ ( sqrt(rhog3*sind(phi)+1i*rhog3*cosd(phi)) *(1-1i*tand(phi./2) ));

n=1;
refn=3;
lambdarhon = lambdarhof(refn, n, rhog3, phi); %Lambda rho for the harmonic of interest
d1 = pd./lambdarhon;

f1pred = (2.*1.*f1^2).*(pd*(real(delfstar2layer(1, d1, phi, pd))./zq));
f3pred = (2.*3.*f1^2).*(pd*(real(delfstar2layer(3, d1, phi, pd))./zq));
f5pred = (2.*5.*f1^2).*(pd*(real(delfstar2layer(5, d1, phi, pd))./zq));
g1pred = (2.*1.*f1^2).*(pd*(imag(delfstar2layer(1, d1, phi, pd))./zq));
g3pred = (2.*3.*f1^2).*(pd*(imag(delfstar2layer(3, d1, phi, pd))./zq));
g5pred = (2.*5.*f1^2).*(pd*(imag(delfstar2layer(5, d1, phi, pd))./zq));

set(handles.delf1,'string',num2str(f1pred));
set(handles.delf3,'string',num2str(f3pred));
set(handles.delf5,'string',num2str(f5pred));
set(handles.delg1,'string',num2str(g1pred));
set(handles.delg3,'string',num2str(g3pred));
set(handles.delg5,'string',num2str(g5pred));


function lrho = lambdarhof(refn, n, grho, phi)
if refn == n
    lrho = 1/(n*5e6)*(grho^0.5)/cosd(phi/2);
else
    lrho = 1/(n*5e6)*((grho*(n^(phi/90))/(refn^(phi/90)))^0.5)/cosd(phi/2);
end

% grhoi = 1.2e8*1000; %Get modulus in the right units
% phii = 20; %Get phi
% refn = 3; %Get refn
% 
% lambdarhon = lambdarhof(refn, n, grhoi, phii); %Lambda rho for the harmonic of interest
% dl = drho./lambdarhon; %Updates dl based on the new drho
% 
% delfstarc = sauerbrey(3,drho)*delfstar(dl, phii);%Calculates the frequency and dissipation
% 
% function F=sauerbrey(n,drho)
% % Calculates the sauerbry shift based on the harmonic and areal density.
% F=2*n*5e6^2*drho/8.84e6;
% 
% function F=delfstar(d,phi)  % input phase angle is in degrees
% % Calculates delfstar (rhs of delf/delfsn equation) with input of d/lambda
% % and phi
% F=-(1/((2*pi*d)*(1-1i*tand(phi/2))))* ...
%     tan(2*pi*d*(1-1i*tand(phi/2)));



    

