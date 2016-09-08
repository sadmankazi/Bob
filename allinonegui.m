function varargout = allinonegui(varargin)

% Copyright (C) 2016 Kazi Sadman (Shull Research Group, Northwestern Uni.)
% 
% This is Version 3.3 of the GUI "BOB."
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>

% ALLINONEGUI MATLAB code for allinonegui.fig
%      ALLINONEGUI, by itself, creates a new ALLINONEGUI or raises the existing
%      singleton*.
%
%      H = ALLINONEGUI returns the handle to a new ALLINONEGUI or the handle to
%      the existing singleton*.
%
%      ALLINONEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALLINONEGUI.M with the given input arguments.
%
%      ALLINONEGUI('Property','Value',...) creates a new ALLINONEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before allinonegui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to allinonegui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help allinonegui

% Last Modified by GUIDE v2.5 07-Sep-2016 22:04:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @allinonegui_OpeningFcn, ...
    'gui_OutputFcn',  @allinonegui_OutputFcn, ...
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

% --- Executes just before allinonegui is made visible.
function allinonegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to allinonegui (see VARARGIN)

set(gcf,'Units','Pixels','Position',get(0,'ScreenSize')) % Make gui full size
% Choose default command line output for allinonegui
handles.output = hObject;
handles.color{1}=[1 0 0]; % red
handles.color{2}=[0 0.5 0]; % dark green
handles.color{3}=[0 0 1]; % blue
handles.marker{1}='+';
handles.marker{2}='o';
handles.marker{3}='s';

% set(0,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0.5 0])
set(0,'defaulttextfontsize',14)
set(0,'defaultlinemarkersize',10);
set(0,'defaultaxeslinewidth',2.25);
set(0,'defaultpatchlinewidth',2.25);
set(0,'defaultaxesfontsize',14)
set(0,'defaultlinelinewidth',1.5)

handles.currentloaded = 0;      % Will be used later to skip code for current data if current isn't loaded
handles.spectrasloaded = 0;     % Will be used later to throw error if spectras are not found
set(hObject,'toolbar','figure');

guidata(hObject,handles)
handles.output = hObject;

% --- Outputs from this function are returned to the command line.
function varargout = allinonegui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadQCM.
function loadQCM_Callback(hObject, eventdata, handles)

% [FileName,PathName] = uigetfile('~/Google Drive/*.mat', 'File Selector');
[FileName,PathName] = uigetfile('File Selector');
handles.filename = FileName;
% If someone cancels out of the uigetfile dialog, filename and pathname will
% both be 0. This checks if that is the case.
if ~FileName
    set(handles.statusupdate, 'String', 'No QCM data loaded!','Foregroundcolor','red');
    return
end

handles.pathname=PathName; %Save this pathname to handle structure so we can open this later for current data

fullpathname = fullfile(PathName,FileName);
[pathstr,name,ext] = fileparts(fullpathname);
handles.filename=name; % only save the name of the file

set(handles.textQCM, 'String', FileName);
guidata(hObject, handles)

m = matfile(fullpathname);
set(handles.dgliqn3,'String',num2str(round(m.freq_shift_ref(2,2))));
handles.dgliq3 = str2num(get(handles.dgliqn3,'String'));      % experimentally calculated delta gamma at third harmonic of liquid alone

b = find(~isnan(m.freq_shift(1:end,1))); %Indices of all elements in time that are not NaN
handles.t = m.freq_shift(b(1:end-1),1);
set(handles.endindex,'String',num2str(length(b)-1));
endidx = str2num(get(handles.endindex,'String'));

handles.delf{1} = m.freq_shift(b(1:end-1),2);
handles.delf{3} = m.freq_shift(b(1:end-1),4);
handles.delg{1} = m.freq_shift(b(1:end-1),3);
handles.delg{3} = m.freq_shift(b(1:end-1),5);
handles.delg{5} = m.freq_shift(b(1:end-1),7);
handles.delf{5} = m.freq_shift(b(1:end-1),6);

% for i = [1 3 5];
%     a = find(isnan(m.freq_shift(1:end,i)));
% 
% end

plot(handles.axes4, handles.t, handles.delf{1},'k+', handles.t, handles.delf{3},'ro',...
    handles.t, handles.delf{5},'bx');
plot(handles.axes6, handles.t, handles.delg{1},'k+', handles.t, handles.delg{3},'ro',...
    handles.t, handles.delg{5},'bx');

ylabel(handles.axes4,'\Deltaf (Hz)','fontweight','bold');
xlabel(handles.axes4,'t (min)','fontweight','bold');
xlim(handles.axes4,[0 handles.t(endidx)]);
legend(handles.axes4,'n=1','n=3','n=5','location','best')
ylabel(handles.axes6,'\Delta\Gamma (Hz)','fontweight','bold');
xlabel(handles.axes6,'t (min)','fontweight','bold');
xlim(handles.axes6,[0 handles.t(endidx)]);

% handles.t = 1.5;
% handles.delf{1} = -11363 ;
% handles.delf{3} = -35052;
% handles.delg{3} = 424;
% handles.delg{1} =  19;

% Save reference frequencies and dissipations
handles.ref(1,1)= m.freq_shift_ref(1,1); % f1
handles.ref(2,1)= m.freq_shift_ref(2,1); % g1
handles.ref(1,2)= m.freq_shift_ref(1,2); % f3
handles.ref(2,2)= m.freq_shift_ref(2,2); % g3
handles.ref(1,3)= m.freq_shift_ref(1,3); % f5
handles.ref(2,3)= m.freq_shift_ref(2,3); % g5

% if get(handles.dontloadspectras,'value')==0;
    if exist([PathName FileName(1:end-4) '_raw_spectras.mat'],'file');
        handles.spectrasloaded = 1;
        rawspectras = load([PathName FileName(1:end-4) '_raw_spectras.mat']);
        
        handles.fieldnames = fields(rawspectras);
        
        for i=1:numel(handles.t);
            data(i).time=num2str(handles.t(i));
        end
        
        for i=1:numel(data);
            dex = find(not(cellfun('isempty', strfind(handles.fieldnames,strrep(data(i).time,'.','dot')))));
            data(i).harmonic1 = rawspectras.(handles.fieldnames{dex(1)});
            data(i).harmonic3 = rawspectras.(handles.fieldnames{dex(2)});
            data(i).harmonic5 = rawspectras.(handles.fieldnames{dex(3)});
        end
        handles.data = data;
    else
                set(handles.statusupdate,'string','Warning: No raw spectras were found.','Foregroundcolor','red');
    end
% else
% end


if exist([PathName '1_03_CA_C01.txt'],'file')
    
    handles.currentloaded = 1;
    fullpathname1 = fullfile(PathName,'1_03_CA_C01.txt');
    
%     inid=fopen(fullpathname1);
%     refelectrode=textscan(inid,'%s %s',25,'delimiter',':');
%     potentialval=textscan(inid,'%s %s',30,'delimiter','\t');
%     handles.numheaderlines = str2num(refelectrode{1, 2}{2, 1}); % Read the number of header lines in text file
%     fclose(inid);
%     
%     set(handles.expinfo, 'String',[{refelectrode{1,2}{25,1}},{potentialval{1,1}{5,1}},{refelectrode{1, 2}{9, 1}}]);
    
    inid=fopen(fullpathname1);
    datavals= textscan(inid,'%f %f  %f  %f  %f %f %f %f %f %f %f %f %*[^\n]','HeaderLines',handles.numheaderlines, 'CollectOutput',true,...
        'delimiter', '\t', 'treatasempty', ',');
    fclose(inid);
    
    datavals= datavals{1};
    corrt= datavals(1:end,8)./60;  % convert to min
    corrt(2:2:end) = [];  % reduce # elements by 2
    corrt=corrt-corrt(1); % rezero first point
    corrv= datavals(:,10);
    corrv(2:2:end) = [];
    corri= datavals(:,11);
    corri(2:2:end) = [];
    
    handles.corrt=corrt;
    handles.corrv=corrv;
    handles.corri=corri;
    handles.corrt=handles.corrt-handles.corrt(1);
    
    handles.corrt(1)= [];  % Delete first data points
    handles.corrv(1)= [];
    handles.corri(1)= [];
else
end

% Update the status depending on what files were found and loaded:
if handles.currentloaded && handles.spectrasloaded;
    set(handles.statusupdate, 'String', 'QCM data, raw spectras and current loaded.','Foregroundcolor',[0 0.5 0]);
elseif ~handles.currentloaded && handles.spectrasloaded;
    set(handles.statusupdate, 'String', 'Warning: No current data was loaded.','Foregroundcolor',[1 0.5 0]);
else
    set(handles.statusupdate, 'String', 'Warning: No spectras were loaded.','Foregroundcolor',[1 0.5 0]);
end

guidata(hObject,handles)

% --- Executes on button press in plotqcmdata.
function plotqcmdata_Callback(hObject, eventdata, handles)

if ~isfield(handles,'filename')
    set(handles.statusupdate, 'String', 'Load QCM data first!','Foregroundcolor','red');
    return
end

handles.edissratio = [];
handles.eharmratio = [];
handles.d1out = [];
handles.drhoout = [];
handles.phiout = [];

f1=5e6;        % Fundamental frequency
zq=8.84e6;     % Modulus of quartz, kg/m^2-s

sauerbrey=@(n,drho) 2*n*f1^2*drho/zq;
delfstarliq=@(n,etarho) (f1^1.5/zq)*(n*etarho/pi)^0.5*(1i-1);

ndrho=3; % harmonic used for thickness calculation

i=3;
nh{1}=[1 3 3];
nh{2}=[3 5 3];
nh{3}=[3 5 5];

if get(handles.onelayer,'value')==0;
    handles.etarhoexpt=pi*handles.dgliq3^2*zq^2/(3*f1^3);
    Rliq=@(n,drho) delfstarliq(n,handles.etarhoexpt)./sauerbrey(n,drho);  % defined in Elizabth's 2-layer paper
elseif get(handles.onelayer,'value')==1;
    handles.etarhoexpt = 0;
    Rliq=@(n,drho) 0+0i;
end

d=@(n,d1,phi) d1.*n.^(1-phi./180);    % d/lambda
Dn=@(n,d1,phi)   2.*pi.*d(n,d1,phi).*(1-1i.*tand(phi./2));    % defined in Elizabth's 2015 Viscoelastic paper
grho=@(n,d1,drho,phi) n.^2.*f1^2.*(cosd(phi./2)).^2.*(1./d(n,d1,phi)).^2.*drho.^2;
delfstardn=@(Dn,Rliq)  -(Dn.^-2+Rliq.^2)./((cot(Dn)./Dn)+Rliq);
delfstar2layer=@(n,d1,phi,drho) delfstardn(Dn(n,d1,phi),Rliq(n,drho));
drhocalc=@(n,df,d1,phi,drho) df.*(zq./(2*n*f1^2))./real(delfstar2layer(n,d1,phi,drho));
rhodelta=@(n,grho,phi) ((grho.^0.5)./(2*pi.*n.*f1.*sind(phi./2)));

soln=[0.05, 45, 0.001];  % [not sure d1/lam/ phase angle/ Drho] These are the initial guesses for d1/lam, phi and drho.
lb=[0.01, 0, 0.0005];    % lower bounds on final solution
ub=[1, 90, 0.5];         % upper boqnds on final solution

inputsoln= soln;

endidx = str2num(get(handles.endindex,'String'));
start = str2num(get(handles.startindex,'String'));
step = str2num(get(handles.stepindex,'String'));

for i = 1:i;
    
    set(handles.statusupdate, 'String', [num2str(round(i/3*100)) ' % complete.'],'Foregroundcolor',[0 0.5 0]);
    drawnow
    
    for k=start:step:endidx
        
        % a, b and c correspond to the n1:n2,n3 calculation
        dfa=handles.delf{nh{i}(1)}(k); dfb=handles.delf{nh{i}(2)}(k); dfc=handles.delf{nh{i}(3)}(k); dgc=handles.delg{nh{i}(3)}(k);
        
        handles.eharmratio(k)=nh{i}(1)*dfb/(nh{i}(2)*dfa);   %Experimental harmonic ratio
        handles.edissratio(k)=-dgc/dfc;      %Experimantal dissipation ratio
        
        % generalzed to treat n1:n2,n3 calculation
        fdissratio=@(d1,phi,drho) -imag(delfstar2layer(nh{i}(3),d1,phi,drho))/(real(delfstar2layer(nh{i}(3),d1,phi,drho)));
        fharmratio=@(d1,phi,drho) real(delfstar2layer(nh{i}(2),d1,phi,drho))/(real(delfstar2layer(nh{i}(1),d1,phi,drho)));
        
        if get(handles.onelayer,'value')==0;
            f3tosolve=@(x) [fharmratio(x(1),x(2),x(3))-handles.eharmratio(k);...
                fdissratio(x(1),x(2),x(3))-handles.edissratio(k);...
                100*(drhocalc(ndrho, handles.delf{ndrho}(k),x(1),x(2),x(3))-x(3))];  %multiply by 100 to get enough accuracy
            inputsoln= soln;                               % Dynamic guesses for solving, i.e., previous solution is new initial guess
        elseif get(handles.onelayer,'value')==1;
            f3tosolve=@(x) [fharmratio(x(1),x(2),0)-handles.eharmratio(k);...
                fdissratio(x(1),x(2),0)-handles.edissratio(k)];
            soln=soln(1:2);  % [not sure d1/lam/ phase angle] These are the initial guesses for d1/lam, phi and drho.
            inputsoln= soln;
        end
        
        options = optimset('display','off','TolFun',10e-8); % Set the solver tolerance
        
        try
            
            if get(handles.onelayer,'value')==0;
                [soln,error]=lsqnonlin(f3tosolve,inputsoln,lb,ub,options);
            elseif get(handles.onelayer,'value')==1;
                [soln,error]=lsqnonlin(f3tosolve,inputsoln(1:2),lb(1:2),ub(1:2),options);
                soln(3) = dfa./real(delfstardn(Dn(ndrho,soln(1),soln(2)),0))*zq./(2*1*f1.^2);
            end
            
        catch Err
            soln = [NaN NaN NaN];
        end
        
        if ~isnan(soln(1))
            inputsoln = soln;
        end
        
        handles.d1out{i}(k) = soln(1);
        handles.phiout{i}(k) = soln(2);
        handles.drhoout{i}(k) = drhocalc(ndrho,handles.delf{ndrho}(k),soln(1),soln(2),soln(3));
        
    end
end
set(handles.statusupdate, 'String', 'Solved!','Foregroundcolor',[0 0.5 0]);

% In case the we don't solve for every point, delete the ones which weren't
% solved for so that the matrix dimensions mtch later for plotting
if step>1
    for i = [1 2 3];
        handles.d1out{i} = handles.d1out{i}(handles.d1out{i}~=0);
        handles.phiout{i} = handles.phiout{i}(handles.phiout{i}~=0);
        handles.drhoout{i} = handles.drhoout{i}(handles.drhoout{i}~=0);
    end
end

handles.d3out = d(3, handles.d1out{2}, handles.phiout{2});
handles.d5out = d(5, handles.d1out{2}, handles.phiout{2});

handles.grho3out{1} = grho(3, handles.d1out{1}, handles.drhoout{1}, handles.phiout{1});
handles.grho3out{2} = grho(3, handles.d1out{2}, handles.drhoout{2}, handles.phiout{2});
handles.grho3out{3} = grho(3, handles.d1out{3}, handles.drhoout{3}, handles.phiout{3});

handles.grho1out = grho(1, handles.d1out{2}, handles.drhoout{2}, handles.phiout{2});
handles.grho5out = grho(5, handles.d1out{2}, handles.drhoout{2}, handles.phiout{2});

handles.sauerbreycorrection1 = real(delfstar2layer(1, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}));
handles.sauerbreycorrection3 = real(delfstar2layer(3, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}));
handles.sauerbreycorrection5 = real(delfstar2layer(5, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}));

handles.rhodel1out = rhodelta(1, handles.grho1out, handles.phiout{2});
handles.rhodel3out = rhodelta(3, handles.grho3out{2}, handles.phiout{2});
handles.rhodel5out = rhodelta(5, handles.grho5out, handles.phiout{2});

% Calculate predicted shifts based on phi, d/lambda, drho that was just
% solved for:
f3pred = (2.*3.*f1^2).*(handles.drhoout{2}.*(real(delfstar2layer(3, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
f1pred = (2.*1.*f1^2).*(handles.drhoout{2}.*(real(delfstar2layer(1, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
g3pred = (2.*3.*f1^2).*(handles.drhoout{2}.*(imag(delfstar2layer(3, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
g1pred = (2.*2.*f1^2).*(handles.drhoout{2}.*(imag(delfstar2layer(1, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
f5pred = (2.*5.*f1^2).*(handles.drhoout{2}.*(real(delfstar2layer(5, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));
g5pred = (2.*5.*f1^2).*(handles.drhoout{2}.*(imag(delfstar2layer(5, handles.d1out{2}, handles.phiout{2}, handles.drhoout{2}))./zq));

% Viscoelastic Plots
plot(handles.axes14, repmat(handles.t(start:step:endidx),1,3), 1e6*cell2mat(handles.drhoout')', '+');
plot(handles.axes10, repmat(handles.t(start:step:endidx),1,3), 0.001*cell2mat(handles.grho3out')', '+'); %Multiply 0.001 to convert from Pa*kg/m^3 to Pa*g/cm^3
plot(handles.axes11, repmat(handles.t(start:step:endidx),1,3), cell2mat(handles.phiout')', '+');
% set(handles.axes11, 'YLim', [0 90])

% Harmonic Plots
cla(handles.axes4)
cla(handles.axes6)

plot(handles.axes4, handles.t(start:step:endidx), handles.delf{1}(start:step:endidx),'k+',...
    handles.t(start:step:endidx), handles.delf{3}(start:step:endidx),'ro',...
    handles.t(start:step:endidx), handles.delf{5}(start:step:endidx),'bx',...
    handles.t(start:step:endidx), f1pred,'c.',...
    handles.t(start:step:endidx), f3pred,'c.',...
    handles.t(start:step:endidx), f5pred,'c.');

plot(handles.axes6, handles.t(start:step:endidx), handles.delg{1}(start:step:endidx),'k+', ...
    handles.t(start:step:endidx), handles.delg{3}(start:step:endidx),'r+', ...
    handles.t(start:step:endidx), handles.delg{5}(start:step:endidx),'bx', handles.t(start:step:endidx), g5pred,'c.',...
    handles.t(start:step:endidx), g1pred,'c.',handles.t(start:step:endidx), g3pred,'c.');

legend(handles.axes4,'n=1','n=3','n=5','Pred')


if handles.currentloaded == 1;
    % Current Data & Plots:
    [~,handles.idx] = min(bsxfun(@(x,y)abs(x-y),handles.corrt',handles.t),[],2); % Find corresponding time points in longer current data
    totalhandles.chargedensity = cumtrapz(60.*handles.corrt,handles.corri)/1.27; % Total charge density delivered to solution until time t(x) in mC/cm^2
    Q = -1e-4.*totalhandles.chargedensity(handles.idx); % Factor 1e-4 changes units to mC/m^2 from mC/cm^2
    
    %Match charge and mass time points
    handles.masstocharge = (1e6.*handles.drhoout)./Q'; %mg/mC
else
end

linkaxes([handles.axes4,handles.axes6,handles.axes14,handles.axes10,...
    handles.axes11],'x');

% Now put all the xlabels and ylabels...
ylabel(handles.axes4,'\Deltaf (Hz)','fontweight','bold');
xlabel(handles.axes4,'t (min)','fontweight','bold');
xlim(handles.axes4,[handles.t(start) handles.t(endidx)])

ylabel(handles.axes6,'\Delta\Gamma (Hz)','fontweight','bold');
xlabel(handles.axes6,'t (min)','fontweight','bold');
xlim(handles.axes6,[handles.t(start) handles.t(endidx)])

xlabel(handles.axes14,'t (min)','fontweight','bold');
ylabel(handles.axes14,'\Delta M_A (mg/m^2)','fontweight','bold');
xlim(handles.axes14,[handles.t(start) handles.t(endidx)])
legend(handles.axes14, '133','353','355','location','best');

ylabel(handles.axes10,'\rho |G*_3| (Pa-g/cm^3)','fontweight','bold');
xlabel(handles.axes10,'t (min)','fontweight','bold');
xlim(handles.axes10,[handles.t(start) handles.t(endidx)])

ylabel(handles.axes11,'\phi (deg.)','fontweight','bold');
xlabel(handles.axes11,'t (min)','fontweight','bold');
xlim(handles.axes11,[handles.t(start) handles.t(endidx)])

zoom on
clear('start','step','endidx')
guidata(hObject, handles);

% --- Executes on button press in makecontours.
function makecontours_Callback(hObject, eventdata, handles)

if ~isfield(handles,'d1o ut')
    set(handles.statusupdate, 'String', 'No solved solutions!','Foregroundcolor','red');
    return
end

kc = str2num(get(handles.kc,'String'));

[val kc] = min(abs(handles.t-kc)); % find the index of the closest time point inputed

disstolerance=0.01;
harmtolerance=0.01;

eharmratio = handles.delf{3}(kc)/(3*handles.delf{1}(kc));
edissratio = -handles.delg{3}(kc)/handles.delf{3}(kc);

dissrange=[(1-disstolerance)*edissratio, (1+disstolerance)*edissratio]; % this is the range of dissipation ratios to consider (n=3)
harmrange=[(1-harmtolerance)*eharmratio, (1+harmtolerance)*eharmratio]; % range of harmonic ratios to consider (delf(3)/delf(1))

f1=5e6;
zq=8.84e6;

n=50; % resolution of map
dplot = linspace(0.0,0.23,n);  phiplot = linspace(0,90,n);
for i = 1:n
    ratio=3^(1-phiplot(i)/180);  % ratio of d/lambda at n=3 to d/lambda at n=1
    harmratio(1:n,1)=1;
    dissratio(1:n,1)=0;
    realdf(1:n,1)=-1;
    imagdf(1:n,1)=0;
    
    % Defining these functions again as Matlab can't view plotqcmdata callback
    sauerbrey=@(n,drho) 2*n*f1^2*drho/zq;
    delfstarliq=@(n,etarho) (f1^1.5/zq)*(n*etarho/pi)^0.5*(1i-1);
    
    if get(handles.onelayer,'value')==0;
        handles.etarhoexpt=pi*handles.dgliq3^2*zq^2/(3*f1^3);
        Rliq=@(n,drho) delfstarliq(n,handles.etarhoexpt)/sauerbrey(n,drho);  % defined in Elizabth's 2-layer paper
    elseif get(handles.onelayer,'value')==1;
        Rliq=@(n,drho) 0;
        handles.etarhoexpt = 0;
    end
    
    d=@(n,d1,phi) d1.*n.^(1-phi./180); %d/lambda
    Dn=@(n,d1,phi)   2*pi*d(n,d1,phi)*(1-1i*tand(phi/2)); %defined in Elizabth's 2-layer paper
    delfstardn=@(Dn,Rliq)  -(Dn^-2+Rliq^2)/((cot(Dn)/Dn)+Rliq);
    delfstar2layer=@(n,d1,phi,drho) delfstardn(Dn(n,d1,phi),Rliq(n,drho));
    
    fdissratio=@(d1,phi,drho) -imag(delfstar2layer(3,d1,phi,drho))/real(delfstar2layer(3,d1,phi,drho));
    fharmratio=@(d1,phi,drho) real(delfstar2layer(3,d1,phi,drho))/real(delfstar2layer(1,d1,phi,drho));
    
    for j = 2:n
        harmratio(i,j) = fharmratio(dplot(j)/ratio, phiplot(i), handles.drhoout{1}(kc));
        dissratio(i,j) = fdissratio(dplot(j)/ratio, phiplot(i), handles.drhoout{1}(kc));
        realdf(i,j) = real(delfstar2layer(3,dplot(j)/ratio,phiplot(i), handles.drhoout{1}(kc)));
        imagdf(i,j) = imag(delfstar2layer(3,dplot(j)/ratio,phiplot(i), handles.drhoout{1}(kc)));
    end
end

handles.d3out1 = d(3,handles.d1out{1},handles.phiout{1});
contourplots=figure;  % plot the dissipation ratio and harmonic ratio

fplot=subplot(2,2,1);
contourf(dplot,phiplot,realdf,linspace(-1.5,0,256),'edgecolor','none');
hold on
colormap(jet(256));
caxis([-1.5,0]);
colorbar
xlabel('d/\lambda_3')
ylabel('\phi','fontweight','bold')
title('(a) \Deltaf_3/\Deltaf_{s3}')

gplot=subplot(2,2,2);
contourf(dplot,phiplot,imagdf,linspace(0,1,256),'edgecolor','none');
hold on
colormap(jet(256));
caxis([0,1]);
colorbar
xlabel('d/\lambda_3')
ylabel('\phi','fontweight','bold')
title('(b) \Delta\Gamma_{3}/\Deltaf_{s3}','fontweight','bold')

plot(fplot, handles.d3out1(kc), handles.phiout{1}(kc),'k+','markersize',16,'color','red')
plot(gplot, handles.d3out1(kc), handles.phiout{1}(kc),'k+','markersize',16,'color','red')

harmplot=subplot(2,2,3);
contourf(dplot,phiplot,harmratio,linspace(0,2,256),'edgecolor','none');
hold on
contour(dplot,phiplot,dissratio,dissrange,'edgecolor','black','linewidth',2,'linestyle','--')
contour(dplot,phiplot,harmratio,harmrange,'edgecolor','black','linewidth',2)
colormap(jet(256));
caxis([0,1.5]);
colorbar
xlabel('d/\lambda_3')
ylabel('\phi','fontweight','bold')
numharmratio = handles.delf{3}(kc)./(3*handles.delf{1}(kc));
title(['(c) r_h= ' num2str(numharmratio)])

dissplot=subplot(2,2,4);
contourf(dplot,phiplot,dissratio,linspace(0,1,256),'edgecolor','none');
hold on
contour(dplot,phiplot,dissratio,dissrange,'edgecolor','black','linewidth',2,'linestyle','--')
contour(dplot,phiplot,harmratio,harmrange,'edgecolor','black','linewidth',2)
colormap(jet(256));
caxis([0,1]);
colorbar
xlabel('d/\lambda_3')
ylabel('\phi','fontweight','bold')
numdissratio = -handles.delg{3}(kc)./(handles.delf{3}(kc));
title(['(d) r_d= ' num2str(numdissratio)])
plot(handles.d3out(kc), handles.phiout{1}(kc),'k+','markersize',16,'color','red')
plot(harmplot, handles.d3out1(kc), handles.phiout{1}(kc),'k+','markersize',16,'color','red')

function dgliqn3_Callback(hObject, eventdata, handles)
% hObject    handle to dgliq32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latex_dgliq3 as text
%        str2double(get(hObject,'String')) returns contents of latex_dgliq3 as a double


% --- Executes on button press in radiobutton2.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2

function kc_Callback(hObject, eventdata, handles)
% hObject    handle to kc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kc as text
%        str2double(get(hObject,'String')) returns contents of kc as a double


% --- Executes during object creation, after setting all properties.
function kc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function dgliqn3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dgliq32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in refresh.
function refresh_Callback(hObject, eventdata, handles)
close(gcbf)
close all
clear all
reset(0)
firstrealgui
clc

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)

if ~isfield(handles,'d1out')
    set(handles.statusupdate, 'String', 'No solved values to plot!','Foregroundcolor','red');
    return
end

% kc = str2num(get(handles.kc,'String'));
% [val kc] = min(abs(handles.t(start:step:endidx)-kc)); % find the index of the closest time point inputed

% tbegin = 0;
% tend = 6;
% w0 = 5e6;

% handles.phiout([2 3]) = 90;
% deletepoints = [1 5 13 17];
% handles.t(deletepoints) = [];
% handles.phiout(deletepoints) = [];
% handles.drhoout(deletepoints) = [];
% handles.grho3out(deletepoints) = [];
% handles.rhodel3out(deletepoints) = [];
% handles.rhodel1out(deletepoints) = [];
% handles.sauerbreycorrection1(deletepoints) = [];
% handles.sauerbreycorrection3(deletepoints) = [];
% handles.t = handles.t-handles.t(1);
calcprops=figure('units','inches','Position', [2.5, 4, 15, 5]);
% Syntax: set(gcf,?position?,[a b W H])
% (a,b) = is the lower left corner
% H = the horizontal length of the window
% W = the vertical width of the window

endidx = str2num(get(handles.endindex,'String'));
start = str2num(get(handles.startindex,'String'));
step = str2num(get(handles.stepindex,'String'));

a = subplot(1,3,1);
plot(a, repmat(handles.t(start:step:endidx),1,3), 1e6*cell2mat(handles.drhoout')', '+');
xlabel('t (min)');
ylabel('\Delta M_A (mg/m^2)');
% xlim([tbegin tend])
% set(a, 'XTick', [tbegin:1:tend]);
title('(a)','fontweight','bold');
legend('133','353','355','location','best')
% hold on
% plot(handles.t(start:step:endidx)(kc),1e6*handles.drhoout(kc),'ro')
% hold off

c = subplot(1,3,3);
plot(c,repmat(handles.t(start:step:endidx),1,3), cell2mat(handles.phiout')', '+');
xlabel('t (min)');
ylabel('\phi (deg.)');
% xlim([tbegin tend])
% set(c, 'XTick', [tbegin:1:tend]);
title('(c)', 'fontweight','bold');
% set(gca,'YAxisLocation','Right');
% hold on
% plot(handles.t(kc-3),handles.phiout(kc-3),'ro')
% hold off

b = subplot(1,3,2);
% yyaxis(b,'left')
plot(b,repmat(handles.t(start:step:endidx),1,3), 0.001*cell2mat(handles.grho3out')', '+'); %Multiply 0.001 to convert from Pa*kg/m^3 to Pa*g/cm^3
ylabel('\rho |G^*_3| (Pa-g/cm^3)')
set(b, 'ycolor','k')
title('(b)')
% yyaxis(b,'right')
% plot(b,handles.t,(0.001*handles.grho3out)./(3*w0),'--','color',[0.4660 0.6740 0.1880])
% ylabel('\rho |\eta^*_3| (Pa*s-g/cm^3)')
% xlabel('t (min)')
% set(b, 'ycolor','k')
% ylim([0 1.65])
% xlim([tbegin tend])
% set(b, 'XTick', [tbegin:1:tend]);
% legend('\rho |G^*_3|','Contour','\rho |\eta^*_3|','location','northwest')

% % calcprops.PaperPosition=[0 0 17 5];
% % calcprops.PaperSize=[17 5];
% % print(calcprops,'viscoelastic plots.eps','-depsc')
% % saveas(calcprops,'viscoelastic plots.svg')

thickness = figure('units','inches','Position', [2.5, 4, 15, 5]);
a = subplot(1,3,1);
plot(a, handles.t(start:step:endidx), handles.d1out{1},'k+',...
    handles.t(start:step:endidx),handles.d3out,'ro',handles.t(start:step:endidx),handles.d5out,'bx');
ylabel('d/\lambda_n');
xlabel('t (min)');
% xlim([tbegin tend])
% set(a, 'XTick', [tbegin:1:tend]);
title('(a)','fontweight','bold');
legend('n=1','n=3','n=5','location','northwest')

b = subplot(1,3,2);
plot(b, handles.t(start:step:endidx), -1./handles.sauerbreycorrection1,'k+', handles.t(start:step:endidx), -1./handles.sauerbreycorrection3,'ro',handles.t(start:step:endidx), -1./handles.sauerbreycorrection5,'bx');
title('(b)','fontweight','bold');
ylabel('\rhod/(\rhod)_{sn}');
xlabel('t (min)');
% ylim(handles.axes16,[0 3])
% xlim([tbegin tend])
% ylim([0.85 1.35])
% set(b, 'XTick', [tbegin:1:tend]);

c = subplot(1,3,3);
plot(c, handles.t(start:step:endidx), handles.rhodel1out*1000.,'k+', handles.t(start:step:endidx), handles.rhodel3out*1000,'ro',handles.t(start:step:endidx), handles.rhodel5out*1000,'bx'); %Multiply by 1000 to convert g/m^2
title('(c)','fontweight','bold');
ylabel('\rho\delta_n (g/m^2)');
xlabel('t (min)');
% xlim([tbegin tend])
% ylim([0 10])
% set(c, 'XTick', [tbegin:1:tend]);

% thickness.PaperPosition=[0 0 12 5];
% thickness.PaperSize=[12 5];
% print(thickness,'thickness plots.eps','-depsc')
% saveas(thickness,'thickness plots.svg')

if handles.currentloaded == 1;
    
    handles.f=figure;
    
    [handles.f,a1,a2] = plotyy(handles.corrt,handles.corri/1.27,handles.corrt(handles.idx),handles.masstocharge,'plot','plot');
    set(handles.f(1),'ycolor','b')
    set(handles.f(2),'ycolor','r')
    set(a2,'Marker','+')
    a1.Color='b';
    a2.Color='r';
    ylabel(handles.f(2),'M/Q (mg/C)','fontweight','bold') % label left y-axis
    ylabel(handles.f(1),'Current density (mA/cm^2)','fontweight','bold') % label right y-axis
    xlabel(handles.f(1),'t (min)','fontweight','bold') % label x-axis
    xlim([0 handles.t(length(handles.t))])
    xlim(handles.f(1),[0 tend])
    xlim(handles.f(2),[0 tend])
    set(handles.f(1),'ytick',[-10:2:0])
    set(handles.f(2),'ytick',[-1*10^-4:2*10^-5:0])
else
end

% --- Executes on button press in phigplot.
function phigplot_Callback(hObject, eventdata, handles)

shiftsimulator;
% calcprops=figure;
% plot(0.001*handles.grho3out,  handles.phiout, 'b+');
% xlabel('\rho |G*_3| (Pa-g/cm^3)','fontweight','bold');
% ylabel('\phi (deg.)','fontweight','bold');
% ylim([0 90]);


% --- Executes on button press in shiftratios.
function shiftratios_Callback(hObject, eventdata, handles)
figure;
subplot(2,2,1);
plot(handles.t, handles.delf{3}./(3*handles.delf{1}), 'm+');
xlabel('t (min)','fontweight','bold');
ylabel('\Deltaf_3/3\Deltaf_1','fontweight','bold');
xlim([0 handles.t(length(handles.t))])
title('r_h','fontweight','bold');

subplot(2,2,2);
plot(handles.t, -handles.delg{3}./handles.delf{3}, 'm+');
xlabel('t (min)','fontweight','bold');
ylabel('-\Delta\Gamma_3/\Deltaf_3','fontweight','bold');
title('r_d','fontweight','bold');
xlim([0 handles.t(length(handles.t))])

subplot(2,2,3);
plot(handles.t, handles.delf{3}/3, 'o','color',handles.color{2});
hold on
plot(handles.t, handles.delf{1}, '+','color',handles.color{1});
hold off
xlabel('t (min)','fontweight','bold');
ylabel('-\Deltaf/n','fontweight','bold');
title('Normalized Shifts','fontweight','bold');
xlim([0 handles.t(length(handles.t))])

subplot(2,2,4);
plot(handles.t, handles.delg{3}/3, 'o','color',handles.color{2});
hold on
plot(handles.t, handles.delg{1}, '+','color',handles.color{1});
hold off
xlabel('t (min)','fontweight','bold');
ylabel('-\Delta\Gamma/n','fontweight','bold');
title('Normalized Shifts','fontweight','bold');
xlim([0 handles.t(length(handles.t))])
legend('n=3','n=1','location','best')


% --- Executes on button press in showhandles.
function showhandles_Callback(hObject, eventdata, handles)
keyboard


% --- Executes on button press in peakfitting.
function peakfitting_Callback(hObject, eventdata, handles)

if ~handles.spectrasloaded
    set(handles.statusupdate, 'String', 'No spectras were loaded to refit!','Foregroundcolor','red');
    return
end

%Call the peakrefitting GUI, whose output is handles.din
handles.din = peakrefitting('y', handles.figure1);

%Save the potentially updated shifts into the handle structure
handles.t = handles.din.t;
handles.delf{1} = handles.din.delf{1};
handles.delf{3} = handles.din.delf{3};
handles.delf{5} = handles.din.delf{5};
handles.delg{1} = handles.din.delg{1};
handles.delg{3} = handles.din.delg{3};
handles.delg{5} = handles.din.delg{5};

set(handles.endindex,'String',num2str(length(handles.t)));

handles.din = [];
set(handles.statusupdate,'string','Spectras successfully refitted!','foregroundcolor',[0 0.5 0])
guidata(hObject,handles)


% --- Executes on button press in dontloadspectras.
function dontloadspectras_Callback(hObject, eventdata, handles)
% hObject    handle to dontloadspectras (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dontloadspectras


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in onelayer.
function onelayer_Callback(hObject, eventdata, handles)
% hObject    handle to onelayer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onelayer


% --- Executes on button press in bulkcalg.
function bulkcalc_Callback(hObject, eventdata, handles)

if ~isfield(handles,'filename')
    set(handles.statusupdate, 'String', 'No QCM data loaded!','Foregroundcolor','red');
    return
end

if handles.delf{1} > handles.delg{1}
       set(handles.statusupdate, 'String', 'Shifts not in the bulk limit!','Foregroundcolor','red'); 
end

calcprops=figure('units','inches','Position', [2.5, 4, 12, 5]);

%Calculate phase angle based on bulk limit equations:
phi1=2.*atan(-handles.delf{1}./handles.delg{1}).*(180/pi);
phi3=2.*atan(-handles.delf{3}./handles.delg{3}).*(180/pi);
phi5=2.*atan(-handles.delf{5}./handles.delg{5}).*(180/pi);

subplot(1,2,1)
plot(handles.t,phi1,'k+',handles.t,phi3,'ro',handles.t,phi5,'bx');
xlabel('Time (min)');
ylabel('\phi (Degrees)');

title('(a) \phi')

% Calculate Complex Modulus from phase angle
pG1=(handles.delf{1}.*pi.*8.84e6./(5e6.*sind(phi1./2))).^2;
pG3=(handles.delf{3}.*pi.*8.84e6./(5e6.*sind(phi3./2))).^2;
pG5=(handles.delf{5}.*pi.*8.84e6./(5e6.*sind(phi5./2))).^2;

subplot(1,2,2)
yyaxis left
plot(handles.t,pG1./1000,'k+',handles.t,pG3./1000,'ro',handles.t,pG5./1000,'bx');
xlabel('Time (min)');
ylabel('\rho|G_n^*| (Pa-g/cm^3)');
title('(b) \rho|G_n^*|')
set(gca, 'ycolor','k')

yyaxis right 
plot(handles.t,pG3./(1000*2*pi*3*5e6),'g--')
legend('n=1','n=3','n=5','\rho|\eta^*|', 'location','best');
ylabel('\rho|\eta_3^*| (Pa*s-g/cm^3)');
set(gca, 'ycolor','k')

function endindex_Callback(hObject, eventdata, handles)
% hObject    handle to endindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endindex as text
%        str2double(get(hObject,'String')) returns contents of endindex as a double


% --- Executes during object creation, after setting all properties.
function endindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stepindex_Callback(hObject, eventdata, handles)
% hObject    handle to stepindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepindex as text
%        str2double(get(hObject,'String')) returns contents of stepindex as a double


% --- Executes during object creation, after setting all properties.
function stepindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startindex_Callback(hObject, eventdata, handles)
% hObject    handle to startindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startindex as text
%        str2double(get(hObject,'String')) returns contents of startindex as a double


% --- Executes during object creation, after setting all properties.
function startindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in thicknessplots.
function thicknessplots_Callback(hObject, eventdata, handles)

if ~isfield(handles,'d1out')
    set(handles.statusupdate, 'String', 'No solved solutions!','Foregroundcolor','red');
    return
end

endidx = str2num(get(handles.endindex,'String'));
start = str2num(get(handles.startindex,'String'));
step = str2num(get(handles.stepindex,'String'));

% Thickness Plots
calcprops=figure('units','inches','Position', [2.5, 4, 15, 5]);
subplot(1,3,1)
plot( handles.t(start:step:endidx), handles.d1out{1},'k+', handles.t(start:step:endidx),handles.d3out,'ro',handles.t(start:step:endidx),handles.d5out,'bx');
ylabel('d/\lambda_n','fontweight','bold');
xlabel('t (min)','fontweight','bold');
xlim([0 handles.t(end)])
title('(a)','fontweight','bold')
legend('n=1','n=3','n=5','Location','best')

subplot(1,3,2)
plot(handles.t(start:step:endidx), -1./handles.sauerbreycorrection1,'k+', handles.t(start:step:endidx), -1./handles.sauerbreycorrection3,'ro',handles.t(start:step:endidx), -1./handles.sauerbreycorrection5,'bx');
ylabel('\rhod/(\rhod)_{sn}','fontweight','bold');
xlabel('t (min)','fontweight','bold');
% ylim(handles.axes16,[-3 3])
xlim([0 handles.t(end)])
title('(b)','fontweight','bold')

subplot(1,3,3)
plot( handles.t(start:step:endidx), handles.rhodel1out*1000.,'k+', handles.t(start:step:endidx), handles.rhodel3out*1000,'ro',handles.t(start:step:endidx), handles.rhodel5out*1000,'bx'); %Multiply by 1000 to convert g/m^2
ylabel('\rho\delta_n (g/m^2)','fontweight','bold');
xlabel('t (min)','fontweight','bold');
% ylim([0 20]);
title('(c)','fontweight','bold')
xlim([0 handles.t(end)])


% --- Executes on button press in currentplot.
function currentplot_Callback(hObject, eventdata, handles)

if handles.currentloaded == 0;
    set(handles.statusupdate,'string','No current data found!','foregroundcolor','red')
    return
end

figure;
yyaxis(handles.axes22,'left')
plot(handles.corrt,handles.corri/1.27,'b-')
ylabel('M/Q (mg/mC)','fontweight','bold') % label left y-axis
    
yyaxis(handles.axes22,'right')
plot(handles.corrt(handles.idx),handles.masstocharge,'r+')
ylabel('Current density (mA/cm^2)','fontweight','bold') % label right y-axis
xlabel('t (min)','fontweight','bold') % label x-axis


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9
