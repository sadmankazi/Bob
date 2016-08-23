function varargout = peakrefitting(varargin)
% PEAKREFITTING MATLAB code for peakrefitting.fig
%      PEAKREFITTING, by itself, creates a new PEAKREFITTING or raises the existing
%      singleton*.
%
%      H = PEAKREFITTING returns the handle to a new PEAKREFITTING or the handle to
%      the existing singleton*.
%
%      PEAKREFITTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEAKREFITTING.M with the given input arguments.
%
%      PEAKREFITTING('Property','Value',...) creates a new PEAKREFITTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before peakrefitting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to peakrefitting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help peakrefitting

% Last Modified by GUIDE v2.5 30-Jun-2016 19:42:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @peakrefitting_OpeningFcn, ...
    'gui_OutputFcn',  @peakrefitting_OutputFcn, ...
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


% --- Executes just before peakrefitting is made visible.
function peakrefitting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to peakrefitting (see VARARGIN)
handles.currenttimeidx =1;
handles.noplotupdate =0;  % Use this to suppress updating plots in in "refitall"
handles.peaklocation{1}(1:1500) = NaN;
handles.peaklocation{3}(1:1500) = NaN;
handles.peaklocation{5}(1:1500) = NaN;
handles.gmax{1}(1:1500) = NaN;
handles.gmax{3}(1:1500) = NaN;
handles.gmax{5}(1:1500) = NaN;

dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'y'));
if (isempty(mainGuiInput)) ...
        || (length(varargin) <= mainGuiInput) ...
        || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.onelayerguiMain = varargin{mainGuiInput+1};
    
    % Obtain handles using GUIDATA with the caller's handle
    mainHandles = guidata(handles.onelayerguiMain);
    
    handles.main = mainHandles;
end

if dontOpen
    set(handles.statusupdate, 'String', 'See command line error!','Foregroundcolor','red');
    disp('-----------------------------------------------------');
    disp('Improper input arguments. Pass a property value pair')
    disp('whose name is "onelayerguiwithcond" and value is the handle')
    disp('to the onelayerguiwithcond figure, e.g:');
    disp('   x = onelayerguiwithcond()');
    disp('   condfig3(''onelayerguiwithcond'', x)');
    disp('-----------------------------------------------------');
else
    
    set(0,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0.5 0])
    set(0,'defaulttextfontsize',14)
    set(0,'defaultlinemarkersize',10);
    set(0,'defaultaxeslinewidth',2.25);
    set(0,'defaultpatchlinewidth',2.25);
    set(0,'defaultaxesfontsize',14)
    set(0,'defaultlinelinewidth',1.5)
end

% Update handles structure
guidata(hObject, handles);
plotshifts(hObject, eventdata, handles)
plotspectra(hObject, eventdata, handles)
set(handles.statusupdate, 'String', 'Spectras succesfully imported.','Foregroundcolor', [0 0.5 0]);

% datacursormode on
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'UpdateFcn',{@myupdatefcn, handles})

% UIWAIT makes peakrefitting wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = peakrefitting_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.main;
delete(hObject);

function plotshifts(hObject, eventdata, handles)
cla(handles.frequencyplot)
cla(handles.dissipationplot)

plot(handles.frequencyplot,handles.main.t, handles.main.delf{1}, 'k+',handles.main.t, handles.main.delf{3}/3, 'ro',handles.main.t, handles.main.delf{5}/5,'bx',...
    handles.main.t(handles.currenttimeidx), handles.main.delf{1}(handles.currenttimeidx),'gx',...
    handles.main.t(handles.currenttimeidx), handles.main.delf{3}(handles.currenttimeidx)/3,'gx',...
    handles.main.t(handles.currenttimeidx), handles.main.delf{5}(handles.currenttimeidx)/5,'gx');

plot(handles.dissipationplot, handles.main.t, handles.main.delg{1}, 'kx',handles.main.t, handles.main.delg{3}/3, 'ro',handles.main.t, handles.main.delg{5}/5,'bx',...
    handles.main.t(handles.currenttimeidx), handles.main.delg{1}(handles.currenttimeidx),'gx',...
    handles.main.t(handles.currenttimeidx), handles.main.delg{3}(handles.currenttimeidx)/3,'gx',...
    handles.main.t(handles.currenttimeidx), handles.main.delg{5}(handles.currenttimeidx)/5,'gx');

ylabel(handles.dissipationplot,'\Delta\Gamma/n (Hz)','fontweight','bold');
xlabel(handles.dissipationplot,'t (min)','fontweight','bold');

ylabel(handles.frequencyplot,'\Deltaf/n (Hz)','fontweight','bold');
xlabel(handles.frequencyplot,'t (min)','fontweight','bold');
legend(handles.frequencyplot,'n=1','n=3','n=5','Active point','location','best')
set(handles.showingtimepoint,'string',num2str(handles.main.t(handles.currenttimeidx)))


function plotspectra(hObject, eventdata, handles)
cla(handles.n1peak)
cla(handles.n3peak)
cla(handles.n5peak)

plot(handles.n1peak, handles.main.data(handles.currenttimeidx).harmonic1(1:end,1), ...
    handles.main.data(handles.currenttimeidx).harmonic1(1:end,2), 'r-',...
    handles.main.data(handles.currenttimeidx).harmonic1(1:end,1), ...
    handles.main.data(handles.currenttimeidx).harmonic1(1:end,4), 'k-',...
    handles.peaklocation{1}(handles.currenttimeidx), handles.gmax{1}(handles.currenttimeidx),'cx');

ylabel(handles.n1peak,'Conductance','fontweight','bold');
xlabel(handles.n1peak,'Frequency (Hz)','fontweight','bold');
legend(handles.n1peak,'Raw','Fitted','location','best')

plot(handles.n3peak, handles.main.data(handles.currenttimeidx).harmonic3(1:end,1),...
    handles.main.data(handles.currenttimeidx).harmonic3(1:end,2), 'r-',...
    handles.main.data(handles.currenttimeidx).harmonic3(1:end,1),...
    handles.main.data(handles.currenttimeidx).harmonic3(1:end,4), 'k-',...
    handles.peaklocation{3}(handles.currenttimeidx), handles.gmax{3}(handles.currenttimeidx),'cx');

ylabel(handles.n3peak,'Conductance','fontweight','bold');
xlabel(handles.n3peak,'Frequency (Hz)','fontweight','bold');

plot(handles.n5peak, handles.main.data(handles.currenttimeidx).harmonic5(1:end,1),...
    handles.main.data(handles.currenttimeidx).harmonic5(1:end,2), 'r-',...
    handles.main.data(handles.currenttimeidx).harmonic5(1:end,1),...
    handles.main.data(handles.currenttimeidx).harmonic5(1:end,4), 'k-',...
    handles.peaklocation{5}(handles.currenttimeidx), handles.gmax{5}(handles.currenttimeidx),'cx');

ylabel(handles.n5peak,'Conductance','fontweight','bold');
xlabel(handles.n5peak,'Frequency (Hz)','fontweight','bold');

set(handles.showingtimepoint,'string',num2str(handles.main.t(handles.currenttimeidx)))

% --- Executes on button press in nextpoint.
function nextpoint_Callback(hObject, eventdata, handles)

if handles.currenttimeidx <= length(handles.main.t) && handles.currenttimeidx >= 1
    set(handles.statusupdate, 'String', 'QCM & raw spectras loaded.','Foregroundcolor',[0 0.5 0]);
end

if handles.currenttimeidx < length(handles.main.t)
    handles.currenttimeidx = handles.currenttimeidx +1;
    set(handles.showingtimepoint,'string',num2str(handles.main.t(handles.currenttimeidx)))
else
    set(handles.statusupdate, 'String', 'It is the last point!','Foregroundcolor','red');
    return
end

cla(handles.n1peak);
cla(handles.n3peak);
cla(handles.n5peak);
guidata(hObject,handles)
plotspectra(hObject,eventdata, handles);
plotshifts(hObject, eventdata, handles);

% --- Executes on button press in previouspoint.
function previouspoint_Callback(hObject, eventdata, handles)

if handles.currenttimeidx <= length(handles.main.t) && handles.currenttimeidx>=1
    set(handles.statusupdate, 'String', 'QCM & raw spectras loaded.','Foregroundcolor',[0 0.5 0]);
end

if handles.currenttimeidx > 1
    handles.currenttimeidx = handles.currenttimeidx -1;
    set(handles.showingtimepoint,'string',num2str(handles.main.t(handles.currenttimeidx)))
else
    set(handles.statusupdate, 'String', 'It is the first point!','Foregroundcolor','red');
    return
end

cla(handles.n1peak);
cla(handles.n3peak);
cla(handles.n5peak);
guidata(hObject,handles)
plotspectra(hObject,eventdata, handles);
plotshifts(hObject, eventdata, handles);

% --- Executes on button press in showhandles.
function showhandles_Callback(hObject, eventdata, handles)
keyboard


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'leftarrow')
    previouspoint_Callback(hObject, eventdata, handles)
elseif strcmp(eventdata.Key, 'rightarrow')
    nextpoint_Callback(hObject, eventdata, handles)
end

function [peaklocation HWHM Gmax fitteddata] = fitspectra(freq, conductance, susceptance, handles, hObject)

[index, peak_detect, locs, width] = findrelevantpeaks(freq, conductance);

numpeaks = length(index);          % Hopefully at this point numpeaks is < 3
phi=0;                             % Assume rotation angle is 0
factor_range_fit= str2num(get(handles.fitfactor,'String'));             % Range of frequency to fit over on either side of Gmax
halfg = 0.5*width(1);              % Value of the half dissipation
coffset = min(conductance);        % Conductance offset
soffset = mean([susceptance(1), susceptance(end)]); %susceptance offset
Gmax = max(conductance);           % Maximum in G relative to the baseline
f0 = freq(index);                  % Frequency at the identified peak

% Determine the range over which to fit:
freq_upperlim = locs(1) + halfg*factor_range_fit;
freq_lowerlim = locs(1) - halfg*factor_range_fit;

% Now find the indexes in freq that match the upper and lower limits:
[c, freqidx_upper] = min(abs(freq-freq_upperlim));
[c, freqidx_lower] = min(abs(freq-freq_lowerlim));

p0= [f0(1) halfg(1) phi  Gmax(1) coffset soffset];            % Array of initial guesses for lorentzcond function
freq1 = freq(freqidx_lower:freqidx_upper);                    % Appropriate freq range from fit factor
cond1 = smooth(conductance(freqidx_lower:freqidx_upper),3);   % Appropriate cond range from fit factor (added some small smoothing)
sus1 = smooth(susceptance(freqidx_lower:freqidx_upper),3);

options=optimset('display','off','tolfun',1e-10,'tolx',1e-10,'MaxIter',5000);
lb(1:length(p0)) = -Inf; %Assigns the lower bound to the parameters to -Inf
ub(1:length(p0)) = Inf;  %Assigns the upper bound to the parameters to Inf
ub(3) = 90;              %Changes the phase angle upper bound to 90

if get(handles.fitcond,'value') == 1;
    p = lsqcurvefit(@lfun4c,p0,freq1,cond1,lb,ub,options);
elseif get(handles.fitsus,'value') == 1;
    p = lsqcurvefit(@lfun4s,p0,freq1,sus1,lb,ub,options);
elseif get(handles.fitcondsus,'value') == 1;
    p = lsqcurvefit(@lfun4cs,p0,freq1,[cond1 sus1],lb,ub,options);
end

fitteddata = [];

fitteddata = p(4).*((((freq.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(freq.^2)).^2)+...
    ((freq.^2).*((2.*p(2)).^2)))).*cosd(p(3))-((((p(1)).^2-freq.^2)).*freq.*(2.*p(2)))./...
    (((((p(1)).^2)-(freq.^2)).^2)+((freq.^2).*((2.*p(2)).^2))).*sind(p(3))) +p(5);

handles.fitteddata = fitteddata;
% [Gmax, I] = max(fitteddata);
% val = (Gmax - min(fitteddata));
% halfG = min(fitteddata) + val/2;
% halfmax = find(fitteddata>halfG); % All indexes of points with fitteddata>Gmax/2
HWHM = p(2);
peaklocation = p(1);
% hpercent = ((HWHM1-HWHM)/HWHM)*100
%  peaklocation = freq(I);
% percent = ((peaklocation1-peaklocation)/peaklocation) *100
% HWHM = 0.5*abs(freq(halfmax(1)) - freq(halfmax(length(halfmax)))); % Half width half max
guidata(hObject,handles)


% --- Executes on button press in refitn1.
function handles = refitn1_Callback(hObject, eventdata, handles)

if get(handles.useprevn1,'value')==0;
    [peaklocation, HWHM, Gmax, fittedfreqdata] = fitspectra(handles.main.data(handles.currenttimeidx).harmonic1(1:end,1),...
        handles.main.data(handles.currenttimeidx).harmonic1(1:end,2),...
        handles.main.data(handles.currenttimeidx).harmonic1(1:end,3),...
        handles, hObject);
else
    [peaklocation, HWHM, Gmax, fittedfreqdata] = fitspectra(handles.main.data(handles.currenttimeidx-1).harmonic1(1:end,1),...
        handles.main.data(handles.currenttimeidx).harmonic1(1:end,2),...
        handles.main.data(handles.currenttimeidx).harmonic1(1:end,3),...
        handles, hObject);
end

if get(handles.fitcond,'value')==1;
    handles.main.data(handles.currenttimeidx).harmonic1(1:end,4) = fittedfreqdata;
elseif get(handles.fitcond,'value')==1;
    handles.main.data(handles.currenttimeidx).harmonic1(1:end,5) = fittedfreqdata;
end

handles.main.delg{1}(handles.currenttimeidx) = HWHM- handles.main.ref(2,1);
handles.main.delf{1}(handles.currenttimeidx) = peaklocation- handles.main.ref(1,1);

handles.gmax{1}(handles.currenttimeidx) = Gmax;
handles.peaklocation{1}(handles.currenttimeidx) = peaklocation;

guidata(hObject,handles)

if handles.noplotupdate ==0;
    plotspectra(hObject,eventdata, handles);
    plotshifts(hObject, eventdata, handles);
else
end



% --- Executes on button press in refitn3.
function handles = refitn3_Callback(hObject, eventdata, handles)

if get(handles.useprevn3,'value')==0;
    [peaklocation, HWHM, Gmax, fittedfreqdata] = fitspectra(handles.main.data(handles.currenttimeidx).harmonic3(1:end,1),...
        handles.main.data(handles.currenttimeidx).harmonic3(1:end,2),...
        handles.main.data(handles.currenttimeidx).harmonic3(1:end,3),...
        handles, hObject);
else
    [peaklocation, HWHM, Gmax, fittedfreqdata] = fitspectra(handles.main.data(handles.currenttimeidx-1).harmonic3(1:end,1),...
        handles.main.data(handles.currenttimeidx).harmonic3(1:end,2),...
        handles.main.data(handles.currenttimeidx).harmonic3(1:end,3),...
        handles, hObject);
end

if get(handles.fitcond,'value')==1;
    handles.main.data(handles.currenttimeidx).harmonic3(1:end,4) = fittedfreqdata;
elseif get(handles.fitsus,'value')==1;
    handles.main.data(handles.currenttimeidx).harmonic3(1:end,5) = fittedfreqdata;
end

handles.main.delg{3}(handles.currenttimeidx) = HWHM- handles.main.ref(2,2);
handles.main.delf{3}(handles.currenttimeidx) = peaklocation- handles.main.ref(1,2);

handles.gmax{3}(handles.currenttimeidx) = Gmax;
handles.peaklocation{3}(handles.currenttimeidx) = peaklocation;

guidata(hObject,handles)
if handles.noplotupdate ==0;
    plotspectra(hObject,eventdata, handles);
    plotshifts(hObject, eventdata, handles);
    
else
end


% --- Executes on button press in refitn5.
function handles = refitn5_Callback(hObject, eventdata, handles)

if get(handles.useprevn5,'value')==0;
    [peaklocation, HWHM, Gmax, fittedfreqdata] = fitspectra(handles.main.data(handles.currenttimeidx).harmonic5(1:end,1),...
        handles.main.data(handles.currenttimeidx).harmonic5(1:end,2),...
        handles.main.data(handles.currenttimeidx).harmonic5(1:end,3),...
        handles, hObject);
else
    [peaklocation, HWHM, Gmax, fittedfreqdata] = fitspectra(handles.main.data(handles.currenttimeidx-1).harmonic5(1:end,1),...
        handles.main.data(handles.currenttimeidx).harmonic5(1:end,2),...
        handles.main.data(handles.currenttimeidx).harmonic5(1:end,3),...
        handles, hObject);
end

if get(handles.fitcond,'value')==1;
    handles.main.data(handles.currenttimeidx).harmonic5(1:end,4) = fittedfreqdata;
elseif get(handles.fitsus,'value')==1;
    handles.main.data(handles.currenttimeidx).harmonic5(1:end,5) = fittedfreqdata;
end

handles.main.delg{5}(handles.currenttimeidx) = HWHM- handles.main.ref(2,3);
handles.main.delf{5}(handles.currenttimeidx) = peaklocation- handles.main.ref(1,3);

handles.gmax{5}(handles.currenttimeidx) = Gmax;
handles.peaklocation{5}(handles.currenttimeidx) = peaklocation;

guidata(hObject,handles)

if handles.noplotupdate ==0;
    plotspectra(hObject,eventdata, handles);
    plotshifts(hObject, eventdata, handles);
else
end


% --- Executes on button press in showpeak.
function showpeak_Callback(hObject, eventdata, handles)

handles = guidata(hObject);
handles.showpeak_status = get(handles.showpeak,'value');
guidata(hObject, handles);


% --- Executes on button press in remove.
function remove_Callback(hObject, eventdata, handles)

handles.main.t(handles.currenttimeidx) = [];
handles.main.delf{1}(handles.currenttimeidx) = [];
handles.main.delf{3}(handles.currenttimeidx) = [];
handles.main.delf{5}(handles.currenttimeidx) = [];

handles.main.delg{1}(handles.currenttimeidx) = [];
handles.main.delg{3}(handles.currenttimeidx) = [];
handles.main.delg{5}(handles.currenttimeidx) = [];


if handles.currenttimeidx < length(handles.main.t)
    
    handles.main.data(handles.currenttimeidx) = [];
    
    nextpoint_Callback(hObject, eventdata, handles)
    guidata(hObject, handles);
    
%     plotspectra(hObject, eventdata, handles)
%     plotshifts(hObject, eventdata, handles);
    
elseif handles.currenttimeidx == length(handles.main.t)
    
    handles.main.data(handles.currenttimeidx) = [];
    
    previouspoint_Callback(hObject, eventdata, handles)
    guidata(hObject, handles);
%     plotspectra(hObject, eventdata, handles)
%     plotshifts(hObject, eventdata, handles);
end


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datacursormode off
uiresume(handles.figure1)



% --- Executes on button press in refitall.
function refitall_Callback(hObject, eventdata, handles)
% hObject    handle to refitall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.useprevn1,'value')==1 || get(handles.useprevn3,'value')==1;
    set(handles.statusupdate, 'String', 'Cannot use previous values to refit all data!','Foregroundcolor','red');
else
    
    set(handles.statusupdate,'String','Refitting...','Foregroundcolor',[0 0.5 0])
    pause(.1)
    handles.noplotupdate =1;

    for i= 1:1:numel(handles.main.t)
        handles.currenttimeidx =i;

        handles =refitn1_Callback(hObject, eventdata, handles);
        handles =refitn3_Callback(hObject, eventdata, handles);
        handles =refitn5_Callback(hObject, eventdata, handles);
    end
end

handles.noplotupdate =0;
    set(handles.statusupdate,'String','Refitting complete!','Foregroundcolor',[0 0.5 0])

guidata(hObject,handles)
plotspectra(hObject,eventdata, handles);
plotshifts(hObject, eventdata, handles);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
freq_shift(:,1) = handles.main.t;
freq_shift(:,2) = handles.main.delf{1};
freq_shift(:,3) = handles.main.delg{1};
freq_shift(:,4) = handles.main.delf{3};
freq_shift(:,5) = handles.main.delg{3};
freq_shift(:,6) = handles.main.delf{5};
freq_shift(:,7) = handles.main.delg{5};

a=length(handles.main.t) +1;
freq_shift(a,1) = NaN;
freq_shift(a,2) = NaN;
freq_shift(a,3) = NaN;
freq_shift(a,4) = NaN;
freq_shift(a,5) = NaN;
freq_shift(a,6) = NaN;
freq_shift(a,7) = NaN;

freq_shift_ref(1,1) = handles.main.ref(1,1);
freq_shift_ref(2,1) = handles.main.ref(2,1);
freq_shift_ref(1,2) = handles.main.ref(1,2);
freq_shift_ref(2,2) = handles.main.ref(2,2);
freq_shift_ref(1,3) = handles.main.ref(1,3);
freq_shift_ref(2,3) = handles.main.ref(2,3);

save([handles.main.pathname handles.main.filename '_refitted.mat'], 'freq_shift', 'freq_shift_ref');
% save([handles.main.pathname 'refitted_spectras.mat'],'handles.main.data');
set(handles.statusupdate, 'String', 'Refitted data has been saved.','Foregroundcolor',[0 0.5 0]);


function [index, peak_detect, locs, width] = findrelevantpeaks(freq,conductance)

smoothcond = smooth(conductance,30);

%Note: peak_deteck and pks are be the same.
[peak_detect,index]=findpeaks(smoothcond,'sortstr','descend');
[pks, locs, width, prominence] = findpeaks(smoothcond,freq,'sortstr','descend'); % Initial guesses

% figure;
% findpeaks(smoothcond)

if isempty(peak_detect)
    disp('No peaks were found!')
    set(handles.statusupdate, 'String', 'No peaks were found!','Foregroundcolor','red');
    peak_detect = [];
    index = [];
    numpeaks = 0;
    return
end

% % Ignore weird peaks identified at the edges
% indices = find((abs(index)>0.80*length(freq)) && (abs(index)<0.15*length(freq)));
% index(indices) = [];
% width(indices) = [];
% locs(indices) = [];
% peak_detect(indices) = [];

% Index of peaks greater than a third of the size of the max peak
peakstofit = find(prominence>prominence(1)./3);
index = index(peakstofit);
locs = locs(peakstofit);
width = width(peakstofit);
peak_detect = peak_detect(peakstofit);

% Now sometimes there are two peaks indentified VERY close to eachother
% This next bit checks for that and picks the ones with idx atleast 20 away
% from main peak
if length(index)>1
    a = find((abs(index(2:end)-index(1))<20)) +1;
    %     index = vertcat(index(1),index(a)');
    index(a) = [];
    peak_detect(a) = [];
    locs(a) = [];
    width(a) = [];
else
end

% Warning if there are more than 3 peaks identified
if length(index) > 3
    set(handles.statusupdate, 'String', 'More than 3 peaks found!','Foregroundcolor','red');
    return
end

% --- Executes on button press in useprevn3.
function useprevn3_Callback(hObject, eventdata, handles)
% hObject    handle to useprevn3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useprevn3


% --- Executes on button press in useprevn1.
function useprevn1_Callback(hObject, eventdata, handles)
% hObject    handle to useprevn1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useprevn1



% function output_txt = myupdatefcn(~,event_obj, handles)
% This is the function that runs when datacursormode is employed. The
% output output-txt is what appears in the box.
%
% Determines output box--this is the usual function of datacursor, modified
% to know what the x axis actually is.
% pos = get(event_obj,'Position');
% output_txt = {['Time: ', num2str(round(pos(1),2)), ' min'],...
%     ['y: ',num2str(round(pos(2),2))]};
%
% handles.currenttimeidx = find(handles.main.t ==pos(1))
%
% plotspectra(hObject,eventdata, handles);
% plotshifts(hObject, eventdata, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg('Do you want to save refitted data & transfer to 2LayerGUI?',...
    'Close Request Function',...
    'Yes','No','Cancel','Yes');
switch selection,
    case 'Yes',
        close_Callback(hObject, eventdata, handles)
    case 'No',
        % Do Nothing
        delete(hObject)
    case 'Cancel'
        % Do Nothing
end



function fitfactor_Callback(hObject, eventdata, handles)
% hObject    handle to fitfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fitfactor as text
%        str2double(get(hObject,'String')) returns contents of fitfactor as a double


% --- Executes during object creation, after setting all properties.
function fitfactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useprevn5.
function useprevn5_Callback(hObject, eventdata, handles)
% hObject    handle to useprevn5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useprevn5

function F_both = lfun4cs(p,x)
F_both = [lfun4c(p,x) lfun4s(p,x)];


function F_cond = lfun4c(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance

F_cond= p(4).*((((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*cosd(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*sind(p(3)))+p(5);

function F_sus = lfun4s(p,x)
%Order of parameters in p is as follows
%p(1): f0 maximum frequency
%p(2): gamma0 dissipation
%p(3): phi phse angle difference
%p(4): Gmax maximum conductance
F_sus = -p(4).*(-(((x.^2).*((2.*p(2)).^2))./(((((p(1)).^2)-(x.^2)).^2)+...
    ((x.^2).*((2.*p(2)).^2)))).*sind(p(3))-((((p(1)).^2-x.^2)).*x.*(2.*p(2)))./...
    (((((p(1)).^2)-(x.^2)).^2)+((x.^2).*((2.*p(2)).^2))).*cosd(p(3)))+p(6);
