function varargout = lamby(varargin)
% LAMBY MATLAB code for lamby.fig
%      LAMBY, by itself, creates a new LAMBY or raises the existing
%      singleton*.
%
%      H = LAMBY returns the handle to a new LAMBY or the handle to
%      the existing singleton*.
%
%      LAMBY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAMBY.M with the given input arguments.
%
%      LAMBY('Property','Value',...) creates a new LAMBY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lamby_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lamby_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lamby

% Last Modified by GUIDE v2.5 15-Feb-2014 23:04:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lamby_OpeningFcn, ...
                   'gui_OutputFcn',  @lamby_OutputFcn, ...
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


% --- Executes just before lamby is made visible.
function lamby_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lamby (see VARARGIN)

% Choose default command line output for lamby
handles.output = hObject;

% Initialize the arrays that hold the dispersion curves values
handles.plate = Lamby.IsotropicPlate();
handles.dispersion.cpa = [];
handles.dispersion.cps = [];
handles.dispersion.cga = [];
handles.dispersion.cgs = [];
handles.dispersion.cpaTol = [];
handles.dispersion.cpsTol = [];
handles.dispersion.cgaTol = [];
handles.dispersion.cgsTol = [];
handles.dispersion.fcpa = [];
handles.dispersion.fcps = [];
handles.dispersion.fcga = [];
handles.dispersion.fcgs = [];
handles.dispersion.f = [];

% Initialize the arrays the hold the imported data
handles.imports = [];

% initialize the properties of the PZT 
handles.mypzt = Lamby.MonolithicPztActuator();
handles.mypzt.plate = handles.plate;

% initialize the response variables
handles.excitation.t = [];
handles.excitation.signal = [];
handles.excitation.f = [];
handles.excitation.signalSpectrum = [];

handles.spectrum.f = [];
handles.spectrum.Hs_strain = [];
handles.spectrum.Hs_disp = [];
handles.spectrum.fsInd = [];
handles.spectrum.Ha_strain = [];
handles.spectrum.Ha_disp = [];
handles.spectrum.faInd = [];

handles.response.Hs_strain = [];
handles.response.Hs_disp = [];
handles.response.fsInd = [];
handles.response.Ha_strain = [];
handles.response.Ha_disp = [];
handles.response.faInd = [];

handles.response.us = [];
handles.response.ua = [];
handles.response.es = [];
handles.response.ea = [];

handles.currentFig = [];
handles.nplot = 0;

% initialize the properties of the plate into the GUI fields
set(handles.edt_plateDensity, 'String', num2str(handles.plate.density))
set(handles.edt_plateThickness, 'String', num2str(handles.plate.thickness))
set(handles.edt_property1, 'String', num2str(handles.plate.Y/1e9));
set(handles.edt_property2, 'String', num2str(handles.plate.nu));

% pzt properties GUI fields
set(handles.edt_pztLength, 'String', num2str(handles.mypzt.length));
set(handles.edt_pztThickness, 'String', num2str(handles.mypzt.thickness));
set(handles.edt_pztY11, 'String', num2str(1e-9/handles.mypzt.S(1,1)));
set(handles.edt_pztd31, 'String', num2str(1e12*handles.mypzt.d(3,1)));
set(handles.edt_pztVoltage, 'String', num2str(handles.mypzt.V));
set(handles.edt_bondThickness, 'String', num2str(handles.mypzt.bondThickness*1e6));
set(handles.edt_bondG, 'String', num2str(handles.mypzt.bondShearModulus*1e-9));
    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lamby wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lamby_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_dispersionShow.
function pb_dispersionShow_Callback(hObject, eventdata, handles)
% hObject    handle to pb_dispersionShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f0 = str2double(get(handles.edt_dispFreqStart, 'String'));
fstep = str2double(get(handles.edt_dispFreqStep, 'String'));
fend = str2double(get(handles.edt_dispFreqEnd, 'String'));

f = (f0:fstep:fend)*1e3;
isNewRange = ~isequal(handles.dispersion.f, f);
handles.dispersion.f = f;    

% compute the values for dispersion velocities
tol = str2double(get(handles.edt_symPhaseTol, 'String'));
tolchanged = ~isequal(handles.dispersion.cpsTol, tol);
if get(handles.cb_dispSymmPhase, 'Value') && ...
        (isNewRange || isempty(handles.dispersion.cps) || tolchanged)

    handles.dispersion.cpsTol = tol;
    [handles.dispersion.cps, handles.dispersion.fcps] = ...
        handles.plate.phaseVelocity(f, 's', tol);
end

tol = str2double(get(handles.edt_asymPhaseTol, 'String'));
tolchanged = ~isequal(handles.dispersion.cpaTol, tol);
if get(handles.cb_dispAsymPhase, 'Value') && ...
        (isNewRange || isempty(handles.dispersion.cpa) || tolchanged)

    handles.dispersion.cpaTol = tol;
    [handles.dispersion.cpa, handles.dispersion.fcpa] = ...
        handles.plate.phaseVelocity(f, 'a', tol);
end    

tol = str2double(get(handles.edt_symGroupTol, 'String'));
tolchanged = ~isequal(handles.dispersion.cgsTol, tol);
if get(handles.cb_dispSymmGroup, 'Value') && ...
        (isNewRange || isempty(handles.dispersion.cgs) || tolchanged)

    handles.dispersion.cgsTol = tol;
    [handles.dispersion.cgs, handles.dispersion.fcgs] = ...
        handles.plate.groupVelocity(f, 's', tol);
end

tol = str2double(get(handles.edt_asymGroupTol, 'String'));
tolchanged = ~isequal(handles.dispersion.cgaTol, tol);
if get(handles.cb_dispAsymGroup, 'Value') && ...
        (isNewRange || isempty(handles.dispersion.cga) || tolchanged)

    handles.dispersion.cgaTol = tol;
    [handles.dispersion.cga, handles.dispersion.fcga] = ...
        handles.plate.groupVelocity(f, 'a', tol);
end

% if range is scalar, just show values in the messages box
if isscalar(f)
    set(handles.txt_messages, 'String', ''); 
    str = {};
    if get(handles.cb_dispSymmPhase, 'Value')
        for i=1:length(handles.dispersion.cps)
            str = [{}; get(handles.txt_messages, 'String')];
            set(handles.txt_messages, 'String', [str; ['Phase Velocity S' num2str(i-1) ...
                ' = ' num2str(handles.dispersion.cps{i}) ' m/s']]);
        end
    end    

    if get(handles.cb_dispAsymPhase, 'Value')
        for i=1:length(handles.dispersion.cpa)
            str = [{}; get(handles.txt_messages, 'String')];
            set(handles.txt_messages, 'String', [str; ['Phase Velocity A' num2str(i-1) ...
                ' = ' num2str(handles.dispersion.cpa{i}) ' m/s']]);
        end
    end
    
    if get(handles.cb_dispSymmGroup, 'Value')
        for i=1:length(handles.dispersion.cgs)
            str = [{}; get(handles.txt_messages, 'String')];
            set(handles.txt_messages, 'String', [str; ['Group Velocity S' num2str(i-1) ...
                ' = ' num2str(handles.dispersion.cgs{i}) ' m/s']]);
        end
    end
    
    if get(handles.cb_dispAsymGroup, 'Value')
        for i=1:length(handles.dispersion.cga)
            str = [{} get(handles.txt_messages, 'String')];
            set(handles.txt_messages, 'String', [str; ['Group Velocity A' num2str(i-1) ...
                ' = ' num2str(handles.dispersion.cga{i}) ' m/s']]);
        end
    end
else % if frequency range is a vector, make a plot
    figure;
    if get(handles.cb_dispSymmPhase, 'Value')
        props.prefix = 'S';       
        props.spec = '-ro';
        Lamby.drawVelocityDispersion(handles.dispersion.cps, ...
            handles.dispersion.fcps, f, props);
    end
    if get(handles.cb_dispAsymPhase, 'Value')
        props.prefix = 'A';        
        props.spec = '-bs';
        Lamby.drawVelocityDispersion(handles.dispersion.cpa, ...
            handles.dispersion.fcpa, f, props);
    end
    
    if get(handles.cb_dispSymmGroup, 'Value')
        props.prefix = 'S';        
        props.spec = '-ro';
        Lamby.drawVelocityDispersion(handles.dispersion.cgs, ...
            handles.dispersion.fcgs, f, props);
    end
    if get(handles.cb_dispAsymGroup, 'Value')
        props.prefix = 'A';
        props.spec = '-bs';
        Lamby.drawVelocityDispersion(handles.dispersion.cga, ...
            handles.dispersion.fcga, f, props);
    end    
end

guidata(hObject, handles);


% --- Executes on selection change in pm_propertiesBy.
function pm_propertiesBy_Callback(hObject, eventdata, handles)
% hObject    handle to pm_propertiesBy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_propertiesBy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_propertiesBy

sel = get(hObject, 'Value');

if sel == 1
    set(handles.txt_property1, 'String', 'Y (GPa)'); 
    set(handles.edt_property1, 'String', num2str(handles.plate.Y/1e9));
    set(handles.txt_property2, 'String', 'nu');
    set(handles.edt_property2, 'String', num2str(handles.plate.nu));
elseif sel == 2
    set(handles.txt_property1, 'String', 'cl (m/s)'); 
    set(handles.edt_property1, 'String', num2str(handles.plate.cl));
    set(handles.txt_property2, 'String', 'ct (m/s)');
    set(handles.edt_property2, 'String', num2str(handles.plate.ct));
elseif sel == 3
    set(handles.txt_property1, 'String', 'mu (GPa)'); 
    set(handles.edt_property1, 'String', num2str(handles.plate.mu/1e9));
    set(handles.txt_property2, 'String', 'lambda (GPa)');
    set(handles.edt_property2, 'String', num2str(handles.plate.lambda/1e9));
end


function edt_plateDensity_Callback(hObject, eventdata, handles)
% hObject    handle to edt_plateDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

density = str2double(get(hObject,'String'));
if ~isPositiveNumber(density)
    set(hObject, 'String', num2str(handles.plate.density));
else
    handles.plate.density = density;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end    


function edt_plateThickness_Callback(hObject, eventdata, handles)
% hObject    handle to edt_plateThickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

thickness = str2double(get(hObject,'String'));
if ~isPositiveNumber(thickness)
    set(hObject, 'String', num2str(handles.plate.thickness));
else
    handles.plate.thickness = thickness;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end


function edt_property1_Callback(hObject, eventdata, handles)
% hObject    handle to edt_property1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prop1 = str2double(get(hObject,'String'));
sel =  get(handles.pm_propertiesBy, 'Value');

if ~isPositiveNumber(prop1)
    if sel == 1
        set(hObject, 'String', num2str(handles.plate.Y/1e9));
    elseif sel == 2
        set(hObject, 'String', num2str(handles.plate.cl));
    elseif sel == 3
        set(hObject, 'String', num2str(handles.plate.mu/1e9));
    end
else
    if sel == 1
        handles.plate.Y = prop1*1e9;
    elseif sel == 2
        handles.plate.cl = prop1;
    elseif sel == 3
        handles.plate.mu = prop1*1e9;
    end
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end


function edt_property2_Callback(hObject, eventdata, handles)
% hObject    handle to edt_property2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prop2 = str2double(get(hObject,'String'));
sel =  get(handles.pm_propertiesBy, 'Value');

if ~isPositiveNumber(prop2)
    if sel == 1
        set(hObject, 'String', num2str(handles.plate.nu));
    elseif sel == 2
        set(hObject, 'String', num2str(handles.plate.ct));
    elseif sel == 3
        set(hObject, 'String', num2str(handles.plate.lambda/1e9));
    end
else
    if sel == 1
        handles.plate.nu = prop2;
    elseif sel == 2
        handles.plate.ct = prop2;
    elseif sel == 3
        handles.plate.lambda = prop2*1e9;
    end
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end


% --- Executes on button press in pb_modeShapeShow.
function pb_modeShapeShow_Callback(hObject, eventdata, handles)
% hObject    handle to pb_modeShapeShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = str2double(get(handles.edt_modeShapeFreq, 'String'));

if ~isPositiveNumber(f) 
    set(handles.txt_messages, 'String', 'Mode shape frequency should be a positive number.');
    return;
end

if get(handles.cb_symmModeShape, 'Value')
    [Usx Usy zs] = handles.plate.symmetricModeshape(f);
    
    for i=1:length(Usx)
        figure;plot(zs, Usx{i}, zs, Usy{i});        
    end
end

if get(handles.cb_asymModeShape, 'Value')
    [Uax Uay za] = handles.plate.antisymmetricModeshape(f);        
    for i=1:length(Uax)
        figure;plot(za, Uax{i}, za, Uay{i});        
    end
end



function edt_pztLength_Callback(hObject, eventdata, handles)
% hObject    handle to edt_pztLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
len = str2double(get(hObject,'String'));
if ~isPositiveNumber(len)
    set(hObject, 'String', num2str(handles.mypzt.length));
else
    handles.mypzt.length = len;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles)
end    


function edt_pztY11_Callback(hObject, eventdata, handles)
% hObject    handle to edt_pztY11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

y11 = str2double(get(hObject,'String'));
if ~isPositiveNumber(y11)
    set(hObject, 'String', num2str(1e-9/handles.mypzt.S(1,1)));
else
    handles.mypzt.S(1,1) = y11*1e9;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles)
end  

function edt_pztd31_Callback(hObject, eventdata, handles)
% hObject    handle to edt_pztd31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d31 = str2double(get(hObject,'String'));
if ~isscalar(d31) || ~isnumeric(d31)
    set(hObject, 'String', num2str(1e12*handles.mypzt.d(3,1)));
else
    handles.mypzt.d(3,1) = d31*1e-12;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end  

function edt_pztVoltage_Callback(hObject, eventdata, handles)
% hObject    handle to edt_pztVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = str2double(get(hObject,'String'));
if ~isscalar(v) || ~isnumeric(v)
    set(hObject, 'String', num2str(handles.mypzt.V));
else
    handles.mypzt.V = v;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end  


function edt_bondThickness_Callback(hObject, eventdata, handles)
% hObject    handle to edt_bondThickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

thick = str2double(get(hObject,'String'));
if ~isPositiveNumber(thick)
    set(hObject, 'String', num2str(handles.mypzt.bondThickness*1e6));
else
    handles.mypzt.bondThickness = thick*1e-6;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end


function edt_bondG_Callback(hObject, eventdata, handles)
% hObject    handle to edt_bondG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

G = str2double(get(hObject,'String'));
if ~isPositiveNumber(G)
    set(hObject, 'String', num2str(1e-9*handles.mypzt.bondShearModulus));
else
    handles.mypzt.bondShearModulus = 1e9*G;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end  

function edt_pztThickness_Callback(hObject, eventdata, handles)
% hObject    handle to edt_pztThickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

thick = str2double(get(hObject,'String'));
if ~isPositiveNumber(thick)
    set(hObject, 'String', num2str(handles.mypzt.thickness));
else
    handles.mypzt.thickness = thick;
    guidata(hObject, handles);
    resetResponseVariables(hObject, handles);
end  


function edt_excitFreqStart_Callback(hObject, eventdata, handles)
% hObject    handle to edt_excitFreqStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prevnum = get(hObject, 'UserData');
num = str2double(get(hObject,'String'));
numEnd = str2double(get(handles.edt_excitFreqEnd, 'String'));

if ~isPositiveNumber(num) || (num > numEnd)
    set(hObject, 'String', num2str(prevnum));
else
    set(hObject, 'UserData', num);
    if get(handles.cb_spectrumShowPlot, 'Value') && prevnum ~= num
        updateSpectrumModes(hObject, handles);
    end
end

function edt_excitFreqStep_Callback(hObject, eventdata, handles)
% hObject    handle to edt_excitFreqStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prevnum = get(hObject, 'UserData');
num = str2double(get(hObject,'String'));

if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    set(hObject, 'UserData', num);
    if get(handles.cb_spectrumShowPlot, 'Value') && prevnum ~= num
        updateSpectrumModes(hObject, handles);
    end
end

function edt_excitFreqEnd_Callback(hObject, eventdata, handles)
% hObject    handle to edt_excitFreqEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prevnum = get(hObject, 'UserData');
num = str2double(get(hObject,'String'));
numStart = str2double(get(handles.edt_excitFreqStart, 'String'));

if ~isPositiveNumber(num) || (numStart > num)
    set(hObject, 'String', num2str(prevnum));
else
    set(hObject, 'UserData', num);
    if get(handles.cb_spectrumShowPlot, 'Value') && prevnum ~= num
        updateSpectrumModes(hObject, handles);
    end
end


function edt_gaussExcitT_Callback(hObject, eventdata, handles)
% hObject    handle to edt_gaussExcitT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prevnum = get(hObject, 'UserData');
num = str2double(get(hObject,'String'));

if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    set(hObject, 'UserData', num);
    if get(handles.cb_gaussShowPlot, 'Value') && prevnum ~= num
        updateGaussModes(hObject, handles);
    end
end


function edt_gaussExcitF0_Callback(hObject, eventdata, handles)
% hObject    handle to edt_gaussExcitF0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prevnum = get(hObject, 'UserData');
num = str2double(get(hObject,'String'));

if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    set(hObject, 'UserData', num);
    if get(handles.cb_gaussShowPlot, 'Value') && prevnum ~= num
        updateGaussModes(hObject, handles);
    end
end


function edt_gaussExcitBW_Callback(hObject, eventdata, handles)
% hObject    handle to edt_gaussExcitBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prevnum = get(hObject, 'UserData');
num = str2double(get(hObject,'String'));

if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    set(hObject, 'UserData', num);
    if get(handles.cb_gaussShowPlot, 'Value') && prevnum ~= num
        updateGaussModes(hObject, handles);
    end
end


% --- Executes on button press in pb_importData.
function pb_importData_Callback(hObject, eventdata, handles)
% hObject    handle to pb_importData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName, PathName, FilterIndex] = uigetfile('*','Select data file(s)', 'MultiSelect', 'on');

str = get(handles.lb_importedData, 'String');
if isempty(str)
   str = {};
end

if FilterIndex
    set(handles.txt_messages, 'String', 'Importing data...');
    pause(0.1);
    userDat = get(hObject, 'UserData');
    if ~iscell(FileName)
        FileName = {FileName};
    end
    for i = 1:length(FileName)
        dat = importdata([PathName '/' FileName{i}], '\t', 1);
        
        % intialize the variables to hold the imported data
        handles.imports(end+1).t = dat.data(:,1);
        handles.imports(end).data = dat.data(:,2:end);
        ndat = length(handles.imports(end).data(1,:));
        handles.imports(end).startInd = ones(ndat, 1);
        handles.imports(end).endInd = length(handles.imports(end).data(:,1))*ones(ndat, 1);
        handles.imports(end).f = [];
        handles.imports(end).spectrum = cell(ndat, 1);
        
        fname = strsplit(FileName{i},'.');
        fname = fname{1};
        handles.imports(end).name = fname;
        for j = 1:size(handles.imports(end).data,2)            
            str = [str; [fname '_' num2str(j)]];
            set(handles.lb_importedData, 'String', str);
        end
        userDat(end+1) = size(handles.imports(end).data,2);
    end
    set(hObject, 'UserData', userDat);
end
set(handles.txt_messages, 'String', 'Data imported successfully.');
guidata(hObject, handles);


% --- Executes on selection change in lb_importedData.
function lb_importedData_Callback(hObject, eventdata, handles)
% hObject    handle to lb_importedData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sel = get(hObject, 'Value');

if isscalar(sel)
    [ind, datInd] = getSelectionIndex(handles, sel);
    
    set(handles.edt_startIndex, 'Enable', 'on');
    set(handles.edt_startTime, 'Enable', 'on');
    set(handles.edt_endIndex, 'Enable', 'on');
    set(handles.edt_endTime, 'Enable', 'on');

    startInd = handles.imports(ind).startInd(datInd);
    endInd = handles.imports(ind).endInd(datInd);
    set(handles.edt_startIndex, 'String', num2str(startInd));
    set(handles.edt_endIndex, 'String', num2str(endInd));
    set(handles.edt_startTime, 'String', num2str(1e6*handles.imports(ind).t(startInd)));
    set(handles.edt_endTime, 'String', num2str(1e6*handles.imports(ind).t(endInd)));
else
    set(handles.edt_startIndex, 'Enable', 'off');
    set(handles.edt_startTime, 'Enable', 'off');
    set(handles.edt_endIndex, 'Enable', 'off');
    set(handles.edt_endTime, 'Enable', 'off');
end


function edt_startTime_Callback(hObject, eventdata, handles)
% hObject    handle to edt_startTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = 1e-6*str2double(get(hObject, 'String'));
endInd = str2double(get(handles.edt_endIndex, 'String'));
sel = get(handles.lb_importedData, 'Value');

if isscalar(sel)
    [ind, datInd] = getSelectionIndex(handles, sel);
    [~, startInd] = min(abs(handles.imports(ind).t - num));        
    if isPositiveNumber(num) && startInd < endInd
        handles.imports(ind).startInd(datInd) = startInd;
    else
        startInd = handles.imports(ind).startInd(datInd);
    end
    set(hObject, 'String', num2str(1e6*handles.imports(ind).t(startInd)));
    set(handles.edt_startIndex, 'String', num2str(startInd));
end
guidata(hObject, handles);


function edt_startIndex_Callback(hObject, eventdata, handles)
% hObject    handle to edt_startIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(hObject, 'String'));
endIndex = str2double(get(handles.edt_endIndex, 'String'));
sel = get(handles.lb_importedData, 'Value');

if isscalar(sel)
    [ind, datInd] = getSelectionIndex(handles, sel);
    if isPositiveNumber(num) && num < endIndex
        handles.imports(ind).startInd(datInd) = num;
    else
        num = handles.imports(ind).startInd(datInd);
        set(hObject, 'String', num2str(num));
    end

    set(handles.edt_startTime, 'String', num2str(1e6*handles.imports(ind).t(num)));
end
guidata(hObject, handles);



function edt_endTime_Callback(hObject, eventdata, handles)
% hObject    handle to edt_endTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = 1e-6*str2double(get(hObject, 'String'));
startInd = str2double(get(handles.edt_startIndex, 'String'));
sel = get(handles.lb_importedData, 'Value');

if isscalar(sel)
    [ind, datInd] = getSelectionIndex(handles, sel);
    maxInd = length(handles.imports(ind).t);
    [~, endInd] = min(abs(handles.imports(ind).t - num));
    
    if isPositiveNumber(num) && endInd > startInd && endInd <= maxInd
        handles.imports(ind).endInd(datInd) = endInd;
    else
        endInd = handles.imports(ind).endInd(datInd);
    end
    set(hObject, 'String', num2str(1e6*handles.imports(ind).t(endInd)));
    set(handles.edt_endIndex, 'String', num2str(endInd));
end
guidata(hObject, handles);


function edt_endIndex_Callback(hObject, eventdata, handles)
% hObject    handle to edt_endIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(hObject, 'String'));
startIndex = str2double(get(handles.edt_startIndex, 'String'));
sel = get(handles.lb_importedData, 'Value');

if isscalar(sel)
    [ind, datInd] = getSelectionIndex(handles, sel);
    if isPositiveNumber(num) && num > startIndex
        handles.imports(ind).endInd(datInd) = num;
    else
        num = handles.imports(ind).endInd(datInd);
        set(hObject, 'String', num2str(num));
    end

    set(handles.edt_endTime, 'String', num2str(1e6*handles.imports(ind).t(num)));
end
guidata(hObject, handles);


% --- Executes on button press in pb_resetAll.
function pb_resetAll_Callback(hObject, eventdata, handles)
% hObject    handle to pb_resetAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nfiles = length(handles.imports);
set(handles.lb_importedData, 'Value', 1);
for i=1:nfiles
    ndat = length(handles.imports(i).startInd);
    for j=1:ndat
        if i==1 && j == 1
            set(handles.edt_startIndex, 'String', num2str(1));
            set(handles.edt_startTime, 'String', num2str(1e6*handles.imports(1).t(1)));
            set(handles.edt_endIndex, 'String', num2str(length(handles.imports(1).t)));
            set(handles.edt_endTime, 'String', num2str(1e6*handles.imports(1).t(end)));
        end
        handles.imports(i).startInd(j) = 1;
        handles.imports(i).endInd(j) = length(handles.imports(i).t);
    end
end
guidata(hObject, handles);


% --- Executes on button press in pb_applyToFile.
function pb_applyToFile_Callback(hObject, eventdata, handles)
% hObject    handle to pb_applyToFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sel = get(handles.lb_importedData, 'Value');
startInd = str2double(get(handles.edt_startIndex, 'String'));
endInd = str2double(get(handles.edt_endIndex, 'String'));
if isscalar(sel) 
    [ind, datInd] = getSelectionIndex(handles, sel);
    
    for i=1:length(handles.imports(ind).startInd)
        handles.imports(ind).startInd(i) = startInd;
        handles.imports(ind).endInd(i) = endInd;        
    end 
end
guidata(hObject, handles);


% --- Executes on button press in pb_plotResponse.
function pb_plotResponse_Callback(hObject, eventdata, handles)
% hObject    handle to pb_plotResponse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
normalize = get(handles.cb_normalizeData, 'Value');
strainDisp = get(handles.rb_strain, 'Value');
freqplot = get(handles.rb_plotFrequency, 'Value');
showLgnd = get(handles.cb_showLegend, 'Value');
addModes = get(handles.cb_addModes, 'Value');
holdplot = get(handles.cb_holdplot, 'Value');

% spectrum plots 
plotSpectrum = get(handles.cb_spectrumShowPlot, 'Value');
if plotSpectrum
    selModes = get(handles.lb_spectrumModes, 'Value');
    nmodes = get(handles.lb_spectrumModes, 'UserData');
    selSModes = selModes(selModes <= nmodes(1));
    selAModes = selModes(selModes > nmodes(1))-nmodes(1);
end

% gaussian pulse plots 
plotAnalytical = get(handles.cb_gaussShowPlot, 'Value');
if plotAnalytical
    gaussSelModes = get(handles.lb_gaussModes, 'Value');
    gaussnmodes = get(handles.lb_gaussModes, 'UserData');
    gaussSModes = gaussSelModes(gaussSelModes <= gaussnmodes(1));
    gaussAModes = gaussSelModes(gaussSelModes > gaussnmodes(1))-gaussnmodes(1);
end

% imports plots
selImport = get(handles.lb_importedData, 'Value');

% calculate the total number of plots for the figure
totplots = 0;
if plotSpectrum && freqplot 
    totplots = totplots + length(selSModes)+length(selAModes);
end
if plotAnalytical 
    totplots = totplots + length(gaussSModes) + length(gaussAModes);
end

totplots = totplots+length(selImport);

lnspecs = {'-o', '-d', '-s', '--o', '--d', '--s', '-.o', '-.d', '-.s'};
lnspecsT = {'-', '--', '-.', '-.', '-..', '--.'};
lgnd = {};

if holdplot 
    if ~isempty(handles.currentFig)
        figure(handles.currentFig);
    end
end

if ~holdplot  || isempty(handles.currentFig)    
    handles.nplot = 0;
    handles.currentFig = figure('Position', [2200 300 1300 800]);
end

box off;
grid off;
hold on;
colrs = lines(handles.nplot+totplots);

set(handles.txt_messages, 'String', 'Calculating response . . .');
pause(0.1);

% ------------------------------------------------------------------------%
% compute and plot in the frequency domain -------------------------------%
if freqplot
    % calculate the spectrum 
    if plotSpectrum
        fstart = 1e3*str2double(get(handles.edt_excitFreqStart, 'String'));
        fstep = 1e3*str2double(get(handles.edt_excitFreqStep, 'String'));
        fend = 1e3*str2double(get(handles.edt_excitFreqEnd, 'String'));
        f = fstart:fstep:fend;
        % symmetric modes 
        if (~isequal(f, handles.spectrum.f) || isempty(handles.spectrum.Hs_strain)) && ....
                ~isempty(selSModes)        
            [handles.spectrum.Hs_strain, handles.spectrum.Hs_disp] = ...
                handles.mypzt.frequencyResponse(f, 's');
            
            if isempty(selAModes)
                handles.spectrum.Ha_strain = [];
                handles.spectrum.Ha_disp = [];
                handles.spectrum.faInd = [];            
            end
        end
        % antisymmetric modes
        if (~isequal(f, handles.spectrum.f) || isempty(handles.spectrum.Ha_strain)) && ...
                ~isempty(selAModes)
            [handles.spectrum.Ha_strain, handles.spectrum.Ha_disp] = ...
                handles.mypzt.frequencyResponse(f,'a');
            
            if isempty(selSModes)
                handles.spectrum.Hs_strain = [];
                handles.spectrum.Hs_disp = [];
                handles.spectrum.fsInd = [];            
            end
        end
        handles.spectrum.f = f;
        % make the plots
        for i=1:length(selSModes)
%             f = handles.spectrum.f(handles.spectrum.fsInd{selSModes(i)});
            if strainDisp
                y = abs(handles.spectrum.Hs_strain{selSModes(i)});
            else
                y = abs(handles.spectrum.Hs_disp{selSModes(i)});
            end
            if normalize 
                y = y/max(y);
            end
            plot(f/1e3, y, lnspecs{mod(handles.nplot,length(lnspecs))+1}, ...
                'Color', colrs(handles.nplot+1,:),'MarkerEdgeColor', colrs(handles.nplot+1,:), ...
                'MarkerFaceColor', colrs(handles.nplot+1,:),'MarkerSize', 12, 'LineWidth', 4);                
            lgnd{handles.nplot+1} = ['Spectrum S' num2str(i-1)];
            handles.nplot = handles.nplot+1;            
        end
        for i=1:length(selAModes)
%             f = handles.spectrum.f(handles.spectrum.faInd{selAModes(i)});
            if strainDisp
                y = abs(handles.spectrum.Ha_strain{selAModes(i)});
            else
                y = abs(handles.spectrum.Ha_disp{selAModes(i)});
            end
            if normalize
                y = y/max(y);
            end
            plot(f/1e3, y, lnspecs{mod(handles.nplot,length(lnspecs))+1}, ...
                'Color', colrs(handles.nplot+1,:), 'MarkerEdgeColor', ...
                colrs(handles.nplot+1,:), 'MarkerFaceColor', colrs(handles.nplot+1,:),...
                'MarkerSize', 12, 'LineWidth', 4);
            lgnd{handles.nplot+1} = ['Spectrum A' num2str(i-1)];
            handles.nplot = handles.nplot+1;
        end
    end
    
    % calculate the frequency domain of the response to arbitrary signal
    if plotAnalytical 
        x = str2double(get(handles.edt_sensorPos, 'String'));
        % symmetric modes
        if isempty(handles.response.Hs_strain) && ~isempty(gaussSModes)
            [handles.response.Hs_strain, handles.response.Hs_disp] = ...
                handles.mypzt.frequencyResponse(handles.excitation.f, 's', ...
                x, handles.excitation.signalSpectrum, -30, 0.08);
        end
        
        % antisymmetric modes
        if isempty(handles.response.Ha_strain) && ~isempty(gaussAModes)
            [handles.response.Ha_strain, handles.response.Ha_disp] = ...
                handles.mypzt.frequencyResponse(handles.excitation.f, 'a',...
                x, handles.excitation.signalSpectrum, -30, 0.08);
        end
        f = handles.excitation.f;
        
        f = fftshift(f);
        indleft = find(f<0);
        f(indleft) = [];        
        
        % make the plots
        for i=1:length(gaussSModes)
%             fs = f(handles.response.fsInd{gaussSModes(i)});
            if strainDisp
                y = fftshift(abs(handles.response.Hs_strain{gaussSModes(i)}));
            else
                y = fftshift(abs(handles.response.Hs_disp{gaussSModes(i)}));               
            end
            y(indleft) = [];
            if normalize 
                y = y/max(y);
            end
            plot(f/1e3, y, lnspecs{mod(handles.nplot,length(lnspecs))+1}, 'Color', ...
                colrs(handles.nplot+1,:), 'MarkerEdgeColor', colrs(handles.nplot+1,:), ...
                'MarkerFaceColor', colrs(handles.nplot+1,:),'MarkerSize', 8, 'LineWidth', 4);
            lgnd{handles.nplot+1} = ['Response S' num2str(i-1)];
            handles.nplot = handles.nplot+1;
        end
        for i=1:length(gaussAModes)
%             fa = f(handles.response.faInd{gaussAModes(i)});
            if strainDisp
                y = fftshift(abs(handles.response.Ha_strain{gaussAModes(i)}));
            else
                y = fftshift(abs(handles.response.Ha_disp{gaussAModes(i)}));
            end
            if normalize
                y = y/max(y);
            end
            y(indleft) = [];
            plot(f/1e3, y, lnspecs{mod(handles.nplot,length(lnspecs))+1}, 'Color', ...
                colrs(handles.nplot+1,:), 'MarkerEdgeColor', colrs(handles.nplot+1,:), ...
                'MarkerFaceColor', colrs(handles.nplot+1,:), 'MarkerSize', 8, 'LineWidth', 4);
            lgnd{handles.nplot+1} = ['Response A' num2str(i-1)];
            handles.nplot = handles.nplot+1;
        end        
    end
    
    % calculate the frequency domain of imported data;
    if ~isempty(selImport)        
        for i=1:length(selImport) 
            [ind, datInd] = getSelectionIndex(handles, selImport(i));
            dt = mean(diff(handles.imports(ind).t));
            fs = 1/dt;
            NFFT = length(handles.imports(ind).t);
            df = fs/NFFT;
            handles.imports(ind).f = (0:(NFFT-1))*df;
            handles.imports(ind).f(handles.imports(ind).f >= fs/2) = ...
                handles.imports(ind).f(handles.imports(ind).f >= fs/2) - fs;
            
            % get the truncated signal according to the window start and end
            sig = handles.imports(ind).data(:,datInd);
%             sig = sig - mean(sig);
            sig(1:handles.imports(ind).startInd(datInd)) = 0;
            sig(handles.imports(ind).endInd(datInd):end) = 0;
            
            % calculate the FFT
            handles.imports(ind).spectrum{datInd} = fft(sig, NFFT)/NFFT;           
    
            % Get the one sided FFT
            negf = find(handles.imports(ind).f  < 0);
            handles.imports(ind).f(negf) = [];
            handles.imports(ind).spectrum{datInd}(negf) = [];
            
            % Find the part of spectrum that is within 30 dB the max
            delind = find(abs(handles.imports(ind).spectrum{datInd}) < ...
                    max(abs(handles.imports(ind).spectrum{datInd}))*1e-3);
            handles.imports(ind).f(delind) = [];
            handles.imports(ind).spectrum{datInd}(delind) = [];

            spc = abs(handles.imports(ind).spectrum{datInd});
            if normalize
                spc = spc/max(spc);
            end
            
            plot(handles.imports(ind).f/1e3, spc,...
                lnspecs{mod(handles.nplot,length(lnspecs))+1}, 'Color', colrs(handles.nplot+1,:),...
                'MarkerEdgeColor', colrs(handles.nplot+1,:), 'MarkerFaceColor',...
                colrs(handles.nplot+1,:), 'MarkerSize', 8, 'LineWidth', 4);
            lgnd{handles.nplot+1} = ['Import ' num2str(ind) '-' num2str(datInd)];
            handles.nplot = handles.nplot+1;
        end        
    end
    xlabel('Frequency (KHz)');
    ylabel('Spectrum');
% ------------------------------------------------------------------------%
else %% plot the time domain signals --------------------------------------    
    % plot the analytical time domain signal
    if plotAnalytical        
        x = str2double(get(handles.edt_sensorPos, 'String'));
        % symmetric modes
        if isempty(handles.response.us) && ~isempty(gaussSModes)
            [handles.response.us, handles.response.es] = handles.mypzt.timeResponse(...
                handles.excitation.f, 's', x, handles.excitation.signalSpectrum, -30, 0.08);
        end
        
        % antisymmetric modes
        if isempty(handles.response.ua) && ~isempty(gaussAModes)
            [handles.response.ua, handles.response.ea] = handles.mypzt.timeResponse(...
                 handles.excitation.f, 'a', x, handles.excitation.signalSpectrum, -30, 0.08);
        end
        
        % make the plots
        y = zeros(1,length(handles.excitation.t));
        for i=1:length(gaussSModes)
            if addModes
                if strainDisp
                    y = y + 1e6*handles.response.es{gaussSModes(i)};
                else
                    y = y + 1e9*handles.response.us{gaussSModes(i)};
                end
            else
                if strainDisp
                    y = 1e6*handles.response.es{gaussSModes(i)};
                else
                    y = 1e9*handles.response.us{gaussSModes(i)};
                end
                if normalize 
                    y = y/max(abs(y));
                end
            
                plot(handles.excitation.t*1e6, y, lnspecsT{mod(handles.nplot,length(lnspecsT))+1},...
                    'Color', colrs(handles.nplot+1,:), 'LineWidth', 4, 'MarkerSize', 6);
                lgnd{handles.nplot+1} = ['Analytical S' num2str(i-1)];
                handles.nplot = handles.nplot+1;
            end
        end        
        
        for i=1:length(gaussAModes)            
            if addModes 
                if strainDisp
                    y = y + 1e6*handles.response.ea{gaussAModes(i)};
                else
                    y = y + 1e9*handles.response.ua{gaussAModes(i)};
                end
            else                
                if strainDisp
                    y = 1e6*handles.response.ea{gaussAModes(i)};
                else
                    y = 1e9*handles.response.ua{gaussAModes(i)};
                end
                if normalize
                    y = y/max(abs(y));
                end
                plot(handles.excitation.t*1e6, y, lnspecsT{mod(handles.nplot,length(lnspecsT))+1},...
                    'Color', colrs(handles.nplot+1,:), 'LineWidth', 4, 'MarkerSize', 6);
                lgnd{handles.nplot+1} = ['Analytical A' num2str(i-1)];
                handles.nplot = handles.nplot+1;
            end
        end
        if addModes 
            if normalize 
                y = y/max(abs(y));
            end
            plot(handles.excitation.t*1e6, y, lnspecsT{mod(handles.nplot,length(lnspecsT))+1},...
                    'Color', colrs(handles.nplot+1,:), 'LineWidth', 4, 'MarkerSize', 6);
            lgnd{handles.nplot+1} = 'Analytical';
            handles.nplot = handles.nplot+1;
        end
    end
    
    % plot the imported time data 
    if ~isempty(selImport)        
        for i=1:length(selImport) 
            [ind, datInd] = getSelectionIndex(handles, selImport(i));

            % get the truncated signal according to the window start and end            
            sig = handles.imports(ind).data(:,datInd);
%             sig = sig - mean(sig);
            sig(1:handles.imports(ind).startInd(datInd)) = 0;
            sig(handles.imports(ind).endInd(datInd):end) = 0;
            if strainDisp 
                sig = 1e6*sig;
            else
                sig = 1e9*sig;
            end
            if normalize
                sig = sig/max(abs(sig));
            end

            plot(handles.imports(ind).t*1e6, sig, lnspecsT{mod(handles.nplot,length(lnspecsT))+1},...
                'Color', colrs(handles.nplot+1,:), 'LineWidth', 4, 'MarkerSize', 6);
            lgnd{handles.nplot+1} = ['Import ' num2str(ind) '-' num2str(datInd)];
            handles.nplot = handles.nplot+1;
        end        
    end
    
    xlabel('Time (\mus)');
    if strainDisp
        if normalize
            ylabel('Normalized strain');
        else            
            ylabel('\mu Strain');
        end
    else
        if normalize
            ylabel('Normalized displacement');
        else
            ylabel('Displacement (nm)');
        end
    end
end
if showLgnd
    legend(lgnd);
end
set(handles.txt_messages, 'String', 'Response calculated.');
guidata(hObject, handles)

% --- Executes on button press in cb_spectrumShowPlot.
function cb_spectrumShowPlot_Callback(hObject, eventdata, handles)
% hObject    handle to cb_spectrumShowPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_spectrumShowPlot
 isSel = get(hObject, 'Value');
 str = get(handles.lb_spectrumModes, 'String');
 
 if isSel && isempty(str) 
     updateSpectrumModes(hObject, handles);
 end


% --- Executes on button press in cb_gaussShowPlot.
function cb_gaussShowPlot_Callback(hObject, eventdata, handles)
% hObject    handle to cb_gaussShowPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of cb_spectrumShowPlot
 isSel = get(hObject, 'Value');
 str = get(handles.lb_gaussModes, 'String');
 
 if isSel && isempty(str) 
     updateGaussModes(hObject, handles);
 end
 
 function edt_sensorPos_Callback(hObject, eventdata, handles)
% hObject    handle to edt_sensorPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(hObject,'String'));
prevnum = get(hObject, 'UserData');
if ~isscalar(num) || ~isnumeric(num)
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetResponseVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end

 
%--- Executes on button press in pb_plotGauss.
function pb_plotGauss_Callback(hObject, eventdata, handles)
% hObject    handle to pb_plotGauss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BW = str2double(get(handles.edt_gaussExcitBW, 'String'));
f0 = 1e3*str2double(get(handles.edt_gaussExcitF0, 'String'));
T = 1e-6*str2double(get(handles.edt_gaussExcitT, 'String'));

[~,  ~, ~, ~] = Lamby.gauspulse(f0, BW, T, 0, 1);


% --- Executes on button press in pb_clearSel.
function pb_clearSel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_clearSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.lb_importedData, 'Value', []);



% ---------------------------------------------------------------------- %%
%%%%%%% user functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateSpectrumModes(hObject, handles)

maxFreq = 1e3*str2double(get(handles.edt_excitFreqEnd, 'String'));

set(handles.lb_spectrumModes, 'String', '');

[cs, ~] = handles.plate.phaseVelocity(maxFreq, 's', 0.01);
[ca, ~] = handles.plate.phaseVelocity(maxFreq, 'a', 0.01);

str = cell(length(cs)+length(ca),1);
for i=1:length(cs)
    str{i} = ['S'  num2str(i-1)];
end

for i=1:length(ca)
    str{i+length(cs)} = ['A'   num2str(i-1)];
end

set(handles.lb_spectrumModes, 'String', str);
set(handles.lb_spectrumModes, 'UserData', [length(cs) length(ca)]);

% update the mode list for the Gaussian Pulse excitation
function updateGaussModes(hObject, handles)

BW = str2double(get(handles.edt_gaussExcitBW, 'String'));
f0 = 1e3*str2double(get(handles.edt_gaussExcitF0, 'String'));
T = 1e-6*str2double(get(handles.edt_gaussExcitT, 'String'));

[handles.excitation.t,  handles.excitation.signal, handles.excitation.f, ...
 handles.excitation.signalSpectrum] = Lamby.gauspulse(f0, BW, T, 0, 0);

f = fftshift(handles.excitation.f);
signalSpectrum = fftshift(handles.excitation.signalSpectrum);
negf = find(f < 0);
f(negf) = [];
signalSpectrum(negf) = [];

[ind1, ind2] = Lamby.findSignalLimit(signalSpectrum, -30);

% handles.excitation.f(ind) = [];
% handles.excitation.signalSpectrum(ind) = [];

% indright = find(handles.excitation.f >= 0);           
% fright = handles.excitation.f(indright);
% sigright = handles.excitation.signalSpectrum(indright);
% inds = find(abs(sigright) > 1e-3*max(abs(sigright)));
% fmax = max(fright(inds));
fmax = max(handles.excitation.f(ind2));

set(handles.lb_gaussModes, 'String', '');

[cs, ~] = handles.plate.phaseVelocity(fmax, 's', 0.01);
[ca, ~] = handles.plate.phaseVelocity(fmax, 'a', 0.01);

str = cell(length(cs)+length(ca),1);
for i=1:length(cs)
    str{i} = ['S'  num2str(i-1)];
end

for i=1:length(ca)
    str{i+length(cs)} = ['A'   num2str(i-1)];
end

set(handles.lb_gaussModes, 'String', str);
set(handles.lb_gaussModes, 'UserData', [length(cs) length(ca)]);

guidata(hObject, handles);
resetResponseVariables(hObject, handles);


% get the indices for the selected item in the imported data listbox
function [ind, datInd] = getSelectionIndex(handles, sel)
% handles   a structure containing all the handles and data in the GUI
% sel       the selected item in the listbox
sizes = get(handles.pb_importData, 'UserData');
sm = cumsum(sizes);    
ind = find(sel <= sm,1);
if ind > 1 
    datInd = sel - sum(sizes(1:ind-1));
else
    datInd = sel;
end    

% empty the response variables if some variable changed 
function resetResponseVariables(hObject, handles)
handles.response.Hs_strain = [];
handles.response.Hs_disp = [];
handles.response.fsInd = [];
handles.response.Ha_strain = [];
handles.response.Ha_disp = [];
handles.response.faInd = [];
handles.response.us = [];
handles.response.es = [];
handles.response.ua = [];
handles.response.ea = [];
handles.spectrum.f = [];
handles.spectrum.Hs_strain = [];
handles.spectrum.Hs_disp = [];
handles.spectrum.fsInd = [];
handles.spectrum.Ha_strain = [];
handles.spectrum.Ha_disp = [];
handles.spectrum.faInd = [];
guidata(hObject, handles);


% empty the response variables if some variable changed 
function resetDispersionVariables(hObject, handles)
handles.dispersion.cpa = [];
handles.dispersion.cps = [];
handles.dispersion.cga = [];
handles.dispersion.cgs = [];
handles.dispersion.cpaTol = [];
handles.dispersion.cpsTol = [];
handles.dispersion.cgaTol = [];
handles.dispersion.cgsTol = [];
handles.dispersion.fcpa = [];
handles.dispersion.fcps = [];
handles.dispersion.fcga = [];
handles.dispersion.fcgs = [];
handles.dispersion.f = [];
guidata(hObject, handles);



% --------------------------------------------------------------------
function mnu_save_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
response = handles.response;
spectrum = handles.spectrum;
excitation = handles.excitation;
save lamby_data response spectrum excitation;



function edt_symPhaseTol_Callback(hObject, eventdata, handles)
% hObject    handle to edt_symPhaseTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(hObject,'String'));
prevnum = get(hObject, 'UserData');
if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetDispersionVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end



function edt_asymPhaseTol_Callback(hObject, eventdata, handles)
% hObject    handle to edt_asymPhaseTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(hObject,'String'));
prevnum = get(hObject, 'UserData');
if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetDispersionVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end



function edt_symGroupTol_Callback(hObject, eventdata, handles)
% hObject    handle to edt_symGroupTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num = str2double(get(hObject,'String'));
prevnum = get(hObject, 'UserData');
if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetDispersionVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end

function edt_asymGroupTol_Callback(hObject, eventdata, handles)
% hObject    handle to edt_asymGroupTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_asymGroupTol as text
%        str2double(get(hObject,'String')) returns contents of edt_asymGroupTol as a double
num = str2double(get(hObject,'String'));
prevnum = get(hObject, 'UserData');
if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetDispersionVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end

function edt_dispFreqEnd_Callback(hObject, eventdata, handles)
% hObject    handle to edt_dispFreqEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_dispFreqEnd as text
%        str2double(get(hObject,'String')) returns contents of edt_dispFreqEnd as a double
num = str2double(get(hObject,'String'));
fstart = str2double(get(handles.edt_dispFreqStart, 'String'));
prevnum = get(hObject, 'UserData');
if ~isPositiveNumber(num) || num < fstart
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetDispersionVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end

function edt_dispFreqStep_Callback(hObject, eventdata, handles)
% hObject    handle to edt_dispFreqStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_dispFreqStep as text
%        str2double(get(hObject,'String')) returns contents of edt_dispFreqStep as a double
num = str2double(get(hObject,'String'));
prevnum = get(hObject, 'UserData');
if ~isPositiveNumber(num)
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetDispersionVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end

function edt_dispFreqStart_Callback(hObject, eventdata, handles)
% hObject    handle to edt_dispFreqStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_dispFreqStart as text
%        str2double(get(hObject,'String')) returns contents of edt_dispFreqStart as a double
num = str2double(get(hObject,'String'));
fend = str2double(get(handles.edt_dispFreqEnd, 'String'));
prevnum = get(hObject, 'UserData');
if ~isPositiveNumber(num) || num > fend
    set(hObject, 'String', num2str(prevnum));
else
    if num ~= prevnum 
        resetDispersionVariables(hObject, handles);
    end
    set(hObject, 'UserData', num);
end
