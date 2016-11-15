function [data,header] = ReadBinarySAFT(fn,imgFlag,saveFlag,saveDir)
% Coverts SAFT binary files to *.mat files
% INPUTS:
%   fn =      *.rf file (full path) -- **SAFT binary format**
%   imgFlag = true/false on whether to save images (default = false)
%   saveFlag = true/false on whether to save data (default = true)
%   saveDir = alternate location for saving outputs
%       ('same' for with the files, default = pwd)
% OUTPUTS (can be supressed):
%   data =    wfm matrix (Usound x Scan x Index)
%   header =  all header information that could be extracted
%   coords = axis-location matricies for each data point
%
% --- [*.mat/*.jpg] files saved ---

%% ESTABLISH INPUTS
if nargin < 1
    fprintf('no filename input.\n');
    return;
elseif nargin < 2
    imgFlag = false;
    saveFlag = false;
elseif nargin < 3
    saveFlag = false;
end

ind = [0 strfind(fn,'\')];
sName = fn(ind(end)+1:end-3);

if nargin < 4
    saveDir = [];
elseif strcmpi(saveDir,'same')
    saveDir = fn(1:ind(end));
else
    if ~strcmpi(saveDir(end),'\'); saveDir(end+1)='\'; end
end

h_NUM = 2^11; %number of bytes in the file-header

%% CHECK FOR *.MAT VERSION
fMat = [fn(1:end-3) '.mat'];
try
    load(fMat,'header','data','coords')
    cscan = coords.C;
    bscan = coords.B;
    dscan = coords.D;
    S = coords.scanLoc;
    I = coords.indexLoc;
    U = coords.usoundLoc;
catch
if strcmp(fn(end-2:end),'mat'); return; end
%% LOAD THE FILE
fid = fopen(fn);
data8bit = fread(fid, 'ubit8')';
fclose(fid);

% Read out the header information (always 8-bit)
Line = data8bit(1:h_NUM);
%% BUILD HEADER
header = buildHeader(Line);

%% RECONSTRUCT THE WAVEFORMS
% Actual data may not be 8-bit
if header.dataType_bit ~= 8
    header.dataType_bit = 16;
    fid = fopen(fn);
    data16bit = fread(fid, 'ubit16')';
    fclose(fid);
    
    rawData = data16bit(h_NUM/2+1:end);
else
    rawData = data8bit(h_NUM+1:end);
end

bits = header.dataType_bit;
dh_NUM = 2^5/(bits/8); %number of bytes in the data-header
AscanLen = header.ascanLen;

totalLen = length(rawData);

% Check to make sure that the file is intact
numAscans = totalLen/(AscanLen+dh_NUM);
% % % % points = header.numSteps;
points = [header.xPoints header.yPoints];
if numAscans ~= points(1)*points(2)
    fprintf('Number of Ascans incorrect.\n')
    return;
end

% Extract waveforms and build the 3D data matrix
data = NaN*ones(AscanLen,points(1),points(2));
for ii = 1:points(2) %index
    for jj = 1:points(1) %scan
        loc = ((ii-1)*points(1)+jj-1)*(dh_NUM+AscanLen);
        Ascan = rawData(loc+dh_NUM+(1:AscanLen))';
        data(:,jj,ii) = Ascan-2^(bits-1);
    end
end

%% IMAGE CONSTRUCTION
% refractAng = header.scanDir_deg;
% 
% % Axes
% span = header.sampStop_in-header.sampStart_in;
% lim.ascan = header.sampStart_in:span/(header.ascanLen-1):header.sampStop_in;
% lim.usound = lim.ascan*cosd(refractAng);
% lim.scan = header.xStart_in(1):header.xStep_in:header.xStop_in;
% lim.scan = lim.scan(1:end-1); 
% lim.index = header.yStart_in:header.yStep_in:header.yStop_in;
% lim.index = lim.index(1:end-1);
% 
% S_temp = NaN*ones(size(data(:,:,1)));
% for jj = 1:size(data,2)
%     sCol = lim.scan(jj)+lim.ascan'*sind(refractAng);
%     S_temp(:,jj) = sCol;
% end
% U = NaN*ones(size(data)); S = U; I = U;
% for kk = 1:size(data,3)
%     U(:,:,kk) = lim.usound'*ones(1,size(data,2));
%     S(:,:,kk) = S_temp;
%     I(:,:,kk) = lim.index(kk);
% end
% res = header.xStep_in(1);
% S = round(S/res)*res;
% lim.scan = unique(S)';
% 
% clearvars S_temp
% % D-scan
% dscan = max(abs(data),[],2);
% dscan = reshape(dscan,numel(lim.usound),numel(lim.index),1);
% 
% % B-scan
% bscan = max(abs(data),[],3);
% 
% % C-scan
% cscan = NaN*ones(numel(lim.index),numel(lim.scan));
% for ii = 1:size(cscan,1) %index
%     temp = data(:,:,ii);
%     for jj = 1:size(cscan,2) %scan
%         TF = S(:,:,ii)==lim.scan(jj); TF = TF(:);
%         cscan(ii,jj) = max(temp(TF));
%     end
% end
% 
% coords = struct('C',cscan,'B',bscan,'D',dscan,...
%     'scanLoc',S,'indexLoc',I,'usoundLoc',U);

end

if imgFlag
    % Generate\save plot image (*.jpg)
    figure(1)
    set(gcf,'units','normalized','position',[.01 .06 .98 .84])
    subplot(2,2,[1 3])
    h(3) = surface(unique(S),unique(I),cscan);
    axis ij; axis image;
    ylabel('Index (in)'); xlabel('Scan (in)');
    subplot(2,2,2)
    h(2) = surface(S(:,:,1),U(:,:,1),bscan);
    axis ij; axis image;
    ylabel('Usound (in)'); xlabel('Scan (in)');
    subplot(2,2,4)
    h(1) = surface(unique(I),unique(U),dscan);
    axis ij; axis image;
    ylabel('Usound (in)'); xlabel('Index (in)');
    set(h,'edgealpha',0);
    
    print(gcf,[saveDir sName '.jpg'],'-djpeg','-r1000');
    close(gcf)
end

%% Save data (*.mat)
if saveFlag
    save([saveDir sName '.mat'],'header','data','coords')
end

end

function header = buildHeader(SAFT_text)

fieldList = {...
    10  ''; % "ASCII"
    81  'title'; % title (...user-text?)
    9   'date'; % date
    9   'time'; % time
    2   ''; % domainFlag (0=time, 1=freq)
    12  ''; % number of dataSets
    7   'dataMin' % dataMin (digi-units)
    7   'dataMax' % dataMax (digi-units)
    17  'dataAve' % dataAve (digi-units)
    4   ''; % projectionFlag (0=No, 1=Yes)
    2   ''; % units (0=inches)
    7   ''; % sDataType (0=8bit, 1=16bit)
    51  ''; % sCalFilename (...user-text?)
    %%
    81  'probeComment'; % transducer-comment (user-text)
    17  'freq_MHz'; % frequency (MHz)
    17  'RxWedgePath_in'; % reciever wedge-path (in)
    17  'TxWedgePath_in'; % transmit wedge-path (in)
    17  'RxWedgeVel_inPERs'; % wedge-velocity - R (in/s)
    17  'TxWedgeVel_inPERs'; % wedge-velocity - T (in/s)
    17  'beamDia_in'; % beam diameter (in)
    17  ''; % ??? (refracted angle(deg)??)
    17  'incidentAng_deg'; % incident angle wedge (deg)
    17  'skew_deg'; % squint (skew) direction (0 = x, 90.0 = y)
    2   ''; % mode (0=PE, 1=TSAFT, 2=TSAFT-2)
    17  ''; % initial x-offset of T to R (in)
    17  ''; % f_number (???)
    17  ''; % x-offset in wedge
    17  ''; % y-offset in wedge
    13  ''; % transducer spares (???)
    %%
    81  'matComment'; % material-comment (user-text)
    17  'matVel_inPERs'; % material velocity (in/s)
    17  'refractAng_deg'; % refracted angle (deg)
    17  'matThickness_in'; % thickness (in)
    17  ''; % diameter of pipe (in) (ID? OD?)
    17  ''; % diameter of track (in)
    7   ''; % material type (0=unknown, 1=plate, 2=pipe, 3=nozzle)
    23  ''; % material spares (???)
    %%
    81  'sampComment'; % sampling-comment (user-text)
    17  ''; % delay increment (10ns)
    17  ''; % initial delay (10ns)
    7   'ascanLen'; % number of points in an Ascan (#)
    17  'sampStart_in'; % sound-path to start of data-window (in)
    17  'sampStop_in'; % depth to end of data-window (in)
    7   'averages'; % number of averages (#)
    17  ''; % min time between pulses (s)
    17  ''; % step along each wave path (in)
    11  'windowStart_ns'; % time to start of sampling (ns)
    11  'windowStop_ns'; % time to end of sampling (ns)
    1   ''; % depth to end of data window
    %%
    81  'scanComment'; % scanning-comment (user-text)
    17  'scanDir_deg'; % scan direction (0=primary axis)
    17  'xStart_in'; % initial reciever primary coordinate (in)
    17  'yStart_in'; % initial reciever secondary coordinate (in)
    17  'xStop_in'; % final reciever primary coordinate (in)
    17  'yStop_in'; % final reciever secondary coordinate (in)
    17  'xStep_in'; % x-increment (in)
    17  'yStep_in'; % y-increment (in)
    7   'xPoints'; % total points in x-direction
    7   'yPoints'; % total points in y-direction
    2   ''; % downstreamFlag (Y=scanner points downstream, N=upstream)
    12  ''; % transmit_half_vees
    12  ''; % receive_half_vees
    17  ''; % number of estimated half V's before signal encounters object plane (arrow)
    17  ''; % init_pos_tertiary_axis
    17  ''; % final_pos_tertiary_axis
    2   ''; % scan_toward_trackFlag (Y=collects moving towards, N=away)
    4   ''; % scanner type (???)
    7   ''; % scan pattern (???)
    17  ''; % z-increment
    %% PROCESSING
    81+9+9+12+12+12+7+2+7+7+7+7+20+7+7+7+17+7+7+12+17+17+18 '';
    %% NOZZLE
    17+17+17+17 '';
    %% OTHER
    17+17+52 '';
    %% TVG
    17+7+1+17+17+7+17+7 '';
    %%
    7   ''; % digi type (0=no digitizer, 1=str864, 2=CS12100)
    7   ''; % TVG type (0=unknown, 1=tek bin, 2=ISA card)
    7   ''; % pulser type (0=unknown, 1=pcpr100)
    17  'Vpp'; % Vpp
    7   ''; % sync mode (0=hardware, 1=software, 2=rearm mode)
    7   ''; % nCalScreenValue
    7   ''; % nCalBaseline
    17  ''; % dCalVolts
    7   ''; % nCalGainSwitch
    7   ''; % nCalGainNob
    7   ''; % nBaseline
    7   ''; % nGainSwitch
    7   ''; % nGainNob
    17  ''; % dCalCountsPerVolt
    17  ''; % dTimeOfFlight1QT
    17  ''; % dTimeOfFlightHalfT
    17  ''; % dTimeOfFlight3QT
    13  ''; % positions spare (???)
    %%
    7+7+7+7+4 '';
    };
fieldSize = cell2mat(fieldList(:,1));

count = 0;
header = struct();
for ii = 1:numel(fieldSize)
    cSet = count+1:sum(fieldSize(1:ii));
    temp = SAFT_text(cSet);
    
    temp = temp(temp~=205);
    while ~isempty(temp) && (temp(end)==0 || temp(end)==32); temp = temp(1:end-1); end
    temp = char(temp);
    if ~isempty(fieldList{ii,2})
% % fprintf('(%02.0f)\t%s:\t%s\n',ii,fieldList{ii,2},temp);
        if ~isnan(str2double(temp)); temp = str2double(temp); end
        header.(fieldList{ii,2}) = temp;
    end
    
    switch ii
        case 12
            temp = str2double(temp);
            header.dataType_bit = 8+8*(temp>0);
        case 39
            temp = str2double(temp);
            header.Ts = (temp*10)/1e9;
            header.Fs = 1/header.Ts;
    end
    
    count = cSet(end);
end
% % fprintf('** %d **\n',count)

end