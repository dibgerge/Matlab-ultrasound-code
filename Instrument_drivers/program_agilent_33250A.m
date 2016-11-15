%% -- Parameter initialization
paddress = 10;
fc = 20e3;
N = 4;

%% --- setup the signal
T = N/fc;
t = linspace(0, T, 2^14);
y = 0.5*(1-cos(2*pi*fc*t/N)).*cos(2*pi*fc*t);

%% --- Open the funciton generator object
h = instrfind('Type', 'gpib', 'BoardIndex', 0, 'PrimaryAddress', paddress, 'Tag', '');
if isempty(h)
    h = gpib('NI', 0, paddress);
else
    fclose(h);
    h = h(1);
end
set(h, 'InputBufferSize', 512000, 'OutputBufferSize', 512000, 'Timeout', 130);
fopen(h);
 
fprintf(h, ['FREQ ' num2str(1/T)]);
fprintf(h, 'APPL:SIN');
fprintf(h, 'VOLT 10');

data = num2str(y, '%1.15f,');       %fixedpoint numbers seperated by commas
data(end) = [];                     %delete final comma
data = strrep(data,' ','');         %in case any spaces snuck in, delete them.

fprintf(h, ['data volatile, ', data]);
fprintf(h, ['data:copy HAN' num2str(N) 'CYC']);
fclose(h);

