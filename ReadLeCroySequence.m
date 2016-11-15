function w = ReadLeCroySequence(pth, fname, range)
%READLECROYSEQUENCE Reads a collection of files in a given directory, name
%       prefix, and given range
%Inputs: 
%  pth: Directory path where the data files are. 
%  fname: A name prefix for the files to be loaded. the filenames are
%         '<fname>xxxxx.trc'
%  range:  A vector indicating the file numbers to be loaded. This is an
%         optional argument. If it is not defined, all files in directory 'pth' with prefix 'fname'
%         are loaded.
%Outputs: 
%  t: The common time vector for all loaded files. 
%  y: A matrix with the number of columns corresponding to the number of
%     files loaded.
%NOTE. All files should have the same number of entries (same vector
%      length, sampling frequency, and time vector), otherwise this
%      function might give inconsistent results

if nargin == 2
    d = dir(pth);
    dnames = {d.name};
    filenames = dnames(strncmp(fname, dnames, length(fname)));
    N = length(filenames);
    range = [];
else
    N = length(range);
end

w = UTlib.utcollection;

for i=1:N
    if ~isempty(range)
        fnamefull = fullfile(pth, sprintf([fname, '%05d.trc'], range(i)));
    else
        fnamefull = fullfile(pth, filenames{i});
    end
    
    d = ReadLeCroyBinaryWaveform(fnamefull);
    
    w.ut{i} = UTlib.utsignal(d.x, d.y);
end

end

