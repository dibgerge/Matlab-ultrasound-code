function device = getVisaObject(manufacturer, serialNum)

% Check arguments
if ~ischar(manufacturer) || ~ischar(serialNum)
    error('manufacturer and serial arguments should be a string')
end
    
manufacturer = upper(manufacturer);
serialNum = upper(serialNum);

hwinfo = instrhwinfo('visa');

hwvendor =[];
for i=1:length(hwinfo.InstalledAdaptors)
    k = findstr(upper(hwinfo.InstalledAdaptors{i}), manufacturer);
    if ~isempty(k)
        hwvendor = hwinfo.InstalledAdaptors{i};
        break;
    end
end

hwinfoVendor = instrhwinfo('visa', hwvendor);
visaAddress = [];
for i=1:length(hwinfoVendor.ObjectConstructorName)
    k = findstr(upper(hwinfoVendor.ObjectConstructorName{i}), serialNum);
    if ~isempty(k)
        ioresourcedescriptor = regexpi(hwinfoVendor.ObjectConstructorName{i}, ...
            '[\w|\.|:|-]*''', 'match'); 
        visaAddress = ioresourcedescriptor{end}(1:end-1);
        break;
    end
end

device = visa('NI', visaAddress);
device.InputBufferSize = 500000;
device.OutputBufferSize = 500000;
end

