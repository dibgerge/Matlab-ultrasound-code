function device = startOscilloscope()

% --- DSO 1004A oscilloscope initialization --- %
serialnum_oscilloscope = 'CN50056834'; %get it from Agilent IO suite connection expert
serialnum_oscilloscope = 'CN53262520'; 
manufacturer = 'agilent';

h = instrfind('Type', 'scope', 'Name', 'scope-agilent_DSO1004A');

if isempty(h)
    visaObjOscilloscope = getVisaObject(manufacturer, serialnum_oscilloscope);
    device = icdevice('agilent_DSO1004A.mdd', visaObjOscilloscope);
else
    disconnect(h);
    device = h(1);
end

connect(device)

end

