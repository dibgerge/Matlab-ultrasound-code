function device = startGenerator()

% --- Agilent 33220A function generator initialization --- %
serialnum = 'MY44048138'; %get it from Agilent IO suite connection expert
%serialnum = 'MY44060125'; %get it from Agilent IO suite connection expert
manufacturer = 'NI';

%h = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0957::0x0407::MY44048138::0::INSTR');

visaObj = getVisaObject(manufacturer, serialnum);
device = icdevice('agilent_33220a.mdd', visaObj);    

disconnect(device)
connect(device)



end

