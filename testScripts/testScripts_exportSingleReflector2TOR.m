clear all

single_reflector = singleReflector(1.2,0.8,0.2);
%exportSingleReflectorToTOR(obj,fullpathName,freqValGhz,prefixName)
single_reflector.exportSingleReflectorToTOR('C:\testing_area\singleReflTestScript.tor',4,'sngle');

%Yes... it is that simple :P