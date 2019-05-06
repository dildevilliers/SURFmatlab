function writeVariable2TOR(fileID,variableName, variableValue)

fprintf(fileID, '%s  real_variable  \n',variableName);
fprintf(fileID, '(\n');
writeTORstandardLine(fileID,'value',variableValue);
fprintf(fileID, ')\n\n');

