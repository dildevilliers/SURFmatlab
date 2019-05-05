function writeFreqInfo2TOR(fileID,freqName,freqVal)
fprintf(fileID, '%s  frequency  \n',freqName);
fprintf(fileID, '(\n');
if checkIfParse(freqVal)
    freqVal = parseThisString(freqVal);
    fprintf(fileID, '  frequency_list   : sequence("%s" GHz)\n',freqVal);
else
    freqVal = num2str(freqVal);
    fprintf(fileID, '  frequency_list   : sequence(%s GHz)\n',freqVal);
end
fprintf(fileID, ')\n\n');