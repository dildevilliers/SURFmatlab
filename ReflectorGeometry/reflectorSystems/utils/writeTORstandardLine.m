function writeTORstandardLine(fileID,propertyName,variableValue,useUnit,useComma,useQuotes)
unitAndComma = '';
quotes = '"';
if nargin == 3
    useUnit = 0;
    useComma = 0;
    useQuotes = 1;
elseif nargin == 4
    useComma = 0;
    useQuotes = 1;
elseif nargin == 5
    useQuotes = 1;
end
if useUnit
    unitAndComma = ' m';
end
if useComma
    unitAndComma = strcat(unitAndComma,',');
end
if ~useQuotes
    quotes = '';
end
if checkIfParse(variableValue)
    variableValue = parseThisString(variableValue);
    fprintf(fileID, '  %-17s: %s%s%s%s\n',propertyName,quotes,variableValue,quotes,unitAndComma);
else
    variableValue = num2str(variableValue);
    fprintf(fileID, '  %-17s: %s%s\n',propertyName,variableValue,unitAndComma);
end