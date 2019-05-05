function writeTOR3PartLine(fileID,propertyName,part1Name,part1Value,part2Name,part2Value,part3Name,part3Value,useUnit,useComma)
comma = '';
unit1 = ", ";
unit2 = "";
if nargin == 8
    useComma = 0;
    useUnit = 0;
elseif nargin == 9
    useComma = 0;
end
if useComma
    comma = ',';
end
if useUnit
    unit1 = " m, ";
    unit2 = " m";
end

if checkIfParse(part1Value)
    part1Value = parseThisString(part1Value);
    part1 = strcat(part1Name,": """,part1Value,"",unit1);
else
    part1Value = num2str(part1Value);
    part1 = strcat(part1Name,": ",part1Value,unit1);
end
if checkIfParse(part2Value)
    part2Value = parseThisString(part2Value);
    part2 = strcat(part2Name,": """,part2Value,"",unit1);
else
    part2Value = num2str(part2Value);
    part2 = strcat(part2Name,": ",part2Value,unit1);
end
if checkIfParse(part3Value)
    part3Value = parseThisString(part3Value);
    part3 = strcat(part3Name,": """,part3Value,"",unit2);
else
    part3Value = num2str(part3Value);
    part3 = strcat(part3Name,": ",part3Value,unit2);
end
catString = strcat(part1,part2,part3);
fprintf(fileID, '  %-17s: struct(%s)%s\n',propertyName,catString,comma);