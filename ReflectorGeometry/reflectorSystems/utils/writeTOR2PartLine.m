function writeTOR2PartLine(fileID,propertyName,part1Name,part1Value,part2Name,part2Value,useUnit,useComma)
comma = '';
unit1 = ", ";
unit2 = "";
if nargin == 6
    useComma = 0;
    useUnit = 0;
elseif nargin == 7
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
    part2 = strcat(part2Name,": """,part2Value,"",unit2);
else
    part2Value = num2str(part2Value);
    part2 = strcat(part2Name,": ",part2Value,unit2);
end
catString = strcat(part1,part2);
fprintf(fileID, '  %-17s: struct(%s)%s\n',propertyName,catString,comma);