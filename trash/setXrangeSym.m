        function obj = setXrangeSym(obj)
            % Attempts to set the x-range (ph, az, or ep) for the angular gridTypes in the FarField
            % object to symmetrical - [0,360] is transformed to [-180,180]
            % with the redundant 360 cut replaced by a redundant -180 cut
            % (if it was found).  The resulting object has the same number
            % of field points as the input object.
            
            if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                i360 = find(obj.x == 2*pi);   % Redundant
                i180 = find(obj.x == pi);      % Will become redundant
                Nredun = min(numel(i360),numel(i180));    % How many we have to remove/add
                % Truncate the indexes to the shortest length
                i360 = i360(1:Nredun);
                i180 = i180(1:Nredun);
                redunFound = Nredun > 0;
                % First remove the ph=360 from the matrix...
                if redunFound
                    obj.x(i360) = [];
                    obj.y(i360) = [];
                    obj.E1(i360,:) = [];
                    obj.E2(i360,:) = [];
                    obj.E3(i360,:) = [];
                end
                % Now shift the x values
                lt0 = find(obj.x > pi);
                obj.x(lt0) = obj.x(lt0) - 2*pi;
                % if the 360 was removed, add a -180 cut
                if redunFound
                    obj.x = [obj.x;ones(numel(i180),1).*-pi];
                    obj.y = [obj.y;obj.y(i180)];
                    obj.E1 = [obj.E1;obj.E1(i180,:)];
                    obj.E2 = [obj.E2;obj.E2(i180,:)];
                    obj.E3 = [obj.E3;obj.E3(i180,:)];
                end
                % Sort
                obj = obj.sortGrid;
            else
                warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
            end
            obj = setRangeTypes(obj);
        end