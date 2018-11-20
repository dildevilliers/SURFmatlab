        function obj = setXrangePos(obj)
            % Attempts to set the x-range (ph, az, or ep) for the angular gridTypes in the FarField
            % object to positive - [-180,180] is transformed to [0,360]
            % with the redundant -180 cut replaced by a redundant 360 cut
            % (if it was found).  The resulting object has the same number
            % of field points as the input object.
            
            if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                i180 = find(obj.x == -pi);   % Redundant
                i0 = find(obj.x == 0);      % Will become redundant
                Nredun = min(numel(i180),numel(i0));    % How many we have to remove/add
                % Truncate the indexes to the shortest length
                i180 = i180(1:Nredun);
                i0 = i0(1:Nredun);
                redunFound = Nredun > 0;
                % First remove the ph=-180 from the matrix...
                if redunFound
                    obj.x(i180) = [];
                    obj.y(i180) = [];
                    obj.E1(i180,:) = [];
                    obj.E2(i180,:) = [];
                    obj.E3(i180,:) = [];
                end
                % Now shift the x values
                lt0 = find(obj.x < 0);
                obj.x(lt0) = obj.x(lt0) + 2*pi;
                % if the -180 was removed, add a 360 cut
                if redunFound
                    obj.x = [obj.x;ones(numel(i0),1).*2.*pi];
                    obj.y = [obj.y;obj.y(i0)];
                    obj.E1 = [obj.E1;obj.E1(i0,:)];
                    obj.E2 = [obj.E2;obj.E2(i0,:)];
                    obj.E3 = [obj.E3;obj.E3(i0,:)];
                end
                % Sort
                obj = obj.sortGrid;
            else
                warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
            end
            obj = setRangeTypes(obj);
        end