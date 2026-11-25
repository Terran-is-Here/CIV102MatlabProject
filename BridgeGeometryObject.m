classdef BridgeGeometryObject
    % Bridge Geometry Constructor Class
    % Contains the data pertaining to the bridge geometry's PROPERTIES within a region. 
    % Maybe I can implement a geometric constructor segment later on
    % pertaining to calculating the properties within a region based on
    % some array of geometric objects
    properties
        length {mustBeNumeric} = 0 %section length
        momentOfInertia {mustBeNumeric} = 0 %value of I
        centroidalAxis {mustBeNumeric}  = 0% value of Y Bar / Centroidal / Neutral Axis, from the bottom. 
        height {mustBeNumeric} = 0% Height from the bottom of the member to the top

        % Start and End Point of corresponding bridge geometry that this
        % points towards. 
        startPoint {mustBeNumeric} = 0
        endPoint {mustBeNumeric} = 0
    end

    methods
        % Object Constructor Method 
        function obj = BridgeGeometryObject(sectionLength, momentOfInertia, centroidalAxis, height, start, end)
            obj.length = sectionLength; 
            obj.momentOfInertia = momentOfInertia; 
            obj.centroidalAxis = centroidalAxis; 
            obj.height = height;
            obj.startPoint = start;
            obj.endPoint = end; 
        end
        % Creates an (???) dimensioned array, where the corresponding
        % values present represent the corresponding fail 
        function arr = BridgeNodeConstructor(obj, MaterialObj, variableDep)  
            arguments
                obj BridgeGeometryObject
                MaterialObj MaterialProperties % expecting a material properties block. 
                variableRatios (1,2) = [0,0] % If variable geometry exists, passes these properties into a polyval expression
            end

            
            
        end
    end
end