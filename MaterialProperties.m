classdef MaterialPropertiesConstructor
    %   Data Container Class pertaining to the material properties of a
    %   substance

    properties
        %using default values for matboard for instantiating an object
        tensileStressMax {mustBeNumeric} = 30 %MPa
        compressiveStressMax {mustBeNumeric}= 6 %MPa
        shearStressMax {mustBeNumeric} = 4 %MPa
        elasticModulus {mustBeNumeric} = 4000 %MPa
        poissonRatio = 0.2
        
    end

    methods
        %Constructor Method 
        function obj = MaterialProperties(tensileStressMax,compressiveStressMax,shearStressMax,elasticModulus,poissonRatio)
            obj.tensileStressMax = tensileStressMax; 
            obj.compressiveStressMax = compressiveStressMax; 
            obj.shearStressMax = shearStressMax; 
            obj.elasticModulus = elasticModulus; 
            obj.poissonRatio = poissonRatio; 
        end
    end
end