clc, clear

%% Data Import
load("bridgeloads.mat") %loads data corresponding to max shear, max bending, and total train weight experienced

%% Bridge Geometry

% Bridge Geometry
classdef BridgeGeometryConstructor
    properties
        length {mustBeNumeric}
        
    end
end

geometry = []