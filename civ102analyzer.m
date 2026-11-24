clc, clear all, clf %clear code parameters as per usual
%% Documentation
% In terms of the bridge, assume it takes the following structure in terms
% of notation: 

% A=========================================B
% Where A is the left-most endpoint, and B is the right-most endpoint. 
% This analyzer solves from traversing a train from left TO right. 

%% Todo
% Rewrite Load Simulation to accept arbitrary assumptions for resolution. 
%% Simulation Parameters
resolution = 1200; %break bridge up into 1200 segments, resolution of 1mm per node 
L = 1200;
n_train = 3; % train locations present, evenly split across span 
global dx
dx = L/resolution; %spatial step

shearForceDiagrami = zeros(n_train,resolution+1);
bendingMomentDiagrami = zeros(n_train, resolution+1);

% Material Characteristics 
materialThickness = 1.27;  %in mm

% Train Loading Characteristics

totalTrainWeight = 400; % in terms of Newtons

cartWeights = ones(1,3).*totalTrainWeight/3; % Assume that the cartWeights is arranged in Cart 2, Cart 1, Locomotive. 

% in terms of millimeters, adjusts based on ratio of resolution and
% distance
trainLoadPositions = [-856 -680 -516	-340	-176	0].*dx; %right-most value is locomotive

trainLoads = ones(1,6).*totalTrainWeight/6*-1; %initial load is fully uniformly distributed

%% Bridge Geometry
%Parameters are [x of cross-section (or cross-sectional change location),
%etc. 
params = [0,0,0,0;
    0,0,0,0;
    0,0,0,0];


%% Functions

% numerical integral
function out = NumericalIntegral(arr)
    global dx;
    out = zeros(1,size(arr,2));
    currentvalue = 0;
    for index = 1:size(arr,2);
        currentvalue  = currentvalue + arr(index)*dx;
        out(index) = currentvalue;
    end
    out
end

function arr = ReplaceNegatives(arr) 
    arr = max(arr,0); %rep
end

function out = NonNegativeToOne(arr)
    out= arr;
    out(out>=0) = 1;
    out(out<0) = 0;
end

% Take in a new set of cart weights, update axle loads and total train
% weight respectively. 
function [newTotalTrainWeight, newTrainLoads] =  updateTrainLoads(newCartWeights)
    newTotalTrainWeight = sum(newCartWeights,2) 
    newTrainLoads = ones(1,6) 
    newTrainLoads = newTrainLoads(1:2).*newCartWeights(1)/2 %cart 2
    newTrainLoads = newTrainLoads(3:4).*newCartWeights(2)/2 % cart 1
    newTrainLoads = newTrainLoads(5:6).*newCartWeights(3)/2 %locomotive
end 

%% Load Simulation
reactions = [0 0]; % Reaction at A, B respectively. 
for iteration = 1:n_train %step through evenly located train positions 
    %distance from the left
    effectiveLoadLocations = (trainLoadPositions + resolution*(iteration/n_train)); %automatically trim negatives, as these are off the bridge
    
    effectiveLoads = NonNegativeToOne(effectiveLoadLocations).*trainLoads;
    
    % Solving for reaction on far side: 
    reactions(2) = -1*sum(effectiveLoads.*effectiveLoadLocations)/(resolution);
    reactions(1) = -1*sum(effectiveLoads) - reactions(2);
    
    % solving for applied loads
    appliedLoads = zeros(1,resolution+1);
    appliedLoads(1) = reactions(1);
    appliedLoads(end) = reactions(2);
    
    for loadLocation = 1:6 % Iterate through each axle
        if effectiveLoadLocations(loadLocation) < 0 %if load location is negative, and thus off the track
            continue % pass this location as it's not needed
        end
        appliedLoads(effectiveLoadLocations(loadLocation)) = effectiveLoads(loadLocation);
    end 
    appliedLoads
    shearForceDiagrami(iteration,:) = NumericalIntegral(appliedLoads);
    bendingMomentDiagrami(iteration,:) = NumericalIntegral(shearForceDiagrami(iteration,:))
end

%plot(linspace(1,1200,1201),bendingMomentDiagrami(1,:))
hold on
plot(linspace(1,1200,1201),shearForceDiagrami(3,:))
yline(0)
hold off
%% Stress / Potential Failure Mode Simulation