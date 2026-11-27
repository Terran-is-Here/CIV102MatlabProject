clc, clear, clf, close all %clear code parameters as per usual
%% Documentation
% In terms of the bridge, assume it takes the following structure in terms
% of notation: 

% A=========================================B
% Where A is the left-most endpoint, and B is the right-most endpoint. 
% This analyzer solves from traversing a train from left TO right. 

%% Todo
% Rewrite Load Simulation to accept arbitrary assumptions for resolution. 
%% Simulation Parameters
resolution = 2400; %break bridge up into 1200 segments, resolution of 1mm per node 
L = 1200; %actual bridge length
n_train = 2*resolution; % train locations present, evenly split across span + bridge unloading 
dx = L/resolution; %Convertion ratio from realdistance / nodes 

shearForceDiagrami = zeros(n_train,resolution+1);
bendingMomentDiagrami = zeros(n_train, resolution+1);

% Material Characteristics 
materialThickness = 1.27;  %in mm

% Train Loading Characteristics

totalTrainWeight = 400; % in terms of Newtons

cartWeights = ones(1,3).*totalTrainWeight/3; % Assume that the cartWeights is arranged in Cart 2, Cart 1, Locomotive. 



% in terms of millimeters, adjusts based on ratio of resolution and
% distance
trainLoadPositions = [-856 -680 -516	-340	-176	0]/dx; %right-most value is locomotive

 trainLoads = ones(1,6).*totalTrainWeight/6*-1; %initial load is fully uniformly distributed


% Custom Weight Section
% [totalTrainWeight, trainLoads] =updateTrainLoads([135,135,482]);
%% Bridge Geometry
%Parameters are [x of cross-section (or cross-sectional change location),
%etc. 
params = [0,0,0,0;
    0,0,0,0;
    0,0,0,0];


%% Functions

% numerical integral
function out = NumericalIntegral(arr, dx)
    arguments
    arr double  %data array to "integrate" over
    dx double = 1 %scaling value 
    end
    
    out = zeros(1,size(arr,2));
    currentvalue = 0;
    for index = 1:size(arr,2)
        currentvalue  = currentvalue + arr(index)*dx; %scale by resolution if required 
        out(index) = currentvalue;
    end
    
end

function arr = ReplaceNegatives(arr) 
    arr = max(arr,0); %rep
end

function mask = TrimLoads(nodePositions, resolution) % Trims loads to not be effective if they are OFF the bridge
    arguments
        nodePositions double
        resolution double = 1200 
    end

    % If nodes are in [0,resolution] -> 1
    % else, not included (i.e. off the bridge) 
    mask = double((nodePositions >=0 ) & (nodePositions <= resolution));
    
   
end

% Take in a new set of cart weights, update axle loads and total train
% weight respectively. 
function [newTotalTrainWeight, newTrainLoads] =  updateTrainLoads(newCartWeights)
    newTotalTrainWeight = sum(newCartWeights,2); 
    newTrainLoads = zeros(1,6);
    newTrainLoads(1:2) = newCartWeights(1) / -2; % cart 2
    newTrainLoads(3:4) = newCartWeights(2) / -2; % cart 1
    newTrainLoads(5:6) = newCartWeights(3) / -2; % locomotive
end 

%% Load Simulation
reactions = [0 0]; % Reaction at A, B respectively. 
for iteration = 1:n_train %step through each node, where iteration is the node at a point
    
    %distance from the left, scaled by resolution size 
    effectiveLoadLocations = (trainLoadPositions + iteration); % move train right by a node
    effectiveLoads = TrimLoads(effectiveLoadLocations,resolution).*trainLoads; %adds loads
    
    % Solving for reaction on far side: 
    reactions(2) = -1*sum(effectiveLoads.*effectiveLoadLocations)/(resolution); %considers resolution as distance basically
    reactions(1) = -1*sum(effectiveLoads) - reactions(2); 
    % solving for applied loads
    appliedLoads = zeros(1,resolution+1);
    appliedLoads(1) = reactions(1);
    appliedLoads(end) = reactions(2);
    
    for loadLocation = 1:6 % Iterate through each axle
        if effectiveLoadLocations(loadLocation) < 0 %if load location is negative, and thus off the track
            continue % pass this location as it's not needed
        end
        if effectiveLoadLocations(loadLocation) > resolution % if load location is off the bridge, do not include
            continue %pass this location as it's not needed. 
        end
        %step through all axle loads and put it on the beam
        loadPosition = round(effectiveLoadLocations(loadLocation))+1;
        appliedLoads(loadPosition) = appliedLoads(loadPosition)+ effectiveLoads(loadLocation);
    end 
    shearForceDiagrami(iteration,:) = NumericalIntegral(appliedLoads);
    
    %scale moment integral by scaling factor (as each step is equivalent to
    %L/resolution in terms of distance )
    bendingMomentDiagrami(iteration,:) = NumericalIntegral(shearForceDiagrami(iteration,:),dx); %scale by factor
end

% note that we only really care about the magnitude of the shear force
% present for our member, not the direction 
maxShearForceData = max(shearForceDiagrami,[],"ComparisonMethod","abs"); 
maxBendingMomentdata = max(bendingMomentDiagrami,[],"ComparisonMethod","abs");
%plot(linspace(1,1200,1201),bendingMomentDiagrami(1,:))


x_space = linspace(1,L,resolution+1);
figure(1)
hold on
plot(x_space,maxShearForceData)
yline(0)
legend("Shear Force")
ylabel("Shear Force (N)")
xlabel("Distance along Beam (mm)")
hold off


figure(2)
hold on
plot(x_space,maxBendingMomentdata)
set(gca, 'YDir', 'reverse')
yline(0)
legend("Bending Moment")
ylabel("Bending Moment (N*mm)")
xlabel("Distance along Beam (mm)")
hold off


save("bridgeloads.mat","maxShearForceData", "maxBendingMomentdata", "totalTrainWeight","resolution")


%% Stress / Potential Failure Mode Simulation