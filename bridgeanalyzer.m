clc, clear, clf
close all
%% Data Import
load("bridgeloads.mat") %loads data corresponding to max shear, max bending, and total train weight experienced

%% SImulation Settings
simConfigs = struct;
simConfigs.L = 1200; 
simConfigs.resolution = 2400;
simConfigs.dx = simConfigs.L/simConfigs.resolution; %scaling factor for resolution / distance

%% Material Property Fields
materialProperties = struct(...
    "maxTensileStress",30,... %MPa
    "maxCompressiveStress",6, ... %MPa
    "maxShearStress",4, ... %MPa
    "elasticModulus",4000, ... %MPa
    "poissonRatio",0.2, ... %ratio as per usual
    "glueStrength", 4); % MPa

%% Bridge Geometry Functions

% Wrapper function which provides an array of Structs pertaining to the
% thickness of a glue (for use in tau blah blah) and the corresponding
% depth in the strcture that it is present in. 
function gluePropertyObject = constructGlueObject(glueThickness, glueDepth)
    gluePropertyObject = struct("glueThickness", glueThickness, "glueDepth", glueDepth);
end 

% Wrapper function which provides named fields for objects.
% Height and width are per usual, initial Y value should be measured from
% the top. 
% TODO: Will extend this towards like, shear failure properties
function bridgeGeometryObject = constructBridgeGeometryObject(width, height, initialYValue)
arguments
    width (1,1) double
    height (1,1) double
    initialYValue (1,1) double
end  

    bridgeGeometryObject = struct( ...
        "width",width, ...
        "height",height, ...
        "initialYValue",initialYValue);
end

% Calculates I, ybar for a given cross-section and returns it's values. 
function arr = BridgeNodeConstructor(bridgeGeometry,bridgeGlueGeometry,L,deckThickness,configs,endFlag)  
    arguments
        bridgeGeometry 
        bridgeGlueGeometry 
        L 
        deckThickness
        configs 
        endFlag = false  
    end
    % first, calculate corresponding centroidal depth 
    totalArea = 0;
    numerator = 0; 
    for number = 1:length(bridgeGeometry)
        geometricObject = bridgeGeometry(number);
        objectArea = geometricObject.width * geometricObject.height;
        totalArea = totalArea + objectArea; %get A_total
    
        % A*(y_bar local + y_actual) 
        numerator = numerator + objectArea*(geometricObject.height*0.5 + geometricObject.initialYValue);
    end 
    nodeCentroid = numerator/totalArea;
    
    % calculate corresponding moment of inertia 
    nodeMomentOfInertia = 0; 
    for number = 1:length(bridgeGeometry)
        geometricObject = bridgeGeometry(number);
        objectArea = geometricObject.width * geometricObject.height; 
        objectCentroidHeight = geometricObject.height*0.5 + geometricObject.initialYValue; %assuming rectangles 

        nodeMomentOfInertia = nodeMomentOfInertia ...
            + (geometricObject.width*(geometricObject.height^3))/12; % iterate static (bh^3)/12
        
        nodeMomentOfInertia = nodeMomentOfInertia...
            + objectArea*((objectCentroidHeight - nodeCentroid)^2); %calculate AD^2
    end 

    % find section height 
    sectionHeight = max([bridgeGeometry.height] + [bridgeGeometry.initialYValue]);


    % return array of structs
    effectiveNodeLength = round(configs.resolution * (L / configs.L)) + endFlag;
    arr = repmat(struct( ...
        "nodeMomentOfInertia",nodeMomentOfInertia, ... % Return corresponding term of I for the node
        "nodeCentroid",nodeCentroid, ... % Return corresponding term of ybar of the node from the top
        "sectionHeight",sectionHeight,... % Return corresponding height of the section overall
        "bridgeGeometryObjects", bridgeGeometry,... % Return corresponding objects that went into constructing this object. 
        "bridgeGlueGeometry", bridgeGlueGeometry, ...
        "deckThickness", deckThickness)... % Return applicable glue geometry at some certain point. 
        ,effectiveNodeLength,1); % scale by sim scaling factor
end

% Defaults to traversing downwards, indicate negative direction to traverse upwards. 
function Q = findQEffective(depth,centroidalAxis,geometryObjects, direction)
    % FINDQEFFECTIVE finds the effective First Moment of Area, Q for some set of geometric objects, at some 
    % depth, and in a certain direction from a shear-free surface
    arguments
        depth (1,1) double
        centroidalAxis (1,1) double
        geometryObjects struct
        direction (1,1) double = 1 
        
    end

    % Concept: At the depth, find the values where our object's end less than that value. I.e. find where v_start is less than our thing. 
    effectiveObjects = geometryObjects([]); % 
    for i  = 1:numel(geometryObjects)%iterate through all Geometric Objects
        geometricObject = geometryObjects(i);
        % case in which direction is positive 
        if (geometricObject.initialYValue < depth) & direction == 1 % if our object starts higher up than our depth, then adjust accordingly
            effectiveObjects(end+1) = geometricObject; % not preallocating but geometry should be simple anyways
        elseif (geometricObject.initialYValue+geometricObject.height) > depth & direction == -1 % if our object starts lower than our depth, then adjust accordingly
            effectiveObjects(end+1) = geometricObject; % not preallocating but geometry should be simple anyways
        end
    end

    % compute effective height and centroidal distance 
    Q = 0;
    for i = 1:numel(effectiveObjects)
        effectiveObject = effectiveObjects(i);
        objectBottom = effectiveObject.height + effectiveObject.initialYValue;
        
        % neeeds to be fixed for more complex geometries
        if direction == 1
            if (objectBottom > depth)
                effectiveHeight = depth-effectiveObject.initialYValue;
            else
                effectiveHeight = effectiveObject.height; 
            end
            effectiveCentroidDistance = abs(effectiveHeight*0.5 + effectiveObject.initialYValue-centroidalAxis);
        elseif direction == -1
            if (effectiveObject.initialYValue < depth) 
                effectiveHeight = objectBottom - depth ;
                effectiveCentroidDistance = abs(depth + effectiveHeight*0.5 - centroidalAxis);
            else
                effectiveHeight = effectiveObject.height;
                effectiveCentroidDistance = abs(effectiveObject.initialYValue + effectiveObject.height*0.5 - centroidalAxis); 
            end
        end

        effectiveArea = effectiveHeight * effectiveObject.width;  %find effective area
        
        Q = Q + effectiveArea*effectiveCentroidDistance; %AD
    end 
end 

% backwards conversion ratio
function V = findMaxShearForce(tauMax,I,B,Q)
    V = (tauMax*I*B)/Q;
end 

% forwards converstion function
function tau = findShearStress(V,Q,I,B)
    tau = (V*Q)/(I*B);
end

function M = findLocalBucklingMax(materialProperties,momentOfInertia,k,depth,thickness,width,webHeight,webSpacing)
    % Returns max safe value of M that can be supplied. 
    arguments
        materialProperties;
        momentOfInertia;
        k;
        depth;
        thickness;
        width;
        webHeight = 0 ;
        webSpacing = 0; 
    end
    if webHeight && webSpacing % if WebHeight and WebSpacing Exist, consider that this is now in Case 4

    else % Else, Calculate max moment. 
        ratio = thickness/width;
        M = ((k*(pi^2)*materialProperties.elasticModulus*momentOfInertia)/(12*depth*(1-materialProperties.poissonRatio^2)))*ratio^2;

    end
end

%% Bridge Geometry Construction
% Define the corresponding bridge geometry here, I'd make a function for it
% but eh. 

bridgeProperties = struct; 
bridgeProperties.webWidth = 80; % Width between main webs, used for Shear Analysis 

% Standard Glue Depths 
mainDeckGlueProperties = constructGlueObject([20 20], [1.57*2 30]);

% design 0 as an example: 
% FIRST MEMBER ALWAYS CORRESPONDS TO VERTICAL BEAM
design0BridgeGeometry = [
    constructBridgeGeometryObject(1.27,73.73,1.27); %vertical beams 
    constructBridgeGeometryObject(100,1.27,0);
    constructBridgeGeometryObject(1.27,73.73,1.27); %vertical beams
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(80,1.27,75);
];



design1BridgeGeometry = [
    constructBridgeGeometryObject(1.27,173.73,1.27); %vertical beams
    constructBridgeGeometryObject(100,1.27,0);
    constructBridgeGeometryObject(1.27,173.73,1.27); %vertical beams 
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(80,1.27,175);
    ];

%w,h,y_init
design1properBridgeGeometry = [
    constructBridgeGeometryObject(1.27,98.73,1.27) %side beam
    constructBridgeGeometryObject(1.27,98.73,1.27) %side beam
    constructBridgeGeometryObject(100,1.27,0) %top span
    
    constructBridgeGeometryObject(65,1.27,1.27) %main beam below
    
    constructBridgeGeometryObject(1.27,50,2.54) %right flap vertical beam
    constructBridgeGeometryObject(1.27,50,2.54) %right flap vertical beam
    constructBridgeGeometryObject(14.96,1.27,1.27) %right flap flat
    constructBridgeGeometryObject(14.96,1.27,1.27) %right flap flat
]


% Process corresponding cross-sections as a set of nodes
%bridgeNodeContainer1 = BridgeNodeConstructor(design0BridgeGeometry,mainDeckGlueProperties,600,1.27,simConfigs);
%bridgeNodeContainer2 = BridgeNodeConstructor(design1BridgeGeometry,mainDeckGlueProperties,600,1.27,simConfigs, true);

% Package up nodes into a desired cross-sectional node profile. 
%bridgeNodeContainer = [bridgeNodeContainer1; bridgeNodeContainer2];
bridgeNodeContainer = BridgeNodeConstructor(design1properBridgeGeometry,mainDeckGlueProperties,1200,2.54,simConfigs, true);

% Bridge Member Locations
bridgeDiaphragmLocations = round(linspace(1,simConfigs.resolution+1,10)) %assume evenly spaced


%% Failure Mode Analysis
for nodeNumber = 1:length(bridgeNodeContainer) %For every node within the bridge..
    
    % For each pixel, get properties of that pixel
    % "Individual Bridge Node"
    IBN= bridgeNodeContainer(nodeNumber);

    % First check for bending at the top
    % Since Y-Axis is measured from the top being 0: 
    if maxBendingMomentdata(nodeNumber) > 0 % Check for if top is in tens/comp, positive -> top in comp. 
        bridgeNodeContainer(nodeNumber).maxCompressiveBendingMoment = ...   
        (materialProperties.maxCompressiveStress*IBN.nodeMomentOfInertia)...
        /IBN.nodeCentroid;

        bridgeNodeContainer(nodeNumber).maxTensileBendingMoment = ...   
            (materialProperties.maxTensileStress*IBN.nodeMomentOfInertia)...
            /(IBN.sectionHeight-IBN.nodeCentroid);
    
    elseif maxBendingMomentdata(nodeNumber) < 0 %reversed sign
        bridgeNodeContainer(nodeNumber).maxTensileBendingMoment = ...   
        (materialProperties.maxCompressiveStress*IBN.nodeMomentOfInertia)...
        /IBN.nodeCentroid;

        bridgeNodeContainer(nodeNumber).maxCompressiveBendingMoment = ...   
            (materialProperties.maxTensileStress*IBN.nodeMomentOfInertia)...
            /(IBN.sectionHeight-IBN.nodeCentroid);
    end % does not account for moment = 0, analysis is irrelevant in these locations. 

    
    % Finding the Max Shear Stress experienced for a certain member 
    % Finding Q
    topQ = findQEffective(IBN.deckThickness,IBN.nodeCentroid,IBN.bridgeGeometryObjects,1);
    nodeMaxShearStress = findShearStress(maxShearForceData(nodeNumber),topQ,IBN.nodeMomentOfInertia,IBN.bridgeGeometryObjects(1).width*2);
    bridgeNodeContainer(nodeNumber).nodeMaxShearStress = nodeMaxShearStress;
    
    
    % Checking for Glue Failure at Depth
    % Bridge Glue Geometry Object Index
    bridgeNodeContainer(nodeNumber).maxGlueShearForce = []; %buffer
    for bridgeGGOIndex = 1:length(IBN.bridgeGlueGeometry.glueThickness) % assume lists are equal length
        glueDepth = IBN.bridgeGlueGeometry.glueDepth(bridgeGGOIndex);
        glueThickness = IBN.bridgeGlueGeometry.glueThickness(bridgeGGOIndex);
        QEffective = findQEffective(glueDepth, IBN.nodeCentroid, IBN.bridgeGeometryObjects); % Find effective Q at this height 

        bridgeNodeContainer(nodeNumber).maxGlueShearForce(end+1) =  findMaxShearForce(materialProperties.glueStrength,IBN.nodeMomentOfInertia, glueThickness,QEffective);
    end

    % Checking for Wall Shear Failure 
    wallQ = findQEffective(IBN.nodeCentroid,IBN.nodeCentroid,IBN.bridgeGeometryObjects,1);
    bridgeNodeContainer(nodeNumber).maxWallShearForce = ...
    findMaxShearForce(materialProperties.maxShearStress,IBN.nodeMomentOfInertia,1.27,wallQ);


    %% Wrong, will edit. 
    % Checking for MainSpan - Failure Mode 1 from shear
    
    bridgeNodeContainer(nodeNumber).maxTopLocalMoment = ...
        findLocalBucklingMax(materialProperties,IBN.nodeMomentOfInertia,4,IBN.nodeCentroid,IBN.deckThickness,75); 
    
    
    % Checking for MainSpan - Failure Mode 2 (edges) from shear 
    
    bridgeNodeContainer(nodeNumber).maxTopEdgeLocalMoment = ...
        findLocalBucklingMax(materialProperties,IBN.nodeMomentOfInertia,0.425,IBN.nodeCentroid,IBN.deckThickness,12.5); 
    
    % Checking for Web Member - Failure Mode 3 (Webs) from Shear 
    % first, find vertical unanchored web member, assume geometricobject 1 represents this member. 
    case3TestObject = IBN.bridgeGeometryObjects(1); 
    % find corresponding effective height
    effectiveHeight = IBN.nodeCentroid-case3TestObject.initialYValue; %assume that object is inherently clipped (of course) 
    bridgeNodeContainer(nodeNumber).maxCase3Moment = ...
        findLocalBucklingMax(materialProperties,IBN.nodeMomentOfInertia,6,effectiveHeight,1.27,12.5); 
    

end 

% Shear Thin Plate Buckling 
% Iterate for every region..
% Shear Thin Plate Buckling Constant
sTPBC = ((5*(pi^2)*materialProperties.elasticModulus)/(12*(1-materialProperties.poissonRatio)))

for region = 2:length(bridgeDiaphragmLocations)
    regionStart = bridgeDiaphragmLocations(region-1);
    regionEnd = bridgeDiaphragmLocations(region);
    


    for nodeIndex = regionStart:regionEnd  %update corresponding fail values per node
        IBN = bridgeNodeContainer(nodeIndex); 
        verticalComponent = IBN.bridgeGeometryObjects(1)
        ratio1 = verticalComponent.width/verticalComponent.height; 
        ratio2 = verticalComponent.width/((regionEnd-regionStart)*simConfigs.dx);
        maxShearStress = sTPBC*(ratio1^2 + ratio2^2); 
        bridgeNodeContainer(nodeIndex).maxShearThinPlateBuckleStress = maxShearStress;


    end 
end


%% Plots
x = linspace(0,simConfigs.L,simConfigs.resolution+1);
% Flexural Stress Plots
figure(1) 
hold on
set(gca, 'YDir', 'reverse')
plot(x,maxBendingMomentdata)
plot(x,[bridgeNodeContainer.maxCompressiveBendingMoment])
plot(x,[bridgeNodeContainer.maxTensileBendingMoment])
yline(0)
ax = gca;
ax.YAxis.Exponent = 0;
title("Bending Moment vs Flexural Failure Modes")
legend("Bending Moment","Failure under Compressive Failure","Failure under Tensile Failure")
ylabel("Bending Moment (N*mm)")
xlabel("Distance along Beam (mm)")

hold off

% Shear Stress Plots
% Glue Stress Plots (per location of glue) 
row1Cells = arrayfun(@(node) node.maxGlueShearForce(:, 1), bridgeNodeContainer, 'UniformOutput', false);
maxGlueShearForce1 = [row1Cells{:}];

figure(2) 
hold on
plot(x,abs(maxShearForceData))
plot(x,maxGlueShearForce1.')
yline(0)
ax = gca;
ax.YAxis.Exponent = 0;
legend("Max Shear Force", "Glue Shear Stress Fail 1")
xlabel("Distance along Beam(mm)")
ylabel("Shear Force (N)")
title("Shear Force vs Glue Failure Modes")
hold off

figure(3) 
hold on
plot(x,abs(maxShearForceData))
plot(x,[bridgeNodeContainer.maxWallShearForce])
yline(0)
ax = gca;
ax.YAxis.Exponent = 0;
title("Shear Force vs Wall Shear Failure")
legend("Shear Stress","Wall Shear Failure")
ylabel("Shear Force (N)")
xlabel("Distance along Beam (mm)")
hold off 

figure(4) 
hold on
set(gca, 'YDir', 'reverse')
plot(x,maxBendingMomentdata)
plot(x,[bridgeNodeContainer.maxTopLocalMoment])
yline(0)
ax = gca;
ax.YAxis.Exponent = 0;
title("Bending Moment vs Main Deck Local Buckling(?)")
legend("Bending Moment","Main Deck Local Buckling Failure")
ylabel("Bending Moment (N*mm)")
xlabel("Distance along Beam (mm)")

hold off

figure(5) 
hold on
set(gca, 'YDir', 'reverse')
plot(x,maxBendingMomentdata)
plot(x,[bridgeNodeContainer.maxTopEdgeLocalMoment])
yline(0)
ax = gca;
ax.YAxis.Exponent = 0;
title("Bending Moment vs Main Deck (Edge) Local Buckling(?)")
legend("Bending Moment","Main Deck (Edge) Local Buckling(?)Edge Local Shear Failure")
ylabel("Bending Moment (N*mm)")
xlabel("Distance along Beam (mm)")

hold off

figure(6) 
hold on
set(gca, 'YDir', 'reverse')
plot(x,maxBendingMomentdata)
plot(x,[bridgeNodeContainer.maxCase3Moment])
yline(0)
ax = gca;
ax.YAxis.Exponent = 0;
title("Bending Moment vs Case 3 Max Moment")
legend("Bending Moment","Web Compression Failure *for one section* (Case 3)")
ylabel("Bending Moment (N*mm)")
xlabel("Distance along Beam (mm)")

hold off
figure(7) 
hold on
plot(x,abs([bridgeNodeContainer.nodeMaxShearStress]))
plot(x,[bridgeNodeContainer.maxShearThinPlateBuckleStress])
yline(0)
ax = gca;
title("Shear Stress vs max Case 4 Shear Stress")
legend("Max Shear Stress","Thin Plate Buckling Max Shear Stress")
ylabel("Shear Stress (MPa)")
xlabel("Distance along Beam (mm)")

hold off


% Stress Data Display
x_space = linspace(1,simConfigs.L,resolution+1);
figure(10)
hold on
plot(x_space,maxShearForceData)
yline(0)
legend("Shear Force")
ylabel("Shear Force (N)")
xlabel("Distance along Beam (mm)")
hold off


figure(11)
hold on
plot(x_space,maxBendingMomentdata)
set(gca, 'YDir', 'reverse')
yline(0)
legend("Bending Moment")
ylabel("Bending Moment (N*mm)")
xlabel("Distance along Beam (mm)")
hold off
%% Factor of Safety Calculations 
%figure(3) 
%hold on

%hold off 
