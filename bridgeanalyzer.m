clc, clear, clf

%% Data Import
load("bridgeloads.mat") %loads data corresponding to max shear, max bending, and total train weight experienced

%% SImulation Settings
L = 1200; 
resolution = 1200;

%% Material Property Fields
materialProperties = struct(...
    "maxTensileStress",30,...
    "maxCompressiveStress",6, ...
    "maxShearStress",4, ...
    "elasticModulus",4000, ...
    "poissonRatio",0.2, ...
    "glueStrength", 4);
% All in terms of MPa. 

%% Bridge Geometry Functions

% Wrapper function which provides named fields for objects.
% Height and width are per usual, initial Y value should be measured from
% the top. 
% TODO: Will extend this towards like, shear failure properties
function bridgeGeometryObject = constructBridgeGeometryObject(width, height, initialYValue,glueYValue,glueArea)
arguments
    width 
    height 
    initialYValue 
    glueYValue  =0
    glueArea  =0
end  

    bridgeGeometryObject = struct( ...
        "width",width, ...
        "height",height, ...
        "initialYValue",initialYValue, ...
        "glueYValue",glueYValue, ...
        "glueArea", glueArea);
end

% Calculates I, ybar for a given cross-section and returns it's values. 
function arr = BridgeNodeConstructor(bridgeGeometry,L)  
    bridgeGeometry
    % first, calculate corresponding centroidal depth 
    totalArea = 0;
    numerator = 0; 
    for number = 1:length(bridgeGeometry)
        geometricObject = bridgeGeometry(number)
        objectArea = geometricObject.width * geometricObject.height;
        totalArea = totalArea + objectArea; %get A_total
    
        % A*(y_bar local + y_actual) 
        numerator = numerator + objectArea*(geometricObject.height*0.5 + geometricObject.initialYValue);
    end 
    nodeCentroid = numerator/totalArea;
    
    % calculate corresponding moment of inertia 
    nodeMomentOfInertia = 0; 
    for number = 1:length(bridgeGeometry)
        geometricObject = bridgeGeometry(number)
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
    arr = repmat(struct( ...
        "nodeMomentOfInertia",nodeMomentOfInertia, ...
        "nodeCentroid",nodeCentroid, ...
        "sectionHeight",sectionHeight),L,1);
end


%% Bridge Geometry Field 
% Define the corresponding bridge geometry here, I'd make a function for it
% but eh. 

% design 0 as an example: 
design0BridgeGeometry = [
    constructBridgeGeometryObject(100,1.27,0);
    constructBridgeGeometryObject(1.27,73.73,1.27);
    constructBridgeGeometryObject(1.27,73.73,1.27); %vertical beams
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(80,1.27,75);
];

bridgeNodeContainer1 = BridgeNodeConstructor(design0BridgeGeometry,600);

design1BridgeGeometry = [
    constructBridgeGeometryObject(100,1.27,0);
    constructBridgeGeometryObject(1.27,173.73,1.27);
    constructBridgeGeometryObject(1.27,173.73,1.27); %vertical beams
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(5,1.27,1.27);
    constructBridgeGeometryObject(80,1.27,175);
    ];

bridgeNodeContainer2 = BridgeNodeConstructor(design1BridgeGeometry,600+1);


bridgeNodeContainer = [bridgeNodeContainer1; bridgeNodeContainer2]

for nodeNumber = 1:length(bridgeNodeContainer);
    nodeNumber
    % For each pixel, get properties
    IBN= bridgeNodeContainer(nodeNumber);

    % First check for bending at the top
    % Since Y-Axis is measured from the top being 0: 
    if maxBendingMomentdata(nodeNumber) > 0 % Check for if top is in tens/comp, positive -> top in comp. 
        bridgeNodeContainer(nodeNumber).maxCompressiveBendingMoment = ...   
        (materialProperties.maxCompressiveStress*IBN.nodeMomentOfInertia)...
        /IBN.nodeCentroid

        bridgeNodeContainer(nodeNumber).maxTensileBendingMoment = ...   
            (materialProperties.maxTensileStress*IBN.nodeMomentOfInertia)...
            /(IBN.sectionHeight-IBN.nodeCentroid)
    end

    % Check for glue failure 
    % First, get subarray of 



end 


hold on
plot(linspace(1,L,resolution+1),maxBendingMomentdata)
plot(linspace(1,L,resolution+1),[bridgeNodeContainer.maxCompressiveBendingMoment])
plot(linspace(1,L,resolution+1),[bridgeNodeContainer.maxTensileBendingMoment])
set(gca, 'YDir', 'reverse')
yline(0)
legend("Bending Moment","Failure under Compressive Failure","Failure under Tensile Failure")
ylabel("Bending Moment (Nmm)")
xlabel("Distance along Beam (mm)")
hold off