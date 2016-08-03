
tic

%% Network creation
%fileName = 'pnm_extraction_results.mat'
voxelLength=2.2e-6;
network=CreateNetworkImageFromPythonData(fileName,voxelEdgeLength);


%% Simulations

linkDiameter=2*network.GetLinkDataList.('CapillaryRadius');
%linkDiameter=sqrt(voxelLength^2*double(network.GetLinkData('RawData_GeometricSurface'))./pi);
network.AddNewLinkData(linkDiameter,'Diameter');


pureContactAngle=[110,130];
linkCassieBaxterContactAngle = LocalScaleComputeCassieBaxterContactAngle(network,pureContactAngle);


%Pressure parameters for PNM
pressureStep=[1400,2200,2800,3900,5300];
pressureCode=[114,122,128,139,153];
myContactAngle=115;

%Pressure parameters for Full Morphology
% gamma=72e-3;
% maxBallRadius = 55;  
% ballRadius = 1:4:maxBallRadius;
% pressureStep =  2*gamma./(ballRadius.*voxelLength) ;
% pressureCode =  100+ballRadius;
% myContactAngle=180;


clusterOptions.ContactAngle=myContactAngle;
clusterOptions.ThroatPressure='LaplaceCylinder';
image_IP_Result=RunDrainage_IP_PSI(network,pressureStep,pressureCode,clusterOptions);

%Post-processing

codeForSolid=255;

[satProfPNM,capPressPNM]=ProcessSaturationCurves(image_IP_Result(end:-1:1,:,:),pressureStep,pressureCode,[1,0,0],codeForSolid,' ');


satProfExp = open('24BAdrainage_satProf_experimental.fig');
axSatProfExp = get(satProfExp, 'Children');
capPressExp = open('24BAdrainage_capPressure_experimental.fig');
axcapPressExp = get(capPressExp, 'Children');


axSatProfPNM = get(satProfPNM, 'Children');
for i = 1 : numel(axSatProfPNM)
    ax2Children = get(axSatProfPNM(i),'Children');
    copyobj(ax2Children, axSatProfExp(i));
end

axCapPressPNM = get(capPressPNM, 'Children');
for i = 1 : numel(axCapPressPNM)
    ax2Children = get(axCapPressPNM(i),'Children');
    copyobj(ax2Children, axcapPressExp(i));
end

GetPNMOutputImage


