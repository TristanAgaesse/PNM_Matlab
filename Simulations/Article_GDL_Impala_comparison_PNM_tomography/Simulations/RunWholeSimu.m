%Network creation
tic
%load pnm extraction results
InputScriptImageBasedPNM
voxelLength=2.2e-6;
inputContainerMap('VoxelEdgeLength')=voxelLength;
network = CreateImageBasedPoreNetwork(inputContainerMap);


%Simulations

linkDiameter=2*network.GetLinkDataList.('CapillaryRadius');
%linkDiameter=sqrt(voxelLength^2*double(network.GetLinkData('RawData_GeometricSurface'))./pi);
network.AddNewLinkData(linkDiameter,'Diameter');

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

% Holder = 200, Top membrane=180 Fibers=50, Bottom membrane=75
PsiFusionImage=ReadTiff('../3DSamples/PSI_FusionImages_2540.tif');
[fusionlabelEnds,fusionOrderLabels,fusionLabelIndices] = PoreNetworkImageBased.ParseLabeledImage(PsiFusionImage);
fullMorphoImage=ReadTiff('../ResultsFullMorphology/24BA_2540_FullMorpho_theta115.tif');

holder=PoreNetworkImageBased.GetVoxelsOfLabel(200,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
fibers=PoreNetworkImageBased.GetVoxelsOfLabel(50,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
topMembrane=PoreNetworkImageBased.GetVoxelsOfLabel(180,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
bottomMembrane=PoreNetworkImageBased.GetVoxelsOfLabel(75,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);

image_IP_Result(topMembrane)=codeForSolid;
image_IP_Result(holder)=codeForSolid;
image_IP_Result(fibers)=codeForSolid;
image_IP_Result(bottomMembrane)=114;

experimentalImage = PsiFusionImage;
experimentalImage(topMembrane)=codeForSolid;
experimentalImage(holder)=codeForSolid;
experimentalImage(fibers)=codeForSolid;
experimentImage(bottomMembrane)=114;

%Voxel match simu/exp
voxelByVoxelMatch=ComputeVoxelbyVoxelMatch(experimentalImage(85:-1:15,:,:),image_IP_Result(85:-1:15,:,:),pressureCode);

figure
plot(pressureStep,voxelByVoxelMatch)
title('Voxel by voxel match water increments')

%Voxel match simu/exp cumulative
voxelByVoxelError_CumulativePSI=ComputeVoxelbyVoxelMatch_CumulativePSI(experimentalImage(85:-1:15,:,:),image_IP_Result(85:-1:15,:,:),pressureCode);

figure
plot(pressureStep,voxelByVoxelError_CumulativePSI)
title('Voxel by voxel match water distributions')

%Voxel match simu/fullMorpho
fullMorphoImage(topMembrane)=codeForSolid;
fullMorphoImage(holder)=codeForSolid;
fullMorphoImage(fibers)=codeForSolid;
fullMorphoImage(bottomMembrane)=114;

voxelByVoxelError_CumulativeFullMorpho=ComputeVoxelbyVoxelMatch_CumulativePSI(fullMorphoImage(85:-1:15,:,:),image_IP_Result(85:-1:15,:,:),pressureCode);

figure
plot(pressureStep,voxelByVoxelError_CumulativeFullMorpho)
title('Voxel by voxel match PNM agains Full Morphology')


%[satProfPNM,capPressPNM]=ProcessSaturationCurves(image_IP_Result(70:-1:24,:,:),pressureStep,pressureCode,[1,0,0],codeForSolid,' ');
[satProfPNM,capPressPNM]=ProcessSaturationCurves(image_IP_Result(85:-1:15,:,:),pressureStep,pressureCode,[1,0,0],codeForSolid,' ');

satProfExp = open('24BAdrainage_satProf_experimental_slices15-85.fig');
axSatProfExp = get(satProfExp, 'Children');


axSatProfPNM = get(satProfPNM, 'Children');
for i = 1 : numel(axSatProfPNM)
    ax2Children = get(axSatProfPNM(i),'Children');
    copyobj(ax2Children, axSatProfExp(i));
end

capPressExp = open('24BAdrainage_capPressure_experimental_slices15-85.fig');
axcapPressExp = get(capPressExp, 'Children');

axCapPressPNM = get(capPressPNM, 'Children');
for i = 1 : numel(axCapPressPNM)
    ax2Children = get(axCapPressPNM(i),'Children');
    copyobj(ax2Children, axcapPressExp(i));
end

%GetPNMOutputImage


toc