tic

%read Full Morpho result image
imageFullMorpho = ReadTiff('../ResultsFullMorphology/24BA_2540_FullMorpho_theta115.tif');

%parameters for Full Morphology
gamma=72e-3;
voxelLength=2.2e-6;
codeForSolid=255;
pressureStep=[1400,2200,2800,3900,5300];
pressureCode=[114,122,128,139,153];
myContactAngle=115;


% maxBallRadius = 55;  
% ballRadius = 1:4:maxBallRadius;
% pressureStep =  2*gamma./(ballRadius.*voxelLength) ;
% pressureCode =  100+ballRadius;
% myContactAngle=180;





%Post-processing

PsiFusionImage=ReadTiff('../3DSamples/PSI_FusionImages_2540.tif');
[fusionlabelEnds,fusionOrderLabels,fusionLabelIndices] = PoreNetworkImageBased.ParseLabeledImage(PsiFusionImage);

% Holder = 200, Top membrane=180 Fibers=50, Bottom membrane=75
holder         = PoreNetworkImageBased.GetVoxelsOfLabel(200,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
fibers         = PoreNetworkImageBased.GetVoxelsOfLabel(50,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
topMembrane    = PoreNetworkImageBased.GetVoxelsOfLabel(180,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
bottomMembrane = PoreNetworkImageBased.GetVoxelsOfLabel(75,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);

imageFullMorpho(topMembrane) = codeForSolid;
imageFullMorpho(holder)      = codeForSolid;
imageFullMorpho(fibers)      = codeForSolid;

experimentalImage = PsiFusionImage;
experimentalImage(topMembrane)  = codeForSolid;
experimentalImage(holder)       = codeForSolid;
experimentalImage(fibers)       = codeForSolid;
experimentImage(bottomMembrane) = 114;

voxelByVoxelMatch=ComputeVoxelbyVoxelMatch(experimentalImage(85:-1:15,:,:),imageFullMorpho(85:-1:15,:,:),pressureCode);

figure
plot(pressureStep,voxelByVoxelMatch)
title('Voxel by voxel match')

%Voxel match simu/exp cumulative
voxelByVoxelError_CumulativePSI=ComputeVoxelbyVoxelMatch_CumulativePSI(experimentalImage(85:-1:15,:,:),imageFullMorpho(85:-1:15,:,:),pressureCode);

figure
plot(pressureStep,voxelByVoxelError_CumulativePSI)
title('Voxel by voxel match water distributions')



[satProfSimu,capPressSimu]=ProcessSaturationCurves(imageFullMorpho(85:-1:15,:,:),pressureStep,pressureCode,[1,0,0],codeForSolid,' ');

satProfExp = open('24BAdrainage_satProf_experimental_slices15-85.fig');
axSatProfExp = get(satProfExp, 'Children');


axSatProfPNM = get(satProfSimu, 'Children');
for i = 1 : numel(axSatProfPNM)
    ax2Children = get(axSatProfPNM(i),'Children');
    copyobj(ax2Children, axSatProfExp(i));
end

capPressExp = open('24BAdrainage_capPressure_experimental_slices15-85.fig');
axcapPressExp = get(capPressExp, 'Children');

axCapPressPNM = get(capPressSimu, 'Children');
for i = 1 : numel(axCapPressPNM)
    ax2Children = get(axCapPressPNM(i),'Children');
    copyobj(ax2Children, axcapPressExp(i));
end



%Get output Image

imageFullMorpho(topMembrane)        = 180;
imageFullMorpho(holder)             = 200;
imageFullMorpho(fibers)             = 50;
imageFullMorpho(bottomMembrane)     = 75;
imageFullMorpho(imageFullMorpho==0) = 255; %void

toc