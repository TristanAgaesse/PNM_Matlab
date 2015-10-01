


experimentImage=ReadTiff('../3DSamples/PSI_FusionImages_2540.tif');

codeForSolid=255;
pressureStep=[1400,2200,2800,3900,5300];
pressureCode=[114,122,128,139,153];

[labelEnds,orderLabels,labelIndices] = PoreNetworkImageBased.ParseLabeledImage(experimentImage);

%Holder = 200, Top membrane=180, Fibers=50
linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(200,labelEnds,orderLabels,labelIndices);
experimentImage(linearIndices)=codeForSolid;

linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(180,labelEnds,orderLabels,labelIndices);
experimentImage(linearIndices)=codeForSolid;

linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(50,labelEnds,orderLabels,labelIndices);
experimentImage(linearIndices)=codeForSolid;

%, Bottom membrane=75, filled with water

linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(75,labelEnds,orderLabels,labelIndices);
experimentImage(linearIndices)=114;


%[satProfPNM,capPressPNM]=ProcessSaturationCurves(image_IP_Result(70:-1:24,:,:),pressureStep,pressureCode,[1,0,0],codeForSolid,' ');
[satProfPNM,capPressPNM]=ProcessSaturationCurves(experimentImage(85:-1:15,:,:),pressureStep,pressureCode,[1,0,0],codeForSolid,' ');