
%Add information regarding the experimental set up to the image which
%contains the water distribution simulated with PNM. This information is
%taken from PsiFusionImage


PsiFusionImage=ReadTiff('../3DSamples/PSI_FusionImages_2540.tif');

[labelEnds,orderLabels,labelIndices] = PoreNetworkImageBased.ParseLabeledImage(PsiFusionImage);


PNMOutputImage = image_IP_Result;



 

%Holder = 200, Top membrane=180

% holderAndTopMembrane= PsiFusionImage>179;
% PNMOutputImage(holderAndTopMembrane)=PsiFusionImage(holderAndTopMembrane);

linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(200,labelEnds,orderLabels,labelIndices);
PNMOutputImage(linearIndices)=200;

linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(180,labelEnds,orderLabels,labelIndices);
PNMOutputImage(linearIndices)=180;


%Fibers=50, Bottom membrane=75

% fibersAndBottomMembrane = PsiFusionImage<76;
% PNMOutputImage(fibersAndBottomMembrane)=PsiFusionImage(fibersAndBottomMembrane);

linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(50,labelEnds,orderLabels,labelIndices);
PNMOutputImage(linearIndices)=50;

linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(75,labelEnds,orderLabels,labelIndices);
PNMOutputImage(linearIndices)=75;

%Void
PNMOutputImage(PNMOutputImage==0)=255;

