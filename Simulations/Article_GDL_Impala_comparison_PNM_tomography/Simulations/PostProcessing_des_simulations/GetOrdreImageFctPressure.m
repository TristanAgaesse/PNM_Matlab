[cluster,breakthroughPressure,invasionPressureList]=ComputeInvasionPercolation(network,inletLink,[],'hydrophobic');
invadedPores=cluster.GetInvadedPores;

nStep=50;
pressureStep=quantile(invasionPressureList,nStep*5);
pressureStep=pressureStep(1:nStep);
invadedPoresFctPressure=zeros(1,network.GetNumberOfPores);

lastfirstPore=1;
for iStep=1:nStep
    firstPore=find(invasionPressureList>pressureStep(iStep),1)-1;
    invadedPoresFctPressure(invadedPores(lastfirstPore:firstPore))=iStep;
    lastfirstPore=firstPore;
end
network.AddNewPoreData(invadedPoresFctPressure,'ordreInvasionFctPressure')
ordreImagePressure=network.GetImagePoreData('ordreInvasionFctPressure');
ordreImagePressure(ordreImagePressure==0)=100;