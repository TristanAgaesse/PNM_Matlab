function image_IP_Result=RunDrainage_IP_PSI(network,pressureStep,pressureCode,clusterOptions)
%Script to run invasion percolation on PSI images. 
%Produces an image containing water distribution as a function of pressure.


    contactAngle=clusterOptions.ContactAngle;
    contactAngleList = contactAngle*(pi/180)*ones(1,network.GetNumberOfLinks);
    network.AddNewLinkData(contactAngleList,'ContactAngle');
%     clusterOptions.ThroatPressure = 'LaplaceCylinder';
    clusterOptions.Coalescence = 'none' ;
    clusterOptions.SurfaceTension=72e-3;
    [cluster,breakthroughPressure,invasionPressureList]=ComputeInvasionPercolation(network,network.GetLinksFrontiere(2),[],'currentWettability',clusterOptions);


    disp('Post-processing : image water distribution as a function of pressure')

    invadedPores=cluster.GetInvadedPores;
    invadedPoresFctPressure=zeros(1,network.GetNumberOfPores);

    nStep=length(pressureStep);
    lastfirstPore=1;
    for iStep=1:nStep
        firstPore=find(invasionPressureList>pressureStep(iStep),1);
        invadedPoresFctPressure(invadedPores(lastfirstPore:firstPore))=pressureCode(iStep);
        lastfirstPore=firstPore;
    end
    network.AddNewPoreData(invadedPoresFctPressure,'ordreInvasionFctPressure');


    disp('Post-processing : building the image water distribution')
    image_IP_Result = uint8(network.GetImagePoreData('ordreInvasionFctPressure'));
%     solid = network.MaterialImage>0 ;
%     image_IP_Result(solid) = 255 ;

end

