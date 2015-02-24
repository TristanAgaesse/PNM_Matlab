function capillaryPressureCurve = ComputeCapillaryPressureCurve(network,inletLink,wettability,clusterOptions)
%COMPUTECAPILLARYPRESSURECURVE Calcule la courbe de pression capillaire
% input : - network
%           - numBoundaryInlet
%           - wettability : 'currentWettability', 'hydrophobic', 'hydrophilic' or 'random'.
%           - clusterOptions (voir ClusterMonophasique)
% output : capillaryPressureCurve (saturation en fonction de la pression)
    
    disp('Computing Capillary Pressure Curve')
    tic;
    
    [cluster,~,invasionPressureList]=ComputeInvasionPercolation(network,inletLink,[],wettability,clusterOptions);
    
    invadedPores=cluster.GetInvadedPores;
    
    minPressure=min(invasionPressureList);
    maxPressure=max(invasionPressureList);
    
    nStep=1000;
    pressureStep=minPressure:(maxPressure-minPressure)/(nStep):maxPressure;
    pressureStep=pressureStep(1:end-1);
    
    capillaryPressureCurve=horzcat(transpose(pressureStep),zeros(nStep,1));
    
    for iStep=1:nStep
       
        firstPore=find(invasionPressureList>pressureStep(iStep),1)-1;
        fooCluster=network.CreateVoidCluster;
        fooCluster.InvadedPores=invadedPores(1:firstPore);
        options.type ='totalSaturation';
        capillaryPressureCurve(iStep,2)=ComputeSaturation(fooCluster,network,options);

    end

    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Computing Capillary Pressure Curve finished. Time spent : %d minutes %f s. \n',minutes,secondes);
end

