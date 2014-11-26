function capillaryPressureCurve = ComputeCapillaryPressureCurve(network,numBoundaryInlet,wettability,varargin)
%COMPUTECAPILLARYPRESSURECURVE Calcule la courbe de pression capillaire
% input : - network
%           - numBoundaryInlet
%           - wettability : 'currentWettability', 'hydrophobic', 'hydrophilic' or 'random'.
%           - varargin (optionnel) : clusterOptions (voir ClusterMonophasique)
% output : capillaryPressureCurve (saturation en fonction de la pression)
    
    
    [cluster,breakthroughPressure,invasionPressureList]=ComputeInvasionPercolation(network,network.GetLinksFrontiere(numBoundaryInlet),[],wettability,varargin);
    
    invadedPores=cluster.GetInvadedPores;
    
    minPressure=min(invasionPressureList);
    maxPressure=max(invasionPressureList);
    
    nStep=500;
    pressureStep=minPressure:(maxPressure-minPressure)/(nStep-1):maxPressure;
    
    capillaryPressureCurve=horzcat(transpose(pressureStep),zeros(nStep,1));
    
    for iStep=1:nStep
       
        firstPore=find(invasionPressureList>pressureStep(iStep),1)-1;
        fooCluster=network.CreateVoidCluster;
        fooCluster.InvadedPores=invadedPores(1:firstPore);
        capillaryPressureCurve(iStep+1,2)=ComputeSaturation(fooCluster,network,'totalSaturation');

    end


end

