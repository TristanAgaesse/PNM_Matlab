function [cluster,outputInformation] = ComputeCondensation( network, options )
%COMPUTECONDENSATION Summary of this function goes here
%   Detailed explanation goes here
%
%Algorithm : - Demarrage de la condensation dans le pore ayant l'humidité
%                   relative la plus élevée
%            - parmis les pores condensables proches du cluster de condensation, 
%                   envahir celui qui peut être envahi avec la pressions 
%                   capillaire la plus faible  
%
%Input : network, options
%       network : pore network on which the algorithm is run
%       options.TemperatureInletLinks : links inlet for heat
%       options.TemperatureOutletLinks : links outlet for heat
%       options.TemperatureInlet : temperature inlet (in K)
%       options.TemperatureOutlet : temperature outlet (in K)
%       options.LiquidWaterOutletLinks : links outlet for water
%       options.AirPressure : uniform pressure imposed to air
%       options.VaporInletLinks
%       options.VaporOutletLinks
%       options.RelativeHumidityInlet
%       options.RelativeHumidityOutlet
%       options.ClusterOptions : options for the beheavior of liquid phase (see ClusterMonophasique)
%                                  
%Output : [cluster,outputInformation]
%       cluster : liquid water cluster resulting from condensation
%       outputInformation : information to analyse the degradation process 

%---------------------------------------------------------------------------------------------    


    %Compute temperature field
    
    [temperature,heatTransferCoefficient] = ComputeTemperatureField(network,options.TemperatureInlet,options.TemperatureOutlet,options.TemperatureInletLinks,options.TemperatureOutletLinks) ;
    
    outputInformation.TemperatureField = temperature;
    outputInformation.HeatTransferCoefficient = heatTransferCoefficient;
    
    
    
    %Compute equilibrum vapor pressure field
    
    equilibriumVaporPressure = ComputeEquilibriumVaporPressure(network,temperature);
    
    outputInformation.EquilibriumVaporPressure = equilibriumVaporPressure;
    
    
    
    %Compute partial pressure field
    voidCluster=network.CreateVoidCluster;
    fullCluster=voidCluster.GetComplementaryCluster;
    inletVaporPressure = options.RelativeHumidityInlet*options.AirPressure;
    outletVaporPressure = options.RelativeHumidityOutlet*options.AirPressure;
    partialVaporPressure = ComputePartialVaporPressure(network,fullCluster,inletVaporPressure,outletVaporPressure,options.VaporInletLinks,options.VaporOutletLinks,options.AirPressure);
    
    
    
    %invade the pore which has the max partial pressure if > equilibrium
    %partial pressure
    
    condensationRatio = partialVaporPressure ./ equilibriumVaporPressure ;
    [maxRatio,indexMaxRatio] = max(condensationRatio);
    
    if maxRatio>1
        firstInvadedPore=indexMaxRatio;
    else
        return
    end
    
    %Initialise cluster
    
    clusterOptions
    cluster = ClusterMonophasique.InitialiseCondensationCluster(network,clusterOptions,firstInvadedPore,options.LiquidWaterOutletLinks);
    
    
    outlet_reached = false;
    outletPores = network.GetPoresFrontiere(outletLink);
    
    %Find accessible pores
    nPoreAccessible=FindNumberOfAccessiblePores(network,inletLink);
    
    
    while not(outlet_reached) && time<nPoreAccessible
        
        %Update partial pressure field
        
        partialVaporPressure = ComputePartialVaporPressure(network,cluster.GetComplementaryCluster,inletVaporPressure,outletVaporPressure,options.VaporInletLinks,options.VaporOutletLinks,options.AirPressure);
        
        %new invasion : des pores condensables proches d'une zone envahie, 
        %choisir celui qui peut être envahi par IP (Pc la plus faible)  
        [minPressure,indexMinPressureLink] = FindNextInvadedLink(cluster,partialVaporPressure,equilibriumPressure);
        
        
        interfaceChangeInformation=cluster.InvadeNewPore(indexMinPressureLink);
        cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink);
        
        %verifier si outlet_reached
        if ismember(invadedPore,outletPores)
            outlet_reached = true;
            breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),outletLink);
            cluster.InvadeOutletLink(breakthroughLinks);
        end
        
        
    end
        
        
end




%---------------------------------------------------------------------------------------------    
function  [temperature,heatTransferCoefficient] = ComputeTemperatureField(network,temperatureInlet,temperatureOutlet,temperatureInletLinks,temperatureOutletLinks) 
    %Temperature field resulting from a temperature difference between 
    %temperature inlet and temperature outlet

    boundaryConditions=struct;
    boundaryConditions.inletLink = temperatureInletLinks;
    boundaryConditions.outletLink = temperatureOutletLinks;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    boundaryConditions.inletValue = temperatureInlet;
    boundaryConditions.outletValue = temperatureOutlet;

    voidCluster=network.CreateVoidCluster;
    fullCluster=voidCluster.GetComplementaryCluster;

    [ temperature, ~, ~, heatTransferCoefficient ]=ComputeDiffusion(network,fullCluster, boundaryConditions);


end




%---------------------------------------------------------------------------------------------
function equilibriumVaporPressure = ComputeEquilibriumVaporPressure(temperature)
    %http://fr.wikipedia.org/wiki/Pression_de_vapeur_saturante

    %Rankine formula
    equilibriumVaporPressure = 1e5*exp(13.7-5120/temperature) ;

end



%---------------------------------------------------------------------------------------------        
function partialVaporPressure = ComputePartialVaporPressure(network,cluster,inletVaporPressure,outletVaporPressure,vaporInletLinks,vaporOutletLinks,airPressure)
    
    %Convert inletVaporPressure,outletVaporPressure to concentration
    R = 8.314 ;
    airConcentration = airPressure/(R*temperature);
    
    inletConcentration = inletVaporPressure*airConcentration ;%airConcentration(vaporInletLinks);
    outletConcentration = outletVaporPressure*airConcentration;%airConcentration(vaporOutletLinks);
    
    %Compute diffusion of vapor
    boundaryConditions=struct;
    boundaryConditions.inletLink = vaporInletLinks;
    boundaryConditions.outletLink = vaporOutletLinks;
    boundaryConditions.inletType = 'Dirichlet' ;
    boundaryConditions.outletType = 'Dirichlet' ;
    boundaryConditions.inletValue = inletConcentration;
    boundaryConditions.outletValue = outletConcentration;
    
    waterConcentration = ComputeDiffusion(network,cluster, boundaryConditions);    %TODO : change conductance to ln for multicomponent diffusion
    
    %Convert back concentrations to vapor pressure
    
    partialVaporPressure = waterConcentration./airConcentration ;
    
end




%---------------------------------------------------------------------------------------------
function [minPressure,indexMinPressureLink] = FindNextInvadedLink(cluster,partialVaporPressure,equilibriumVaporPressure)
    %new invasion : des pores condensables proches d'une zone envahie, 
    %choisir celui qui peut être envahi par IP (Pc la plus faible)  

    
    %List of condensable pores next to cluster
    poresCondensables = partialVaporPressure>equilibriumVaporPressure ;
    
    clusterLink = cluster.GetInterfaceLinks;
    outwardPores = cluster.GetOutwardPore(1:length(clusterLink));
    voisinCondensable = poresCondensables(outwardPores);
    clusterLinkCondensable = clusterLink(voisinCondensable);

    %Chose condensable pore accessible with min capillary pressure
    criticalPressures  =  cluster.GetCriticalPressures;
    boundaryCriticalPressures = criticalPressures(cluster.GetInterfaceLinks);
    [minPressure,indexMinPressureLink]  =  min(boundaryCriticalPressures(clusterLinkCondensable));

    indicesPoresCondensables = find(poresCondensables);
    indexMinPressureLink=indicesPoresCondensables(indexMinPressureLink);

end



%---------------------------------------------------------------------------------------------
function nPoreAccessible=FindNumberOfAccessiblePores(network,inletLink)

    fooCluster=network.CreateVoidCluster;
    totalFloodCluster=fooCluster.GetComplementaryCluster;
    percoPath=totalFloodCluster.FindPercolationPath(inletLink,1:network.GetNumberOfLinks);
    nPoreAccessible=0;
    for i=1:length(percoPath)
        nPoreAccessible=nPoreAccessible+length(percoPath{i}.GetInvadedPores);
    end

end



