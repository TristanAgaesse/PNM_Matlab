function outputInformation = ComputeCondensation( network, options )
%COMPUTECONDENSATION Summary of this function goes here
%   Detailed explanation goes here
%
%Algorithm : - 
%
%Input : network, inletLink,outletLink,options
%       network : pore network on which the algorithm is run
%       options.nIterations= number of iterations
%       options.TemperatureInletLinks : links inlet for heat
%       options.TemperatureOutletLinks : links outlet for heat
%       options.DeltaTemperature : temperature difference between inlet and outlet
%       options.LiquidWaterOutletLinks : links outlet for water
%       options.VaporInletLinks
%       options.VaporOutletLinks
%       options.RelativeHumidityInlet
%       options.RelativeHumidityOutlet
%       options.ClusterOptions : options for the beheavior of liquid phase (see ClusterMonophasique)
%                                  
%Output : outputInformation : information to analyse the degradation process 

%---------------------------------------------------------------------------------------------    


    %Compute temperature field
    
    [temperature,heatTransferCoefficient] = ComputeTemperatureField(network,options.TemperatureInletLinks,options.TemperatureOutletLinks,options.DeltaTemperature)  ;
    
    outputInformation.TemperatureField = temperature;
    outputInformation.HeatTransferCoefficient = heatTransferCoefficient;
    
    
    
    %Compute equilibrum vapor pressure field
    
    equilibriumVaporPressure = ComputeEquilibriumVaporPressure(network,temperature);
    
    outputInformation.EquilibriumVaporPressure = equilibriumVaporPressure;
    
    
    
    %Compute partial pressure field
    
    partialVaporPressure = ComputePartialVaporPressure(network,temperature);
    
    
    
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
     
    for iIteration=1:nIterations

        %Update partial pressure field
    
    
    
        %new invasion : des pores condensables proches d'une zone envahie, 
        %choisir celui qui peut être envahi par IP (Pc la plus faible)  
        [minPressure,indexMinPressureLink] = FindNextInvadedLink
        
        interfaceChangeInformation=cluster.InvadeNewPore(indexMinPressureLink);
        cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink);
    
        
    end
        
        
    
    
%---------------------------------------------------------------------------------------------    
    function  [temperature,heatTransferCoefficient] = ComputeTemperatureField(network,temperatureInlet,temperatureOutlet,deltaTemperature) 
        %Temperature field resulting from a temperature difference between 
        %temperature inlet and temperature outlet
        
        
        boundaryConditions=struct;
        boundaryConditions.inletLink = temperatureInlet;
        boundaryConditions.outletLink = temperatureOutlet;
        boundaryConditions.inletType = 'Dirichlet' ;
        boundaryConditions.outletType = 'Dirichlet' ;
        boundaryConditions.inletValue = deltaTemperature;
        boundaryConditions.outletValue = 0;
        
        voidCluster=network.CreateVoidCluster;
        fullCluster=voidCluster.GetComplementaryCluster;

        [ temperature, ~, ~, heatTransferCoefficient ]=ComputeDiffusion(network,fullCluster, boundaryConditions);
        
        
    end
    



%---------------------------------------------------------------------------------------------
    function equilibriumVaporPressure = ComputeEquilibriumVaporPressure(temperature)
        
        
        
    end


%---------------------------------------------------------------------------------------------        
    function partialVaporPressure = ComputePartialVaporPressure(temperature)

        
        
    end


%---------------------------------------------------------------------------------------------
    function [minPressure,indexMinPressureLink] = FindNextInvadedLink
        %new invasion : des pores condensables proches d'une zone envahie, 
        %choisir celui qui peut être envahi par IP (Pc la plus faible)  
        
        poresCondensables = partialPressure>equilibriumPressure ;
        
        clusterLinkCondensable = clusterLink(voisinCondensable);
        
        criticalPressures  =  cluster.GetCriticalPressures;
        boundaryCriticalPressures = criticalPressures(cluster.GetInterfaceLinks);
        [minPressure,indexMinPressureLink]  =  min(boundaryCriticalPressures(clusterLinkCondensable));
        
        indicesPoresCondensables = find(poresCondensables);
        indexMinPressureLink=indicesPoresCondensables(indexMinPressureLink);
        
    end





end

