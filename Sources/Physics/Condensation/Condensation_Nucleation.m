function [ nucleationClusters, nucleationInfos ] = Condensation_Nucleation(network,...
                                        equilibriumVaporPressure,options,diffusionConductances,temperature)
%CONDENSATION_NUCLEATION Nucleation step of the condensation algorithm.
%   Condensation occurs when RH>1 in at least one pore. In that case, the 
%   nucleation pores are found with the following algorithm.
%   The first condensation pore is the one with the maximum RH. We invade 
%   this pore and impose RH=1 at his links. Then we update the RH field,
%   solving the vapor diffusion equation with these new boundary
%   conditions. If somes pores have RH>1, we invade the pores with the
%   maximum RH, then update RH and so on. The nucleation ends when RH<1 in 
%   all remaining pores. 
    
    %% Initialisation of the algorithm
    superSaturation = 1; % critical value of RH for condensation
    
    nucleation = true;
    nPore = network.GetNumberOfPores;
    invadedPore = false(nPore,1);
    
    while nucleation
        %% Compute relative humidity in each gaz pore
        
        %boundary conditions for diffusion
        [gasTransportPores,inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks] = Condensation_GetBoundaryConditionsForDiffusion(network,...
                                        options,invadedPore,equilibriumVaporPressure);
        
        
        partialVaporPressure = Condensation_ComputePartialVaporPressure(network,...
                gasTransportPores,diffusionConductances,...
                inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks,options.AirPressure,temperature);
        
        nucleationInfos.PartialVaporPressure{1} = partialVaporPressure;
        
        
        
        condensationRatio = partialVaporPressure ./ equilibriumVaporPressure ;
        
        %% invade the pore which has the max partial pressure if > equilibrium
        
        [maxRatio,indexMaxRatio] = max(condensationRatio);
        
        if maxRatio>superSaturation
            nucleation = true;
            invadedPore(indexMaxRatio)=true;
        else
            nucleation = false;
        end
        
        %partial pressure
        
    end

    %% Define nucleationClusters from invadedPore
    conComp = FindComposantesConnexes(network,find(invadedPore));
    
    nCluster = length(conComp);
    nucleationClusters = cell(nCluster,1);
    for i=1:nCluster
        clusterPores=conComp{i};
        nucleationClusters{i} = ClusterMonophasique.InitialiseCondensationCluster(network,...
                options.ClusterOptions,...
                clusterPores,options.LiquidWaterOutletLinks);
    end
    
end

