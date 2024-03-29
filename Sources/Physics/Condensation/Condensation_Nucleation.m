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
    superSaturation = 1.1; % critical value of RH for condensation
    
    nucleation = true;
    nPore = network.GetNumberOfPores;
    invadedPore = false(nPore,1);
    step = 0;
    nucleationInfos.PartialVaporPressure={};
    nucleationInfos.MaxRH={};
    nucleationInfos.InvadedPore={};
    nucleationInfos.EffectiveDiffusion={};
    nucleationInfos.WaterConcentration={};
    nucleationInfos.WaterDiffusionFlux={};
    
    while nucleation
        %% Compute relative humidity in each gaz pore
        step = step+1;
        
        %boundary conditions for diffusion
        [gasTransportPores,inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks] = Condensation_GetBoundaryConditionsForDiffusion(network,...
                                        options,invadedPore,equilibriumVaporPressure);
        
        
        [partialVaporPressure,waterConcentration,diffusionFlux,effDiffusion] = Condensation_ComputePartialVaporPressure(network,...
                gasTransportPores,diffusionConductances,...
                inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks,options.AirPressure,temperature);
        
        condensationRatio = partialVaporPressure ./ equilibriumVaporPressure ;
        
        nucleationInfos.EffectiveDiffusion{step} = effDiffusion;
        nucleationInfos.WaterDiffusionFlux{step}=diffusionFlux;
        invadedPores = horzcat(nucleationInfos.InvadedPore{:});
        partialVaporPressure(invadedPores)=equilibriumVaporPressure(invadedPores); %put RH=1 in invaded pores
        nucleationInfos.PartialVaporPressure{step} = partialVaporPressure;
        waterConcentration(invadedPores)=partialVaporPressure(invadedPores)./(8.314*temperature(invadedPores)); %put RH=1 in invaded pores
        nucleationInfos.WaterConcentration{step}=waterConcentration;
        
        
        
        %% invade the pore which has the max partial pressure if > equilibrium
        
        [maxRH,indexMaxRH] = max(condensationRatio);
        
        gasTransportPoresBoolean=zeros(1,nPore);
        gasTransportPoresBoolean(gasTransportPores)=1;
        
        if maxRH>superSaturation
            nucleation = true;
            assert(gasTransportPoresBoolean(indexMaxRH)==1);
            invadedPore(indexMaxRH)=true;
            nucleationInfos.MaxRH{step} = maxRH;
            nucleationInfos.InvadedPore{step}=indexMaxRH;
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
    
    nucleationInfos.MaxRH = cell2mat(nucleationInfos.MaxRH);
    nucleationInfos.InvadedPore = cell2mat(nucleationInfos.InvadedPore);
    nucleationInfos.EffectiveDiffusion = cell2mat(nucleationInfos.EffectiveDiffusion);
end

