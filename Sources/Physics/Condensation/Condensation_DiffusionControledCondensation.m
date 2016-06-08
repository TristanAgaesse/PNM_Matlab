function [condensationClusters, condensationInfos] = Condensation_DiffusionControledCondensation(...
                network,nucleationClusters,options,diffusionConductances)
%CONDENSATION_DIFFUSIONCONTROLEDCONDENSATION Deuxieme partie de l'algorithme 
% de condensation de Marc Prat. 
%   
%   Lors de cette seconde etape, la condensation se poursuit sur les 
%   bords des clusters d'eau nuclees car de la vapeur d'eau arrive par 
%   diffusion. Le  rythme de croissance des clusters est controle par le 
%   flux de vapeur qui arrive par diffusion. On procede donc a un 
%   envahissement selon un algorithme temporel tenant compte d'un bilan de 
%   matiere d'eau. La direction de croissance des clusters est controlee 
%   par les forces capillaires : les liens ayant la pression capillaire 
%   critique la plus faible sont envahis en premier. Plus de détails dans 
%   la fonction DiffusionControledCondensation.
    
    
    %% Initialise algorithm
    %nPoreAccessible=FindNumberOfAccessiblePores(network,1:nPore);
    outputInformation.PartialVaporPressure={};
    outputInformation.EffectiveDiffusion={};
    
    nPore = network.GetNumberOfPores;
    allInvadedPore = zeros(1,nPore);
    nCluster = length(nucleationClusters);
    condensationClusters = cell(1,nCluster);
    for iCluster=1:nCluster
        cluster = nucleationClusters{iCluster}.CopyCluster;
        condensationClusters{iCluster} = cluster;
        allInvadedPore(cluster.GetInvadedPores) = 1;
    end
    
    outletPores = network.GetPoresFrontiere(options.LiquidWaterOutletLinks);
    
    
    
    %% Begin invasion loop
    iteration = 0;
    outlet_reached = false;
    while not(outlet_reached) %%&& iteration<nPoreAccessible
        iteration = iteration+1;
        
        %% Compute water flux on the boundary of each cluster
        
        % Compute water vapor diffusion in the GDL
        [gasTransportPores,inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks] = GetBoundaryConditionsDiffusion(network,options,allInvadedPore);
        
        gasTransportPores = cluster.GetInvadedPoresComplementary; %TODO:check redondance with BC above
        
        [partialVaporPressure,diffusionFlux,effectiveDiffusion ] = Condensation_ComputePartialVaporPressure(network,...
                gasTransportPores,diffusionConductances,...
                inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks,options.AirPressure,temperature);
        
        outputInformation.EffectiveDiffusion{end+1} = effectiveDiffusion;
        outputInformation.PartialVaporPressure{end+1} = partialVaporPressure;
        
        % sum the diffusion flux on the boundary of each cluster
        nCluster = length(condensationClusters);
        totalFlux = zeros(1,nCluster);
        for iCluster=1:nCluster
            cluster = condensationClusters{iCluster};
            totalFlux(iCluster) = diffusionFlux(cluster.GetInterfaceLinks);
        end
            
        %% Find next invasion time and link
        
        %Trouver quel pore est en train d'etre envahi pour chaque cluster
        poreBeingInvaded = zeros(1,nCluster);    
        for iCluster=1:nCluster
            cluster = condensationClusters{iCluster};
            poreBeingInvaded(iCluster) = cluster.GetMinimalPressureLink;
        end    
            
            
        %Trouver le temps d'envahissement de ce pore pour chaque cluster
        for iCluster=1:nCluster
            cluster = condensationClusters{iCluster};
            currentSaturation
            invasionTime
            
            poreBeingInvaded(iCluster) = cluster.GetMinimalPressureLink;
        end     
            
            
        %% Updater chaque cluster : envahissement T=min(Tcluster)
            
            
            
            
        [~,indexMinPressureLink] = Condensation_FindNextInvadedLink(cluster,...
                            partialVaporPressure,equilibriumVaporPressure);
        
        if indexMinPressureLink>0
            invadedPore = cluster.GetOutwardPore(indexMinPressureLink);
            outputInformation.InvadedPore{end+1} = invadedPore;
            
            %TODO : Gérer l'évolution du cluster
            %Condensation_UpdateClustersCinetic
                %TODO : envahir les pores actifs des autres clusters en fonction du temps et des débits
                %TODO : update allInvadedPore
            interfaceChangeInformation=cluster.InvadeNewPore(indexMinPressureLink);
            cluster.UpdateCriticalPressure(interfaceChangeInformation,[],options.LiquidWaterOutletLinks);
            
            
            %verifier si outlet_reached
            if ismember(invadedPore,outletPores)
                outlet_reached = true;
                breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),...
                                              options.LiquidWaterOutletLinks);
                cluster.InvadeOutletLink(breakthroughLinks);
            end
        else
            %no new pore condensable near cluster
            return
        end
    end
    
end

