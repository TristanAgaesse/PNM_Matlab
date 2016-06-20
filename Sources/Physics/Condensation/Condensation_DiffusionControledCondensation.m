function [condensationClusters, condensationInfos] = Condensation_DiffusionControledCondensation(...
                network,nucleationClusters,options,diffusionConductances,equilibriumVaporPressure,temperature)
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
    condensationInfos.PartialVaporPressure={};
    condensationInfos.EffectiveDiffusion={};
    condensationInfos.InvadedPore={};
    condensationInfos.InvasionTime={};
    
    nPore = network.GetNumberOfPores;
    allInvadedPore = zeros(1,nPore);
    nCluster = length(nucleationClusters);
    condensationClusters = cell(1,nCluster);
    for iCluster=1:nCluster
        cluster = nucleationClusters{iCluster}.CopyCluster;
        condensationClusters{iCluster} = cluster;
        allInvadedPore(cluster.GetInvadedPores) = 1;
    end
    poreSaturation=allInvadedPore;
    
    outletLink = options.LiquidWaterOutletLinks;
    %outletPores = network.GetPoresFrontiere(options.LiquidWaterOutletLinks);
    
    poreVolume = network.GetPoreData('Volume');
    
    if nCluster==0
       return 
    end
    
    
    %% Begin invasion loop
    iteration = 0;
    outlet_reached = false;
    invasionTime = 0;
    while invasionTime<100 && iteration<100
        iteration = iteration+1;
        
        %% Compute water flux on the boundary of each cluster
        
        % Compute water vapor diffusion in the GDL
        [gasTransportPores,inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks] = Condensation_GetBoundaryConditionsForDiffusion(network,options,allInvadedPore,equilibriumVaporPressure);
        
        %gasTransportPores = cluster.GetInvadedPoresComplementary; %TODO:check redondance with BC above
        
        [partialVaporPressure,diffusionFlux,effectiveDiffusion ] = Condensation_ComputePartialVaporPressure(network,...
                gasTransportPores,diffusionConductances,...
                inletVaporPressure,outletVaporPressure,vaporInletLinks,...
                vaporOutletLinks,options.AirPressure,temperature);
        
        condensationInfos.EffectiveDiffusion{end+1} = effectiveDiffusion;
        condensationInfos.PartialVaporPressure{end+1} = partialVaporPressure;
        
        % sum the diffusion flux on the boundary of each cluster
        nCluster = length(condensationClusters);
        totalFlux = zeros(1,nCluster);
        for iCluster=1:nCluster
            cluster = condensationClusters{iCluster};
            totalFlux(iCluster) = sum(diffusionFlux(cluster.GetInterfaceLinks));
        end
        assert(all(totalFlux>0),'Evaporation at one cluster !') 
        
        
        %% Find next invasion time and link
        
        %Trouver quel pore est en train d'etre envahi pour chaque cluster
        poreBeingInvaded = zeros(1,nCluster);    
        for iCluster=1:nCluster
            cluster = condensationClusters{iCluster};
            indexInvadedLink = cluster.GetMinimalPressureLink;
            invadedPore = cluster.GetOutwardPore(indexInvadedLink);
            while invadedPore<=0 %Check for outlet link invasion
                disp('Outlet link invaded')
                linkAbsoluteIndex = cluster.GetInterfaceLinkAbsoluteNumber(indexInvadedLink);
                assert(ismember(linkAbsoluteIndex,outletLink))
                cluster.InvadeOutletLink(linkAbsoluteIndex);
                indexInvadedLink = cluster.GetMinimalPressureLink;
                invadedPore = cluster.GetOutwardPore(indexInvadedLink);
            end
            poreBeingInvaded(iCluster) = invadedPore;
        end    
        
        
        %Trouver le temps d'envahissement de ce pore pour chaque cluster
        clusterInvasionTime = zeros(1,nCluster); 
        for iCluster=1:nCluster
            currentSaturation = poreSaturation(poreBeingInvaded(iCluster));
            volume=poreVolume(poreBeingInvaded(iCluster));
            clusterInvasionTime(iCluster) = volume*(1-currentSaturation)/totalFlux(iCluster);
        end  
        
        assert(all(clusterInvasionTime>0),'Negative invasion time at one cluster !')    
        
        
        %% Updater chaque cluster : envahissement T=min(Tcluster)
        
        [invasionTime,numCluster]=min(clusterInvasionTime);    
        condensationInfos.InvasionTime{iteration}=invasionTime;
            
        %Mettre a jour le cluster envahi
        cluster = condensationClusters{iCluster};
        relativeIndexPore = cluster.GetMinimalPressureLink;

%             %verifier si outlet_reached
%             if ismember(poreBeingInvaded(numCluster),outletPores)
%                 outlet_reached = true;
%                 breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),...
%                                               options.LiquidWaterOutletLinks);
%                 cluster.InvadeOutletLink(breakthroughLinks);
%             end

        interfaceChangeInformation = cluster.InvadeNewPore(relativeIndexPore);
        cluster.UpdateCriticalPressure(interfaceChangeInformation,[],options.LiquidWaterOutletLinks);
        poreSaturation(poreBeingInvaded(numCluster)) = 1;
        condensationInfos.InvadedPore{end+1} = poreBeingInvaded(numCluster);

        %Mettre a jour les autres cluster
        for iCluster = setdiff(1:nCluster,numCluster) 
            currentSaturation = poreSaturation(poreBeingInvaded(iCluster));
            volume = poreVolume(poreBeingInvaded(iCluster));
            saturationIncrease = totalFlux(iCluster)*invasionTime/volume;                
            poreSaturation(poreBeingInvaded(iCluster)) = currentSaturation+saturationIncrease;
            assert(poreSaturation(poreBeingInvaded(iCluster))<1)
        end
            
    end
    
    condensationInfos.EffectiveDiffusion=cell2mat(condensationInfos.EffectiveDiffusion);
    condensationInfos.InvadedPore=cell2mat(condensationInfos.InvadedPore);
    condensationInfos.InvasionTime=cell2mat(condensationInfos.InvasionTime);    
    
end

