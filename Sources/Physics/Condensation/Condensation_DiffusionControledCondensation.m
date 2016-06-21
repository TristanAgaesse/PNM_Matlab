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
        evaporationBool=not(all(totalFlux>0));
        if evaporationBool
            disp('Evaporation at one cluster !')
            break
        end
        
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
            iPore=poreBeingInvaded(iCluster);
            currentSaturation = poreSaturation(iPore);
            volume=poreVolume(iPore);
            clusterInvasionTime(iCluster) = (1-currentSaturation)*(volume/totalFlux(iCluster));
        end  
        
        assert(all(clusterInvasionTime>0),'Negative invasion time at one cluster !')    
        
        [invasionTime,iCriticalCluster]=min(clusterInvasionTime);  
        
        % Gerer le cas d'un envahissement cooperatif: plusieurs clusters envahissent le meme pore
        [uniquePoreInvaded,ia,ic]= unique(poreBeingInvaded); %[C,ia,ic]=unique(A); C=A(ia) and A=C(ic).
        nCooperativeCluster = length(uniquePoreInvaded);
        cooperativeBool = not(nCooperativeCluster==nCluster);
        cooperativeInvasionTime = zeros(1,nCooperativeCluster);
        if cooperativeBool
            %Compute cooperative invasion time for all clusters
            for iCooperativeCluster=1:nCooperativeCluster
                iPore=poreBeingInvaded(iCooperativeCluster);
                currentSaturation = poreSaturation(iPore);
                volume=poreVolume(iPore);
                cooperativeflux = sum(totalFlux(ia(iCooperativeCluster))); %Check formula 
                cooperativeInvasionTime(iCooperativeCluster) = (1-currentSaturation)*(volume/cooperativeflux);
            end
            
            % update effective invasion time if cooperation is responsible for the next pore invasion 
            [cooperativeInvasionTime,iCriticalCooperativeCluster]=min(cooperativeInvasionTime);
            if cooperativeInvasionTime<invasionTime
                cooperativeBool = true;
                invasionTime = cooperativeInvasionTime;
            else
                cooperativeBool = false;
            end
            
        end
        
        condensationInfos.InvasionTime{iteration}=invasionTime;
        
        
        %% Updater chaque cluster 
        
        % Mettre a jour les clusters ou il y a invasion complete d un pore
        if not(cooperativeBool)
            
            iPore = poreBeingInvaded(iCriticalCluster);
            poreSaturation(iPore) = 1;
            allInvadedPore(iPore) = 1;
            condensationInfos.InvadedPore{end+1} = iPore;
            
            cluster = condensationClusters{iCluster};
            relativeIndexPore = cluster.GetMinimalPressureLink;
            assert(iPore==cluster.GetOutwardPore(relativeIndexPore))
            interfaceChangeInformation = cluster.InvadeNewPore(relativeIndexPore);
            cluster.UpdateCriticalPressure(interfaceChangeInformation,[],options.LiquidWaterOutletLinks);
            
            simpleFillingClustersList = setdiff(1:nCluster,iCriticalCluster);
            
        elseif cooperativeBool
            
            criticalClusters = ia(iCriticalCooperativeCluster);
            iPore = unique(poreBeingInvaded(criticalClusters));
            assert(length(iPore)==1)
            
            poreSaturation(iPore) = 1;
            allInvadedPore(iPore) = 1; 
            
            % fusion de deux clusters
            iCriticalCluster=criticalClusters(1);
            cluster = condensationClusters{iCriticalCluster};
            relativeIndexPore = cluster.GetMinimalPressureLink;
            assert(iPore==cluster.GetOutwardPore(relativeIndexPore))
            interfaceChangeInformation = cluster.InvadeNewPore(relativeIndexPore);
            cluster.UpdateCriticalPressure(interfaceChangeInformation,[],options.LiquidWaterOutletLinks);
            
            newCluster=cluster.CopyCluster;
            condensationClusters{criticalClusters(1)}.delete;
            for i=criticalClusters(2:end)
                newCluster = newCluster.FuseClusters(condensationClusters{i}).CopyCluster;
                condensationClusters{i}.delete;
            end
            
            simpleFillingClustersList = setdiff(1:nCluster,criticalClusters);
            nCluster = length(simpleFillingClustersList)+1;
            newCondensationClusters=cell(1,nCluster);
            for i=1:nCluster-1
                newCondensationClusters{i}=condensationClusters{simpleFillingClustersList(i)}.CopyCluster;
            end
            newCondensationClusters{nCluster}=newCluster;
            condensationClusters=newCondensationClusters;
        end
        
        % Mettre a jour les saturations des pores partiellement envahis
        for iCluster = simpleFillingClustersList
            iPore = poreBeingInvaded(iCluster);
            currentSaturation = poreSaturation(iPore);
            volume = poreVolume(iPore);
            saturationIncrease = invasionTime*(totalFlux(iCluster)/volume);                
            poreSaturation(iPore) = currentSaturation+saturationIncrease;
            assert(poreSaturation(iPore)<=1)
        end
        
        
        %%verifier si outlet_reached
%             if ismember(poreBeingInvaded(numCluster),outletPores)
%                 outlet_reached = true;
%                 breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),...
%                                               options.LiquidWaterOutletLinks);
%                 cluster.InvadeOutletLink(breakthroughLinks);
%             end
        
    end
    
    condensationInfos.EffectiveDiffusion=cell2mat(condensationInfos.EffectiveDiffusion);
    condensationInfos.InvadedPore=cell2mat(condensationInfos.InvadedPore);
    condensationInfos.InvasionTime=cell2mat(condensationInfos.InvasionTime);    
    
end

