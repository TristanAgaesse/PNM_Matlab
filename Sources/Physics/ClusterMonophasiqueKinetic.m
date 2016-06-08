classdef ClusterMonophasiqueKinetic < ClusterMonophasique
    % ClusterMonophasiqueKinetic Subclass of ClusterMonophasique, used for
    % kinetic invasions
    
    properties
        PoreSaturation
    end
    
    methods
        function cluster = ClusterMonophasiqueKinetic(invadedPores,interfaceLinks,interfacePoresOutward,...
                                booleanInvadedLinks,network,criticalPressures,clusterOptions,poreSaturation)
            %Constructeur
            
            cluster = cluster@ClusterMonophasique(invadedPores,interfaceLinks,interfacePoresOutward,...
                                booleanInvadedLinks,network,criticalPressures,clusterOptions);
            
            nPore = network.GetNumberOfPores;
            assert(length(poreSaturation)==nPore)
            assert(all(poreSaturation<=1) && all(poreSaturation>=0))
            
            cluster.PoreSaturation = poreSaturation;               
        end
        
        function interfaceChangeInformation = InvadeNewPore(cluster,indexInvadedLink)
            %Gere l'envahissement d'un nouveau pore 
            %input : -cluster
            %        -indexInvadedLink (indice du lien envahi dans la liste
            %        cluster.InterfaceLinks)
            %output : interfaceChangeInformation : 
            %       interfaceChangeInformation{iPoreOutaward}  =  { {numPore} ,{[indices des liens interface donnant sur ce pore]}
            
            assert(0<indexInvadedLink && indexInvadedLink<=length(cluster.InterfaceLinks),...
                        'indexInvadedLink must be the reference of a link in cluster.InterfaceLinks')
            
            %envahir le pore associe
            time  =  length(find(cluster.InvadedPores))+1;
            
            if time<length(cluster.InvadedPores)+1
                newInvadedPore  =  cluster.InterfacePoresOutward(indexInvadedLink);
                if newInvadedPore<=0
                    disp('Erreur : Impossible d''envahir l''exterieur du domaine !')
                    return
                end
                cluster.InvadedPores(time)  =  newInvadedPore;
                
                %mettre a jour les listes reperant la position de la frontiere
                %Beware, the invaded pore can be a neighbour of the cluster through several links, not just indexInvadedLink
                invadedPoreLinks  =  cluster.Network.GetLinksOfPore(newInvadedPore);
                
                [linkAllreadyInInterface, loc ]= ismember(invadedPoreLinks,cluster.InterfaceLinks); %loop over interface links
                linksToBeInvaded = invadedPoreLinks(linkAllreadyInInterface); 
                
                cluster.BooleanInvadedLinks(linksToBeInvaded) = 1*ones(1,length(linksToBeInvaded));
                
                newInterfaceLinks  =  invadedPoreLinks(~linkAllreadyInInterface);  
                
                linksToDelete = zeros(1,length(cluster.InterfaceLinks));            
                linksToDelete(loc(linkAllreadyInInterface)) = ones(1,length(loc(linkAllreadyInInterface)));  
                
                cluster.InterfaceLinks  =  [cluster.InterfaceLinks(~linksToDelete),newInterfaceLinks];
                
                pore_outward_of_new_links  =  zeros(1,length(newInterfaceLinks));
                allNeighbourPores=cluster.Network.GetPoresOfLink(newInterfaceLinks);
                for iLink  =  1:length(newInterfaceLinks)
                    neighbourPores = allNeighbourPores(iLink,:);
                    pore_outward_of_new_links(iLink)  =  neighbourPores(neighbourPores~=newInvadedPore);  
                end
                cluster.InterfacePoresOutward  =  [cluster.InterfacePoresOutward(~linksToDelete),pore_outward_of_new_links];

                
                %liste des liens dont il faut updater la pression
                unique_new_pore_outward  =  unique(pore_outward_of_new_links);
                interfaceChangeInformation  =  cell(1,length(unique_new_pore_outward));
                for iPore = 1:length(unique_new_pore_outward)
                    interfaceChangeInformation{iPore}{1} = unique_new_pore_outward(iPore);
                    interfaceChangeInformation{iPore}{2} = cluster.InterfaceLinks(...
                                cluster.InterfacePoresOutward==unique_new_pore_outward(iPore)); %loop over interface pores outward
                end
            end
        end
        
        
        function interfaceChangeInformation = GetInterfaceChangeInformation(cluster,indexLinksToUpdate)
            %Donne la structure interfaceChangeInformation utilisee par
            %UpdateCriticalPressure pour mettre a jour les pressions
            %capillaires
            %input: cluster, indexLinksToUpdate (indices dans la liste des liens d'interface du cluster)
            %output : -interfaceChangeInformation : 
            %          interfaceChangeInformation{iPoreOutward} = { {numPore} ,{[indices des liens interface donnant sur ce pore]}
            
            poreOutward = cluster.GetOutwardPore(indexLinksToUpdate);
            
            [unique_pore_outward,~,iU] = unique(poreOutward);
            interfaceChangeInformation = cell(1,length(unique_pore_outward));
            for i = 1:length(unique_pore_outward)
                interfaceChangeInformation{i} = cell(1,2);
                interfaceChangeInformation{i}{1} = unique_pore_outward(i);
            end
            
            for i = 1:length(indexLinksToUpdate)
                assert(interfaceChangeInformation{iU(i)}{1} == poreOutward(i));
                interfaceChangeInformation{iU(i)}{2} = [interfaceChangeInformation{iU(i)}{2},indexLinksToUpdate(i)];
            end

        end
        
        
        function UpdateCriticalPressure(cluster,interfaceChangeInformation,linkInlet,linkOutlet)         
            %Met a jour les pressions critiques a la frontiere d'un amas
            %liquide lorsqu'un nouveau pore est envahi.
            %input : -cluster
            %        -interfaceChangeInformation = cell array, {{num_pore_outward},{[num_liens]}}
            %        -linkInlet
            %        -linkOutlet
            
            options = cluster.GetClusterOptions;
                        
            theta = cluster.Network.LinkDataList.LinkDatas.ContactAngle;
            
            nInterfacePore = length(interfaceChangeInformation);
            for iInterfacePore = 1:nInterfacePore
                pore = interfaceChangeInformation{iInterfacePore}{1};
                liens = interfaceChangeInformation{iInterfacePore}{2};
                nLien = length(liens);
                
                if pore == -1 %liens sur une frontiere envahis depuis l'interieur du domaine
                    
                    isLinkInlet=ismember(liens,linkInlet);%inlet
                    cluster.CriticalPressures(liens(isLinkInlet))=Inf;
                    
                    isLinkOutlet=ismember(liens,linkOutlet);%outlet
                    
                    cluster.CriticalPressures(liens(isLinkOutlet)) = ...
                        LocalScaleComputeCriticalPressureWithoutCoalescence(cluster.Network,liens(isLinkOutlet),options);
                    
                    cluster.CriticalPressures(liens(not(or(isLinkInlet,isLinkOutlet)))) = Inf ;%wall
                    
                    
                else %liens internes ou inlet envahis depuis l'exterieur
                    
                    Pc = LocalScaleComputeCriticalPressureWithoutCoalescence(cluster.Network,liens,options);
                    
                    if strcmp(options.Coalescence,'none')
                        cluster.CriticalPressures(liens) = Pc;
                    
                    elseif strcmp(options.Coalescence,'numberOfInvadedNeighbours')
                        coordinance = length(cluster.Network.GetLinksOfPore(pore));
                        for iLien = liens
                            if theta(iLien)<pi/2 %cas hydrophile
                                coalescenceFactor = (1-(nLien-1)/coordinance);
                            else %cas hydrophobe
                                coalescenceFactor = 1;
                            end
                        end
                        
                        cluster.CriticalPressures(liens) = coalescenceFactor*Pc;
                    end
                end
            end
            
        end

        
        function poreSaturation = GetPoreSaturation(cluster)
            poreSaturation = cluster.PoreSaturation;
        end 
        
        function pores = GetInvadedPores(cluster)
            pores = cluster.InvadedPores(cluster.InvadedPores>0);
        end        
                

        function booleanArray = GetInvadedPoresBooleans(cluster)
            %input: cluster
            %output : booleanArray : 1*nPore logical array : is this pore invaded ?
            booleanArray = zeros(1,cluster.Network.GetNumberOfPores);
            for iPore = cluster.GetInvadedPores
                booleanArray(iPore) = 1;
            end
        end
        
        
        function compCluster = GetComplementaryCluster(cluster)
            %Renvoie le cluster complementaire correspondant aux pores non
            %envahis d'un cluster donne.
            
            compInvadedPores = setdiff((1:cluster.Network.GetNumberOfPores),cluster.GetInvadedPores);
            compBooleanInvadedLinks = not(cluster.BooleanInvadedLinks);
            compInterfaceLinks = cluster.GetInterfaceLinks;
            clusterOptions = cluster.GetClusterOptions;
            
            poreSaturation=zeros(1,cluster.Network.GetNumberOfPores);
            poreSaturation(compInvadedPores)=1;     % TO DO
            
            compInterfacePoresOutward = [];             % TO DO
            criticalPressures = zeros(1,length(compInterfaceLinks)); % TO DO
            
            compCluster = ClusterMonophasiqueKinetic(compInvadedPores,compInterfaceLinks,...
                        compInterfacePoresOutward,compBooleanInvadedLinks,...
                        cluster.Network,criticalPressures,clusterOptions,poreSaturation);
        end
        
        
        function fusionCluster = FuseClusters(cluster1,cluster2)
            %Renvoit le cluster forme de l'union de deux clusters
            %Input: cluster1,cluster2
            %Output : fusionCluster
            
            clusterOptions = cluster1.GetClusterOptions;
            
            fusionInvadedPores = union(cluster1.GetInvadedPores,cluster2.GetInvadedPores);
            fusionBooleanInvadedLinks = or(cluster1.BooleanInvadedLinks,cluster2.BooleanInvadedLinks);
            
            interfaceIntersection = intersect(cluster1.GetInterfaceLinks,cluster2.GetInterfaceLinks);
            [xorInterface1,xorIndices1]=setdiff(cluster1.GetInterfaceLinks,interfaceIntersection);
            [xorInterface2,xorIndices2]=setdiff(cluster2.GetInterfaceLinks,interfaceIntersection);
            
            fusionInterfaceLinks = [xorInterface1,xorInterface2];
            
            poreSaturation = zeros(1,cluster1.Network.GetNumberOfPores);
            poreSaturation(fusionInvadedPores)=1 ; %TODO : put saturation <1 for pores not completely invaded
            
            compInterfacePoresOutward = [cluster1.InterfacePoresOutward(xorIndices1),cluster2.InterfacePoresOutward(xorIndices2)]; 
            criticalPressures = zeros(1,length(cluster1.Network.GetNumberOfLinks));
            criticalPressures(xorInterface1)=cluster1.CriticalPressures(xorInterface1);
            criticalPressures(xorInterface2)=cluster2.CriticalPressures(xorInterface2);
            
            fusionCluster = ClusterMonophasiqueKinetic(fusionInvadedPores,fusionInterfaceLinks,...
                compInterfacePoresOutward,fusionBooleanInvadedLinks,...
                cluster1.Network,criticalPressures,clusterOptions,poreSaturation);
        end
        
        function copiedCluster = CopyCluster(cluster)
            %CopyCluster : fait une copie du cluster (utile car cluster est
            %un handle donc n'accepte pas facilement d'etre copie).
            copiedCluster = ClusterMonophasiqueKinetic(cluster.InvadedPores,cluster.InterfaceLinks,...
                cluster.InterfacePoresOutward ,cluster.BooleanInvadedLinks,...
                cluster.Network,cluster.CriticalPressures,cluster.GetClusterOptions,cluster.GetPoreSaturation);
        end
        
    end
    
    methods (Static = true)
        function cluster = CreateVoidCluster(network)
            %input : network
            %output: void cluster
            invadedPores = [];
            poreSaturation = zeros(1,network.GetNumberOfPores);
            interfaceLinks = [];
            poresFrontiereOutward = [];
            booleanInvadedLinks = zeros(1,network.GetNumberOfLinks);
            criticalPressures = zeros(1,network.GetNumberOfLinks);
            clusterOptions=struct;
            
            cluster = ClusterMonophasiqueKinetic(invadedPores,interfaceLinks,...
                poresFrontiereOutward,booleanInvadedLinks,...
                network,criticalPressures,clusterOptions,poreSaturation);
            
        end
        
        function cluster = InitialiseCondensationCluster(network,clusterOptions,firstInvadedPores,outletLink)
            %input: network,clusterOptions,firstInvadedPore,inletLink,outletLink
            %output : cluster paramétré de façon à commencer
            %       une condensation depuis firstInvadedPores
            
            %listes qui reperent les pores et les liens envahis
            nPore = network.GetNumberOfPores;
            poreSaturation = zeros(1,nPore);
            poreSaturation(firstInvadedPores)=1;
                        
            nLink = network.GetNumberOfLinks;
            booleanInvadedLinks = zeros(1,nLink);
            
            %listes qui reperent la position de la frontiere
            interfaceLinks = network.GetPoreRegionBoundaryLinks(firstInvadedPores);
            allInterfacePore = network.GetPoresOfLink(interfaceLinks);
            
            interfacePoreOutward = zeros(1,length(interfaceLinks)); %si exterieur : 0 pour outlet, -1 pour wall
            foo1=not(ismember(allInterfacePore(:,1),firstInvadedPores));
            foo2=not(ismember(allInterfacePore(:,2),firstInvadedPores));
            assert(all((foo1+foo2)==1),'un link frontiere a 0 ou 2 pores voisins dans le cluster !')
            interfacePoreOutward(foo2) = allInterfacePore(foo2,2);
            interfacePoreOutward(foo1) = allInterfacePore(foo1,1);
            interfacePoreOutward(allInterfacePore(foo1,1)==-1)=-1;
            interfacePoreOutward(ismember(interfaceLinks,outletLink)) = 0;
            
            criticalPressures = zeros(1,network.GetNumberOfLinks);
            
            cluster = ClusterMonophasiqueKinetic(firstInvadedPores,interfaceLinks,...
                interfacePoreOutward,booleanInvadedLinks,...
                network,criticalPressures,clusterOptions,poreSaturation);
            
            %initialisation des pressions critiques
            [poreOutwardUnique,~,m] = unique(interfacePoreOutward);
            nporeOutwardUnique = length(poreOutwardUnique);
            
            linksToInitialise = cell(1,nporeOutwardUnique);
            
            for iporeOutwardUnique = 1:nporeOutwardUnique
                linksToInitialise{iporeOutwardUnique}{1} = poreOutwardUnique(iporeOutwardUnique);%numero du pore
                linksToInitialise{iporeOutwardUnique}{2} = [];
            end
            for iInterfaceLink = 1:length(interfaceLinks)
                iporeOutwardUnique = m(iInterfaceLink);
                linksToInitialise{iporeOutwardUnique}{2} = ...
                    [linksToInitialise{iporeOutwardUnique}{2},interfaceLinks(iInterfaceLink)]; %numeros des liens associe a ce pore
            end
            
            
            
            cluster.UpdateCriticalPressure(linksToInitialise,[],outletLink);
            
        end
        
        
    end
    
end

