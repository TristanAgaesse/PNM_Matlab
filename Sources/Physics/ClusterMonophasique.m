classdef ClusterMonophasique < handle
    %CLUSTERMONOPHASIQUE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        InvadedPores
        BooleanInvadedLinks
        InterfaceLinks
        InterfacePoresOutward %si exterieur : 0 pour outlet, -1 pour wall
        CriticalPressures    %tableau 1*network.GetNumberOfLink contenant les pressions critiques d'invasion
        Network
        ClusterOptions  %(optionnel) clusterOptions.Coalescence = 'none' or 'numberOfInvadedNeighbours'
                           %     clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder','PurcellToroid'
                           %    clusterOptions.SurfaceTension
    end
    
    methods
        function cluster = ClusterMonophasique(invadedPores,interfaceLinks,interfacePoresOutward,...
                                booleanInvadedLinks,network,criticalPressures,clusterOptions)
            %Constructeur
            cluster.InvadedPores  =  invadedPores;
            cluster.InterfaceLinks  =  interfaceLinks;
            cluster.InterfacePoresOutward  =  interfacePoresOutward;
            cluster.BooleanInvadedLinks  =  booleanInvadedLinks;
            cluster.Network  =  network;
            cluster.CriticalPressures  =  criticalPressures;
            
            if isfield(clusterOptions,'CapillaryPressureLaw')
                if strcmp(clusterOptions.CapillaryPressureLaw,'PurcellToroid')
                    clusterOptions.CapillaryPressureLaw  = 'PurcellToroid' ;
                end
            else
                clusterOptions.CapillaryPressureLaw  = 'LaplaceCylinder' ;
            end
            
            if isfield(clusterOptions,'Coalescence')
                if strcmp(clusterOptions.Coalescence,'none')
                    clusterOptions.Coalescence  =  'none';
                end
            else
                clusterOptions.Coalescence  =  'numberOfInvadedNeighbours';
            end
            
            cluster.ClusterOptions  =  clusterOptions;
        end
        
        
        function [indexMinPressureLink,minPressure]  =  GetMinimalPressureLink(cluster)
            %Donne le lien de plus petite pression capillaire parmi les
            %liens a l'interface du cluster
            %input : cluster
            %output : - indexMinPressureLink : indice du lien dans la liste
            %                   interne au cluster des liens d'interfaces 
            %         - minPressure : valeur de la pression capillaire min
            criticalPressures  =  cluster.GetCriticalPressures;
            [minPressure,indexMinPressureLink]  =  min(criticalPressures(cluster.GetInterfaceLinks));
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
        
        
        function InvadeOutletLink(cluster,invadedOutletLink)
            %Gere l'envahissement d'un lien situe à l'outlet
            %input: cluster, invadedOutletLink
                        
            cluster.BooleanInvadedLinks(invadedOutletLink) = 1*ones(1,length(invadedOutletLink));
            
            [cluster.InterfaceLinks,index] = setdiff(cluster.InterfaceLinks,invadedOutletLink);
            
            cluster.InterfacePoresOutward = cluster.InterfacePoresOutward(index);
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

        
        function links = GetInvadedLinks(cluster)
            links = find(cluster.BooleanInvadedLinks);
        end
        
        function linksBoolean = GetInvadedLinksBoolean(cluster)
            linksBoolean = cluster.BooleanInvadedLinks;
        end
        
        function pores = GetInvadedPores(cluster)
            pores = cluster.InvadedPores(cluster.InvadedPores>0);
        end        
                
        function pores=GetInvadedPoresComplementary(cluster)
            allPores=1:cluster.Network.GetNumberOfPores;
            pores=setdiff(allPores,cluster.GetInvadedPores);
        end
        
        function pore = GetOutwardPore(cluster,indexLink)
            pore = cluster.InterfacePoresOutward(indexLink);
        end
        
        function links = GetInterfaceLinks(cluster)
            links = cluster.InterfaceLinks;
        end
        
        function interfaceLength = GetInterfaceLength(cluster)
            interfaceLength = length(cluster.InterfaceLinks);
        end
        
        function options = GetClusterOptions(cluster)
            options = cluster.ClusterOptions;
        end
        
        function SetClusterOptions(cluster,clusterOptions)
            cluster.ClusterOptions=clusterOptions;
        end
        
        
        function criticalPressures = GetCriticalPressures(cluster)
            criticalPressures = cluster.CriticalPressures;
        end
        
        
        function setCriticalPressures(cluster,criticalPressures)
            assert(length(criticalPressures) == cluster.Network.GetNumberOfLinks,...
                'critical pressures doit etre un tableau 1*network.GetNumberOfLink')
            cluster.CriticalPressures = criticalPressures;
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
            
            compInterfacePoresOutward = [];                                                            % TO DO
            criticalPressures = zeros(1,length(compInterfaceLinks));                              % TO DO
            
            compCluster = ClusterMonophasique(compInvadedPores,compInterfaceLinks,...
                        compInterfacePoresOutward,compBooleanInvadedLinks,...
                        cluster.Network,criticalPressures,clusterOptions);
        end
        
        
        function fusionCluster = FuseClusters(cluster1,cluster2)
            %Renvoit le cluster forme de l'union de deux clusters
            %Input: cluster1,cluster2
            %Output : fusionCluster
            
            clusterOptions = cluster1.GetClusterOptions;
            
            fusionInvadedPores = union(cluster1.GetInvadedPores,cluster2.GetInvadedPores);
            fusionBooleanInvadedLinks = or(cluster1.BooleanInvadedLinks,cluster2.BooleanInvadedLinks);
            fusionInterfaceLinks = setxor(cluster1.GetInterfaceLinks,cluster2.GetInterfaceLinks);
            
            compInterfacePoresOutward = [];                                                            % TO DO
            criticalPressures = zeros(1,length(fusionInterfaceLinks));                              % TO DO
            
            fusionCluster = ClusterMonophasique(fusionInvadedPores,fusionInterfaceLinks,...
                compInterfacePoresOutward,fusionBooleanInvadedLinks,...
                cluster1.Network,criticalPressures,clusterOptions);
        end
        
        
        
        
        function poresPercolants = FindPercolationPath(cluster,linkInlet,linkOutlet)
            %Trouver les chemins de percolation : composantes connexes de 
            %pores envahis d�bouchant sur l'inlet et l'outlet
            %input : cluster,linkInlet,linkOutlet,poreNetwork
            %output: poresPercolants
            invadedPoresList = cluster.GetInvadedPores;
            pores_outlet = cluster.Network.GetPoresFrontiere(linkOutlet);
            pores_inlet = cluster.Network.GetPoresFrontiere(linkInlet);
            
            pore_envahis_outlet = intersect(invadedPoresList,pores_outlet);
            pore_envahis_inlet = intersect(invadedPoresList,pores_inlet);
            
            network = cluster.Network;
            
            num_chemin = 1;
            poresPercolants = {};
            composantes_connexes_du_fluide = network.FindComposantesConnexes(invadedPoresList);
            for num_composante = 1:length(composantes_connexes_du_fluide)
                if ~isempty(intersect(composantes_connexes_du_fluide{num_composante},pore_envahis_outlet)) ...
                        && ~isempty(find(ismember(composantes_connexes_du_fluide{num_composante},pore_envahis_inlet),1))
                    
                    %find invaded links
                    invadedPores = composantes_connexes_du_fluide{num_composante};
                    linkCount = zeros(1,network.GetNumberOfLinks);
                    for iPore = invadedPores
                        thoseLinks = network.Pores{iPore};
                        linkCount(thoseLinks) = linkCount(thoseLinks)+1;
                    end
                    
                    booleanInvadedLinks = linkCount>0;
                    interfaceLinks = find(linkCount == 1);
                    poresFrontiereOutward = [];
                    criticalPressures = zeros(1,network.GetNumberOfLinks);
                    clusterOptions=cluster.GetClusterOptions;
                    
                    percolatingCluster = ClusterMonophasique(invadedPores,interfaceLinks,...
                        poresFrontiereOutward,booleanInvadedLinks,...
                        network,criticalPressures,clusterOptions);
                    poresPercolants{num_chemin} = percolatingCluster;
                    num_chemin = num_chemin+1;
                end
            end

        end
        
        
        function copiedCluster = CopyCluster(cluster)
            %CopyCluster : fait une copie du cluster (utile car cluster est
            %un handle donc n'accepte pas facilement d'etre copie).
            copiedCluster = ClusterMonophasique(cluster.InvadedPores,cluster.InterfaceLinks,...
                cluster.InterfacePoresOutward ,cluster.BooleanInvadedLinks,...
                cluster.Network,cluster.CriticalPressures,cluster.GetClusterOptions);
        end
        
    end
    
    methods (Static = true)
        function cluster = CreateVoidCluster(network)
            %input : network
            %output: void cluster
            invadedPores = zeros(1,network.GetNumberOfPores);
            interfaceLinks = [];
            poresFrontiereOutward = [];
            booleanInvadedLinks = zeros(1,network.GetNumberOfLinks);
            criticalPressures = zeros(1,network.GetNumberOfLinks);
            clusterOptions=struct;
            
            cluster = ClusterMonophasique(invadedPores,interfaceLinks,...
                poresFrontiereOutward,booleanInvadedLinks,...
                network,criticalPressures,clusterOptions);
            
        end
        
        function cluster = InitialiseInvasionCluster(linkInlet,linkOutlet,poreNetwork,varargin)
            %input: linkInlet,linkOutlet,poreNetwork, varargin(optionnel : clusterOptions)
            %output : cluster parametre de façon à commencer
            %       une invasion percolation depuis les linkInlet
            
            %listes qui rep�rent les pores et les liens envahis
            nPore = poreNetwork.GetNumberOfPores;
            invadedPores = zeros(1,nPore);
            
            nLink = poreNetwork.GetNumberOfLinks;
            booleanInvadedLinks = zeros(1,nLink);
            
            %listes qui reperent la position de la frontiere
            interfaceLinks = linkInlet;
            inletPoreOutward = poreNetwork.GetPoresFrontiere(linkInlet);   
            
            criticalPressures = zeros(1,poreNetwork.GetNumberOfLinks);
            
            clusterOptions = struct;
            if not(isempty(varargin))
                clusterOptions = varargin{1};
            end
            
            cluster = ClusterMonophasique(invadedPores,interfaceLinks,...
                inletPoreOutward,booleanInvadedLinks,...
                poreNetwork,criticalPressures,clusterOptions);
            
            %initialisation des pressions critiques
            [poreOutwardUnique,~,m] = unique(inletPoreOutward);
            nporeOutwardUnique = length(poreOutwardUnique);
            
            linksToInitialise = cell(1,nporeOutwardUnique);
            
            for iporeOutwardUnique = 1:nporeOutwardUnique
                linksToInitialise{iporeOutwardUnique}{1} = poreOutwardUnique(iporeOutwardUnique);%numero du pore
                linksToInitialise{iporeOutwardUnique}{2} = [];
            end
            for iInletLink = 1:length(linkInlet)
                iporeOutwardUnique = m(iInletLink);
                linksToInitialise{iporeOutwardUnique}{2} = ...
                    [linksToInitialise{iporeOutwardUnique}{2},linkInlet(iInletLink)]; %numeros des liens associe a ce pore
            end
            
            
            cluster.UpdateCriticalPressure(linksToInitialise,linkInlet,linkOutlet);
            
        end
        

        
        function cluster = InitialiseCondensationCluster(network,clusterOptions,firstInvadedPore,outletLink)
            %input: network,clusterOptions,firstInvadedPore,inletLink,outletLink
            %output : cluster paramétré de façon à commencer
            %       une condensation depuis firstInvadedPore
            
            %listes qui reperent les pores et les liens envahis
            nPore = network.GetNumberOfPores;
            invadedPores = zeros(1,nPore);
            invadedPores(firstInvadedPore)=1;
            
            nLink = network.GetNumberOfLinks;
            booleanInvadedLinks = zeros(1,nLink);
            
            %listes qui rep�rent la position de la frontiere
            interfaceLinks = network.GetLinksOfPore(firstInvadedPore);
            interfacePoreOutward = network.GetPoresVoisinsOfPore(firstInvadedPore);   
            
            criticalPressures = zeros(1,network.GetNumberOfLinks);
            
           
            cluster = ClusterMonophasique(invadedPores,interfaceLinks,...
                interfacePoreOutward,booleanInvadedLinks,...
                network,criticalPressures,clusterOptions);
            
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

