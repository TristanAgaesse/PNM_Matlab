classdef  PoreNetwork
    %PoreNetwork : structure representant un reseau de pores fondamental,
    %contenant de relations entre pores et liens de type graph
    
    
    properties %(SetAccess  =  protected, GetAccess  =  protected)
        LinkOwners  %%tableau 1*NombreLiens, LinkOwners(i) = num�ro d'un pore voisin
        LinkNeighbours  %tableau 1*NombreLiens, LinkNeighbours(i) = -1 si face sur une fronti�re
        Pores  %structure, Pores{i} = tableau avec numero des liens de la cellule i
        Boundaries  %structure Boundaries.Boundary{i} = infos structur�es sur une fronti�re
        Dimension
        NombreLiens
        NombrePores
        PoreDataList
        LinkDataList
        MacroscopicGeometry
    end
       
    methods
        
        %Constructeur
        function network = PoreNetwork(dimension,pores,owners,neighbours,boundaries,myGeometry)
            %Constructeur
            %input : dimension,pores,owners,neighbours,boundaries
            %output : poreNetwork
            network.Dimension = dimension;
            network.NombreLiens = length(owners);
            network.Pores = pores;
            network.NombrePores = length(pores);
            network.LinkOwners = owners;
            network.LinkNeighbours = neighbours;
            network.Boundaries = boundaries;
            network.MacroscopicGeometry=myGeometry;
            
            data_pore_list = DataPoreList(network.NombrePores);
            network.PoreDataList = data_pore_list;
            data_link_list = DataLinkList(network.NombreLiens);
            network.LinkDataList = data_link_list;
            
        end
        
        %Some Getters
        
        function number = GetNumberOfPores(poreNetwork)
            %input : poreNetwork
            %output : number
            number = poreNetwork.NombrePores;
        end
        
        function number = GetNumberOfLinks(poreNetwork)
            %input : poreNetwork
            %output : number
            number = poreNetwork.NombreLiens;
        end
        
        function dimension = GetDimension(poreNetwork)
            %input : poreNetwork
            %output : dimension
            dimension = poreNetwork.Dimension;
        end
        
        function nBoundary=GetNumberOfBoundaries(network)
            nBoundary=length(network.Boundaries.Boundary);
        end
        
        
        %Topology
        
        function liste_pores_voisins = GetPoresVoisinsOfPore(poreNetwork,num_pore)
            %input : poreNetwork,num_pore
            %output :liste_pores_voisins
            
            assert(length(num_pore)==1);
            liste_liens = poreNetwork.Pores{num_pore};
            liste_pores_voisins = zeros(1,length(liste_liens));
            for i = 1:length(liste_liens)
                num_owner = poreNetwork.LinkOwners(liste_liens(i));
                num_neighbour = poreNetwork.LinkNeighbours(liste_liens(i));
                if num_owner == num_pore
                    liste_pores_voisins(i) = num_neighbour;
                else
                    liste_pores_voisins(i) = num_owner;
                end
            end
            
        end
        
        function indices = GetPoresOfLink(network,num_link)
            %input : network,num_link
            %output : indices : indices(1)=neighbourPore (-1 if
            %                   boundaryLink), indices(2)=ownerPore
            
            if size(num_link,1)~=1
                assert(size(num_link,2)==1 || isempty(num_link));
                num_link=transpose(num_link);
            end
            
            indices=transpose(vertcat(network.LinkNeighbours(num_link),network.LinkOwners(num_link)));
        end
        
        function indices = GetLinksOfPore(poreNetwork,num_pore)
            %input : poreNetwork,num_pore
            %output : indices
            if isempty(num_pore)
                indices=[];
            else
            	indices = poreNetwork.Pores{num_pore};
            end
        end

        function face = GetCommonFace(poreNetwork,pore1,pore2) 
            %input : poreNetwork,pore1,pore2
            %output : face
            
            face = poreNetwork.Pores{pore1}(ismember(poreNetwork.Pores{pore1},poreNetwork.Pores{pore2}));
            assert(length(face) == 1,'Pb GetCommonFace');
        end
        
        
        %Data lists handling
        
        function data = GetPoreData(network,dataName)
            data = network.GetPoreDataList.(dataName);
        end
        
        function data_struct = GetPoreDataList(poreNetwork)
            %input : poreNetwork
            %output : data_struct
            data_struct = poreNetwork.PoreDataList.PoreDatas;
        end       
        
        function AddNewPoreData(network,data,name)
            %input : poreNetwork,data,name
            network.PoreDataList.AddData(data,name);
        end
               
        function RemovePoreData(poreNetwork,name)
            %input : poreNetwork,name
            poreNetwork.PoreDataList.RemoveData(name);
        end        

        function data = GetLinkData(network,dataName)
            data = network.GetLinkDataList.(dataName);
        end
        
        function data_struct = GetLinkDataList(poreNetwork)
            %input : poreNetwork
            %output : data_struct
            data_struct = poreNetwork.LinkDataList.LinkDatas;
        end       
        
        function AddNewLinkData(poreNetwork,data,name)
            %input : poreNetwork,data,name
            poreNetwork.LinkDataList.AddData(data,name);
        end
        
        function RemoveLinkData(poreNetwork,name)
            %input : poreNetwork,name
            poreNetwork.LinkDataList.RemoveData(name);
        end
        
        function diameter = GetLinkDiameter(poreNetwork,iLink)
            %input : poreNetwork,iLink
            %output : diameter
            
            if size(iLink,1)~=1
                assert(size(iLink,2)==1 || isempty(iLink));
                iLink=transpose(iLink);
            end
            
            if not(isfield(poreNetwork.GetLinkDataList,'Diameter'))
                disp('Calcul du diam�tre des liens...');
                tic;
                nLink = poreNetwork.GetNumberOfLinks;
                diameters = zeros(1,nLink);
               for iLink = 1:nLink
                    diameters(iLink) = poreNetwork.ComputeLinkDiameter(iLink);
               end
               poreNetwork.AddNewLinkData(diameters,'Diameter');
               duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
               disp(sprintf('Calcul du diam�tre des liens termin�. Dur�e : %d minutes %f s.',minutes,secondes));    
            end
            
            diameter = poreNetwork.GetLinkDataList.Diameter(iLink);
        end
        
        function diameter = GetPoreDiameter(poreNetwork,iPore)
            %input : poreNetwork,iPore
            %output : diameter
            
            if size(iPore,1)~=1
                assert(size(iPore,2)==1 || isempty(iPore));
                iPore=transpose(iPore);
            end
            
            if not(isfield(poreNetwork.GetPoreDataList,'Diameter'))
            
                if not(isfield(poreNetwork.GetPoreDataList,'Volume'))
                    disp('Calcul du volume des pores...');tic;
                    volumes = ComputeAllPoreVolume(poreNetwork);
                    toc;
                    poreNetwork.AddNewPoreData(volumes,'Volume');
                    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
                    disp(sprintf('Calcul du volume des pores termin�. Dur�e : %d minutes %f s.',minutes,secondes));
                end
                
                volumes = poreNetwork.GetPoreDataList.Volume;
                diameters = (24*volumes/(4*pi)).^(1/3);
                poreNetwork.AddNewPoreData(diameters,'Diameter');
            end
                
            diameter = poreNetwork.GetPoreDataList.Diameter(iPore);
            
        end
        
        
        %Boundaries
        
        
        function liens_frontiere = GetLinksFrontiere(poreNetwork,num_frontieres)
            %renvoie la liste des liens attach�s � une liste de frontieres
            %macroscopiques. Si num_frontiere = 0, renvoie la liste des liens
            %internes.
            %input : poreNetwork,num_frontieres
            %output : liens_frontiere
            
            if size(num_frontieres,1)~=1
                assert(size(num_frontieres,2)==1 );
                num_frontieres=transpose(num_frontieres);
            end
            
            faces_frontieres = cell(1,length(num_frontieres));
            indice = 1;
            for i = num_frontieres
                if i == 0
                    start_face = poreNetwork.Boundaries.Boundary(end).ATTRIBUTE.StartFace...
                        +poreNetwork.Boundaries.Boundary(end).ATTRIBUTE.NombreFaces;
                    faces_frontieres{indice} = start_face:poreNetwork.GetNumberOfLinks;
                elseif i>0 && i <= length(poreNetwork.Boundaries.Boundary)
                    start_face = poreNetwork.Boundaries.Boundary(i).ATTRIBUTE.StartFace;
                    end_face = start_face+poreNetwork.Boundaries.Boundary(i).ATTRIBUTE.NombreFaces-1;
                    faces_frontieres{indice} = start_face:end_face;
                else
                    error('Trying to access links of non existing boundary');
                end
                indice = indice+1;
            end
            liens_frontiere = cell2mat(faces_frontieres);            
        end
        
        function numFrontiere = GetFrontiereOfLink(poreNetwork,numLink)
            %input : poreNetwork,numLink
            %output : numFrontiere
            numFrontiere = 0;
            nBoundary = length(poreNetwork.Boundaries.Boundary); 
            
            assert(length(numLink)==1);
            
            start = poreNetwork.Boundaries.Boundary(nBoundary).ATTRIBUTE.StartFace;
            finish = start+poreNetwork.Boundaries.Boundary(nBoundary).ATTRIBUTE.NombreFaces-1;
            if numLink>finish
                return
            end
            
            for iBoundary = 1:nBoundary
                start = poreNetwork.Boundaries.Boundary(iBoundary).ATTRIBUTE.StartFace;
                finish = start+poreNetwork.Boundaries.Boundary(iBoundary).ATTRIBUTE.NombreFaces-1;
                if and( numLink>=start, numLink<=finish )
                   numFrontiere = iBoundary ;
                   return
                end
            end
        end
        
        function pores_frontiere = GetPoresFrontiere(poreNetwork,linkFrontiere)
            %renvoie la liste des pores attaches a une liste de liens situes sur une frontière. 
            %Un pore peut apparaitre plusieurs fois dans la liste.
            %input : poreNetwork,linkFrontiere
            %output :pores_frontiere
            
            if size(linkFrontiere,1)~=1
                assert(size(linkFrontiere,2)==1 || isempty(linkFrontiere));
                linkFrontiere=transpose(linkFrontiere);
            end
            
            pores_frontiere = poreNetwork.LinkOwners(linkFrontiere);

        end
        
        function [boundaryLinks,innerLinks]=GetPoreRegionBoundaryLinks(network,poreList)
            %Input: network, poreList
            %Output: boundaryLinks,innerLinks
            
            if size(poreList,1)~=1
                assert(size(poreList,2)==1 || isempty(poreList));
                poreList=transpose(poreList);
            end
            
            linkCount = zeros(1,network.GetNumberOfLinks);
            for iPore = poreList
                thoseLinks = network.GetLinksOfPore(iPore); 
                linkCount(thoseLinks) = linkCount(thoseLinks)+1;
            end
            innerLinks = find(linkCount==2);
            boundaryLinks = find(linkCount == 1);
            assert(not(any(linkCount>2)));
        end
        
        
        
        %Utilities
        
        function cluster=CreateVoidCluster(network)
            cluster=ClusterMonophasique.CreateVoidCluster(network);
        end
        
        function cluster=CreateFullCluster(network)
            voidCluster=ClusterMonophasique.CreateVoidCluster(network);
            cluster=voidCluster.GetComplementaryCluster;
        end
        
        function linkData = InterpolatePoreDataToLink(network,poreData)
            %Returns a linkData obtained by interpolation of a poreData
            % For innerLinks, linkData = mean of poreData on the 2
            % neighboor pores. For boundaryLinks, linkData= poreData on the
            % 1 inner neighboor pore.
            %Input : network,poreData
            %Output : linkData       
            
            nPore = network.GetNumberOfPores;
            nLink = network.GetNumberOfLinks;
            assert(length(poreData)==nPore,'poreData must have a length network.GetNumberOfPores')
            linkData = zeros(nLink,1);
            
            innerLinks = network.GetLinksFrontiere(0);
            neighboorPores = network.GetPoresOfLink(innerLinks);
            linkData(innerLinks) = 0.5*(poreData(neighboorPores(:,1))+poreData(neighboorPores(:,2)));
            
            boundaryLinks = network.GetLinksFrontiere(1:network.GetNumberOfBoundaries);
            neighboorPores = network.GetPoresOfLink(boundaryLinks);
            linkData(boundaryLinks) = poreData(neighboorPores(:,2));
        end
        
        
        function adjacencyMatrix = CreateAdjacencyMatrix(network,invadedPoreList)
            %Returns the adjacency matrix of the sub-network given by a
            %list of invaded pores.
            %Input : network,invadedPoreList
            %Output : adjacencyMatrix
                
            if size(invadedPoreList,1)~=1
                assert(size(invadedPoreList,2)==1 || isempty(invadedPoreList));
                invadedPoreList=transpose(invadedPoreList);
            end
            
            nPore = network.GetNumberOfPores;

            internalLink=network.GetLinksFrontiere(0);
            theirPore=network.GetPoresOfLink(internalLink);

            firstIndice=theirPore(:,1);
            secondIndice=theirPore(:,2);

            newFirstIndice=vertcat(firstIndice,secondIndice); 
            newSecondIndice=vertcat(secondIndice,firstIndice); 
            value = ones(length(newFirstIndice),1);

            adjacencyMatrix = sparse(newFirstIndice,newSecondIndice,value,nPore,nPore);

            adjacencyMatrix=adjacencyMatrix(invadedPoreList,invadedPoreList);
        end
        
        
        function conComp = FindComposantesConnexes(network,liste_pores_envahis,varargin)
            %Retourne les composantes connexes des zones envahies par une
            %phase. 
            %input : poreNetwork,liste_pores_envahis
            %output : conComp : cell, conComp{i}=liste des pores 
            %               de la composante i
            
            option='HomeMadeTarjan';
            if ~isempty(varargin)
                option=varargin{1};
            end
            
            if size(liste_pores_envahis,1)~=1
                assert(size(liste_pores_envahis,2)==1);
                liste_pores_envahis=transpose(liste_pores_envahis);
            end
            
            if not(verLessThan('matlab','8.5')) 
                adjacencyMatrix = network.CreateAdjacencyMatrix(liste_pores_envahis);
                myGraph=graph(adjacencyMatrix);
                bins = conncomp(myGraph);
                nConnectedComponents=max(bins);
                conComp = cell(1,nConnectedComponents);
                for iComponent = 1:nConnectedComponents
                    conComp{iComponent} = liste_pores_envahis(bins == iComponent);
                end
                
                
            elseif strcmp(option,'graphconncomp')
                %Uses a free C++ implementation of the Matlab function
                %graphconncomp
                adjacencyMatrix = network.CreateAdjacencyMatrix(liste_pores_envahis);
                
                [labels,nConnectedComponents] = graph_conn_comp(adjacencyMatrix);
                
                conComp = cell(1,nConnectedComponents);
                for iComponent = 1:nConnectedComponents
                    conComp{iComponent} = liste_pores_envahis(labels == iComponent);
                end
                
                
            elseif strcmp(option,'HomeMadeTarjan')
                %algo : parcours en profondeur d'un amas
                %impl�mentation avec pile LIFO contenant les pores � explorer :  
                %on d�pile un pore et on empile ses voisins non encore explor�s.
                
                iComponent = 1;
                pores_envahis_explores = zeros(1,length(liste_pores_envahis));
                %pores_envahis_explores(i) = 0 si liste_pores_envahis(i) n'est
                %pas encore explor�,  = num_composante sinon

                %boucle sur les composantes connexes
                stack_pores_a_explorer = java.util.Stack();
                while not(length(find(pores_envahis_explores)) == length(liste_pores_envahis))                    

                    foo = find(pores_envahis_explores == 0);
                    stack_pores_a_explorer.push(foo(1)); 

                    while ~stack_pores_a_explorer.empty
                        ref_pore = stack_pores_a_explorer.pop;
                        %ref_pore : position dans la liste_pore_envahis   
                        if pores_envahis_explores(ref_pore) == 0
                            ref_voisins_envahis = find(ismember(liste_pores_envahis,...
                                network.GetPoresVoisinsOfPore(liste_pores_envahis(ref_pore))));
                            %les voisins sont les pores ayant un lien commun,
                            %on suppose que les liens sont automatiquement
                            %envahis
                            ref_voisins_envahis_non_explores = ref_voisins_envahis(...
                                ismember(ref_voisins_envahis,find(~pores_envahis_explores)));   
                            for i = 1:length(ref_voisins_envahis_non_explores)
                                stack_pores_a_explorer.push(ref_voisins_envahis_non_explores(i));
                            end
                            pores_envahis_explores(ref_pore) = iComponent;
                        end
                    end
                    iComponent = iComponent+1;
                end

                conComp = cell(1,max(pores_envahis_explores));
                for iComponent = 1:max(pores_envahis_explores)
                    conComp{iComponent} = liste_pores_envahis(pores_envahis_explores == iComponent);
                end
            
            end
        end
        
        function output_struct = PrivateInternalOutputStruct(poreNetwork)
            %cr�ation de la structure output qui �crite dans le
            %fichier r�seau de pores ou visualis�e
            %input : poreNetwork
            %output : output_struct
            
            output_struct = struct;
            
            output_struct.LinkOwners = poreNetwork.LinkOwners;
            output_struct.LinkNeighbours = poreNetwork.LinkNeighbours;               
            output_struct.Boundaries = poreNetwork.Boundaries;            
            output_struct.Cells.Cell = poreNetwork.Pores;    
            
            attribute.Dimension = poreNetwork.Dimension;
            attribute.NombreLiens = poreNetwork.NombreLiens;
            attribute.NombrePores = poreNetwork.NombrePores;
            attribute.Type = 'PoreNetwork';
            output_struct.ATTRIBUTE = attribute;
            
        end
        
    end

end
