classdef PoreNetworkMesh < PoreNetworkEuclidien
    %PoreNetworkMesh Sous-classe de PoreNetworkEuclidien
    %   Reseau de pore associe a un maillage. Les pores et liens sont
    %   formes par les mailles d'un maillage. On a donc des informations sur
    %   les sommets du maillage.
    
    
    properties %(SetAccess  =  protected, GetAccess  =  protected)
        Vertices
        Faces %structure, Faces{i} = tableau des vertices de la face i dans l'ordre des sommets
        NombreVertices
        CellsToVertices
        VerticeDataList
    end
    
    methods

        function pore_network_mesh = PoreNetworkMesh(dimension,faces,pores,cells_to_vertices,owners,neighbours,boundaries,vertices,myGeometry)
            %constructeur
            %input : dimension,faces,pores,cells_to_vertices,owners,neighbours,boundaries,vertices
            %output : pore_network_mesh
            
            nPore = length(pores);
            poreCenter = zeros(nPore,dimension);
            for iPore = 1:nPore
                poreCenter(iPore,:) = mean(vertices(cells_to_vertices{iPore},:));
            end
            
            nLink = length(faces);
            linkCenter = zeros(nLink,dimension);           
            for iLink = 1:nLink
                if dimension == 2
                    linkCenter(iLink,:) = mean(vertices([faces{iLink}(1),faces{iLink}(2)],:));
                else
                    linkCenter(iLink,:) = mean(vertices(faces{iLink},:));
                end
            end
             
            pore_network_mesh = pore_network_mesh@PoreNetworkEuclidien(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry);
                        
            pore_network_mesh.Vertices = vertices;
            pore_network_mesh.NombreVertices = length(vertices(:,1));
            pore_network_mesh.Faces = faces;
            pore_network_mesh.CellsToVertices = cells_to_vertices;
            dataVerticeList = DataVerticeList(pore_network_mesh.NombreVertices);
            pore_network_mesh.VerticeDataList = dataVerticeList;
        end
        
        
        
        function vertices = GetVerticesOfLink(poreNetwork,num_lien)
            %input : poreNetwork,num_lien
            %output : vertices
            if poreNetwork.GetDimension == 2
                vertices = poreNetwork.Vertices([poreNetwork.Faces{num_lien}(1),poreNetwork.Faces{num_lien}(2)],:);
            else
                vertices = poreNetwork.Vertices(poreNetwork.Faces{num_lien},:);
            end
        end
        
        function vertices = GetVerticesOfPore(poreNetwork,num_pore)
            %input : poreNetwork,num_pore
            %output : vertices
            vertices = poreNetwork.Vertices(poreNetwork.CellsToVertices{num_pore},:);
        end
        
        function verticeNumber = GetVerticesOfLinkNumber(poreNetwork,numLink)
            %input : poreNetwork,numLink
            %output : verticeNumber
            verticeNumber = poreNetwork.Faces{numLink};
        end         
        
        function verticeNumber = GetVerticesOfPoreNumber(poreNetwork,numPore)
            %input : poreNetwork,numPore
            %output : verticeNumber
            verticeNumber = poreNetwork.CellsToVertices{numPore};
        end   
         
        function vertice = GetVertice(poreNetwork,numVertice)
            %input : poreNetwork,numVertice
            %output : vertice
            vertice = poreNetwork.Vertices(numVertice,:);
        end
        
        function number = GetNumberOfVertices(poreNetwork)
            %input : poreNetwork
            %output : number
            number = poreNetwork.NombreVertices;
        end
        
        
        
        function dataStruct = GetVerticeDataList(poreNetwork)
            %input : poreNetwork
            %output : dataStruct
            dataStruct = poreNetwork.VerticeDataList.VerticeDatas;
        end       
        
        function AddNewVerticeData(poreNetwork,data,name)
            %input : poreNetwork,data,name
            poreNetwork.VerticeDataList.AddData(data,name);
        end
        
        function RemoveVerticeData(poreNetwork,name)
            %input : poreNetwork,name
            poreNetwork.VerticeDataList.RemoveData(name);
        end  
        
        
        function coord = GetAllVerticesCoordinates(poreNetwork)
            %input : poreNetwork
            %output : coord
            coord = poreNetwork.Vertices;
        end
        
        function output_struct = PrivateInternalOutputStruct(pore_network_mesh)
            %r�cup�ration du l'output_struct par la m�thode de la superclasse g�om�trie
            %input : pore_network_mesh
            %output : output_struct
            
            output_struct = pore_network_mesh.PrivateInternalOutputStruct@PoreNetwork;
            %ajouter champs manquants et changer ATTRIBUTE.Type
            output_struct.Vertices = pore_network_mesh.Vertices;
            output_struct.Faces.Face = pore_network_mesh.Faces;
            output_struct.ATTRIBUTE.Type = 'PoreNetworkUnstructuredMesh';
            output_struct.ATTRIBUTE.NombreVertices = pore_network_mesh.NombreVertices;
            output_struct.CellsToVertices = pore_network_mesh.CellsToVertices;
            
        end
        
    end
    
end
