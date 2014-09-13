classdef PoreNetworkMeshFibrous < PoreNetworkMesh
    %PoreNetworkMeshFibrous Sous-classe de PoreNetworkMesh.
    %   Reseau de pores genere par un maillage, avec les parties solides
    %   constitues par des fibres cylindriques sur les aretes. On a des informations sur
    %   les aretes, et des outils pour calculer les proprietes geometriques
    %   des pores et liens en fonction des fibres.
    
    
    
    properties %(SetAccess = protected, GetAccess = protected)
        Edges %Tableau 2*NombreEdges contenant les num�ros des vertices d�finissant les edges.
        NombreEdges
        VerticesToEdges
        FacesToEdges
        EdgeDataList
        EdgesToFaces
    end
    
    
    
    methods                       
        
        function network=PoreNetworkMeshFibrous(dimension,faces,pores,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,epaisseur_edges,faces_to_edges,edges_to_faces,myGeometry)
            %constructeur � partir de l'input geometrie_macroscopique
            %doit renseigner toutes les properties de la geometrie consid�r�e
            %input : dimension,faces,pores,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,epaisseur_edges,faces_to_edges,edges_to_faces
            %output : network
            
            network = network@PoreNetworkMesh(dimension,faces,pores,cells_to_vertices,owners,neighbours,boundaries,vertices,myGeometry);

            network.Edges = edges;
            network.NombreEdges = length(edges(:,1));
            dataEdgeList = DataEdgeList(network.NombreEdges);
            network.EdgeDataList = dataEdgeList;
            dataEdgeList.AddData(epaisseur_edges,'FiberDiameter');
            network.VerticesToEdges = vertices_to_edges;
            network.FacesToEdges = faces_to_edges;
            network.EdgesToFaces = edges_to_faces;
        end
        
        
        
        function vertices = GetVerticesOfEdge(network,numEdge)
            %input : network,numEdge
            %output : vertices
            vertices = network.Vertices(network.Edges(numEdge,:),:);
        end
        
        
        
        function vertNum = GetVerticesOfEdgeNumber(network,numEdge)
            %input : network,numEdge
            %output : vertNum
            vertNum = network.Edges(numEdge,:);
        end
           
        
        
        function edges = GetEdgesOfLink(network,numLink)
            %input : network,numLink
            %output : edges
            edges = network.FacesToEdges{numLink};
        end
        
        
        
        function edges = GetEdgesOfPore(network,numPore)
            %input : network,numPore
            %output : edges
            links = network.GetLinksOfPore(numPore);
            edges = [];
            for iLink = links
                edges = [edges,network.GetEdgesOfLink(iLink)];
            end
            edges = unique(edges);
        end
        
        
        
        function edges = GetEdgesOfVertice(network,numVertice)
            %input : network,numVertice
            %output : edges
            edges = network.VerticesToEdges{numVertice};
        end
        
        
        
        function links = GetLinksOfEdge(network,numEdge)
            %input : network,numEdg
            %output : links
            
            links = network.EdgesToFaces{numEdge};
        end
        
        
        
        function number = GetNumberOfEdges(network)
            %input : network
            %output : number
            number = network.NombreEdges;
        end
        
        
        
        function neighbourLinks = GetNeighbourLinksOfEdgeInPore(network,numPore,numEdge)
            %input : network,numPore,numEdge
            %output : neighbourLinks
            poreLinks = network.GetLinksOfPore(numPore);
            
            linksOfEdge = network.GetLinksOfEdge(numEdge);
            
            neighbourLinks = linksOfEdge(ismember(linksOfEdge,poreLinks));
            assert(length(neighbourLinks) == 2);
        end
        
        
        
        function neighbourEdges = GetNeighbourEdgesOfVerticeInPore(network,numPore,numVertice)
            %input : network,numPore,numVertice
            %output : neighbourEdges
            neighbourEdges = [];
            edgesOfVertice = network.GetEdgesOfVertice(numVertice);
            verticesOfPore = network.GetVerticesOfPore(numPore);
            
            for iEdge = edgesOfVertice
                edgeVertices = network.GetVerticesOfEdgeNumber(iEdge);
                otherVertice = edgeVertices(edgeVertices~=numVertice);
                if ismember(otherVertice,verticesOfPore)
                    neighbourEdges = [neighbourEdges,iEdge];
                end
                
            end
        end
        
        
        
        function dataStruct = GetEdgeDataList(network)
            %input : network
            %output : dataStruct
            dataStruct = network.EdgeDataList.EdgeDatas;
        end       
        
        
        
        function AddNewEdgeData(network,data,name)
            %input : network,data,name
            network.EdgeDataList.AddData(data,name);
        end
        
        
        
        function RemoveEdgeData(network,name)
            %input : network,name
            network.EdgeDataList.RemoveData(name);
        end  
        
        
        
        function volume = ComputePoreVolume(network, numPore)
            %input : network, numPore
            %output : volume
            methode = 'approximate_intersections';
            
            fibreDiameter = network.GetEdgeDataList.('FiberDiameter');
            vertices = network.GetVerticesOfPore(numPore);
            
            switch network.Dimension
                
                case 2
                    volume = abs(polygonArea(vertices));

                case 3
                    centrePore = mean(vertices);
                    edges = network.GetEdgesOfPore(numPore);
                    
                    if strcmp(methode,'hardcore_clipping')
                        nodes = vertices;
                        links = network.GetLinksOfPore(numPore);
                        nFace = length(links);
                        faces = minConvexHull(nodes);

                        for iEdge = edges
                            neighbourFaces = network.GetNeighbourLinksOfEdgeInPore(numPore,iEdge);

                            foo = network.GetVerticesOfLink(neighbourFaces(1)) ;
                            pointsPlan1 = foo([1 2 3 3],:);
                            vect1 = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(pointsPlan1,centrePore);

                            foo = network.GetVerticesOfLink(neighbourFaces(2)) ;
                            pointsPlan2 = foo([1 2 3 3],:);
                            vect2 = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(pointsPlan2,centrePore);

                            dihedralAngle = pi-acos(dot(vect1,vect2));
                            planeNormal = -(vect1+vect2)/norm(vect1+vect2);
                            
                            distance = fibreDiameter(iEdge)/2*sqrt(dihedralAngle/(2*tan(dihedralAngle/2))) ;
                            
                            foo = network.GetVerticesOfEdge(iEdge);
                            P0 = foo(1,:)+distance*planeNormal;
                            
                            plane = createPlane(P0,-planeNormal);
                            try 
                                [nodes, faces]  =  clipConvexPolyhedronHP(nodes, faces, plane);
                            catch
                                    nodes = 0;
                                    disp(sprintf('Probleme en calculant volume du pore %d',numPore));
                            end
                        end
                        if size(nodes,1)<3
                            
                            volume = 0;
                            
                        else
                            
                            [~, volume] = convhulln(nodes);
                            
                        end
                        
                    elseif strcmp(methode,'approximate_intersections')
                        
                        [~, volume] = convhulln(vertices);
                        
                        nEdge = length(edges);
                        dihedralAngle = zeros(1,nEdge);
                        
                        for i=1:nEdge
                            %enlever le volume de chaque fibre
                            neighbourFaces = network.GetNeighbourLinksOfEdgeInPore(numPore,edges(i));

                            foo = network.GetVerticesOfLink(neighbourFaces(1)) ;
                            pointsPlan1 = foo([1 2 3 3],:);
                            vect1 = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(pointsPlan1,centrePore);

                            foo = network.GetVerticesOfLink(neighbourFaces(2)) ;
                            pointsPlan2 = foo([1 2 3 3],:);
                            vect2 = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(pointsPlan2,centrePore);
                        
                            dihedralAngle(i) = pi-acos(dot(vect1,vect2));
                            
                            edgeVertice = network.GetVerticesOfEdge(edges(i));
                            edgeLength = norm(edgeVertice(1,:)-edgeVertice(2,:));
                            
                            volumeFibre = pi*(fibreDiameter(edges(i))/2)^2*edgeLength*dihedralAngle(i)/(2*pi);
                            
                            volume = volume-volumeFibre;
                            
                        end
                        
                        verticesNumber = network.GetVerticesOfPoreNumber(numPore);
                        for iVertice = verticesNumber
                            %rajouter le volume approximatif des 
                            %intersections entre fibres au niveau des
                            %sommets
                            
                            neighbourEdges = network.GetNeighbourEdgesOfVerticeInPore(numPore,iVertice);
                            nNeighbourEdge = length(neighbourEdges);
                            [fdMax,iMax] = max(fibreDiameter(neighbourEdges));
                            
                            for iEdge = neighbourEdges(setdiff((1:nNeighbourEdge),iMax))
                                angle = dihedralAngle(edges==iEdge);
                                volume = volume+pi*(fibreDiameter(iEdge)/2)^2*(fdMax/2)*angle/(2*pi);
                            end
                            
                            if volume<0
                               disp(sprintf('Probleme en calculant volume du pore %d',numPore)); 
                               volume = 0; 
                            end
                            
                        end
                    end
            end
        end
        
        
        
        function volume = ComputeAllPoreVolume(network)
            %input : network
            %output : volume
            nPore = network.GetNumberOfPores;
            volume = zeros(1,nPore);
            for iPore=1:nPore
                volume(iPore) = network.ComputePoreVolume(iPore);                
            end
        end        
        
        
        
        function diameter = ComputeAllLinkDiameter(network)
            %input : network
            %output : diameter
            nLink = network.GetNumberOfLinks;
            diameter = zeros(1,nLink);
            
            for iLink=1:nLink
                diameter(iLink) = network.ComputeLinkDiameter(iLink);                
            end
        end  

        function linkDiameter = ComputeLinkDiameter(network, numLink)
            %Calcule le diametre equivalent d'un lien en extrudant les
            %fibres
            %input : network, numLink
            %output : linkDiameter
            
            dimension = network.Dimension;
            
            numVertices = network.Faces{numLink};
            numEdges = network.FacesToEdges{numLink};
            vertices = network.Vertices(numVertices,:);
            diametreFibres=network.GetEdgeDataList.FiberDiameter(numEdges);
            

            switch network.Dimension
                case 2
                    verticesExtrudes = PoreNetworkMeshFibrous.ComputeExtrudePolygonParFibre(dimension, vertices, diametreFibres);
                    if ~isempty(verticesExtrudes)
                        linkDiameter = norm(verticesExtrudes(1,:)-verticesExtrudes(2,:));
                    else
                        linkDiameter = 0;
                    end
                case 3
                    verticesExtrudes = PoreNetworkMeshFibrous.ComputeExtrudePolygonParFibre(dimension, vertices, diametreFibres);
                    surface = abs(polygonArea(verticesExtrudes));%librairie mathgeom
                    linkDiameter = sqrt(4*surface/pi); %diam�tre �quivalent
            end
        end
        
        
        
        
        
        function outputStruct = PrivateInternalOutputStruct(network)
            %cr�ation en m�moire de la structure output qui pourra �tre 
            %visualis�e par le visualisateur interne
            %input : network
            %output : outputStruct
            outputStruct = network.PrivateInternalOutputStruct@PoreNetworkMesh;
            %ajouter champs manquants et changer ATTRIBUTE.Type
            outputStruct.Edges = network.Edges;
            outputStruct.EdgeData.FiberDiameter = network.EdgeDataList.EdgeDatas.FiberDiameter;
            outputStruct.ATTRIBUTE.Type = 'PoreNetworkUnstructuredMeshFibrous';
            outputStruct.ATTRIBUTE.NombreEdges = network.NombreEdges;
        end
        
        
        function vtk_struct = PrivateVTKOutputStructPolydataMesh(network)
            %cr�ation en m�moire de la structure d'un fichier VTK POLYDATA 
            %pour afficher le r�seau de pores dans paraview.
            %input : network
            %output :vtk_struct
            vertices = network.Vertices;
            edges = network.Edges;
            faces = network.Faces;
            cells = network.Pores;
            pore_data = network.GetPoreDataList;
            link_data = network.GetLinkDataList;
            edge_data = network.GetEdgeDataList;
            vertice_data = network.GetVerticeDataList;
            dimension = network.Dimension;

            %Partie POINTS du fichier
            nCell = length(cells);
            nombre_points_pour_cells = 0;
            for iCell = 1:nCell
                liste_faces = cells{iCell};
                for iFace = liste_faces
                   nombre_points_pour_cells = nombre_points_pour_cells+length(faces{iFace}); 
                end
            end
            nEdge = length(edges(:,1));
            nVertice = length(vertices(:,1));
            nPoints = 2*nEdge+sum(cellfun('length',faces))+nombre_points_pour_cells+nVertice;
            points_output = zeros(nPoints,3);
                        
            edge_to_points = cell(1,nEdge);
            x_extension = max(vertices(:,1))-min(vertices(:,1));
            y_extension = max(vertices(:,2))-min(vertices(:,2));
            z_extension = (x_extension+y_extension)/20;
            for iEdge=1:nEdge
                if dimension==3
                    points_output((2*iEdge-1),:) = vertices(edges(iEdge,1),:);
                    points_output((2*iEdge),:) = vertices(edges(iEdge,2),:);
                else
                    points_output((2*iEdge-1),:) = [vertices(edges(iEdge,1),:),0];
                    points_output((2*iEdge),:) = [vertices(edges(iEdge,2),:),-z_extension];
                end
                edge_to_points{iEdge} = [(2*iEdge-1),(2*iEdge)];
            end
          
            compteur_debut = 2*nEdge+1;

            nFace = length(faces);
            face_to_points = cell(1,nFace);
            for iFace=1:nFace
                compteur_fin = compteur_debut+length(faces{iFace})-1;
                if dimension==3
                    points_output((compteur_debut:compteur_fin),:) = vertices(faces{iFace},:);
                else
                    points_output((compteur_debut:compteur_fin),:) = horzcat(vertices(faces{iFace},:),zeros(length(faces{iFace}),1));
                end
                face_to_points{iFace} = (compteur_debut):(compteur_fin);
                compteur_debut = compteur_fin+1;
            end
            
            cell_to_points = cell(1,nCell);
            for iCell=1:nCell
                liste_faces = cells{iCell};
                cell_to_points{iCell} = cell(1,length(liste_faces));
                centre = mean(vertices(network.CellsToVertices{iCell},:));
                i = 0;
                for iFace=liste_faces
                    i = i+1;
                    compteur_fin = compteur_debut+length(faces{iFace})-1;
                    %new_vertices = 0.999*(vertices(faces{iFace},:)-ones(length(vertices(faces{iFace})),1)*centre)+ones(length(vertices(faces{iFace})),1)*centre;
                    if dimension==3
                        %points_output((compteur_debut:compteur_fin),:) = new_vertices;
                        points_output((compteur_debut:compteur_fin),:) = vertices(faces{iFace},:);
                    else
                        %points_output((compteur_debut:compteur_fin),:) = horzcat(new_vertices,zeros(length(faces{iFace}),1));
                        points_output((compteur_debut:compteur_fin),:) = horzcat(vertices(faces{iFace},:),zeros(length(faces{iFace}),1));
                    end                   
                    cell_to_points{iCell}{i} = (compteur_debut):(compteur_fin);
                    compteur_debut = compteur_fin+1;
                end
            end
            
            compteur_fin = compteur_debut+nVertice-1;
            if dimension==3
                points_output((compteur_debut:compteur_fin),:) = vertices;
            else
                points_output((compteur_debut:compteur_fin),:) = horzcat(vertices,zeros(nVertice,1));
            end
            vertices_to_points = cell(1,nVertice);
            for i=1:nVertice
                vertices_to_points{i} = compteur_debut+i-1;
            end
            compteur_debut = compteur_fin+1;
            
            assert(compteur_debut==nPoints+1);
            
            %Partie LINES du fichier
            nLine = nEdge;
            lines_output = zeros(nLine,3);
            for iLine=1:nLine
                lines_output(iLine,1) = 2;
                lines_output(iLine,2) = edge_to_points{iLine}(1)-1;
                lines_output(iLine,3) = edge_to_points{iLine}(2)-1;
            end
            
            %Partie POLYGONS du fichier
            nCellPolygone = sum(cellfun('length',cells));
            nPolygone = nFace+nCellPolygone;
            polygons_output = cell(nPolygone,1);
            face_to_polygon = cell(1,nFace);
            for iFace=1:nFace
                foo = face_to_points{iFace}-1;
                polygons_output{iFace} = [length(foo),foo];
                face_to_polygon{iFace} = iFace;
            end
            compteur = nFace+1;
            cell_to_polygon = cell(1,nCell);
            for iCell=1:nCell
                cell_to_polygon{iCell} = compteur:(compteur+length(cells{iCell})-1);
                for i=1:length(cells{iCell})
                    foo = cell_to_points{iCell}{i}-1;
                    polygons_output{compteur} = [length(foo),foo];
                    compteur = compteur+1;
                end
            end
            assert((compteur-1)==length(polygons_output));
            
            %Partie VERTICES du fichier
            vertices_output = horzcat(ones(nVertice,1),transpose((nPoints-nVertice+1):nPoints));
            
            
            %Partie POINT_DATA du fichier
            point_data_output = struct;
                %transcription des edge data
            names = fieldnames(edge_data);
            for i=1:length(names)
                data_output = zeros(nPoints,1);
                data = edge_data.(names{i});
                for iEdge=1:nEdge
                    data_output(edge_to_points{iEdge}(1)) = data(iEdge);
                    data_output(edge_to_points{iEdge}(2)) = data(iEdge);
                end
                point_data_output.(strcat('Edge_',names{i})) = data_output;
            end
                %transcription des vertices data
            names = fieldnames(vertice_data);
            for i=1:length(names)
                data_output = zeros(nPoints,1);
                data = vertice_data.(names{i});
                for iVertice=1:nVertice
                    data_output(vertices_to_points{iVertice}) = data(iVertice);
                end
                point_data_output.(strcat('Vertice_',names{i})) = data_output;
            end    
                
              
            
            %Partie CELL_DATA du fichier
            nCellData = nLine+nPolygone+nVertice;
            cell_data_output = struct;
                %transcription des link data
            names = fieldnames(link_data);
            for i=1:length(names)
                data_output = zeros(nCellData,1);
                data = link_data.(names{i});
                for iFace=1:nFace
                    data_output(face_to_polygon{iFace}+nLine+nVertice) = data(iFace);
                end
                cell_data_output.(strcat('Link_',names{i})) = data_output;
            end
                %transcription des pore data
            names = fieldnames(pore_data);
            for i=1:length(names)
                data_output = zeros(nCellData,1);
                data = pore_data.(names{i});
                for iCell=1:nCell
                    data_output(cell_to_polygon{iCell}+nLine+nVertice) = data(iCell)*ones(length(cell_to_polygon{iCell}),1);
                end
                cell_data_output.(strcat('Pore_',names{i})) = data_output;
            end

            vtk_struct.Points = points_output;
            vtk_struct.Lines = lines_output;
            vtk_struct.Polygons = polygons_output;
            vtk_struct.Vertices = vertices_output;
            vtk_struct.PointData = point_data_output;
            vtk_struct.CellData = cell_data_output;
        end       
        
        
        function ExportToParaview(network,filename,varargin)
            %Cree un fichier de visualisation pour Paraview
            %input : - network
            %        - filename 
            %        - varargin : string optionnelle contenant les options
            %        'BallAndStick' ou 'Mesh'. Par défaut c'est 'Mesh' pour
            %        un PoreNetworkMeshFibrous
            
            writer = FileWriterVTK(filename);
            tic;
            disp('G�n�ration du fichier VTK...');
            if strcmp(varargin,'BallAndStick')
                donnees_output = network.PrivateVTKOutputStructBallAndStick;
            else
                donnees_output = network.PrivateVTKOutputStructPolydataMesh;
            end
            writer.Write(donnees_output);
            duree=toc;minutes=floor(duree/60);secondes=duree-60*minutes;
            disp(sprintf('Fichier VTK g�n�r�. Dur�e : %d minutes %f s.',minutes,secondes));
            writer.delete;clear('writer','file_name','file_content','donnees_output','duree','minutes','secondes');
        end
                       
    end
    
    
    methods (Static = true)
    
    
        function verticesExtrudes=ComputeExtrudePolygonParFibre(dimension, vertices, diametreFibres)
            %renvoie les sommets du polygone obtenu en extrudant une face
            %par des fibres.
            %input : dimension, vertices, diametreFibres
            %output : verticesExtrudes
            
            
            
            switch dimension
                case 2
                    foo = vertices(2,:)-vertices(1,:);
                    largeur = norm(foo);
                    foo = foo/largeur;
                    
                    r1 = diametreFibres(1)/2;
                    r2 = diametreFibres(2)/2;
                    if largeur>r1+r2
                        verticesExtrudes = zeros(2,2);
                        verticesExtrudes(1,:) = vertices(1,:)+r1*foo;
                        verticesExtrudes(2,:) = vertices(2,:)-r2*foo;
                    else
                        verticesExtrudes = [];
                    end
                    
                case 3
                    methode = 'hardcore_clipping';
                    
                    vertices = PolygoneProjetterDansPlan(vertices);
                    
                    ordreTrigonometrique=PolygoneOrdonnerSommetsDansPlan(vertices);
                    vertices = vertices(ordreTrigonometrique,:);
                    diametreFibres=diametreFibres(ordreTrigonometrique);
                    
                    polygon=vertices;
                    
                    nVertice=length(vertices(:,1));
                    if strcmp(methode,'hardcore_clipping')
                        for i=1:nVertice
                            thisVertice = vertices(i,:);
                            if i<nVertice
                                otherVertice = vertices(i+1,:);
                            else
                                otherVertice = vertices(1,:);
                            end
                            rayonFibre = diametreFibres(i)/2;
                            
                            fibreAxis=createLine(thisVertice, otherVertice);
                            line = parallelLine(fibreAxis,-rayonFibre);
                            polygon = clipPolygonHP(polygon, line);%librairie mathgeom
                        end 
                        verticesExtrudes = polygon;
                        
%                     elseif strcmp(methode,'smart_clipping')
%                         %Deplacer les sommets de la face au niveau de
%                         %l'intersection entre les deux fibres voisines
%                         verticesExtrudes = polygon;
%                         
%                         for i=1:length(numVertices)
%                             
%                             thisVerticeNumber = numVertices(i);
%                             thisVertice = network.Vertices(thisVerticeNumber,:);
%                             if i<length(numVertices)
%                                 nextVerticeNumber = numVertices(i+1);
%                             else
%                                 nextVerticeNumber = numVertices(1);
%                             end
%                             nextVertice = network.Vertices(nextVerticeNumber,:);
%                             if i>1
%                                 lastVerticeNumber = numVertices(i-1);
%                             else
%                                 lastVerticeNumber = numVertices(end);
%                             end
%                             lastVertice = network.Vertices(lastVerticeNumber,:);
%                             
%                             numNextFibre = numEdges(i);
%                             diametreNextFibre = network.GetEdgeDataList.FiberDiameter(numFibre);
%                             if i==1
%                                 numLastFibre = length(numVertices);
%                             else
%                                 numLastFibre = i-1;
%                             end
%                             diametreLastFibre = network.GetEdgeDataList.FiberDiameter(numLastFibre);
%                             
%                             angle = acos(nextVertice-thisVertice,lastVertice-thisVertice);
%                             L = sqrt(((diametreNextFibre-diametreLastFibre*cos(pi-angle))/sin(pi-angle))^2+diametreLastFibre^2);
%                             %verticesExtrudes(i,:);

                        %Verifier s'il n'y a pas d'arête complètement
                        %recouverte, auquel cas l'étape précédente doit
                        %être modifiée                       
                    end
            end
        end
        
        
        
    end
    
    
end
