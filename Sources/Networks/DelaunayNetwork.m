classdef DelaunayNetwork
    %DelaunayNetwork Reseau de pores construit a partir d'une triangulation de
    %Delaunay. Les pores sont les points de base et les liens sont les
    %aretes de la triangulation
    
    properties
        Edges %Tableau 2*NombreEdges contenant les num�ros des vertices d�finissant les edges.
        NombreEdges
        VerticesToEdges
        FacesToEdges
        EdgeDataList
        Cells  %structure, Cells{i}=tableau avec numero des faces de la cellule i
        Boundaries  %fronti�res ext�rieures, structure Boundaries.Boundary(i)=infos structur�es pour l'output r�seau de pores
        FaceOwners
        FaceNeighbours
        Dimension
        NombreFaces
        NombreCells
        Vertices
        Faces %structure, Faces{i}=tableau des vertices de la face i dans l'ordre des sommets
        NombreVertices
        CellsToVertices
        VerticeDataList
        EdgesToFaces
        MacroscopicGeometry
    end
    
    methods
        %Constructeur
        function delaunayNetwork=DelaunayNetwork(dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,epaisseur_edges,faces_to_edges,edges_to_faces,myGeometry)
            delaunayNetwork.Dimension=dimension;
            
            delaunayNetwork.Vertices=vertices;
            delaunayNetwork.NombreVertices=length(vertices(:,1));
            
            delaunayNetwork.Cells=cells;
            delaunayNetwork.NombreCells=length(cells);
            delaunayNetwork.CellsToVertices=cells_to_vertices;
            
            delaunayNetwork.Faces=faces;
            delaunayNetwork.NombreFaces=length(owners);
            delaunayNetwork.FaceOwners=owners;
            delaunayNetwork.FaceNeighbours=neighbours;
            
            delaunayNetwork.Boundaries=boundaries;

            delaunayNetwork.Edges=edges;
            delaunayNetwork.NombreEdges=length(edges(:,1));
            delaunayNetwork.VerticesToEdges=vertices_to_edges;
            delaunayNetwork.FacesToEdges=faces_to_edges;
            delaunayNetwork.EdgesToFaces=edges_to_faces;
            delaunayNetwork.MacroscopicGeometry=myGeometry;
            
            dataVerticeList=DataVerticeList(delaunayNetwork.NombreVertices);
            delaunayNetwork.VerticeDataList=dataVerticeList;
            dataEdgeList=DataEdgeList(delaunayNetwork.NombreEdges);
            delaunayNetwork.EdgeDataList=dataEdgeList;
            dataEdgeList.AddData('DiametreCanal',epaisseur_edges);
            
        end

        %Fonctions pour acc�der aux propri�t�s du maillage
        
        function number=GetNumberOfEdges(network)
            number=network.NombreEdges;
        end
        
        function number=GetNumberOfVertices(network)
            number=network.NombreVertices;
        end
        
        function vertices=GetVerticesOfEdge(network,num_edge)
            %renvoie les coordonnees des deux sommets de l'ar�te, une ar�te
            %par ligne
            vertices=network.Vertices(network.Edges(num_edge,:),:);
        end
        
        
        %Fonctions pour manipuler les diam�tres de canaux (edge data)
        
        function dataStruct=GetEdgeDataList(network)
            dataStruct=network.EdgeDataList.EdgeDatas;
        end       
        
        function NewEdgeData(network,data,name)
            network.EdgeDataList.AddData(name,data);
        end
        
        %Fonctions pour manipuler les diam�tres des pores (vertice data)
        function dataStruct=GetVerticeDataList(network)
            dataStruct=network.VerticeDataList.VerticeDatas;
        end       
        
        function NewVerticeData(network,data,name)
            network.VerticeDataList.AddData(name,data);
        end
        
        %Fonctions output pour visualisation
        
        function outputStruct=InternalOutputStruct(network)
            %cr�ation en m�moire de la structure output qui pourra �tre 
            %visualis�e par le visualisateur interne

            outputStruct=struct;
            
            outputStruct.FaceOwners=network.FaceOwners;
            outputStruct.FaceNeighbours=network.FaceNeighbours;               
            outputStruct.Boundaries=network.Boundaries;            
            outputStruct.Cells.Cell=network.Cells;    
            
            attribute.Dimension=network.Dimension;
            attribute.NombreFaces=network.NombreFaces;
            attribute.NombreCells=network.NombreCells;
            attribute.NombreVertices=network.NombreVertices;
            attribute.NombreEdges=network.NombreEdges;
            attribute.Type='DelaunayNetwork';
            outputStruct.ATTRIBUTE=attribute;
            
            outputStruct.Vertices=network.Vertices;
            outputStruct.Faces.Face=network.Faces;
            outputStruct.CellsToVertices=network.CellsToVertices;
            outputStruct.Edges=network.Edges;
            
            names=fieldnames(network.EdgeDataList.EdgeDatas);
            for i=1:length(names)
                outputStruct.EdgeData.(names{i})=network.EdgeDataList.EdgeDatas.(names{i});
            end
            
            names=fieldnames(network.VerticeDataList.VerticeDatas);
            for i=1:length(names)
                outputStruct.VerticeData.(names{i})=network.VerticeDataList.VerticeDatas.(names{i});
            end
        end
        
        function vtk_struct=VTKOutputStruct(pore_network)
            %cr�ation en m�moire de la structure d'un fichier VTK POLYDATA 
            %pour afficher le r�seau de pores dans paraview.
            vertices=pore_network.Vertices;
            edges=pore_network.Edges;
            faces=pore_network.Faces;
            cells=pore_network.Cells;
            pore_data=struct;
            link_data=struct;
            edge_data=pore_network.EdgeDataList.EdgeDatas;
            vertice_data=pore_network.VerticeDataList.VerticeDatas;
            dimension=pore_network.Dimension;

            %Partie POINTS du fichier
            nCell=length(cells);
            nombre_points_pour_cells=0;
            for iCell=1:nCell
                liste_faces=cells{iCell};
                for iFace=liste_faces
                   nombre_points_pour_cells=nombre_points_pour_cells+length(faces{iFace}); 
                end
            end
            nEdge=length(edges(:,1));
            nVertice=length(vertices(:,1));
            nPoints=2*nEdge+sum(cellfun('length',faces))+nombre_points_pour_cells+nVertice;
            points_output=zeros(nPoints,3);
                        
            edge_to_points=cell(1,nEdge);
            x_extension=max(vertices(:,1))-min(vertices(:,1));
            y_extension=max(vertices(:,2))-min(vertices(:,2));
            for iEdge=1:nEdge
                if dimension==3
                    points_output((2*iEdge-1),:)=vertices(edges(iEdge,1),:);
                    points_output((2*iEdge),:)=vertices(edges(iEdge,2),:);
                else
                    points_output((2*iEdge-1),:)=[vertices(edges(iEdge,1),:),0];
                    points_output((2*iEdge),:)=[vertices(edges(iEdge,2),:),0];
                end
                edge_to_points{iEdge}=[(2*iEdge-1),(2*iEdge)];
            end
          
            compteur_debut=2*nEdge+1;

            nFace=length(faces);
            face_to_points=cell(1,nFace);
            for iFace=1:nFace
                compteur_fin=compteur_debut+length(faces{iFace})-1;
                if dimension==3
                    points_output((compteur_debut:compteur_fin),:)=vertices(faces{iFace},:);
                else
                    points_output((compteur_debut:compteur_fin),:)=horzcat(vertices(faces{iFace},:),zeros(length(faces{iFace}),1));
                end
                face_to_points{iFace}=(compteur_debut):(compteur_fin);
                compteur_debut=compteur_fin+1;
            end
            
            cell_to_points=cell(1,nCell);
            for iCell=1:nCell
                liste_faces=cells{iCell};
                cell_to_points{iCell}=cell(1,length(liste_faces));
                centre=mean(vertices(pore_network.CellsToVertices{iCell},:));
                i=0;
                for iFace=liste_faces
                    i=i+1;
                    compteur_fin=compteur_debut+length(faces{iFace})-1;
                    %new_vertices=0.999*(vertices(faces{iFace},:)-ones(length(vertices(faces{iFace})),1)*centre)+ones(length(vertices(faces{iFace})),1)*centre;
                    if dimension==3
                        %points_output((compteur_debut:compteur_fin),:)=new_vertices;
                        points_output((compteur_debut:compteur_fin),:)=vertices(faces{iFace},:);
                    else
                        %points_output((compteur_debut:compteur_fin),:)=horzcat(new_vertices,zeros(length(faces{iFace}),1));
                        points_output((compteur_debut:compteur_fin),:)=horzcat(vertices(faces{iFace},:),zeros(length(faces{iFace}),1));
                    end                   
                    cell_to_points{iCell}{i}=(compteur_debut):(compteur_fin);
                    compteur_debut=compteur_fin+1;
                end
            end
            
            compteur_fin=compteur_debut+nVertice-1;
            if dimension==3
                points_output((compteur_debut:compteur_fin),:)=vertices;
            else
                points_output((compteur_debut:compteur_fin),:)=horzcat(vertices,zeros(nVertice,1));
            end
            vertices_to_points=cell(1,nVertice);
            for i=1:nVertice
                vertices_to_points{i}=compteur_debut+i-1;
            end
            compteur_debut=compteur_fin+1;
            
            assert(compteur_debut==nPoints+1);
            
            %Partie LINES du fichier
            nLine=nEdge;
            lines_output=zeros(nLine,3);
            for iLine=1:nLine
                lines_output(iLine,1)=2;
                lines_output(iLine,2)=edge_to_points{iLine}(1)-1;
                lines_output(iLine,3)=edge_to_points{iLine}(2)-1;
            end
            
            %Partie POLYGONS du fichier
            nCellPolygone=sum(cellfun('length',cells));
            nPolygone=nFace+nCellPolygone;
            polygons_output=cell(nPolygone,1);
            face_to_polygon=cell(1,nFace);
            for iFace=1:nFace
                foo=face_to_points{iFace}-1;
                polygons_output{iFace}=[length(foo),foo];
                face_to_polygon{iFace}=iFace;
            end
            compteur=nFace+1;
            cell_to_polygon=cell(1,nCell);
            for iCell=1:nCell
                cell_to_polygon{iCell}=compteur:(compteur+length(cells{iCell})-1);
                for i=1:length(cells{iCell})
                    foo=cell_to_points{iCell}{i}-1;
                    polygons_output{compteur}=[length(foo),foo];
                    compteur=compteur+1;
                end
            end
            assert((compteur-1)==length(polygons_output));
            
            %Partie VERTICES du fichier
            vertices_output=horzcat(ones(nVertice,1),transpose((nPoints-nVertice+1):nPoints));
            
            
            %Partie POINT_DATA du fichier
            point_data_output=struct;
                %transcription des edge data
            names=fieldnames(edge_data);
            for i=1:length(names)
                data_output=zeros(nPoints,1);
                data=edge_data.(names{i});
                for iEdge=1:nEdge
                    data_output(edge_to_points{iEdge}(1))=data(iEdge);
                    data_output(edge_to_points{iEdge}(2))=data(iEdge);
                end
                point_data_output.(strcat('Edge_',names{i}))=data_output;
            end
                %transcription des vertices data
            names=fieldnames(vertice_data);
            for i=1:length(names)
                data_output=zeros(nPoints,1);
                data=vertice_data.(names{i});
                for iVertice=1:nVertice
                    data_output(vertices_to_points{iVertice})=data(iVertice);
                end
                point_data_output.(strcat('Vertice_',names{i}))=data_output;
            end    
                
                
            
            %Partie CELL_DATA du fichier
            nCellData=nLine+nPolygone;
            cell_data_output=struct;
                %transcription des link data
            names=fieldnames(link_data);
            for i=1:length(names)
                data_output=zeros(nCellData,1);
                data=link_data.(names{i});
                for iFace=1:nFace
                    data_output(face_to_polygon{iFace}+nLine)=data(iFace);
                end
                cell_data_output.(strcat('Link_',names{i}))=data_output;
            end
                %transcription des pore data
            names=fieldnames(pore_data);
            for i=1:length(names)
                data_output=zeros(nCellData,1);
                data=pore_data.(names{i});
                for iCell=1:nCell
                    data_output(cell_to_polygon{iCell}+nLine)=data(iCell)*ones(length(cell_to_polygon{iCell}),1);
                end
                cell_data_output.(strcat('Pore_',names{i}))=data_output;
            end

            vtk_struct.Points=points_output;
            vtk_struct.Lines=lines_output;
            vtk_struct.Polygons=polygons_output;
            vtk_struct.Vertices=vertices_output;
            vtk_struct.PointData=point_data_output;
            vtk_struct.CellData=cell_data_output;
        end       
        
    end
       
end

