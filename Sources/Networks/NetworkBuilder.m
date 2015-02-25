classdef NetworkBuilder 
    %NetworkBuilder Outil de construction des reseaux
    %Generation de reseau de pore a partir d'une geometrie macroscopique.
    %Contient notamment de nombreuses fonctions pour construire un reseau
    %de pores a partir d'un maillage aleatoire
    
    properties (SetAccess = immutable)
        MacroscopicGeometry
    end
    
    methods
        %Constructeur
        function networkBuilder = NetworkBuilder(myGeometry)
            %constructeur a partir de output_struct provenant d'un reader
            %initialisation des proprietes
            myGeometry.ConvertScale;
            %output_struct.MacroscopicGeometry.Vertices = output_struct.MacroscopicGeometry.ATTRIBUTE.ConvertToMeters*output_struct.MacroscopicGeometry.Vertices;
            networkBuilder.MacroscopicGeometry = myGeometry;
        end
        
        %BuildNetwork 
        function network = BuildNetwork(networkBuilder)
            %construction d'un network a partir de l'input geometrie_macroscopique
            %qui doit renseigner toutes les properties de la geometrie consideree
            assert(isa(networkBuilder,'NetworkBuilder'),'BuildNetwork prend un NetworkBuilder en argument');
            
            disp('Generation du reseau...');
            tic;
            
            myGeometry = networkBuilder.MacroscopicGeometry;
            
            %Initialise the random generator
            rng(myGeometry.GetRandomSeed);
            
            type_reseau = myGeometry.GetNetworkType;
            switch type_reseau
                case 'PoreNetworkMesh'
                    [dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,~,~,~,~] = NetworkBuilder.GenerateMesh(myGeometry);
                    
                    network = PoreNetworkMesh(dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,myGeometry);
                
                case 'PoreNetworkMeshFibrous'
                    [dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,faces_to_edges,edges_to_faces] = NetworkBuilder.GenerateMesh(myGeometry);
                    epaisseur_edges = NetworkBuilder.GenerateEdgeThickness(edges,vertices,myGeometry);    
                    
                    network = PoreNetworkMeshFibrous(dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,epaisseur_edges,faces_to_edges,edges_to_faces,myGeometry);
                    
                    diameter = network.ComputeAllLinkDiameter;
                    network.AddNewLinkData(diameter,'Diameter');
                    surface =  network.ComputeAllLinkSurface;
                    network.AddNewLinkData(surface,'Surface');
                    
                case 'DelaunayNetwork'
                    [dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,faces_to_edges,edges_to_faces] = NetworkBuilder.GenerateMesh(myGeometry);
                    epaisseur_edges = NetworkBuilder.GenerateEdgeThickness(edges,vertices,myGeometry); 
                    
                    network = DelaunayNetwork(dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,epaisseur_edges,faces_to_edges,edges_to_faces,myGeometry);
            
                case 'StructuredNetwork'
                    [dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,poreVolume,linkDiameter]=NetworkBuilder.GenerateStructuredNetwork(myGeometry);
                    
                    network=PoreNetworkEuclidien(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry.CopyGeometry);
                    network.AddNewLinkData(linkDiameter,'Diameter');
                    network.AddNewPoreData(poreVolume,'Volume');
            end
            
            duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
            fprintf('Reseau genere. Duree : %d minutes %f s. \n',minutes,secondes);
        end %BuildNetwork 
        
    end   
    
    methods (Static = true)      
        
        %GenerateMesh
        function [dimension,faces,cells,cells_to_vertices,owners,neighbours,boundaries,vertices,edges,vertices_to_edges,faces_to_edges,edges_to_faces] = GenerateMesh(myGeometry)
            %genere un maillage en fonction de la geometrie macroscopique.
            %input : geometrie macroscopique
            %output : 
            %   - vertices : tableau dont les lignes contiennent les coordonnees des vertices
            %   - faces : 1*NombreFaces cell array
            %   - cellules : 1*NombreCells cell array
            %   - owners, neighbours : meme chose que proprietes de
            %     geometrie, avec faces non renumerotees
            %   - faces_des_frontieres_exterieures : cell array taille
            %   NombreFrontieresExterieures ( = toutes les boundaries sauf les
            %   interfaces internes)
            %   faces_des_frontieres_exterieures{i} = {numero de la boundary associee, [numero des faces liees a la frontiere exterieure i]}   
            
            dimension = myGeometry.GetDimension;
            
            %Anisotropie : dilation d'une zone du maillage
            NetworkBuilder.GereAnisotropieDebut(myGeometry);
            
            %tirage aleatoire des points de base de voronoi dans les blocks
            nBlock = myGeometry.GetNumberOfBlocks;
            points_in_blocks = cell(1,nBlock); %points_in_blocks{i} = tableau contenant les coordonnees des points dans le bloc i
            
            for iBlock = 1:nBlock
                block_vertices = myGeometry.GetBlockVertices(iBlock);
                fillingParameters = myGeometry.GetBlockFillingParameters(iBlock) ;
                points_in_blocks{iBlock} = NetworkBuilder.GenerateRandomPoints(block_vertices,fillingParameters);                
            end
            
            %Remodelage des frontieres : le maillage doit etre modifie au
            %voisinage des interfaces entre blocks et des frontieres exterieures
            [points_in_blocks,indices_cellules_esclaves] = NetworkBuilder.RemodelerFrontieres(points_in_blocks,myGeometry);
            
            %Generation du pavage de voronoi a partir des points de
            %tous les blocks mis en commun. Cele permet d'avoir un maillage 
            %coherent aux interfaces entre blocks. 
            %Mise a jour des references vers les cellules esclaves 
            %des frontieres pour tenir compte de la nouvelle numerotation
            %due a la mise en commun des points/cellules.
            [sommets_voronoi,cellules_voronoi,indices_cellules_esclaves] = NetworkBuilder.ProcedurePavage(points_in_blocks,indices_cellules_esclaves,myGeometry);
            
            %Nettoyage du reseau : selection des cellules de voronoi et des
            %sommets adequats. Il faut faire le menage car on a utilise des
            %points en dehors du bloc pour former les frontieres exterieures.
            cellules_inutiles = NetworkBuilder.ListerCellulesInutiles(cellules_voronoi,sommets_voronoi,indices_cellules_esclaves,myGeometry);
            
            %Extraction des faces et de leurs pores voisins : creation de
            %face, owners, neighbour,faces_des_frontieres_exterieures
            [dirty_faces,dirty_owners,dirty_neighbours,faces_des_frontieres_exterieures] = NetworkBuilder.FaceExtraction(cellules_voronoi,cellules_inutiles,indices_cellules_esclaves,dimension,myGeometry,sommets_voronoi);
            
            %Nettoyage des listes cellules et vertices
            [cells_to_vertices,vertices,new_numeros_cellules,new_numeros_sommets,nombre_cellules] = NetworkBuilder.MiseAJourCellulesSommets(cellules_voronoi,sommets_voronoi,cellules_inutiles);                                                                                                          
                        
            %Nettoyage de la liste des faces, renumerotation des faces
            %aux frontieres et mise en ordre des sommets des faces
            [faces,owners,neighbours,boundaries] = NetworkBuilder.RenumerotationFaces(dirty_faces,dirty_owners,dirty_neighbours,faces_des_frontieres_exterieures,myGeometry,new_numeros_cellules,new_numeros_sommets,vertices);
            
            %Creation de la liste cells : cellules to faces
            cells = NetworkBuilder.BuildPoresToLinks(owners,neighbours,nombre_cellules);
            
            %Anisotropie : dilatation inverse
            vertices = NetworkBuilder.GereAnisotropieFin(myGeometry,vertices);
            
            %Extraction des edges
            [edges,vertices_to_edges,faces_to_edges,edges_to_faces] = NetworkBuilder.ExtractEdges(vertices,faces,myGeometry.GetNetworkType);
        end %GenerateMesh
        
        %GenerateStructuredNetwork
        function [dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,poreVolume,linkDiameter]=GenerateStructuredNetwork(myGeometry)
            
            dimension = myGeometry.GetDimension;
            
            assert(myGeometry.GetNumberOfBlocks==1,'Only one block supported for regular network generation. Try a VoronoiMesh network.');
            
            %Regularly spaced points generation
            fillingParameters = myGeometry.GetBlockFillingParameters(1);  
            domain_vertices = myGeometry.GetAllVertices;
            [points,Nxyz] = NetworkBuilder.GenerateRegularlySpacedPoints(domain_vertices,fillingParameters);
            
            %Pores definition
            nPore = length(points(:,1));
            poreCenter = points;
            
            %Link definition
            switch dimension
                case 2
                    nx=Nxyz(1);
                    ny=Nxyz(2);
                    nLink = ny*(nx+1)+nx*(ny+1);
                    owners = zeros(1,nLink);
                    neighbours = zeros(1,nLink);
                    linkCenter=zeros(nLink,2);
                    infos_liens_frontieres=cell(1,4);
                    
                    %Direction X
                    infos_liens_frontieres{1}=zeros(1,ny);
                    infos_liens_frontieres{2}=zeros(1,ny);
                    deltaX=norm(poreCenter(1)-poreCenter(2))/2;
                    for iy=1:ny
                        
                        ilx=1;
                        iLink = ilx+(nx+1)*(iy-1);
                        owners(iLink) = ilx+nx*(iy-1);
                        neighbours(iLink) = -1;
                        linkCenter(iLink,:) = poreCenter(owners(iLink),:)-[deltaX,0]; %TO CHECK
                        infos_liens_frontieres{1}(iy)=iLink;
                        
                        for ilx=2:nx
                            iLink = ilx+(nx+1)*(iy-1);
                            owners(iLink) = ilx+nx*(iy-1)-1;
                            neighbours(iLink) = owners(iLink)+1;
                            linkCenter(iLink,:) = ( poreCenter(owners(iLink),:) + poreCenter(neighbours(iLink),:) )/2;
                        end
                        
                        ilx=nx+1;
                        iLink = ilx+(nx+1)*(iy-1);
                        owners(iLink) = ilx+nx*(iy-1)-1;
                        neighbours(iLink) = -1;
                        linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[deltaX,0] ; %TO CHECK
                        infos_liens_frontieres{2}(iy)=iLink;
                    end
                    %Direction Y
                    shift=(nx+1)*ny;
                    infos_liens_frontieres{3}=zeros(1,nx);
                    infos_liens_frontieres{4}=zeros(1,nx);
                    deltaY=deltaX;
                    for ix=1:nx
                        ily=1;
                        iLink = shift+ily+(ny+1)*(ix-1);
                        owners(iLink) = ily+ny*(ix-1);
                        neighbours(iLink) = -1;
                        linkCenter(iLink,:) = poreCenter(owners(iLink),:)-[0,deltaY]; %TO CHECK
                        infos_liens_frontieres{3}(ix)=iLink;
                        
                        for ily=2:ny
                            iLink = shift+ily+(ny+1)*(ix-1);
                            owners(iLink) = ix+nx*(ily-2);
                            neighbours(iLink) = owners(iLink)+nx;
                            linkCenter(iLink,:) = ( poreCenter(owners(iLink),:) + poreCenter(neighbours(iLink),:) )/2;
                        end
                        
                        ily=ny+1;
                        iLink = shift+ ily+(ny+1)*(ix-1);
                        owners(iLink) = ix+nx*(ily-2);
                        neighbours(iLink) = -1;
                        linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[0,deltaY] ; %TO CHECK
                        infos_liens_frontieres{4}(ix)=iLink;
                    end
                    
                case 3
                    nx=Nxyz(1);
                    ny=Nxyz(2);
                    nz=Nxyz(3);
                    nLink = (nx+1)*ny*nz+nx*(ny+1)*nz+nx*ny*(nz+1);
                    owners = zeros(1,nLink);
                    neighbours = zeros(1,nLink);
                    linkCenter=zeros(nLink,3);
                    infos_liens_frontieres=cell(1,6);
                    
                    %Direction X
                    infos_liens_frontieres{1}=zeros(1,ny*nz);
                    infos_liens_frontieres{2}=zeros(1,ny*nz);
                    deltaX=norm(poreCenter(1)-poreCenter(2))/2;
                    for iz=1:nz
                        for iy=1:ny

                            ilx=1;
                            iLink = ilx+(nx+1)*(iy-1)+(nx+1)*ny*(iz-1);
                            owners(iLink) = ilx+nx*(iy-1)+nx*ny*(iz-1);
                            neighbours(iLink) = -1;
                            linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[-deltaX,0,0]; %TO CHECK
                            infos_liens_frontieres{1}(iy+ny*(iz-1))=iLink;

                            for ilx=2:nx
                                iLink = ilx+(nx+1)*(iy-1)+(nx+1)*ny*(iz-1);
                                owners(iLink) = ilx+nx*(iy-1)-1+nx*ny*(iz-1);
                                neighbours(iLink) = owners(iLink)+1;
                                linkCenter(iLink,:) = ( poreCenter(owners(iLink),:) + poreCenter(neighbours(iLink),:) )/2;
                            end

                            ilx=nx+1;
                            iLink = ilx+(nx+1)*(iy-1)+(nx+1)*ny*(iz-1);
                            owners(iLink) = ilx+nx*(iy-1)-1+nx*ny*(iz-1);
                            neighbours(iLink) = -1;
                            linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[deltaX,0,0] ; %TO CHECK
                            infos_liens_frontieres{2}(iy+ny*(iz-1))=iLink;
                        end
                    end
                    %Direction Y
                    shift=(nx+1)*ny*nz;
                    infos_liens_frontieres{3}=zeros(1,nx*nz);
                    infos_liens_frontieres{4}=zeros(1,nx*nz);
                    deltaY=deltaX;
                    for iz=1:nz
                        for ix=1:nx
                            ily=1;
                            iLink = shift+ily+(ny+1)*(ix-1)+nx*(ny+1)*(iz-1);
                            owners(iLink) = ily+ny*(ix-1)+nx*ny*(iz-1);
                            neighbours(iLink) = -1;
                            linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[0,-deltaY,0]; %TO CHECK
                            infos_liens_frontieres{3}(ix+nx*(iz-1))=iLink;

                            for ily=2:ny
                                iLink = shift+ily+(ny+1)*(ix-1)+nx*(ny+1)*(iz-1);
                                owners(iLink) = ix+nx*(ily-2)+nx*ny*(iz-1);
                                neighbours(iLink) = owners(iLink)+nx;
                                linkCenter(iLink,:) = ( poreCenter(owners(iLink),:) + poreCenter(neighbours(iLink),:) )/2;
                            end

                            ily=ny+1;
                            iLink = shift+ ily+(ny+1)*(ix-1)+nx*(ny+1)*(iz-1);
                            owners(iLink) = ix+nx*(ily-2)+nx*ny*(iz-1);
                            neighbours(iLink) = -1;
                            linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[0,deltaY,0] ; %TO CHECK
                            infos_liens_frontieres{4}(ix+nx*(iz-1))=iLink;
                        end
                    end
                    
                    %Direction Z
                    shift=(nx+1)*ny*nz+nx*(ny+1)*nz;
                    infos_liens_frontieres{5}=zeros(1,nx*ny);
                    infos_liens_frontieres{6}=zeros(1,nx*ny);
                    deltaZ=deltaX;
                    for ix=1:nx
                        for iy=1:ny
                            ilz=1;
                            iLink = shift+ilz+(nz+1)*(ix-1)+(iy-1)*(nz+1)*nx;
                            owners(iLink) = ix+nx*(iy-1);
                            neighbours(iLink) = -1;
                            linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[0,0,-deltaZ]; %TO CHECK
                            infos_liens_frontieres{5}(iy+ny*(ix-1))=iLink;

                            for ilz=2:nz
                                iLink = shift+ilz+(nz+1)*(ix-1)+(iy-1)*(nz+1)*nx;
                                owners(iLink) = ix+nx*(iy-1)+nx*ny*(ilz-2);
                                neighbours(iLink) = owners(iLink)+nx*ny;
                                linkCenter(iLink,:) = ( poreCenter(owners(iLink),:) + poreCenter(neighbours(iLink),:) )/2;
                            end

                            ilz=nz+1;
                            iLink = shift+ ilz+(nz+1)*(ix-1)+(iy-1)*(nz+1)*nx;
                            owners(iLink) = ix+nx*(iy-1)+nx*ny*(ilz-2);
                            neighbours(iLink) = -1;
                            linkCenter(iLink,:) = poreCenter(owners(iLink),:)+[0,0,deltaZ] ; %TO CHECK
                            infos_liens_frontieres{6}(iy+ny*(ix-1))=iLink;
                        end
                    end
                    
                    
                    
                    
            end
            
            %infos_liens_frontieres{indice_frontiere}=listeLiens;
            [boundaries,owners,neighbours] = NetworkBuilder.RenumerotationLiensFrontieres(infos_liens_frontieres,owners,neighbours);
            
            %Pores to links look up table generation
            pores = NetworkBuilder.BuildPoresToLinks(owners,neighbours,nPore);
            
                        

            
            %Volumes pores et diametres liens
            scaleFactor=0.3*norm(poreCenter(1,:) - poreCenter(2));
            poreVolume=rand(1,nPore)*scaleFactor^3;               %TO DO
            linkDiameter=rand(1,nLink)*scaleFactor^2;             %TO DO
            
            
            
            
            
            
        end %GenerateStructuredNetwork       
        
        %ProcedurePavage   
        function [sommets_pavage,cellules_pavage,indices_cellules_esclaves] = ProcedurePavage(points_in_blocks,indices_cellules_esclaves,myGeometry)
            
            %Creation liste commune pour tous les points
            dimension = myGeometry.GetDimension;
            nombre_points_total = sum(cellfun('length',points_in_blocks));
            points = zeros(nombre_points_total,dimension);
            indice_debut_block_courant = 1;
            nBlock = myGeometry.GetNumberOfBlocks;
            indices_debut_blocs = zeros(1,nBlock);
            for iBlock = 1:nBlock
                indices_debut_blocs(iBlock) = indice_debut_block_courant;
                points(indice_debut_block_courant:indice_debut_block_courant+length(points_in_blocks{iBlock})-1,:) = points_in_blocks{iBlock}; 
                indice_debut_block_courant = indice_debut_block_courant+length(points_in_blocks{iBlock});
            end
            assert(indice_debut_block_courant == nombre_points_total+1,'Pb lors de la mise en commun des points generateurs du pavage');
            
            %Mise a jour des references vers les cellules esclaves pour
            %prendre en compte les changements de numeros dus a
            %la mise en commun des points de tous les blocs
            for iBoundary = 1:myGeometry.GetNumberOfBoundaries
                indice_block = myGeometry.FindBlockOfThisBoundary(iBoundary);
                indices_cellules_esclaves{iBoundary} = indices_cellules_esclaves{iBoundary}+indices_debut_blocs(indice_block)-1;
            end
            
            %Generation maillage
            type_pavage = myGeometry.GetPavageType;
            switch type_pavage
                case 'RandomVoronoi'
                    [sommets_pavage,cellules_pavage] = voronoin(points);
                case 'RandomDelaunay'
                    sommets_pavage = points;
                    cellules = delaunay(points(:,1),points(:,2));
                    cellules_pavage = num2cell(cellules,2);
                otherwise
                    error('Type de pavage non reconnu');
            end
            
            
        end %ProcedurePavage
        
        %GenerateRegularlySpacedPoints
        function [points,Nxyz ] = GenerateRegularlySpacedPoints(domain_vertices,fillingParameters)
            
            dimension = length(domain_vertices(1,:));
            nPoints = fillingParameters.PoreNumber;
            
            tol=0.05;
            switch dimension
                case 2
                    xmin = min(domain_vertices(:,1));
                    xmax = max(domain_vertices(:,1));
                    ymin = min(domain_vertices(:,2));
                    ymax = max(domain_vertices(:,2));
                    nombre_points_in = 0;
                    latticeSpacing=((xmax-xmin)*(ymax-ymin)/nPoints)^(1/2);
                    for iteration = 1:5
                        if abs(nombre_points_in-nPoints)>tol*nPoints
                            X = linspace(xmin,xmax,ceil((xmax-xmin)/latticeSpacing));
                            Y = linspace(ymin,ymax,ceil((ymax-ymin)/latticeSpacing));
                            
                            nx = length(X);
                            ny = length(Y);
                            poi = zeros(nx*ny,2);
                            for ix=1:nx
                               for iy=1:ny
                                   poi((iy-1)*nx+ix,:) = [X(ix),Y(iy)];
                               end                              
                            end
                            
                            points = poi(transpose(find(inhull(poi,domain_vertices))),:);
                            nombre_points_in = length(points);
                            latticeSpacing = latticeSpacing*(nombre_points_in/nPoints)^(1/2);
                        end
                    end
                    Nxyz=[nx,ny];
                case 3
                    xmin = min(domain_vertices(:,1));
                    xmax = max(domain_vertices(:,1));
                    ymin = min(domain_vertices(:,2));
                    ymax = max(domain_vertices(:,2));
                    zmin = min(domain_vertices(:,3));
                    zmax = max(domain_vertices(:,3));
                    nombre_points_in = 0;
                    latticeSpacing = ((xmax-xmin)*(ymax-ymin)*(zmax-zmin)/nPoints)^(1/3);

                    for iteration = 1:5
                        if abs(nombre_points_in-nPoints)>tol*nPoints
                            X = linspace(xmin,xmax,ceil((xmax-xmin)/latticeSpacing));
                            Y = linspace(ymin,ymax,ceil((ymax-ymin)/latticeSpacing));
                            Z = linspace(zmin,zmax,ceil((zmax-zmin)/latticeSpacing));
                            
                            nx = length(X);
                            ny = length(Y);
                            nz = length(Z);
                            poi = zeros(nx*ny*nz,3);
                            for ix=1:nx
                               for iy=1:ny
                                   for iz=1:nz
                                    poi(ix+nx*(iy-1)+nx*ny*(iz-1),:) = [X(ix),Y(iy),Z(iz)];
                                   end
                               end                              
                            end
                            
                            points = poi(transpose(find(inhull(poi,domain_vertices))),:);
                            nombre_points_in = length(points);
                            latticeSpacing = latticeSpacing*(nombre_points_in/nPoints)^(1/3);
                        end
                    end
                    Nxyz=[nx,ny,nz];
            end         
            
            
        end %GenerateRegularlySpacedPoints
        
        %GenerateRandomPoints
        function points = GenerateRandomPoints(domain_vertices,fillingParameters)
            %Tire des points aleatoirement dans un domaine polyedrique.
            %input:   - domain_vertices : coordonnees des sommets du
            %       domaine polyedrique dans lequel les points aleatoire
            %       sont tires
            %         - parametres_tirage_aleatoire : parametres de la loi
            %        de probabilite : geometrie_macroscopique.Blocks.Block(indice_block).Remplissage;
            %output : - points : tableau contenant les coordonnees des points [x y (z) ; x y (z) ; ....]  
            
            dimension = length(domain_vertices(1,:));
            nombre_points = fillingParameters.PoreNumber;
            
            %cas oe le bloc est pave aligne avec les axes :
            
            switch dimension
                case 2
                    xmin = min(domain_vertices(:,1));
                    xmax = max(domain_vertices(:,1));
                    ymin = min(domain_vertices(:,2));
                    ymax = max(domain_vertices(:,2));
                    nombre_points_in = 0;
                    for iteration = 1:5
                        if nombre_points_in<nombre_points
                            X = xmin+(xmax-xmin)*rand(2^(iteration-1)*nombre_points,1);
                            Y = ymin+(ymax-ymin)*rand(2^(iteration-1)*nombre_points,1);
                            poi = [X,Y];
                            points = poi(transpose(find(inhull(poi,domain_vertices))),:);
                            nombre_points_in = length(points);
                        end
                        if nombre_points_in>nombre_points
                            points = points(1:nombre_points,:);
                            nombre_points_in = length(points);
                        end
                    end
                case 3
                    xmin = min(domain_vertices(:,1));
                    xmax = max(domain_vertices(:,1));
                    ymin = min(domain_vertices(:,2));
                    ymax = max(domain_vertices(:,2));
                    zmin = min(domain_vertices(:,3));
                    zmax = max(domain_vertices(:,3));
                    nombre_points_in = 0;
                    for iteration = 1:5
                        if nombre_points_in<nombre_points
                            X = xmin+(xmax-xmin)*rand(2^(iteration-1)*nombre_points,1);
                            Y = ymin+(ymax-ymin)*rand(2^(iteration-1)*nombre_points,1);
                            Z = zmin+(zmax-zmin)*rand(2^(iteration-1)*nombre_points,1);
                            poi = [X,Y,Z];
                            points = poi(transpose(find(inhull(poi,domain_vertices))),:);
                            nombre_points_in = length(points);
                        end
                        if nombre_points_in>nombre_points
                            points = points(1:nombre_points,:);
                            nombre_points_in = length(points);
                        end
                    end
            end            
            
            if isfield(fillingParameters,'Density')
                switch fillingParameters.Density
                    case 'Uniforme'
                        points = points;
                        
                    case 'Puissance'
                        direction = fillingParameters.Direction;
                        factor = fillingParameters.Factor;
                        switch direction
                            case 'x'
                                indice = 1;
                            case 'y'
                                indice = 2;
                            case 'z'
                                indice = 3;
                        end
                        Vmin = min(domain_vertices(:,indice));
                        Vmax = max(domain_vertices(:,indice));
                        points(:,indice) = Vmin+(1/(Vmax-Vmin)^(factor-1))*(points(:,indice)-Vmin).^factor;
                        
                end
                
            end
            
            
        end %GenerateRandomPoints
        
        %RemodelerFrontieres
        function [points_in_blocks,indices_cellules_esclaves] = RemodelerFrontieres(points_in_blocks,myGeometry)
            %Remodelage des frontieres : aux frontieres des blocks le
            %maillage change. On veut traiter :
            %   - interfaces 'raide' ou '3e couche' entre deux blocks
            %   - surfaces exterieures 'plat' ou 'rough'
            %   - conditions limites periodiques ('cyclic') entre deux
            %   faces opposees
            %Algorithme : 
            %   - modifier la repartition des points de base voronoi
            %   au voisinage de la frontiere en ajoutant ou enlevant des 
            %points judicieusement choisis. L'ajout de points a ceux tires
            %   dans le corps du bloc se fait  en respectant l'ordre de 
            %   gestion des frontiere suivant : interfaces internes puis 
            %   surfaces rough puis surfaces plan puis limites periodiques
            %   - on cree le maillage de voronoi a partir de tous les
            %   points definis precedemment (interieur blocks + frontieres)
            %   - on nettoye ensuite le maillage des cellules en trop
            %   introduites pres de chaque frontiere.
            
            nBoundary = myGeometry.GetNumberOfBoundaries ;
            indices_cellules_esclaves = cell(1,nBoundary); 
            %indices_cellules_esclaves{i} : cellules potentiellement en trop pres de la boundary i 
            
            %Remodelage frontieres : interfaces entre les blocks
            for iBoundary = 1:nBoundary
%                 if strcmp(geometrie_macroscopique.Boundaries.Boundary(indice_boundary).ATTRIBUTE.Type,'interface')           
%                   
%                     name_neighbour = geometrie_macroscopique.Boundaries.Boundary(indice_boundary).ATTRIBUTE.NeighbourPatch;
% 
%                     indice_block =myGeometry.FindBlockOfThisBoundary(iBoundary);
%                     random_points_in_block = points_in_blocks{indice_block};
%                     infos_block = geometrie_macroscopique.Blocks.Block(indice_block);
%                     infos_frontiere = geometrie_macroscopique.Boundaries.Boundary(indice_boundary);
%                     vertices_geo_macro = geometrie_macroscopique.Vertices;
%                     
%                     [nouveaux_points,ind_cells_esclaves] = NetworkBuilder.GererUneFrontiere(random_points_in_block,infos_block,infos_frontiere,vertices_geo_macro);
%                     points_in_blocks{indice_block} = nouveaux_points;
%                     indices_cellules_esclaves{indice_boundary} = ind_cells_esclaves;
%                 end
            end
            
            %Remodelage frontieres : frontieres exterieures rugueuses ou planes                                       
            for iBoundary = 1:nBoundary
               if strcmp(myGeometry.GetBoundaryType(iBoundary) ,'surface')    
                    iBlock = myGeometry.FindBlockOfThisBoundary(iBoundary);
                    random_points_in_block = points_in_blocks{iBlock};                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %infos_block = myGeometry.Blocks.Block(iBlock);                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
                    %infos_frontiere = myGeometry.Boundaries.Boundary(iBoundary);
                    
                    %[newPoints,ind_cells_esclaves] = NetworkBuilder.GererUneFrontiere(random_points_in_block,infos_block,infos_frontiere,myGeometry.GetAllVertices);
                    [newPoints,ind_cells_esclaves] = NetworkBuilder.GererUneFrontiere(random_points_in_block,myGeometry,iBoundary);
                    points_in_blocks{iBlock} = newPoints;
                    indices_cellules_esclaves{iBoundary} = ind_cells_esclaves; %les numeros des cellules changeront apres mise en commun voronoi
               end
            end
            
            %Remodelage frontieres : limites periodiques
            liste_patch_periodiques_trouves = [];
            for iBoundary = 1:nBoundary
               if strcmp(myGeometry.GetBoundaryType(iBoundary),'cyclic')
                   %on traite la condition limite periodique lorqu'on a
                   %trouve la deuxieme frontiere de la paire
                    if false %la premiere frontiere est deja trouvee                                                           %TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        iBlock = myGeometry.FindBlockOfThisBoundary(iBoundary);
                        random_points_in_block = points_in_blocks{iBlock};
                        %infos_block = geometrie_macroscopique.Blocks.Block(iBlock);
                        %infos_frontiere = geometrie_macroscopique.Boundaries.Boundary(indice_boundary);
                        
                        %[newPoints,ind_cells_esclaves] = NetworkBuilder.GererUneFrontiere(random_points_in_block,infos_block,infos_frontiere,myGeometry.GetAllVertices);
                        [newPoints,ind_cells_esclaves] = NetworkBuilder.GererUneFrontiere(random_points_in_block,myGeometry,iBoundary);
                        points_in_blocks{iBlock} = newPoints;
                        indices_cellules_esclaves{indice_boundary} = ind_cells_esclaves;      
                        
                        %supprimer le nom de la premiere frontiere de 
                        %liste_patch_periodiques_trouves                                                                %TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    else
                        %ajouter le nom de cette frontiere a la liste                                                    %TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                    end
               end
            end
            %verifier que liste_patch_periodiques_trouves est bien vide                                                 %TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            
        end %RemodelerFrontieres
        
        %GererUneFrontiere
        function [newPointsInBlock,indices_cellules_esclaves] = GererUneFrontiere(random_points_in_block,myGeometry,iBoundary)
            %Modifie la repartition des points de base voronoi au voisinage
            %d'une frontiere et indique les cellules esclaves de cette 
            %frontiere, ie celles que cette frontiere pourra supprimer 
            %lors de la phase de nettoyage du maillage
            %input :  - infos_block = geometrie_macroscopique.Blocks.Block(indice_block)
            %           sauf pour conditions cyclic ou interface, auquel cas c'est une
            %           2 cell contenant les infos_block des deux blocs concernes
            %         - infos_frontiere = geometrie_macroscopique.Boundaries.Boundary(indice_boundary)
            %           sauf si conditions cyclic, auquel cas c'est une
            %           2 cell contenant les infos_frontiere des deux blocs concernes
            %
            %output : - newPointsInBlock  : liste des points du
            %           bloc mise a jour apres ajout et retrait de points
            %           pres de la frontiere
            %         - indices_cellules_esclaves
            
            %homogeneisation des input car ils ont une forme differente
            %selon le type de boundary
            
            
            

            
%             if iscell(infos_frontiere)
%                 type = 'cyclic';
%             else
%                 type = infos_frontiere.ATTRIBUTE.Type;
%             end
%             switch type
%                 case 'cyclic'
%                     nbre_frontieres = 2;
%                     inf_frontiere = infos_frontiere;
%                     if iscell(infos_block)
%                         assert(length(infos_block) == 2,'cyclic entre 2 blocks maximum')
%                         nbre_blocks = 2;
%                         initial_points = random_points_in_block;
%                         inf_block = infos_block;
% %                         test_coherence_frontiere_block_1 = isequal(ismember(inf_frontiere{1}.face,inf_block{1}.VertexNumbers),ones(1,length(inf_frontiere{1}.face)));
% %                         test_coherence_frontiere_block_2 = isequal(ismember(inf_frontiere{2}.face,inf_block{2}.VertexNumbers),ones(1,length(inf_frontiere{2}.face)));
% %                         assert(test_coherence_frontiere_block_1&&test_coherence_frontiere_block_2,'Pb infos frontieres bloc cyclic');
%                     else
%                         nbre_blocks = 1;
%                         initial_points{1} = random_points_in_block;         
%                         inf_block{1} = infos_block;
%                     end
%                 case 'interface'
%                     nbre_frontieres = 1;
%                     inf_frontiere{1} = infos_frontiere;
%                     nbre_blocks = 2;
%                     initial_points{1} = random_points_in_block;
%                     initial_points{2} = random_points_in_block;
%                     inf_block = infos_block;
%                 case 'surface'
%                     nbre_frontieres = 1;
%                     inf_frontiere{1} = infos_frontiere;
%                     nbre_blocks = 1;
%                     initial_points{1} = random_points_in_block;
%                     inf_block{1} = infos_block;
%             end
            
            

            %verification que la surface est plane
%             for i = 1:nbre_frontieres
%                 indices_vertex = inf_frontiere{i}.Face;
%                 assert(isequal(unique(indices_vertex),sort(indices_vertex)),'Points face identiques');
%                 vertices_frontiere = vertices_geo_macro(indices_vertex,:);
%                 switch dimension
%                     case 2
%                         assert(length(indices_vertex) == 2,'Frontiere non plane');
%                     case 3
%                         %assert(length(indices_vertex) == 4&&det(vertices_frontiere([1,2,3],:))<10^(-10)*epsilon^3&&det(vertices_frontiere([1,2,4],:))<10^(-10)*epsilon^3,'Frontiere non plane');
%                         
%                                                                     %A FAIRE
%                 end
%             end
            type=myGeometry.GetBoundaryType(iBoundary);
            macroVertices = myGeometry.GetAllVertices;
            dimension = myGeometry.GetDimension;
            epsilon = max(macroVertices(:,1))-min(macroVertices(:,1));
            
            iBlock = myGeometry.FindBlockOfThisBoundary(iBoundary);   

            nbre_frontieres = 1;
             nbre_blocks = 1;
             initial_points{1} = random_points_in_block;

            %calcul des distance des point voronoi a la frontiere pour
            %trouver les points proches de la frontiere
            for num_frontiere = 1:nbre_frontieres
                distances_a_la_frontiere = cell(1,nbre_frontieres);
                for num_block = 1:nbre_blocks
                    distances_a_la_frontiere{num_block} = zeros(length(initial_points{num_block}(:,1)),1);
                    boundaryVertices = myGeometry.GetBoundaryVertices(iBoundary);
                    for num_point = 1:length(initial_points{num_block}(:,1))
                        distances_a_la_frontiere{num_block}(num_point) = NetworkBuilder.DistancePointFrontierePlane(initial_points{num_block}(num_point,:),boundaryVertices);
                    end
                end
            end        
            
             
            
                    
            %modification des points voronoi pres de la frontiere et creation de indices_cellules_esclaves        
            switch type
                
                case 'interface'
                    x_pourcent = 0.3 ;
                    
                    
                    
                    
                    newPointsInBlock = random_points_in_block;
                
                %surface exterieure rugueuse ou plane
                case 'surface'
                    rugosity = myGeometry.GetBoundaryRugosity(iBoundary);
                    [dist_sorted,indices_sorted] = sort(distances_a_la_frontiere{1});
                    nbre_points_initiaux = length(initial_points{1}(:,1));
                    switch rugosity
                        case 'rough'
                            %on inspecte les x% de points plus proches de la frontiere, 
                            x_pourcent = 0.3;
                            %et on tire de nouveaux points de l'autre cote
                            %de la frontiere sur une epaisseur egale a
                            %celle occupee par ces x% de points.
                            nbre_points_proches_in = ceil(x_pourcent*length(indices_sorted));
                            indices_points_proches_in = transpose(indices_sorted(1:nbre_points_proches_in));
                            %definition de la zone ou tirer les points
                            dist_max_in = dist_sorted(nbre_points_proches_in);
                            blockCenter = mean( myGeometry.GetBlockVertices(iBlock));
                            vect_perp = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(boundaryVertices,blockCenter);
                            vertices_decales = boundaryVertices+dist_max_in*ones(length(boundaryVertices(:,1)),1)*vect_perp;
                            zone_vertices = [boundaryVertices;vertices_decales];
                            
                            fooFillingParameter.PoreNumber = nbre_points_proches_in;
                            %tirer les points 
                            points_out_proches_surface = NetworkBuilder.GenerateRandomPoints(zone_vertices,fooFillingParameter);
                            %les ajouter aux points du bloc
                            newPointsInBlock = vertcat(initial_points{1},points_out_proches_surface);
                            %indices_cellules_esclaves : points_in et points_out
                            indices_points_out = (nbre_points_initiaux+1):(nbre_points_initiaux+nbre_points_proches_in);
                            indices_cellules_esclaves = horzcat(indices_points_proches_in,indices_points_out);
                            
                            
                        case 'flat'
                           %on inspecte les x% de points plus proches de la frontiere, 
                            x_pourcent = 0.4;
                            %et on mirore ces points par rapport a la
                            %frontiere
                            nbre_points_proches_in = ceil(x_pourcent*length(indices_sorted));
                            %indices_points_proches_in = transpose(indices_sorted(1:nbre_points_proches_in));
                            
                            blockCenter = mean( myGeometry.GetBlockVertices(iBlock));
                            vect_perp = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(boundaryVertices,blockCenter);
                           
                            points_out_proches_surface = zeros(nbre_points_proches_in,dimension);
                            for num_point = 1:nbre_points_proches_in
                                points_out_proches_surface(num_point,:) = initial_points{1}(indices_sorted(num_point),:)-2*dot(vect_perp,(initial_points{1}(indices_sorted(num_point),:)-boundaryVertices(1,:)))*vect_perp;
                            end
                            %les ajouter aux points du bloc
                            newPointsInBlock = vertcat(initial_points{1},points_out_proches_surface);
                            %indices_cellules_esclaves : points_out
                            indices_points_out = (nbre_points_initiaux+1):(nbre_points_initiaux+nbre_points_proches_in);
                            
                            nbre_points_proches_in = ceil(x_pourcent*length(indices_sorted));
                            indices_points_proches_in = transpose(indices_sorted(1:nbre_points_proches_in));
                            indices_cellules_esclaves = horzcat(indices_points_proches_in,indices_points_out);
                        otherwise
                            newPointsInBlock = random_points_in_block;
                            indices_cellules_esclaves = [];
                            
                    end
                %periodic boundary condition    
                case 'cyclic'
                    %verifier que les sommets des deux patchs sont
                    %geometriquement identiques 
                    type=myGeometry.GetBoundaryType(iBoundary);
                    %translater les points proches des deux frontieres 
                    %indices_cellules_esclaves : points_in proches des deux
                    %       frontieres et points translates
                                
                    newPointsInBlock = random_points_in_block;                                    %A FAIRE

                    indices_cellules_esclaves = [];
            

            
            end

        end %GererUneFrontiere
      
        
        
        
        %TrouverCellulesInutilesFrontiere
        function useless_cells = TrouverCellulesInutilesFrontiere(indices_cellules_esclaves,cellules_inutiles,sommets_voronoi,cellules_voronoi,myGeometry,iBoundary)
            %Trouve les cellules inutiles au maillage parmis les cellules 
            %esclaves d'une frontiere donnee
            %output: - useless_cells : tableau des indices des cellules inutiles
            %        pour cette frontiere.
            
            %homogeneisation des input car ils ont une forme differente
            %selon le type de boundary
            
            
            
%             if iscell(infos_frontiere)
%                 type = 'cyclic';
%             else
%                 type = infos_frontiere.ATTRIBUTE.Type;
%             end
%             switch type
%                 case 'cyclic'
%                     nbre_frontieres = 2;
%                     inf_frontiere = infos_frontiere;
%                     if iscell(infos_block)
%                         assert(length(infos_block) == 2,'cyclic entre 2 blocks maximum')
%                         nbre_blocks = 2;
%                         inf_block = infos_block;
%                         test_coherence_frontiere_block_1 = isequal(ismember(inf_frontiere{1}.face,inf_block{1}.VertexNumbers),ones(1,length(inf_frontiere{1}.face)));
%                         test_coherence_frontiere_block_2 = isequal(ismember(inf_frontiere{2}.face,inf_block{2}.VertexNumbers),ones(1,length(inf_frontiere{2}.face)));
%                         assert(test_coherence_frontiere_block_1&&test_coherence_frontiere_block_2,'Pb infos frontieres bloc cyclic');
%                     else
%                         nbre_blocks = 1;
%                         inf_block{1} = infos_block;
%                     end
%                 case 'interface'
%                     nbre_frontieres = 1;
%                     inf_frontiere{1} = infos_frontiere;
%                     nbre_blocks = 2;
%                     inf_block = infos_block;
%                 case 'surface'
%                     nbre_frontieres = 1;
%                     inf_frontiere{1} = infos_frontiere;
%                     nbre_blocks = 1;
%                     inf_block{1} = infos_block;
%             end

            type=myGeometry.GetBoundaryType(iBoundary);
            iBlock = myGeometry.FindBlockOfThisBoundary(iBoundary);   
            
            switch type
                case 'surface'
                    rugosity = myGeometry.GetBoundaryRugosity(iBoundary) ;
                    criteria  = 'barycenter';
                    switch rugosity
                    
                        case 'rough'
                            %enlever les cellules qui sont liees au sommet 
                            %a l'infini ou qui ne verifient pas le critere
                            useless_cells = [];
                            
                            for num_cell = indices_cellules_esclaves
                                if  cellules_inutiles(num_cell)
                                    useless_cells = [useless_cells,num_cell];
                                else
                                    
                                    switch criteria
%                                         case 'volume'
%                                             %au moins la moitie du volume
%                                             %dans le bloc
%                                             sommets_cell = sommets_voronoi(transpose(cellules_voronoi{num_cell}),:);
%                                             sommets_block = vertices_geo_macro(transpose(inf_block{1}.VertexNumbers),:);
%                                             proportion_volume_in = NetworkBuilder.ProportionVolumeIn(sommets_cell,sommets_block);
%                                             if proportion_volume_in<0.5
%                                                 useless_cells = [useless_cells,num_cell];
%                                             end
                                        case 'barycenter'
                                            %barycentre du bon côté de la
                                            %frontiere
                                            sommets_cell = sommets_voronoi(transpose(cellules_voronoi{num_cell}),:);
                                            barycenter = mean(sommets_cell);
                                            
                                            boundaryVertices = myGeometry.GetBoundaryVertices(iBoundary);
                                            blockCenter = mean( myGeometry.GetBlockVertices(iBlock) );
                                            vect_perp = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(boundaryVertices,blockCenter);
                                            
                                            if dot( (barycenter-boundaryVertices(1,:)),vect_perp)>0
                                                useless_cells = [useless_cells,num_cell];
                                            end  
                                    end
                                end
                            end
                            
                        case 'flat'
                            %enlever les cellules qui sont liees au sommet 
                            %a l'infini ou qui ne verifient pas le critere
                            useless_cells = [];
                            
                            for num_cell = indices_cellules_esclaves
                                if  cellules_inutiles(num_cell)
                                    useless_cells = [useless_cells,num_cell];
                                else
                                    switch criteria
%                                         case 'volume'
%                                             %au moins la moitie du volume
%                                             %dans le bloc
%                                             sommets_cell = sommets_voronoi(transpose(cellules_voronoi{num_cell}),:);
%                                             sommets_block = vertices_geo_macro(transpose(inf_block{1}.VertexNumbers),:);
%                                             proportion_volume_in = NetworkBuilder.ProportionVolumeIn(sommets_cell,sommets_block);
%                                             if proportion_volume_in<0.5
%                                                 useless_cells = [useless_cells,num_cell];
%                                             end
                                        case 'barycenter'
                                            %barycentre du bon côté de la
                                            %frontiere
                                            sommets_cell = sommets_voronoi(transpose(cellules_voronoi{num_cell}),:);
                                            barycenter = mean(sommets_cell);
                                            
                                            boundaryVertices = myGeometry.GetBoundaryVertices(iBoundary);
                                            blockCenter = mean( myGeometry.GetBlockVertices(iBlock) );
                                            vect_perp = NetworkBuilder.ComputeVecteurNormalOrienteExterieurBlock(boundaryVertices,blockCenter);
                                            
                                            if dot( (barycenter-boundaryVertices(1,:)),vect_perp)>0
                                                useless_cells = [useless_cells,num_cell];
                                            end  
                                            
                                    end
                                 end
                             end
                        otherwise    
                            useless_cells = [];
                    
                    end
                    
                case 'cyclic'    
                   useless_cells = [];
                otherwise
                   useless_cells = [] ;                                                       %A FAIRE

            end
            
        end %TrouverCellulesInutilesFrontiere       
        
        %ListerCellulesInutiles
        function cellules_inutiles = ListerCellulesInutiles(cellules_voronoi,sommets_voronoi,indices_cellules_esclaves,myGeometry)
            %Dresse la liste des cellules a supprimer lors du nettoyage.

            cellules_inutiles = zeros(1,length(cellules_voronoi));  %cellules_inutiles(i) = 1 si la cellule i est a enlever, 0 sinon

            %Voronoin genere un sommet a l'infini et parfois des sommets
            %tres eloignes des blocks. Il faudra les enlever.
            tol = 1.3;
            sommets_proches = zeros(1,length(sommets_voronoi));
            for iBlock = 1:myGeometry.GetNumberOfBlocks
                vertex_block = myGeometry.GetBlockVertices(iBlock);
                zone_proche_du_bloc = tol*vertex_block-(tol-1)*ones(length(vertex_block),1)*mean(vertex_block);
                sommets_proches(transpose(inhull(sommets_voronoi,zone_proche_du_bloc))) = 1;
            end
            sommets_infinis = 1-sommets_proches;
            
            foo = find(sommets_infinis);
            for num_cell = 1:length(cellules_voronoi)
                if ~isempty(find(ismember(cellules_voronoi{num_cell},foo),1))
                    cellules_inutiles(num_cell) = 1;
                end
            end 
            
            %Chaque boundary peut demander a supprimer des cellules parmi ses
            %cellules esclaves, par exemple si cette cellule servait
            %uniquement d'artifice pour remodeler la frontiere.
            for iBoundary = 1:myGeometry.GetNumberOfBoundaries 
                indices_useless_cells_frontiere = NetworkBuilder.TrouverCellulesInutilesFrontiere(indices_cellules_esclaves{iBoundary},cellules_inutiles,sommets_voronoi,cellules_voronoi,myGeometry,iBoundary);
                cellules_inutiles(indices_useless_cells_frontiere) = 1;
            end
                        
        end %ListerCellulesInutiles
        
        %FaceExtraction
        function [dirty_faces,dirty_owners,dirty_neighbours,faces_des_frontieres_exterieures] = FaceExtraction(cellules_voronoi,cellules_inutiles,indices_cellules_esclaves,dimension,myGeometry,sommets_voronoi)
            %Trouve les faces entre les cellules et renvoie les indices de
            %leurs sommets et leurs deux cellules voisines. 
            %Algorithme : parcourir les cellules en verifiant parmi leurs 
            %voisines si elles partagent assez de sommets pour avoir une 
            %face en commun.
            %Gere les faces exterieures en mettant -1 en neighbour et en 
            %les associant a une frontiere exterieure si elles
            %ont une cellule a l'interieur, sinon ne les compte pas. 
            
            
            %trouver les faces communes des cellules, ie plus de 2
            %sommets communs en 2D et 3 en 3D (cellules convexes)
            cellules_en_contact = cell(1,length(cellules_voronoi)); 
            %cellules_en_contact{num_cell}{num_face} = {num_autre_cell,[sommets_partages]}
            %cellules_en_contact{num_cell} = {} si pas de face avec une
            %cellule utile.
            %On remplit seulement si num_autre_cell>num_cell donc on
            %ne compte qu'une fois chaque face

            
            %vertices_to_cells{i} = cell contenant les num des cellules ayant
            %le sommet i. Utile pour trouver les cellules adjacente a une
            %cellule donnee
            vertices_to_cells = cell(1,length(sommets_voronoi));
            for iCellule = 1:length(cellules_voronoi)
                for iVertice = cellules_voronoi{iCellule}
                    if iVertice ~= 1
                        vertices_to_cells{iVertice}{length(vertices_to_cells{iVertice})+1} = iCellule;
                    end
                end
            end
            for iVertice = 1:length(sommets_voronoi)
                vertices_to_cells{iVertice} = cell2mat(vertices_to_cells{iVertice});
            end
            
            for num_cell = 1:length(cellules_voronoi)-1
                
                sommets_cell = cellules_voronoi{num_cell};
                numero_face = 1;
                faces_cell = cell(1,0);
                cell_useless = cellules_inutiles(num_cell);
                
                %recherche de face commune avec les cellules partageant au
                %moins un sommet
                cellules_adjacentes = cell(1,length(sommets_cell));
                for i = 1:length(sommets_cell)
                    cellules_adjacentes{i} = vertices_to_cells{sommets_cell(i)};
                end
                cellules_adjacentes = cell2mat(cellules_adjacentes);
                liste_cellules_adjacentes_d_indice_plus_grand = unique(cellules_adjacentes(cellules_adjacentes>num_cell));
                
                for num_autre_cell = liste_cellules_adjacentes_d_indice_plus_grand
                    if(not(isempty(num_autre_cell)))    
                        autre_cell_useless = cellules_inutiles(num_autre_cell);
                        if not(cell_useless && autre_cell_useless)
                            sommets_autre_cell = cellules_voronoi{num_autre_cell};
                            bool_sommets_commun = ismember(sommets_cell,sommets_autre_cell);   
                            if sum(bool_sommets_commun) >= dimension;
                                sommets_communs = sommets_cell(logical(bool_sommets_commun));
                                faces_cell{numero_face} = {num_autre_cell,sommets_communs};
                                numero_face = numero_face+1;
                            end
                        end
                    end
                end
                cellules_en_contact{num_cell} = faces_cell; %cellules_en_contact{num_cell}{num_face} = {num_autre_cell,[sommets_partages]}
            end
            
            %construction de dirty_faces, dirty_owners, dirty_neighbours 
            nbre_dirty_faces = sum(cellfun('length',cellules_en_contact)); %les faces ne sont comptees qu'une fois dans cellules_en_contact
            dirty_faces = cell(1,nbre_dirty_faces);     % dirty_faces : faces interieures avant nettoyage
            dirty_owners = zeros(1,nbre_dirty_faces);                                  
            dirty_neighbours = cell(1,nbre_dirty_faces);    
            indice_courant_face = 1;
            for num_cell = 1:length(cellules_voronoi)
                for num_face = 1:length(cellules_en_contact{num_cell})
                    %cellules_en_contact{num_cell}{num_face} = {num_autre_cell,[sommets_partages]}
                    num_autre_cell = cellules_en_contact{num_cell}{num_face}{1};
                    sommets_face = cellules_en_contact{num_cell}{num_face}{2};
                    dirty_faces{indice_courant_face} = sommets_face;
                    
                    cell_useless = cellules_inutiles(num_cell);
                    autre_cell_useless = cellules_inutiles(num_autre_cell);
                    
                    if cell_useless
                        assert(not(autre_cell_useless),'Pb reperage faces')
                        dirty_owners(indice_courant_face) = num_autre_cell; 
                        %trouver frontiere dont la cellule est esclave
                        dirty_neighbours{indice_courant_face} = []; %code pour pas encore de frontiere associee
                        for num_frontiere = 1:length(indices_cellules_esclaves)
                            if ismember(num_cell,indices_cellules_esclaves{num_frontiere}) || ismember(num_autre_cell,indices_cellules_esclaves{num_frontiere})
                                dirty_neighbours{indice_courant_face} = [-num_frontiere,dirty_neighbours{indice_courant_face}]; %-i : code pour frontiere i associee
                            end
                        end
                    elseif autre_cell_useless
                        dirty_owners(indice_courant_face) = num_cell;
                        %trouver frontiere dont la cellule est esclave
                        dirty_neighbours{indice_courant_face} = []; %code pour pas encore de frontiere associee
                        for num_frontiere = 1:length(indices_cellules_esclaves)
                            if ismember(num_cell,indices_cellules_esclaves{num_frontiere}) || ismember(num_autre_cell,indices_cellules_esclaves{num_frontiere})
                                dirty_neighbours{indice_courant_face} = [-num_frontiere,dirty_neighbours{indice_courant_face}]; %-i : code pour frontiere i associee
                            end
                        end
                    else
                        dirty_owners(indice_courant_face) = num_cell; 
                        dirty_neighbours{indice_courant_face} = num_autre_cell; %ici num_autre_cell>num_cell
                    end
                                        
                    indice_courant_face = indice_courant_face+1;
                end
            end
            assert(indice_courant_face == nbre_dirty_faces+1,'Pb reperage des faces');
 
                        
            %Construction de faces_des_frontieres_exterieures permettant de
            %passer d'une frontiere aux faces frontieres associees.
            nombre_frontieres_exterieures = 0;
            for iBoundary = 1:myGeometry.GetNumberOfBoundaries
                if strcmp(myGeometry.GetBoundaryType(iBoundary),'surface')
                    nombre_frontieres_exterieures = nombre_frontieres_exterieures+1;
                    
                end
            end
            faces_des_frontieres_exterieures = cell(1,nombre_frontieres_exterieures);                             
            %faces_des_frontieres_exterieures{i} = {numero de la boundary associee, [numero des faces liees a la frontiere i]}
            
            a = zeros(1,nbre_dirty_faces);
            for i = 1:length(a)
                if length(dirty_neighbours{i})>0
                    if dirty_neighbours{i}(1) <= 0
                        a(i) = 1;
                    end
                end
            end
            
            indice_frontiere_exterieure = 1;
            
            b = ismember(dirty_owners,(1-cellules_inutiles).*(1:length(cellules_inutiles)));
            assert(sum(and(a,not(b))) == 0,'Pb reperage faces exterieures');
            for iBoundary = 1:myGeometry.GetNumberOfBoundaries
                if strcmp(myGeometry.GetBoundaryType(iBoundary),'surface')
                    
                    faces_des_frontieres_exterieures{indice_frontiere_exterieure}{1} = iBoundary;
                    %trouver numeros_faces_frontiere : faces entre
                    %une cellule utile esclave de la frontiere
                    %et une cellule inutile 
                    
                    d = zeros(1,nbre_dirty_faces);
                    for i = 1:length(d)
                      if ismember(-iBoundary,dirty_neighbours{i})
                        d(i) = 1;
                      end
                    end
                                            
                    numeros_faces_frontiere = find(b &d);
                    %critere : neighbour = -indice_boundary, owner est utile
                    faces_des_frontieres_exterieures{indice_frontiere_exterieure}{2} = numeros_faces_frontiere;

                    indice_frontiere_exterieure = indice_frontiere_exterieure+1;
                end
            end     
            nbre_faces_frontieres_repertoriees = 0;
            for i = 1:nombre_frontieres_exterieures
                nbre_faces_frontieres_repertoriees = nbre_faces_frontieres_repertoriees+length(faces_des_frontieres_exterieures{i}{2});
            end
            %cas oe des faces exterieures appartiennent a 2+ frontieres
            face_to_front = zeros(1,nbre_dirty_faces);
            indice_doublon = 1;
            doublons = {};
            if sum(a) ~= nbre_faces_frontieres_repertoriees 
                %reperage des faces appartenant a 2 ou plus frontieres
                for indice_frontiere_exterieure = 1:nombre_frontieres_exterieures
                    for num_face = faces_des_frontieres_exterieures{indice_frontiere_exterieure}{2}
                        if face_to_front(num_face) == 0
                            face_to_front(num_face) = indice_frontiere_exterieure;
                        elseif face_to_front(num_face)<0 %cas doublon deje repere (face appartenant a 3+ frontieres)
                            doublons{-face_to_front(num_face)}{2} = [doublons{-face_to_front(num_face)}{2},indice_frontiere_exterieure];
                        else %cas nouveau doublon
                            doublons{indice_doublon}{1} = num_face; 
                            doublons{indice_doublon}{2} = [face_to_front(num_face),indice_frontiere_exterieure];
                            face_to_front(num_face) = -indice_doublon;
                            indice_doublon = indice_doublon+1;
                        end
                    end
                end
                %rerepartition des doublons dans les frontieres : la face
                %appartient a la  frontiere la plus proche de son centre
                for indice_doublon = 1:length(doublons)
                    num_face = doublons{indice_doublon}{1};
                    sommets = dirty_faces{num_face};
                    centre_face = mean(sommets_voronoi(sommets,:));
                    distances_frontieres = zeros(1,length(doublons{indice_doublon}{2}));
                    for j = 1:length(doublons{indice_doublon}{2})
                        numBoundary=doublons{indice_doublon}{2}(j);
                        vertices_frontiere = myGeometry.GetBoundaryVertices(numBoundary);
                        distances_frontieres(j) = NetworkBuilder.DistancePointFrontierePlane(centre_face,vertices_frontiere);
                    end
                    [~,position_min] = min(distances_frontieres);
                    num_new_frontiere = doublons{indice_doublon}{2}(position_min);
                    %enlever la face des frontieres trop eloignees
                    for k = doublons{indice_doublon}{2}
                        if k ~= num_new_frontiere
                            old_faces = faces_des_frontieres_exterieures{k}{2};
                            new_faces = old_faces(old_faces ~= num_face);
                            faces_des_frontieres_exterieures{k}{2} = new_faces;
                        end
                    end
                end
            end
            %Remettre -1 en neighbour des faces frontieres
            dirty_neighb = dirty_neighbours;
            dirty_neighbours = zeros(1,nbre_dirty_faces);
            for i = 1:nbre_dirty_faces
                if length(dirty_neighb{i})>0 && dirty_neighb{i}(1)>0
                    dirty_neighbours(i) =  dirty_neighb{i}(1);
                else
                    dirty_neighbours(i) =  -1;
                end
              
            end
            
        end %FaceExtraction
        
        %MiseAJourCellulesSommets
        function [cells_to_vertices,vertices,new_numeros_cellules,new_numeros_sommets,nombre_cellules_utiles] = MiseAJourCellulesSommets(cellules_voronoi,sommets_voronoi,cellules_inutiles)                                                                                                         

                %Determiner les sommets utiles a partir des cellules utiles
                %sommet utile = utile pour au moins une cellule
            cellules_utiles = 1-cellules_inutiles;
            
            sommets_utiles = zeros(1,length(sommets_voronoi));   
            for indice_cellule = find(cellules_utiles)
                for indice_sommet = cellules_voronoi{indice_cellule}  
                    sommets_utiles(indice_sommet) = 1;
                end             
            end
                %renumeroter les sommets et reorganiser les cellules
                %determiner les nouveaux numeros des sommets
            new_numeros_sommets = zeros(1,length(sommets_voronoi));  % new_numeros_sommets(ancien_numero) = nouveau_numero
            numero_courant_sommet = 1;
            for i = 1:length(sommets_voronoi)
                if sommets_utiles(i) == 1
                    new_numeros_sommets(i) = numero_courant_sommet;     %vectoriser ces boucles
                    numero_courant_sommet = numero_courant_sommet+1;  
                end
            end
            nombre_sommets_utiles = numero_courant_sommet-1;             
            
                %creer la nouvelle liste des sommets
            dimension = length(sommets_voronoi(1,:));
            vertices = zeros(nombre_sommets_utiles,dimension);
            numero_courant_sommet = 1;
            for i = 1:length(sommets_voronoi)
                if not(new_numeros_sommets(i) == 0)
                    vertices(numero_courant_sommet,:) = sommets_voronoi(i,:);     %vectoriser ces boucles
                    numero_courant_sommet = numero_courant_sommet+1;  
                end
            end

                %determiner nouveau numero des cellules
            new_numeros_cellules = zeros(1,length(cellules_voronoi)); % new_numeros_cellules(ancien_numero) = nouveau_numero                
            numero_courant_cell = 1;
            for i = 1:length(cellules_voronoi)
                if cellules_utiles(i) == 1
                    new_numeros_cellules(i) = numero_courant_cell;     %vectoriser ces boucles
                    numero_courant_cell = numero_courant_cell+1;
                end
            end
            nombre_cellules_utiles = numero_courant_cell-1;
            
                %creer la liste cells_to_vertices
            cells_to_vertices = cell(1,nombre_cellules_utiles);
            numero_courant_cellule = 1;
            for i = find(new_numeros_cellules)
                    new_cell = cellules_voronoi(i);
                    new_content_cell = zeros(1,length(cellules_voronoi{i}));
                    for j = 1:length(cellules_voronoi{i})
                        new_content_cell(j) = new_numeros_sommets(cellules_voronoi{i}(j));
                    end
                    new_cell{1} = new_content_cell;
                    cells_to_vertices(numero_courant_cellule) = new_cell;    
                    numero_courant_cellule = numero_courant_cellule+1;  
            end
            
        end %MiseAJourCellulesSommets
        
        %%BuildPoresToLinks
        function pores = BuildPoresToLinks(owners,neighbours,nPore)
            %construit cellules : cellules{i} = liste des faces de la cellule
            %i
            pores = cell(1,nPore);
            for iLink = 1:length(owners)
                pores{owners(iLink)} = [pores{owners(iLink)},iLink];
                if neighbours(iLink)>0
                    pores{neighbours(iLink)} = [pores{neighbours(iLink)},iLink];
                end
            end
        end %BuildPoresToLinks
        
        %%RenumerotationFaces
        function [faces,owners,neighbours,boundaries] = RenumerotationFaces(dirty_faces,dirty_owners,dirty_neighbours,faces_des_frontieres_exterieures,myGeometry,new_numeros_cellules,new_numeros_sommets,vertices)
            %Mise a jour des listes faces, owners, neighbours.
            %renumerotation des faces pour acces rapide aux faces d'une
            %frontiere et codage des infos frontieres du reseau de pores
            %Input : - faces, cells, owners, neighbours, faces_des_frontieres_exterieures
            %          issus de GenerateMesh, avec les numeros des faces avant renumerotation
            %        - structure boundaries de la geometrie macroscopique : geometrie_macroscopique.Boundaries.Boundary
            %Output : - newfaces, newcells, newowners, newneighbours avec nouveaux
            %           numeros des faces, structurees comme les proprietes
            %           de la geometrie
            %         - newboundaries : infos structurees pour l'output reseau
            %           de pores
            
                        
            %Mise a jour des listes faces, owners, neighbours pour etre
            %coherent avec les nouveaux numeros des sommets et cellules.
            nbre_faces = length(dirty_owners);
            new_faces = cell(1,nbre_faces);
            new_owners = zeros(1,nbre_faces);
            new_neighbours = zeros(1,nbre_faces);
            for indice_face = 1:nbre_faces
                new_faces{indice_face} = new_numeros_sommets(dirty_faces{indice_face});
                new_owners(indice_face) = new_numeros_cellules(dirty_owners(indice_face));
                if not(dirty_neighbours(indice_face) == -1)
                    new_neighbours(indice_face) = new_numeros_cellules(dirty_neighbours(indice_face));
                else
                    new_neighbours(indice_face) = -1;
                end
            end
            
            
            
            %Etablissement de la renumerotation des faces frontieres
            nFrontiere = length(faces_des_frontieres_exterieures); 
            infos_liens_frontieres = cell(1,nFrontiere);
            for i = 1:nFrontiere
                infos_liens_frontieres{i} = faces_des_frontieres_exterieures{i}{2};
            end
            
            [boundaries,owners,neighbours,~,faces] = NetworkBuilder.RenumerotationLiensFrontieres(infos_liens_frontieres,new_owners,new_neighbours,new_faces);
            
            %Ordonner sommmets des faces selon l'ordre trigonometrique ou
            %anti-trigonometrique, de faeon a avoir les aretes des faces
            dimension = length(vertices(1,:));
            switch dimension
                case 2
                    %rien a faire, voronoin renvoie deje les sommets ordonnes en 2D
                case 3
                    nFaces = length(faces);
                    for iFace = 1:nFaces
                        this_face = faces{iFace};
                        coordonnees_polygone = vertices(this_face,:);                          
                        %coordonnees_planes = PolygoneProjetterDansPlan(coordonnees_polygone);
                        %faces{iFace} = this_face(PolygoneOrdonnerSommetsDansPlan(coordonnees_planes));
                        [~,I] = angleSort3d(coordonnees_polygone);
                        faces{iFace} = this_face(I);
                        
                    end
            end
            
            %Ajout d'informations de MacroscopicGeometry dans la structure boundaries
                        
            for indice_frontiere = 1:nFrontiere
                numBoundary=faces_des_frontieres_exterieures{indice_frontiere}{1};
                boundaries.Boundary(indice_frontiere).CONTENT = [];
                boundaries.Boundary(indice_frontiere).ATTRIBUTE.Name = myGeometry.GetBoundaryName( numBoundary );
                boundaries.Boundary(indice_frontiere).ATTRIBUTE.Type = myGeometry.GetBoundaryType( numBoundary );
                
                boundaryVertices = myGeometry.GetBoundaryVertices( numBoundary );
                dimension = myGeometry.GetDimension;
                switch dimension
                    case 3
                        boundaries.Boundary(indice_frontiere).ATTRIBUTE.Area = abs(polygonArea3d(boundaryVertices));
                    case 2
                        boundaries.Boundary(indice_frontiere).ATTRIBUTE.Area = abs(polygonArea(boundaryVertices));
                end
            end
            
            
        end %RenumerotationFaces
        
        %RenumerotationFacesFrontieres
        function [boundaries,owners,neighbours,newOrder,faces] = RenumerotationLiensFrontieres(infos_liens_frontieres,old_owners,old_neighbours,varargin)
            %renumerote les liens pour respecter le codage des liens
            %frontieres. 
            %input : infos_liens_frontieres,old_owners,old_neighbours,
            %(varargin : cells to renumerotate)
            %output : [boundaries,owners,neighbours,newOrder] +(faces), où
            %   newOrder est tel que owners = old_owners(newOrder)
            
            nLink = length(old_owners);
            nExternalBoundary = length(infos_liens_frontieres);      
            boundary(nExternalBoundary) = struct;

            linkNumberLookUp = zeros(1,nLink);
            firstBoundaryLink = 1 ; 
            
            for iBoundary = 1:nExternalBoundary  
                
                nBoundaryLink = length(infos_liens_frontieres{iBoundary});
                boundary(iBoundary).ATTRIBUTE.NombreFaces = nBoundaryLink;
                boundary(iBoundary).ATTRIBUTE.StartFace = firstBoundaryLink;
                
                for iBoundaryLink = 1:nBoundaryLink
                    iLink=infos_liens_frontieres{iBoundary}(iBoundaryLink);
                    assert(linkNumberLookUp(iLink) == 0);
                    linkNumberLookUp(iLink) = firstBoundaryLink+iBoundaryLink-1;
                end
                firstBoundaryLink = firstBoundaryLink+nBoundaryLink;
            end
            boundaries.Boundary = boundary;

            nTotalBoundaryLink = 0;
            for i = 1:nExternalBoundary
                nTotalBoundaryLink = nTotalBoundaryLink+length(infos_liens_frontieres{i});
            end
            
            linkNumberLookUp(linkNumberLookUp == 0) = (nTotalBoundaryLink+1):nLink;
%            assert(isequal(ismember(correspondance_numeros_faces,(1:nombre_faces)),(1:nombre_faces)),'Pb renumerotation faces')

            %Renumerotating boundary links
            
            owners = zeros(1,nLink);
            neighbours = zeros(1,nLink);
            newOrder = zeros(1,nLink);
            for oldIndices = 1:nLink
                newIndices = linkNumberLookUp(oldIndices);
                newOrder(newIndices) = oldIndices;
                owners(newIndices) = old_owners(oldIndices);
                neighbours(newIndices) = old_neighbours(oldIndices);
            end
            
            faces=cell(1,nLink);
            if ~isempty(varargin)
                for oldIndices = 1:nLink
                    newIndices = linkNumberLookUp(oldIndices);
                    faces{newIndices} = varargin{1}{oldIndices};
                end
            end
            
            
        end   %RenumerotationFacesFrontieres
        

        
        %ExtractEdges
        function [edges,vertices_to_edges,faces_to_edges,edges_to_faces] = ExtractEdges(vertices,faces,networkType)
            %Extrait les aretes du maillage.
            %output : - edges : edges(i,:) = [un sommets de l'edge i, l'autre sommet de l'edge i]
            %         - vertices_to_edges{i} = [num des edges du sommet i] :
            %         pour passer des sommets aux aretes
            
            dimension = length(vertices(1,:));
            
            if dimension == 2 && strcmp(networkType,'PoreNetworkUnstructuredMeshFibrous')
                %dans ce cas les edges sont les vertices
                edges = horzcat((1:length(vertices)).',(1:length(vertices)).');
                vertices_to_edges = cell(1,length(vertices));
                for i = 1:length(vertices)
                    vertices_to_edges{i} = i;
                end
                
                nFace = length(faces);
                faces_to_edges = cell(nFace,1);
                nEdge = length(edges);
                edges_to_faces = cell(1,nEdge);
                for iFace = 1:nFace
                    faces_to_edges{iFace} = faces{iFace};
                    
                    for iEdge = faces{iFace}
                        edges_to_faces{iEdge} = [edges_to_faces{iEdge},iFace];
                    end
                end
                
            else
                %Creation d'une liste d'aretes non uniques a partir des
                %sommets ordonnes des polygones. Les sommets d'une arete
                %sont donnes dans l'ordre indice_sommet_1<indice_sommet_2
                nFace = length(faces);
                edges_non_uniques = cell(nFace,1);
                for iFace = 1:nFace
                    nSommet = length(faces{iFace});
                    edges_non_uniques{iFace} = zeros(nSommet,2);
                    for iSommet = 1:(nSommet-1)
                       edges_non_uniques{iFace}(iSommet,:) = sort([faces{iFace}(iSommet),faces{iFace}(iSommet+1)]);
                    end
                    edges_non_uniques{iFace}(nSommet,:) = sort([faces{iFace}(nSommet),faces{iFace}(1)]);
                end
                foo = cell2mat(edges_non_uniques);
                %Filtration des aretes pour obtenir l'unicite :                   
                [edges,~,n] = unique(foo,'rows');

                %Creation de faces_to_edges et edges_to_faces
                nEdge = length(edges);
                faces_to_edges = cell(nFace,1);
                edges_to_faces = cell(1,nEdge);
                ind_courant = 1;
                for iFace = 1:nFace
                    nOldEdge = length(edges_non_uniques{iFace});
                    faces_to_edges{iFace} = zeros(1,nOldEdge);
                    for iOldEdge = 1:nOldEdge
                        newNumEdge = n(ind_courant);
                        faces_to_edges{iFace}(iOldEdge) = newNumEdge;
                        ind_courant = ind_courant+1;
                        edges_to_faces{newNumEdge} = [edges_to_faces{newNumEdge},iFace];
                    end
                end
                assert(ind_courant == length(foo)+1);

                %Creation de vertices_to_edges                                           
                vertices_to_edges = cell(1,length(vertices));   
                for iEdges = 1:length(edges)
                    sommet_1 = edges(iEdges,1);
                    vertices_to_edges{sommet_1}{length(vertices_to_edges{sommet_1})+1} = iEdges;
                    sommet_2 = edges(iEdges,2);
                    vertices_to_edges{sommet_2}{length(vertices_to_edges{sommet_2})+1} = iEdges;                        
                end
                for iVertice = 1:length(vertices)
                    vertices_to_edges{iVertice} = cell2mat(vertices_to_edges{iVertice});
                end
                    
            end
                       
        end %ExtractEdges
        
        %GenerateEdgeThickness
        function edgeThickness = GenerateEdgeThickness(edges,vertices,myGeometry)
            %Attribue un epaisseur a chaque arete en fonction des parametres
            %de la geometrie macroscopique
            %output: epaisseur_edges(i) = epaisseur de l'arete i
            %Algorithme : classer les sommets par block. Attribuer aux
            %aretes dont les deux sommets sont dans le meme block
            %l'epaisseur voulue dans ce block. Pour les aretes ayant leurs
            %sommets dans deux blocks differents, gerer en fonction des
            %conditions frontieres.
            
            %classement des vertices par block
            nVertices = length(vertices(:,1));
            vertices_to_blocks = zeros(1,nVertices);
            nBlock = myGeometry.GetNumberOfBlocks;
            for iBlock = 1:nBlock
                block_vertices = myGeometry.GetBlockVertices(iBlock);
                vertices_to_blocks(inhull(vertices,block_vertices)) = iBlock;
            end
            
            vertices_out = find(vertices_to_blocks == 0);
            tol = 0.01;
            while not(isempty(vertices_out))    
                tol = tol*2;
                for iBlock = 1:nBlock
                    block_vertices = myGeometry.GetBlockVertices(iBlock);
                    zone_proche_du_bloc = (1+tol)*block_vertices-(tol)*ones(length(block_vertices),1)*mean(block_vertices);
                    vertices_to_blocks(vertices_out(inhull(vertices(vertices_out,:),zone_proche_du_bloc))) = iBlock;
                end
                vertices_out = find(vertices_to_blocks == 0);
            end
            
            %Repartition des aretes par block :
            %Les aretes dont les deux sommets sont dans le meme block vont
            %dans ce block le. Celles qui sont a cheval sur deux blocks
            %sont attribues au block de plus petit edge epaisseur.
            blocks_to_edges = cell(1,nBlock);
            nEdges = length(edges(:,1));
            for iEdge = 1:nEdges
                if edges(iEdge,1) == edges(iEdge,2)
                    num_block = vertices_to_blocks(edges(iEdge,1));
                    blocks_to_edges{num_block}{length(blocks_to_edges{num_block})+1} = iEdge;
                else
                    num_block_1 = vertices_to_blocks(edges(iEdge,1));
                    num_block_2 = vertices_to_blocks(edges(iEdge,1));
                    epaisseur1 = myGeometry.GetBlockFillingParameters(num_block_1).FiberThickness;
                    epaisseur2 = myGeometry.GetBlockFillingParameters(num_block_2).FiberThickness;
                    if epaisseur1<epaisseur2
                        blocks_to_edges{num_block_1}{length(blocks_to_edges{num_block_1})+1} = iEdge;
                    else
                        blocks_to_edges{num_block_2}{length(blocks_to_edges{num_block_2})+1} = iEdge;
                    end
                end
            end
            for iBlock = 1:nBlock
                blocks_to_edges{iBlock} = cell2mat(blocks_to_edges{iBlock});
            end
            
            %Attribution des epaisseurs des aretes en fonction des
            %parametres specifies pour chaque block dans la geometrie
            %macroscopique.
            edgeThickness = zeros(1,nEdges);
            for iBlock = 1:nBlock
                liste_sommets_1 = edges(blocks_to_edges{iBlock},1);
                liste_sommets_2 = edges(blocks_to_edges{iBlock},2);
                coordonnes_sommets_edges = horzcat(vertices(liste_sommets_1,:),vertices(liste_sommets_2,:));
                edgeThickness(blocks_to_edges{iBlock}) = NetworkBuilder.GenerateEdgeThicknessOneBlock(coordonnes_sommets_edges,myGeometry,iBlock);
            end
            
        end %GenerateEdgeThickness
        
        %GenerateEdgeThicknessOneBlock
        function epaisseurs = GenerateEdgeThicknessOneBlock(coordonnes_sommets_edges,myGeometry,iBlock)
            %input :
            %  -  coordonnes_sommets_edges  =  
            %  tableau(nbre_edges_dans_le_block,4 ou 6) : sur une ligne les
            %  coordonnees des sommets de l'arete [X1 Y1 (Z1) X2 Y2 (Z2)]
            %output: 
            %  - epaisseurs : tableau(nbre_edges_dans_le_block,1)
            nEdges = length(coordonnes_sommets_edges(:,1));
            %epaisseurs = myGeometry.GetConvertToMeters*myGeometry.GetBlockFillingParameters(iBlock).FiberThickness*ones(1,nEdges);
            epaisseurs = myGeometry.GetBlockFillingParameters(iBlock).FiberThickness*ones(1,nEdges);
            
        end %GenerateEdgeThicknessOneBlock

        %GereAnisotropieDebut
        function GereAnisotropieDebut(myGeometry)
            nBand = myGeometry.GetNumberOfAnisotropyBand ;
            if nBand>0
                switch myGeometry.GetAnisotropyBandDirection(1)
                    case 'x'
                        indice = 1;
                    case 'y'
                        indice = 2;
                    case  'z'
                        indice = 3;
                end
                for iBand = 1:nBand %ATTENTION, BUG POUR NZONE>1
                    location=myGeometry.GetAnisotropyBandLocation(iBand);
                    bas = min(location);
                    haut = max(location);
                    intensity = myGeometry.GetAnisotropyBandIntensity(iBand);
                    nMacroVertice = myGeometry.GetNumberOfVertices;
                    
                    for iMacroVertice = 1:nMacroVertice
                        vertice=myGeometry.GetVertice(iMacroVertice);
                        if vertice(indice)>haut
                            newCoordinates=vertice;
                            newCoordinates(indice) =  newCoordinates(indice)+(haut-bas)*(intensity-1);
                            myGeometry.ChangeVertice(iMacroVertice,newCoordinates);
                        
                        elseif vertice(indice)>bas
                            newCoordinates=vertice;
                            newCoordinates(indice) =  bas+(newCoordinates(indice)-bas)*intensity;
                            myGeometry.ChangeVertice(iMacroVertice,newCoordinates);
                        end
                    end
                end
            end
        end %GereAnisotropieDebut
        
        %GereAnisotropieFin
        function new_vertices = GereAnisotropieFin(myGeometry,vertices)
            nBand = myGeometry.GetNumberOfAnisotropyBand ;
            if  nBand>0
                switch myGeometry.GetAnisotropyBandDirection(1)
                    case 'x'
                        indice = 1;
                    case 'y'
                        indice = 2;
                    case  'z'
                        indice = 3;
                end
                new_vertices = vertices;
                
                for iBand = 1:nBand %ATTENTION, BUG POUR NZONE>1
                    location=myGeometry.GetAnisotropyBandLocation(iBand);
                    bas = min(location);
                    haut = max(location);
                    intensity = myGeometry.GetAnisotropyBandIntensity(iBand);
                    nVerticeMacro = myGeometry.GetNumberOfVertices;
                    
                    for iVerticeMacro = 1:nVerticeMacro
                        vertice=myGeometry.GetVertice(iVerticeMacro);
                        if vertice(indice)>bas
                            newCoordinates=vertice;
                            newCoordinates(indice) =  newCoordinates(indice)-(haut-bas)*(intensity-1);
                            myGeometry.ChangeVertice(iVerticeMacro,newCoordinates);
                        end
                    end
                    nVertice = length(vertices(:,1));
                    for iVertice = 1:nVertice
                        hauteur = new_vertices(iVertice,indice);
                        if hauteur >= bas+(haut-bas)*intensity
                            new_vertices(iVertice,indice) = new_vertices(iVertice,indice)-((haut-bas)*(intensity-1));
                        elseif hauteur>bas
                            new_vertices(iVertice,indice) = bas+(new_vertices(iVertice,indice)-bas)/intensity;
                        end
                    end
                    
                end
            else
                new_vertices = vertices;
            end
            
        end %GereAnisotropieFin
               
        %DistancePointFrontierePlane
        function distance = DistancePointFrontierePlane(point,boundaryVertices)
            %retoune la distance d'un point a une frontiere plane
            
            centreFrontiere = mean(boundaryVertices);
            dimension = length(point);
            assert(dimension == length(boundaryVertices(1,:)),'Pb coherence dimensions');
            normale = NetworkBuilder.ComputeVecteurNormal(boundaryVertices);
            distance = abs(dot(normale,point-centreFrontiere));
            
        end %DistancePointFrontierePlane
        
        %ProportionVolumeIn
        function proportion = ProportionVolumeIn(sommets_cell,blockVertices)
            %retourne la proportion du volume de la cellule qui est a
            %l'interieur du bloc.
            
            if ismember(Inf,sommets_cell)
                proportion = 0;
            else
                [~,volume_cell] = convhulln(sommets_cell);

                [A_cell,b_cell,Aeq_cell,beq_cell] = vert2lcon(sommets_cell);
                [A_block,b_block,Aeq_block,beq_block] = vert2lcon(blockVertices);
                A = vertcat(A_cell,A_block);
                b = vertcat(b_cell,b_block);
                Aeq = vertcat(Aeq_cell,Aeq_block);
                beq = vertcat(beq_cell,beq_block);

                sommets_intersection = lcon2vert(A,b,Aeq,beq,1e-4);
                if isempty(sommets_intersection)
                    volume_intersection = 0;
                else
                   [~,volume_intersection] = convhulln(sommets_intersection);
                end

                proportion = volume_intersection/volume_cell;
            end
        end %ProportionVolumeIn
        
        function vect_perp = ComputeVecteurNormal(planePoints)
            %renvoie les coordonnees d'un vecteur normal unitaire a un plan
            %en 3D, a un segment en 2D.
            %input : points_plan = coordonnees de 2 (resp. 4) points
            %         differents definissant un hyperplan en 2D (resp. 3D).
            dimension = length(planePoints(1,:));
            switch dimension
                case 2
                    assert(length(planePoints(:,1)) > 1,'au moins 2 points requis');
                    base_vect = planePoints(2,:)-planePoints(1,:);
                    foo = norm(base_vect);
                    base_vect = base_vect/foo;
                    vect_perp = [base_vect(2),-base_vect(1)];
                    
                case 3
                    assert(size(planePoints,1) > 2,'au moins 3 points requis');

                    base_vect_1 = (planePoints(2,:)-planePoints(1,:));
                    base_vect_1 = base_vect_1/norm(base_vect_1);
                    base_vect_2 = (planePoints(3,:)-planePoints(1,:));
                    base_vect_2 = base_vect_2/norm(base_vect_2);
                    
                    if dot(base_vect_1,base_vect_2)<1e-8 && size(planePoints,1) > 3
                        base_vect_2 = (planePoints(4,:)-planePoints(1,:));
                        base_vect_2 = base_vect_2/norm(base_vect_2);
                    end
                    foo = dot(base_vect_1,base_vect_2);
                    if abs(foo)<1e-6
                        disp('pb computing normal vector')
                    end
                    vect_perp = cross(base_vect_1,base_vect_2);
                    vect_perp = vect_perp/norm(vect_perp);
            end
        end
        
        function vect_perp = ComputeVecteurNormalOrienteExterieurBlock(points_plan,blockCenter)
            vect_perp = NetworkBuilder.ComputeVecteurNormal(points_plan);
            dir_int = blockCenter-points_plan(1,:);
            dir_int = dir_int/norm(dir_int);
            if dot(vect_perp,dir_int)>0
                vect_perp = -vect_perp;
            end
        end
    end 
end

