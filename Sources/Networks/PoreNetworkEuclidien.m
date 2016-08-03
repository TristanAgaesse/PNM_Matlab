classdef PoreNetworkEuclidien < PoreNetwork
    %PoreNetworkEuclidien : Subclass of PoreNetwork
    %   L'element en plus par rapport à PoreNetwork est que les pores et
    %   les liens ont des centres (points de l'espace)
    
    properties 
        PoreCenter %tableau (NombrePores,Dimension) contenant les coordonnées des centres des pores
        LinkCenter %tableau (NombreLiens,Dimension) contenant les coordonnées des centres des liens
    end
    
    methods

        function pore_network_euclidien=PoreNetworkEuclidien(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry)
            %constructeur
            %input : dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter
            %output : pore_network_euclidien
            pore_network_euclidien=pore_network_euclidien@PoreNetwork(dimension,pores,owners,neighbours,boundaries,myGeometry);
            
            pore_network_euclidien.PoreCenter=poreCenter;
            pore_network_euclidien.LinkCenter=linkCenter;
            
        end
        
        
        function center=GetPoreCenter(poreNetwork,iPore)
            %input : poreNetwork,iPore
            %output : center
            center=poreNetwork.PoreCenter(iPore,:);
        end      
        
        function center=GetLinkCenter(poreNetwork,iLink)
            %input : poreNetwork,iLink
            %output : center
            center=poreNetwork.LinkCenter(iLink,:);
        end
        


        function [boundaryXYZ,lengthXYZ,surfaceXYZ] = GetCuboidDomain(network)
            % If the network is in a cuboid domain, returns the
            % caracteristics of the cuboid domain. Otherwise it may
            % return an error.
            % Input : network
            % Output : boundaryXYZ{1 (resp.2,3,4,5,6) } = boundary code for Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
            %          lengthXYZ = [Lx,Ly,Lz]
            %          surfaceXYZ = [Ly*Lz,Lx*Lz,Lx*Ly]

            %Compute domain dimensions

            centers = network.GetLinkCenter(1:network.GetNumberOfLinks);
            dimension=network.GetDimension;

            Xmin = min(centers(:,1));
            Xmax = max(centers(:,1));
            Ymin = min(centers(:,2));
            Ymax = max(centers(:,2));
            lengthXYZ = [Xmax-Xmin,Ymax-Ymin];
            surfaceXYZ = [lengthXYZ(2),lengthXYZ(1)];
            if dimension==3
                Zmin = min(centers(:,3));
                Zmax = max(centers(:,3));
                lengthXYZ = [Xmax-Xmin,Ymax-Ymin,Zmax-Zmin];
                surfaceXYZ = [lengthXYZ(2)*lengthXYZ(3),lengthXYZ(1)*lengthXYZ(3),lengthXYZ(1)*lengthXYZ(2)];
            end


            % Assign a domain face to each network boundary

            nBoundary = network.GetNumberOfBoundaries;
            boundaryMean = zeros(nBoundary,network.GetDimension);
            for iBoundary = 1:nBoundary
                boundaryMean(iBoundary,:) = mean(network.GetLinkCenter(network.GetLinksFrontiere(iBoundary)),1);
            end

            boundaryMean(:,1) = boundaryMean(:,1)-Xmin;
            boundaryMean(:,1) = boundaryMean(:,1)./lengthXYZ(1);
            boundaryMean(:,2) = boundaryMean(:,2)-Ymin;
            boundaryMean(:,2) = boundaryMean(:,2)./lengthXYZ(2);
            if dimension==3
                boundaryMean(:,3) = boundaryMean(:,3)-Zmin;
                boundaryMean(:,3) = boundaryMean(:,3)./lengthXYZ(3);
            end

            boundaryXYZ = cell(1,2*dimension);   

            for iBoundary = 1:nBoundary

                if boundaryMean(iBoundary,1) < 0.1
                    boundaryXYZ{1}=[boundaryXYZ{1},iBoundary];
                elseif boundaryMean(iBoundary,1) > 0.9
                    boundaryXYZ{2}=[boundaryXYZ{2},iBoundary];
                end

                if boundaryMean(iBoundary,2) < 0.1
                    boundaryXYZ{3}=[boundaryXYZ{3},iBoundary];
                elseif boundaryMean(iBoundary,2) > 0.9
                    boundaryXYZ{4}=[boundaryXYZ{4},iBoundary];
                end

                if dimension==3
                    if boundaryMean(iBoundary,3) < 0.1
                        boundaryXYZ{5}=[boundaryXYZ{5},iBoundary];
                    elseif boundaryMean(iBoundary,3) > 0.9
                        boundaryXYZ{6}=[boundaryXYZ{6},iBoundary];
                    end
                end

            end   

            checkList=zeros(1,nBoundary);
            for i = 1:length(boundaryXYZ)
                checkList(boundaryXYZ{i})=checkList(boundaryXYZ{i})+1;
            end
            assert(~any(checkList-ones(1,nBoundary)),'some boundaries are unassigned to a domain face') 
        end
        
        function linkLength = ComputeAllLinkLength(network)
            %input : network
            %output : linkLength

            %Get geometric data from the network
            nLink = network.GetNumberOfLinks;
            nPore = network.GetNumberOfPores;
            dimension = network.Dimension;
            poreCenter=network.GetPoreCenter(1:nPore);
            linkCenter=network.GetLinkCenter(1:nLink);

            allLinks=1:nLink;
            internalLinks = network.GetLinksFrontiere(0);
            boundaryLinks = GetLinksFrontiere(network,1:network.GetNumberOfBoundaries);

            function myNorm=FastNorm(myVect,dimension)
                %Vectorial version of the norm function
                if dimension==2
                    myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2) ;
                elseif dimension==3
                    myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2+myVect(:,3).^2) ;
                end
            end 
            
            a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
            distance1=FastNorm(a,dimension);

            b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
            distance2=FastNorm(b,dimension);

            %Compute length
            linkLength = zeros(1,nLink);
            linkLength(internalLinks)=transpose(distance1(internalLinks)+distance2);
            linkLength(boundaryLinks)=2*distance1(boundaryLinks);
        end
                
        
        function ExportToBabe(network,numInlet,numOutlet)
            assert(network.GetDimension==3,'ExportToBabe requires a dimension 3 network');
            
            %input : network,numInlet,numOutlet
            file_name='Geometry';
            writer=FileWriterBabe(file_name);
            tic;
            disp('G�n�ration du fichier Geometry pour Babe...');
            donnees_output=network.PrivateBabeOutputStruct(numInlet,numOutlet);
            writer.Write(donnees_output);
            duree=toc;minutes=floor(duree/60);secondes=duree-60*minutes;
            fprintf('Fichier Geometry genere. Dur�e : %d minutes %f s. \n',minutes,secondes);
            writer.delete;clear('writer','file_name','file_content','donnees_output','duree','minutes','secondes');
        end
        
        function ExportToParaview(network,filename)
            %input :network, filename
            writer=FileWriterVTK(filename);
            tic;
            disp('G�n�ration du fichier VTK...');
            donnees_output=network.PrivateVTKOutputStructBallAndStick;
            writer.Write(donnees_output);
            duree=toc;minutes=floor(duree/60);secondes=duree-60*minutes;
            disp(sprintf('Fichier VTK g�n�r�. Dur�e : %d minutes %f s.',minutes,secondes));
            writer.delete;clear('writer','file_name','file_content','donnees_output','duree','minutes','secondes');
        end
        
        
        function vtk_struct=PrivateVTKOutputStructBallAndStick(network)
            %cr�ation en m�moire de la structure d'un fichier VTK POLYDATA 
            %pour afficher le r�seau de pores dans paraview sous forme ball
            %and stick.
            %input : network
            %output :vtk_struct

            dimension=network.Dimension;
            nPore=network.GetNumberOfPores;
            nLink=network.GetNumberOfLinks;
            pore_data=network.GetPoreDataList;
            link_data=network.GetLinkDataList;

            %Partie POINTS du fichier

            nPoints=nPore+2*nLink;
            points_output=zeros(nPoints,3);
            
            pores_to_points=cell(1,nPore); 
            for iPore=1:nPore
                points_output(iPore,1:dimension)=network.GetPoreCenter(iPore);
                pores_to_points{iPore}=iPore;
            end
            
            link_to_points=cell(1,nLink);
            if dimension==3
                for iLink=1:nLink
                    neighbourPore=network.GetPoresOfLink(iLink);
                    points_output(2*(iLink-1)+nPore+1,:)=network.GetPoreCenter(neighbourPore(2));
                    if neighbourPore(1)~=-1
                        points_output(2*(iLink-1)+1+nPore+1,:)=network.GetPoreCenter(neighbourPore(1));
                    else
                        points_output(2*(iLink-1)+1+nPore+1,:)=network.GetPoreCenter(neighbourPore(2));
                    end
                    
                    link_to_points{iLink}=[2*(iLink-1)+nPore+1,2*(iLink-1)+1+nPore+1];
                end
            else
                for iLink=1:nLink
                    neighbourPore=network.GetPoresOfLink(iLink);
                    points_output((2*iLink-1)+nPore,:)=[network.GetPoreCenter(neighbourPore(2)),0];
                    if neighbourPore(1)~=-1
                        points_output((2*iLink)+nPore,:)=[network.GetPoreCenter(neighbourPore(1)),0];
                    else
                        points_output((2*iLink)+nPore,:)=[network.GetPoreCenter(neighbourPore(2)),0];
                    end

                    link_to_points{iLink}=[(2*iLink-1)+nPore,(2*iLink)+nPore];
                end
            end
            
            %Partie VERTICES du fichier
            nVertice=nPore;
            vertices_output=horzcat(ones(nVertice,1),transpose(0:nPore-1));
            
            %Partie LINES du fichier
            nLine=nLink;
            lines_output=zeros(nLine,3);
            for iLine=1:nLine
                lines_output(iLine,1)=2;
                lines_output(iLine,2)=link_to_points{iLine}(1)-1;
                lines_output(iLine,3)=link_to_points{iLine}(2)-1;
            end
            
            %Partie POLYGON du fichier : vide
            polygons_output=cell(0,1);
            
            
            %Partie POINT_DATA du fichier
            point_data_output=struct;
                %transcription des link data
            names=fieldnames(link_data);
            for i=1:length(names)
                data_output=zeros(nPoints,1);
                data=link_data.(names{i});
                for iLink=1:nLink
                    data_output(link_to_points{iLink}(1))=data(iLink);
                    data_output(link_to_points{iLink}(2))=data(iLink);
                end
                point_data_output.(strcat('Link_',names{i}))=data_output;
            end
                %transcription des pore data
            names=fieldnames(pore_data);
            for i=1:length(names)
                data_output=zeros(nPoints,1);
                data=pore_data.(names{i});
                for iPore=1:nPore
                    data_output(pores_to_points{iPore})=data(iPore);
                end
                point_data_output.(strcat('Pore_',names{i}))=data_output;
            end    
                
              
            
            %Partie CELL_DATA du fichier
            nCellData=nLine+nVertice;
            cell_data_output=struct;
                %transcription des link data
            names=fieldnames(link_data);
            for i=1:length(names)
                data_output=zeros(nCellData,1);
                data=link_data.(names{i});
                for iLink=1:nLink
                    data_output(nVertice+iLink)=data(iLink);
                end
                cell_data_output.(strcat('Link_',names{i}))=data_output;
            end
                %transcription des pore data
            names=fieldnames(pore_data);
            for i=1:length(names)
                data_output=zeros(nCellData,1);
                data=pore_data.(names{i});
                for iPore=1:nPore
                    data_output(iPore)=data(iPore);
                end
                cell_data_output.(strcat('Pore_',names{i}))=data_output;
            end

            vtk_struct.Points=points_output;
            vtk_struct.Lines=lines_output;
            vtk_struct.Vertices=vertices_output;
            vtk_struct.Polygons=polygons_output;
            vtk_struct.PointData=point_data_output;
            vtk_struct.CellData=cell_data_output;
            
        end       
        
        
        
        
        function babe_struct=PrivateBabeOutputStruct(network,numInlet,numOutlet)
            %cr�ation en m�moire de la structure d'un fichier Geometry du
            %code Babe (code PNM C++ du CEA)
            %input : network,numInlet,numOutlet
            %output :  babe_struct           
            
            nPore=network.GetNumberOfPores;
            nLink=network.GetNumberOfLinks;
            nAgglomerate=0;
            nConnexion=0;
            
            %Partie du fichier NetworkSize
            
            babe_struct.NetworkPosition=[0 0 0];
            
            babe_struct.NetworkLength=[ 1 1 1  ];   %TO DO !!!!!!!!!!!!!!!!!!!!!!!!!
            
            
            %Partie du fichier Pores
            
            babe_struct.NumberOfPores=nPore;
            domain='CCL';
            shape='Spherical';
            numericData=cell(1,nPore);
            stringData=cell(1,nPore);
            
            for iPore=1:nPore
                diameter=network.GetPoreDiameter(iPore);
                center=network.GetPoreCenter(iPore);
                numericData{iPore}=[iPore-1,diameter,center(1),center(2),center(3)];
                stringData{iPore}{1}=domain;
                stringData{iPore}{2}=shape;
            end
            
            
            babe_struct.Pores.NumericData=numericData;
            babe_struct.Pores.StringData=stringData;
            
            
            %Partie du fichier Throats
            
            babe_struct.NumberOfThroats=nLink;
            domain='CCL';
            shape='Spherical';
            numericData=cell(1,nLink);
            stringData=cell(1,nLink);
            
            for iLink=1:nLink
                diameter=network.GetLinkDiameter(iLink);
                center=network.GetLinkCenter(iLink);
                numericData{iLink}=[iLink-1,diameter,center(1),center(2),center(3)];
                stringData{iLink}{1}=domain;
                
                numFrontiere=network.GetFrontiereOfLink(iLink);
                
                if ismember(numFrontiere,numOutlet)
                    boundaryType='Exit';
                    
                elseif numFrontiere==0
                    boundaryType='Initialized';
                    
                else
                    boundaryType='BoundaryInitialized';
                end
                
                stringData{iLink}{2}=boundaryType;
                stringData{iLink}{3}=shape;
            end

            babe_struct.Throats.NumericData=numericData;
            babe_struct.Throats.StringData=stringData;
            
            
            %Partie du fichier Agglomerates
            
            babe_struct.NumberOfAgglomerates=nAgglomerate;
            
            numericData=cell(1,nAgglomerate);
            stringData=cell(1,nAgglomerate);
            
            babe_struct.Agglomerates.NumericData=numericData;
            babe_struct.Agglomerates.StringData=stringData;
            
            
            %Partie du fichier Connexions
            
            babe_struct.NumberOfConnexions=nConnexion;
            
            numericData=cell(1,nConnexion);
            stringData=cell(1,nConnexion);
            
            babe_struct.Connexions.NumericData=numericData;
            babe_struct.Connexions.StringData=stringData;
            
            
            %Partie du fichier PoreLinking
            
            poreLinking=cell(1,nPore);
            for iPore=1:nPore
                poreLinking{iPore}=[iPore-1,network.GetLinksOfPore(iPore)-1];
            end
            babe_struct.PoreLinking=poreLinking;
            
            
            %Partie du fichier ThroatLinking
            
            throatLinking=cell(1,nPore);
            for iLink=1:nLink
                poreOfLink=network.GetPoresOfLink(iLink);
                throatLinking{iLink}=[iLink-1,poreOfLink(poreOfLink>0)-1];
            end
            babe_struct.ThroatLinking=throatLinking;
            
            
            %Partie du fichier AgglomerateToPoreLinking
            
            agglomeratesToPoreLinking=cell(1,nAgglomerate);
            babe_struct.AgglomeratesToPoreLinking=agglomeratesToPoreLinking;
            
            
            %Partie du fichier ConnexionLinking           
            
            connexionToPoreLinking=cell(1,nConnexion);
            babe_struct.ConnexionLinking=connexionToPoreLinking;

        end       
        

        
    end

end