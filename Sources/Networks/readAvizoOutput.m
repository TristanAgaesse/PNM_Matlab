
%Lit les output de Avizo pour générer un reseau de pores
%Les noms des fichiers générés par Avizo utilisés sont indiqués dans le
%fichier source de cette fonction.


%-------------------------------------------------------------------------------------
%% Indiquer dans cette section les noms des fichiers et des variables .mat
%--------------------------------------------------------------------------------------

%Voxel edge length (taille d'un voxel)
voxelEdgeLength=1e-6;                                             %Parametre modifiable

%Nombre de frontieres
nombreFrontieres=6;                                               %Parametre modifiable

% Pores (numérotés)
Avizo_pore = load('OUTPUT-PoresSegmentes.mat');                     %Parametre modifiable

% Liens internes : interfaces entre les pores (dilatées d'un pixel et numérotées)
Avizo_interface = load('OUTPUT-LiensDilates.mat');                  %Parametre modifiable

% Liens frontieres : ortho-slices des pores au niveau de toutes les frontieres
Avizo_frontiere=cell(1,nombreFrontieres);
Avizo_frontiere{1}=load('OUTPUT-PoresSegmentes-Ortho-Slice.mat');   %Parametre modifiable
Avizo_frontiere{2}=load('OUTPUT-PoresSegmentes-Ortho-Slice2.mat');  %Parametre modifiable
Avizo_frontiere{3}=load('OUTPUT-PoresSegmentes-Ortho-Slice3.mat');  %Parametre modifiable
Avizo_frontiere{4}=load('OUTPUT-PoresSegmentes-Ortho-Slice4.mat');  %Parametre modifiable
Avizo_frontiere{5}=load('OUTPUT-PoresSegmentes-Ortho-Slice5.mat');  %Parametre modifiable
Avizo_frontiere{6}=load('OUTPUT-PoresSegmentes-Ortho-Slice6.mat');  %Parametre modifiable


%Fichiers contenant les informations géométriques calculées par Avizo
%Ils doivent être au format .xlsx . Pour convertir les fichier .xml générés
%par Avizo au format .xlsx, utiliser la commande bash LibreOffice suivante :
% soffice --headless -convert-to xlsx:"Calc MS Excel 2007 XML" *.xml
%cette commande peut être lancée depuis Matlab en faisant system(...)

tablesXLSX=struct;
tablesXLSX.Liens='OUTPUT-LiensStatistiques.xlsx';                    %Parametre modifiable
tablesXLSX.Pores='OUTPUT-PoresSegmentesStatistiques.xlsx';           %Parametre modifiable

tablesXLSX.Frontieres=cell(1,nombreFrontieres);
tablesXLSX.Frontieres{1}='OUTPUT-PoresSegmentes-Ortho-Slice.Components.xlsx';   %Parametre modifiable
tablesXLSX.Frontieres{2}='OUTPUT-PoresSegmentes-Ortho-Slice2.Components.xlsx';  %Parametre modifiable
tablesXLSX.Frontieres{3}='OUTPUT-PoresSegmentes-Ortho-Slice3.Components.xlsx';  %Parametre modifiable
tablesXLSX.Frontieres{4}='OUTPUT-PoresSegmentes-Ortho-Slice4.Components.xlsx';  %Parametre modifiable
tablesXLSX.Frontieres{5}='OUTPUT-PoresSegmentes-Ortho-Slice5.Components.xlsx';  %Parametre modifiable
tablesXLSX.Frontieres{6}='OUTPUT-PoresSegmentes-Ortho-Slice6.Components.xlsx';  %Parametre modifiable






%--------------------------------------------------------------------------
%%    Script construisant le réseau : ne rien modifier a partir d'ici
%--------------------------------------------------------------------------

disp('Construction du reseau de pores extrait de l''image 3D');

%% Récupération des informations géométriques dans les tables .xlsx et les fichiers .mat
geometricalData=struct;
geometricalData.Liens=read_XLSX(tablesXLSX.Liens);
geometricalData.Pores=read_XLSX(tablesXLSX.Pores);
geometricalData.Frontieres=cell(1,nombreFrontieres);
for i=1:nombreFrontieres
    geometricalData.Frontieres{i}=read_XLSX(tablesXLSX.Frontieres{i});
end

foo=fieldnames(Avizo_pore);
Avizo_pore=Avizo_pore.(foo{1});
foo=fieldnames(Avizo_interface);
Avizo_interface=Avizo_interface.(foo{1});
for i=1:nombreFrontieres
    foo=fieldnames(Avizo_frontiere{i});
    Avizo_frontiere{i}=Avizo_frontiere{i}.(foo{1});
end

    


    
    %% Construction du reseau de pores
    
    disp('Construction du reseau de pores');
    tic;
    
    %Construction des tables de voisinage :
    %algorithme : recherche d'intersections entre les voxels des pores et les voxels des interfaces entre pores
    [interfaceToPore,i2p,p2i,parsedPore,orderPore,parsedInterface,orderInterface]=BuildConnectivityTables(poresImage,internalLinkImage);
    
    
    
    %Rescaling des infos geometriques
    poreVolume=voxelEdgeLength^3*double(porePropertyVolume);
    
    poreCenter=voxelEdgeLength*double(porePropertyCenter);
    internalLinkDiameter=voxelEdgeLength*double(internalLinkDiameter);
    internalLinkCenter=voxelEdgeLength*double(internalLinkCenter);
    
    for iBoundary=1:6
        boundaryLinkCenter{iBoundary}=voxelEdgeLength*double(boundaryLinkCenter{iBoundary});
        boundaryLinkPropertyDiameter{iBoundary}=voxelEdgeLength*double(boundaryLinkPropertyDiameter{iBoundary});
    end
    
    %Construction de la liste des liens internes
    [linksOwners,linksNeighbours,linkCenter,linkDiameter,raw_data_link]=ConstructInternalLinkList(interfaceToPore,i2p,internalLinkCenter,internalLinkDiameter,otherInternalLinkProperties);
    

    %Construction de la liste des liens frontiere
    [linksOwners,linksNeighbours,linkDiameter,linkCenter,raw_data_link,infos_liens_frontieres]=AddBoundaryLinks(linksOwners,linksNeighbours,poresImage,linkDiameter,linkCenter,boundaryLinkPropertyDiameter,boundaryLinkCenter,raw_data_link);

    
    %Renumerotation des liens (necessaire pour coder les liens frontieres dans PNM_Matlab)

    [boundaries,owners,neighbours,newOrder]=NetworkBuilder.RenumerotationLiensFrontieres(infos_liens_frontieres,linksOwners,linksNeighbours);

    pores=NetworkBuilder.BuildPoresToLinks(owners,neighbours,nPore);

    linkCenter=linkCenter(newOrder,:);
    linkDiameter=linkDiameter(newOrder);
    
    datanames=fieldnames(otherInternalLinkProperties);
    for iData=1:length(datanames)
        raw_data_link{iData}=raw_data_link{iData}(newOrder);
    end

    
    % Création du réseau de pores
    dimension=3;
    poreNetwork=PoreNetworkImageBased(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,[],voxelEdgeLength,materialImage,parsedPore,orderPore);

    
    % Ajout au reseau des informations geometrique sur les pores et liens
    
    poreNetwork.AddNewPoreData(poreVolume,'Volume')
    
    poreDiameter=(6*poreVolume/pi).^(1/3);
    poreNetwork.AddNewPoreData(poreDiameter,'Diameter')
    
    poreNetwork.AddNewLinkData(linkDiameter,'Diameter');

    names=fieldnames(otherPoreProperties);
    for iData=1:length(names)
        data=transpose(otherPoreProperties.(names{iData}));
        poreNetwork.AddNewPoreData(data,strcat('RawData_',names{iData}))
    end

    names=fieldnames(otherInternalLinkProperties);
    for iData=1:length(names)
        poreNetwork.AddNewLinkData(raw_data_link{iData},strcat('RawData_',names{iData}))
    end

    
    
    %% Fin de la creation du reseau !
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    disp(sprintf('Reseau de pores construit. Time spent : %d minutes %f s.',minutes,secondes));
    
    
    
    
    
    
    
    %% Utilities : fonctions utilisees  
     
    function [interfaceToPore,i2p,p2i,parsedPore,orderPore,parsedInterface,orderInterface]=BuildConnectivityTables(poresImage,internalLinkImage)
        
        disp('Debut de la construction des tables de connectivite');    
        nInterface=max(max(max(internalLinkImage)));
        nPores=max(max(max(poresImage)));

        nVoxelsInterface=numel(internalLinkImage);
        reshapedInterface=reshape(internalLinkImage,[1,nVoxelsInterface]);
        [sortedInterface,orderInterface]=sort(reshapedInterface);
        parsedInterface=parse_sorted_vector(sortedInterface);
        parsedInterface=parsedInterface(2:end);
        interface=cell(1,nInterface);
        for i=1:nInterface
            assert( parsedInterface{i}(1)==i )
            interface{i}= orderInterface(parsedInterface{i}(2):parsedInterface{i}(3));
        end

        nVoxelsPore=numel(internalLinkImage);
        assert(nVoxelsPore==nVoxelsInterface);
        reshapedPore=reshape(poresImage,[1,nVoxelsPore]);
        [sortedPore,orderPore]=sort(reshapedPore);
        parsedPore=parse_sorted_vector(sortedPore);
        parsedPore=parsedPore(2:end);
        pore=cell(1,nPores);
        for i=1:nPores
            assert( parsedPore{i}(1)==i )
            pore{i}=orderPore(parsedPore{i}(2):parsedPore{i}(3));
        end

        i2p=cell(1,nInterface);
        for j=1:nInterface
            i2p{j}=zeros(1,nPores);
            assert(reshapedInterface(interface{j}(1))==j)
            C= unique(reshapedPore(interface{j}));

            for i=C(C~=0)
                i2p{j}(i)=length(find(reshapedPore(interface{j})==i));
            end
        end

        p2i=cell(1,nPores);
        for j=1:nPores
            p2i{j}=zeros(1,nInterface);
            for i=1:nInterface
                p2i{j}(i)=i2p{i}(j);  
            end
        end

        interfaceToPore=cell(1,nInterface);
        for j=1:nInterface
            interfaceToPore{j}=find(i2p{j});
        end

        disp('Tables de connectivite construites');
    end
    
    



    function [linksOwners,linksNeighbours,linkCenter,linkDiameter,raw_data_link]=ConstructInternalLinkList(interfaceToPore,i2p,linkCenter,linkDiameter,otherLinkData)
                
        nInternalLink=length(interfaceToPore);

        linksOwners=zeros(1,nInternalLink);
        linksNeighbours=zeros(1,nInternalLink);
        badLink=zeros(1,nInternalLink);

        for iLink=1:nInternalLink

            %gestion des cas ou il n'y a pas exactement deux pores dans interfaceToPore{iLink}
            poreVoisin=interfaceToPore{iLink};
            nVoisin=length(poreVoisin);

            if nVoisin==2
                linksOwners(iLink)=poreVoisin(1);
                linksNeighbours(iLink)=poreVoisin(2);

            elseif nVoisin<2
                if nVoisin==0
                    disp('interface avec 0 pores voisins !!!');
                end
                badLink(iLink)=1; %suppression de ce lien

            elseif nVoisin>2
                %ne compter que les 2 voisins principaux (les autres sont des erreurs)

                [~,ordreDecroissant]=sort(i2p{iLink}(poreVoisin),'descend');
                linksOwners(iLink)=poreVoisin(ordreDecroissant(1));
                linksNeighbours(iLink)=poreVoisin(ordreDecroissant(2));
            end    

        end

        linksOwners=linksOwners(not(badLink));
        linksNeighbours=linksNeighbours(not(badLink));
        linkDiameter=linkDiameter(not(badLink));
        linkCenter=linkCenter(not(badLink),:);

        datanames=fieldnames(otherLinkData);
        raw_data_link=cell(1,length(datanames));
        for iDatas=1:length(datanames)
            foofoo=otherLinkData.(datanames{iDatas});
            raw_data_link{iDatas}=transpose(foofoo(not(badLink)));
        end
    end
    
    
    
    function [linksOwners,linksNeighbours,linkDiameter,linkCenter,raw_data_link,infos_liens_frontieres]=AddBoundaryLinks(linksOwners,linksNeighbours,poresImage,linkDiameter,linkCenter,boundaryLinkPropertyDiameter,boundaryLinkCenter,raw_data_link)
        
        %geometricalDataFrontieres=boundaryLinkProperties
        
        iCurrentLink=length(linksOwners)+1;
        infos_liens_frontieres=cell(1,6);

        for iFrontiere=1:6

            if iFrontiere==1
                frontiere=poresImage(1,:,:);
                boundaryLinkCenter{iFrontiere}=horzcat(voxelEdgeLength*ones(nPore,1),boundaryLinkCenter{iFrontiere});
            elseif iFrontiere==2
                frontiere=poresImage(end,:,:);
                boundaryLinkCenter{iFrontiere}=horzcat(voxelEdgeLength*imageSize(1)*ones(nPore,1),boundaryLinkCenter{iFrontiere});
            elseif iFrontiere==3
                frontiere=poresImage(:,1,:);
                boundaryLinkCenter{iFrontiere}=horzcat(boundaryLinkCenter{iFrontiere}(:,1),voxelEdgeLength*ones(nPore,1),boundaryLinkCenter{iFrontiere}(:,2));
            elseif iFrontiere==4
                frontiere=poresImage(:,end,:);
                boundaryLinkCenter{iFrontiere}=horzcat(boundaryLinkCenter{iFrontiere}(:,1),voxelEdgeLength*imageSize(2)*ones(nPore,1),boundaryLinkCenter{iFrontiere}(:,2));
            elseif iFrontiere==5
                frontiere=poresImage(:,:,1);  
                boundaryLinkCenter{iFrontiere}=horzcat(boundaryLinkCenter{iFrontiere},voxelEdgeLength*ones(nPore,1));
            elseif iFrontiere==6
                frontiere=poresImage(:,:,end);
                boundaryLinkCenter{iFrontiere}=horzcat(boundaryLinkCenter{iFrontiere},voxelEdgeLength*imageSize(3)*ones(nPore,1));
            end

            poreList=unique(frontiere);
            poreList=setdiff(poreList,0);
            nLinkFrontiere=length(poreList);
            infos_liens_frontieres{iFrontiere}=iCurrentLink:(iCurrentLink+nLinkFrontiere-1);

            linksOwners=[linksOwners,zeros(1,nLinkFrontiere)];
            linksNeighbours=[linksNeighbours,zeros(1,nLinkFrontiere)];
            linkDiameter=[linkDiameter,zeros(1,nLinkFrontiere)];
            linkCenter=vertcat(linkCenter,zeros(nLinkFrontiere,3));

            for i=1:length(raw_data_link)
                raw_data_link{i}=[raw_data_link{i},zeros(1,nLinkFrontiere)];
            end

            for iLinkFrontiere=1:nLinkFrontiere

                numPoreOwner=poreList(iLinkFrontiere);

                linksOwners(iCurrentLink)=numPoreOwner;
                linksNeighbours(iCurrentLink)=-1;

                linkDiameter(iCurrentLink)=boundaryLinkPropertyDiameter{iFrontiere}(numPoreOwner);

                linkCenter(iCurrentLink,:)=boundaryLinkCenter{iFrontiere}(numPoreOwner,:); 


                iCurrentLink=iCurrentLink+1;
            end
        end
    end

