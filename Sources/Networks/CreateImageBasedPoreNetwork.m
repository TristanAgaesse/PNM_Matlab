function poreNetwork = CreateImageBasedPoreNetwork(inputContainerMap)
%CreateImageBasedPoreNetwork 
%Crée un reseau de pores dans Matlab a partir du resultat d'une segmentation watershed des pores d'une image 3D.
%
% Input : - inputContainerMap : objet de type container.Map contenant les input. 
%
%   Description de inputContainerMap :
%
%       Proprietes de l'image 3D du materiau
%       - inputContainerMap('MaterialImage') : image 3D du materiau (tableau de voxels)
%       - inputContainerMap('VoxelEdgeLength') : taille d'un voxel en mètres
%
%       Proprietes des pores (pores = espaces vides séparés par des watershed lines)
%       - inputContainerMap('PoresImage') : image des pores labelises
%       - inputContainerMap('PorePropertyVolume') : tableau nPore contenant le volume de chaque pore
%       - inputContainerMap('PorePropertyCenter') : tableau (nPore,3) contenant le barycentre de chaque pore
%
%       Proprietes des internal links (internal link = interface entre deux pores : watershed lines)
%   	- inputContainerMap('InternalLinkImage') : image des internal links labelises puis dilates d'un pixel
%       - inputContainerMap('InterfaceToPore') 
%       - inputContainerMap('InternalLinkCapillaryRadius') : tableau contenant le rayon déduit de la distance map de chaque internal links
%       - inputContainerMap('InternalLinkPropertyCenter') : tableau (nInternalLink,3) contenant le barycentre de chaque internal links
%
%       Proprietes des liens frontieres (liens frontieres = slices des pores sur les bords de l'image)
%       - inputContainerMap('BoundaryLinkPropertyCenter') : cell (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) contenant les centres des liens frontieres (tableaux (nPore,3) avec NaN si le pore i n'intersecte pas la slice)
%       - inputContainerMap('BoundaryLinkPropertyCapillaryRadius') : cell (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) contenant les diametres des liens frontieres (tableaux (1,nPore) avec NaN si le pore i n'intersecte pas la slice)
%
%       (Optionnel) Pour mettre d'autre infos geometriques scalaires dans le reseau :
%       - inputContainerMap('OtherPoreProperties') : struct contenant des tableaux (1,nPore)
%       - inputContainerMap('OtherInternalLinkProperties') : struct contenant des tableaux (1,nInternalLink)
%       - inputContainerMap('OtherBoundaryLinkProperties') : struct contenant des cells (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) qui contiennent des tableaux (1,nPore)
%
% Output : - poreNetwork : reseau de pores PoreNetworkImageBased
    
    
    
    
    %% Lecture et controle des donnes d'entree
    
    voxelEdgeLength=inputContainerMap('VoxelEdgeLength');
    assert(length(voxelEdgeLength)==1,'VoxelEdgeLength')
    
    materialImage=inputContainerMap('MaterialImage');
    imageSize=size(materialImage);
    
        
    poresImage=inputContainerMap('PoresImage');
    assert(isequal(imageSize,size(poresImage)),'wrong size of PoresImages')
    nPore=max(max(max(poresImage)));
        
    porePropertyVolume=inputContainerMap('PorePropertyVolume');
    assert(length(porePropertyVolume)==nPore,'length(PorePropertyVolume) must equal nPore')
    
    porePropertyCenter=inputContainerMap('PorePropertyCenter');
    assert(isequal(size(porePropertyCenter),[nPore,3]),'size(PorePropertyCenter) must equal [nPore,3]')
        

    interfaceToPore=inputContainerMap('InterfaceToPore');
       

    internalLinkCapillaryRadius=inputContainerMap('InternalLinkCapillaryRadius');
    nInterfaceLink=length(internalLinkCapillaryRadius);
    
    internalLinkCenter=inputContainerMap('InternalLinkPropertyCenter');
    assert(isequal(size(internalLinkCenter),[nInterfaceLink,3]),'size(InternalLinkPropertyCenter) must equal [nInterfaceLink,3]')
       
    
    boundaryLinkCenter=inputContainerMap('BoundaryLinkPropertyCenter');
    assert( iscell(boundaryLinkCenter) && length(boundaryLinkCenter)==6 , 'wrong size of BoundaryLinkPropertyCenter')
    assert( isequal(size(boundaryLinkCenter{1}),[nPore, 2]), 'wrong size of BoundaryLinkPropertyCenter')
    
    boundaryLinkPropertyCapillaryRadius=inputContainerMap('BoundaryLinkPropertyCapillaryRadius');
    assert( iscell(boundaryLinkPropertyCapillaryRadius) && length(boundaryLinkPropertyCapillaryRadius)==6 , 'wrong size of BoundaryLinkPropertyCapillaryRadius')
    assert( length(boundaryLinkPropertyCapillaryRadius{1})==nPore , 'wrong size of BoundaryLinkPropertyCapillaryRadius')
    
    
    otherPoreProperties=struct;
    if inputContainerMap.isKey('OtherPoreProperties')
        otherPoreProperties=inputContainerMap('OtherPoreProperties');
        assert(isstruct(otherPoreProperties),'inputContainerMap(''OtherPoreProperties'') must be a struct')
    end
    
    otherInternalLinkProperties=struct;
    if inputContainerMap.isKey('OtherInternalLinkProperties')
        otherInternalLinkProperties=inputContainerMap('OtherInternalLinkProperties');
        assert(isstruct(otherInternalLinkProperties),'inputContainerMap(''OtherInternalLinkProperties'') must be a struct')
    end

    otherBoundaryLinkProperties=struct;
    if inputContainerMap.isKey('OtherBoundaryLinkProperties')
        otherBoundaryLinkProperties=inputContainerMap('OtherBoundaryLinkProperties');
        assert(isstruct(otherBoundaryLinkProperties),'inputContainerMap(''OtherBoundaryLinkProperties'') must be a struct')
    end
    
    
    
    
    
    %% Construction du reseau de pores
    
    disp('Construction du reseau de pores');
    tic;
    
    %%Construction des tables de voisinage :
    %%algorithme : recherche d'intersections entre les voxels des pores et les voxels des interfaces entre pores
    %interfaceToPore=BuildConnectivityTables(poresImage,internalLinkImage);
    
    
    
    %Rescaling des infos geometriques
    poreVolume=voxelEdgeLength^3*double(porePropertyVolume);
    
    poreCenter=voxelEdgeLength*double(porePropertyCenter);
    internalLinkCapillaryRadius=voxelEdgeLength*double(internalLinkCapillaryRadius);
    internalLinkCenter=voxelEdgeLength*double(internalLinkCenter);
    
    for iBoundary=1:6
        boundaryLinkCenter{iBoundary}=voxelEdgeLength*double(boundaryLinkCenter{iBoundary});
        boundaryLinkPropertyCapillaryRadius{iBoundary}=voxelEdgeLength*double(boundaryLinkCapillaryRadius{iBoundary});
    end
    
    %Construction de la liste des liens internes
    [linksOwners,linksNeighbours,linkCenter,linkDiameter,raw_data_link]=ConstructInternalLinkList(interfaceToPore,internalLinkCenter,internalLinkCapillaryRadius,otherInternalLinkProperties);
    

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
    poreNetwork=PoreNetworkImageBased(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,[],voxelEdgeLength,materialImage,poresImage);
    
    
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
    fprintf('Reseau de pores construit. Time spent : %d minutes %f s.',minutes,secondes);
    
    
    

    
    %% Utilities : fonctions utilisees  
     
    function [linksOwners,linksNeighbours,linkCenter,linkDiameter,raw_data_link]=ConstructInternalLinkList(interfaceToPore,linkCenter,linkDiameter,otherLinkData)
                
        linksOwners=transpose(interfaceToPore(:,1));
        linksNeighbours=transpose(interfaceToPore(:,2));

        datanames=fieldnames(otherLinkData);
        raw_data_link=cell(1,length(datanames));
        for iDatas=1:length(datanames)
            foofoo=otherLinkData.(datanames{iDatas});
            raw_data_link{iDatas}=transpose(foofoo(not(badLink)));
        end
    end
    
    
    
    function [linksOwners,linksNeighbours,linkDiameter,linkCenter,raw_data_link,infos_liens_frontieres]=AddBoundaryLinks(linksOwners,linksNeighbours,poresImage,linkDiameter,linkCenter,boundaryLinkPropertyDiameter,boundaryLinkCenter,raw_data_link)
        
        
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


end

