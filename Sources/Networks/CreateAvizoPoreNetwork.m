function poreNetwork = CreateAvizoPoreNetwork()
%CreateAvizoPoreNetwork Crée un reseau de pores dans Matlab 
%a partir d'une segmentation watershed calculée par Avizo sur une image 3D.
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


    %% Construction des tables de voisinage 
    %utilise des recherches dans les tableaux de voxels .mat 
    %donnant les  images 3D des pores et des interfaces entre pores

    disp('Construction des tables de connectivite');

    nInterface=max(max(max(Avizo_interface(1,:,:,:))));
    nPore=max(max(max(Avizo_pore(1,:,:,:))));

    nVoxelsInterface=numel(Avizo_interface);
    reshapedInterface=reshape(squeeze(uint32(Avizo_interface(1,:,:,:))),[1,nVoxelsInterface]);
    [sortedInterface,orderInterface]=sort(reshapedInterface);
    parsedInterface=parse_sorted_vector(sortedInterface);
    parsedInterface=parsedInterface(2:end);
    interface=cell(1,nInterface);
    for i=1:nInterface
        assert( parsedInterface{i}(1)==i )
        interface{i}= orderInterface(parsedInterface{i}(2):parsedInterface{i}(3));
    end

    nVoxelsPore=numel(Avizo_interface);
    assert(nVoxelsPore==nVoxelsInterface);
    reshapedPore=reshape(squeeze(uint32(Avizo_pore(1,:,:,:))),[1,nVoxelsPore]);
    [sortedPore,orderPore]=sort(reshapedPore);
    parsedPore=parse_sorted_vector(sortedPore);
    parsedPore=parsedPore(2:end);
    pore=cell(1,nPore);
    for i=1:nPore
        assert( parsedPore{i}(1)==i )
        pore{i}=orderPore(parsedPore{i}(2):parsedPore{i}(3));
    end



    i2p=cell(1,nInterface);
    for j=1:nInterface
        i2p{j}=zeros(1,nPore);
        assert(reshapedInterface(interface{j}(1))==j)
        C= unique(reshapedPore(interface{j}));

        for i=C(C~=0)
            i2p{j}(i)=length(find(reshapedPore(interface{j})==i));
        end
    end

    p2i=cell(1,nPore);
    for j=1:nPore
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

    %% Construction de la liste des liens internes

    nInternalLink=nInterface;

    linksOwners=zeros(1,nInternalLink);
    linksNeighbours=zeros(1,nInternalLink);
    linkArea=zeros(1,nInternalLink);
    linkCenter=zeros(nInternalLink,3);
    badLink=zeros(1,nInternalLink);


    for iLink=1:nInternalLink

        linkArea(iLink)=geometricalData.Liens.Area(iLink);
        linkCenter(iLink,1)=geometricalData.Liens.BaryCenterX(iLink);
        linkCenter(iLink,2)=geometricalData.Liens.BaryCenterY(iLink);
        linkCenter(iLink,3)=geometricalData.Liens.BaryCenterZ(iLink);

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
    linkArea=linkArea(not(badLink));
    linkCenter=linkCenter(not(badLink),:);

    names=fieldnames(geometricalData.Liens);
    raw_data_link=cell(1,length(names));
    for i=1:length(names)
        foo=geometricalData.Liens.(names{i});
        raw_data_link{i}=transpose(foo(not(badLink)));
    end

    %% Construction de la liste des liens frontiere

    iCurrentLink=length(linksOwners)+1;
    infos_liens_frontieres=cell(1,nombreFrontieres);

    for iFrontiere=1:nombreFrontieres

        poreList=unique(Avizo_frontiere{iFrontiere});
        poreList=setdiff(poreList,0);
        nLinkFrontiere=length(poreList);
        infos_liens_frontieres{iFrontiere}=iCurrentLink:(iCurrentLink+nLinkFrontiere-1);

        linksOwners=[linksOwners,zeros(1,nLinkFrontiere)];
        linksNeighbours=[linksNeighbours,zeros(1,nLinkFrontiere)];
        linkArea=[linkArea,zeros(1,nLinkFrontiere)];
        linkCenter=vertcat(linkCenter,zeros(nLinkFrontiere,3));
        
        for i=1:length(fieldnames(geometricalData.Liens))
            raw_data_link{i}=[raw_data_link{i},zeros(1,nLinkFrontiere)];
        end
        
        for iLinkFrontiere=1:nLinkFrontiere

            numPoreOwner=poreList(iLinkFrontiere);

            linksOwners(iCurrentLink)=numPoreOwner;
            linksNeighbours(iCurrentLink)=-1;
            
            %%TO DO :check formula
            linkArea(iCurrentLink)=length(Avizo_frontiere{iFrontiere}==numPoreOwner);%Modifier l'aire pour avoir les memes unites que pour les liens internes
            
            foo=0;
            linkCenter(iCurrentLink,1)=foo; %TO DO :check formula
            linkCenter(iCurrentLink,2)=foo;  %TO DO :check formula
            linkCenter(iCurrentLink,3)=foo;	%TO DO :check formula
            
    %        pores{numPoreOwner}=[pores{numPoreOwner},iCurrentLink];

            iCurrentLink=iCurrentLink+1;
        end

    end

    nLink=length(linksOwners);


    


    %% Construction des listes de pores
    poreCenter=zeros(nPore,3);
    for iPore=1:nPore
        poreCenter(iPore,1)=voxelEdgeLength*geometricalData.Pores.BaryCenterX(iPore); %TO DO : check formula
        poreCenter(iPore,2)=voxelEdgeLength*geometricalData.Pores.BaryCenterY(iPore); %TO DO : check formula
        poreCenter(iPore,3)=voxelEdgeLength*geometricalData.Pores.BaryCenterZ(iPore); %TO DO : check formula
    end
    
    poreVolume=zeros(1,nPore);
    voxelVolume=voxelEdgeLength^3;
    for iPore=1:nPore
        poreVolume(iPore)=voxelVolume*geometricalData.Pores.Volume3d(iPore); %TO DO : check formula
    end
    
    
    
    %% Renumerotation des liens

    [boundaries,owners,neighbours,newOrder]=NetworkBuilder.RenumerotationLiensFrontieres(infos_liens_frontieres,linksOwners,linksNeighbours);

    pores=NetworkBuilder.BuildPoresToLinks(owners,neighbours,nPore);

    linkCenter=linkCenter(newOrder,:);
    linkArea=linkArea(newOrder);


    %% Création du réseau de pores

    dimension=3;

    poreNetwork=PoreNetworkEuclidien(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter);

    
    
    poreNetwork.AddNewPoreData(poreVolume,'Volume')
    
    poreDiameter=(6*poreVolume/pi).^(1/3);
    poreNetwork.AddNewPoreData(poreDiameter,'Diameter')
    
    linkDiameter=sqrt(4*linkArea/pi);
    poreNetwork.AddNewLinkData(linkDiameter,'Diameter');

    names=fieldnames(geometricalData.Pores);
    for i=1:length(names)
        data=transpose(geometricalData.Pores.(names{i}));
        poreNetwork.AddNewPoreData(data,strcat('RawData_',names{i}))
    end

    names=fieldnames(geometricalData.Liens);
    for i=1:length(names)
        data=raw_data_link{i}(newOrder);
        poreNetwork.AddNewLinkData(data,strcat('RawData_',names{i}))
    end

end

