function [ totalSaturation, saturationProfile ] = ComputeSaturation( cluster, poreNetwork, options )
%COMPUTESATURATION Calcule la saturation totale ou la courbe de saturation
%d'un reseau de pore envahi par un cluster.
% Input : -cluster
%         -network
%         -options : struct with the following (facultative) fields :
%        - options.type = 'saturationProfile' ou 'totalSaturation'
%        - options.axe = [x,y,z] axe (axe selon lequel faire les tranches de saturation)
%        - options.nPointProfile = number of slices for profile along the axe 
%        - options.codeForLiquid = if using image, code for the liquid phase
%        - options.codeForSolid = if using image, code for the solid phase
%
% Output : [ totalSaturation, saturationProfile ]


    %Checking initial state
    [mode,axe,nPointCurve,codeForLiquid,codeForSolid]=ReadCheckInputs( cluster, poreNetwork, options );
    CheckPoreVolume(poreNetwork);
    
    
    %Calcul de la saturation totale
    volumePore=poreNetwork.GetPoreDataList.('Volume');
    totalVolume=sum(volumePore);
    invadedVolume=sum(volumePore(cluster.GetInvadedPores));
    totalSaturation=invadedVolume/totalVolume;
    if strcmp(mode,'totalSaturation')
        saturationProfile=[];
        return
    end
    
    %Calcul de la courbe de saturation
    if isa(network,'PoreNetworkMesh') || isa(network,'PoreNetworkMeshFibrous') 
        saturationProfile=ComputeSaturationProfileMesh(network,cluster,nPointCurve,axe);
    
    elseif isa(network,'PoreNetworkImageBased')
        network.AddNewPoreData(cluster.GetInvadedPoresBooleans,'fooSatCalc_invadedPores')
        image=network.GetImagePoreData('fooSatCalc_invadedPores');
        
        [~,foo]=max(abs(axe));
        axe=[0 0 0];
        axe(foo(1))=1;
        
        [~, saturationProfile]=ComputeSaturationOnImage(image,nPointCurve,axe,codeForLiquid,codeForSolid);
    end

end


function saturationProfile=ComputeSaturationProfileMesh(network,cluster,nPointCurve,axe)
    %Calcul de Cmin et Cmax des sommets pour chaque cellule (C pour
    %coordinate suivant l'axe=axialCoordinates)

    poreVolume=network.GetPoreDataList.('Volume');
    totalPoreVolume=sum(poreVolume);

    axialCoordinates=network.GetAllVerticesCoordinates*transpose(axe);
    [cBas,indexMin]=min(axialCoordinates);
    cHaut=max(axialCoordinates);

    nPore=network.GetNumberOfPores;
    CmaxCmin=zeros(2,nPore);
    for iPore=1:nPore
        foo=axialCoordinates(network.GetVerticesOfPoreNumber(iPore));
        Cmax=max(foo);
        Cmin=min(foo);
        CmaxCmin(:,iPore)=[Cmax Cmin];
    end

    %Pour chaque point de la courbe, rep�rer les cellules situ�es totalement au
    %dessus ou au dessous du plan de coupe, et pour les cellules
    %intersectant le plan, les cliper.
    increasingInvadedVolume=zeros(1,nPointCurve);
    saturationProfile=zeros(nPointCurve,2);
    for iPointCurve=1:nPointCurve

        CPointCurve=cBas+iPointCurve/nPointCurve*(cHaut-cBas);
        saturationProfile(iPointCurve,1)=CPointCurve;

        invadedVolumeBeneath=0;
        P0=network.GetVertice(indexMin)+(CPointCurve-cBas)*axe  ;


        if network.Dimension==3;
            clippingPlane=createPlane(P0,axe);
        else
            line=createLine(P0(1), P0(2), axe(1), axe(2));
            clippingLine=orthogonalLine(line, P0);
        end



        for iPore=cluster.GetInvadedPores
            signe=sign(CmaxCmin(:,iPore)-CPointCurve);
            if signe(1)<=0 %pore en dessous 
                invadedVolumeBeneath=invadedVolumeBeneath+poreVolume(iPore);

            elseif signe(2)==1 %pore au dessus 

            else %pore intersect�
                assert(signe(1)==1 && signe(2)==-1)
                NODES=network.GetVerticesOfPore(iPore);
                %centrePore=mean(NODES);


                if network.Dimension==3;
                    [~,volumeSansFibres]=convhulln(NODES);

                    %links=poreNetwork.GetLinksOfPore(iPore);
                    %poreVertices=poreNetwork.GetVerticesOfPoreNumber(iPore);
    %                 nFace=length(links);
    %                 FACES=cell(1,nFace);
    %                 for iFace=1:nFace                   
    %                     linkVertices=poreNetwork.GetVerticesOfLinkNumber(links(iFace));
    %                     [~,~,ic]=intersect(linkVertices,poreVertices);
    %                     FACES{iFace}=ic;
    %                      %Ordonner les sommets pour que la normale pointe à
    %                     %l'extérieur.
    %                     foo=cross((NODES(ic(2),:)-NODES(ic(1),:)),(centrePore-NODES(ic(1),:)));
    %                     if dot(foo,(NODES(ic(3),:)-NODES(ic(2),:)))>0
    %                         %inverser l'ordre
    %                         FACES{iFace}=ic(length(ic)+1-(1:length(ic)));
    %                     end
    %                 end
                    FACES = minConvexHull(NODES);
                    [NODES2, ~] = clipConvexPolyhedronHP(NODES, FACES, clippingPlane);

                    if size(NODES2,1)<4
                        volumeBeneathSansFibres=0;
                    else
                        [~,volumeBeneathSansFibres]=convhulln(NODES2);
                    end
                else
                    volumeSansFibres=abs(polygonArea(NODES));
                    POLY2 = clipPolygonHP(NODES, clippingLine);
                    volumeBeneathSansFibres=abs(polygonArea(POLY2));
                end

                %approximation sur le volume clipper pour tenir compte des fibres : on suppose
                %que les fibres sont r�parties uniform�ment sur les
                %ar�tes au dessus et en dessous

                vol=poreVolume(iPore)*volumeBeneathSansFibres/volumeSansFibres;
                invadedVolumeBeneath=invadedVolumeBeneath+vol; 
            end
        end
        increasingInvadedVolume(iPointCurve)=invadedVolumeBeneath;
    end

    volumeTranche=totalPoreVolume/nPointCurve;%Approximation valable si le r�seau est un pav� � fronti�res planes 
    for iPointCurve=1:nPointCurve
        if iPointCurve>1
            saturationProfile(iPointCurve,2)=(increasingInvadedVolume(iPointCurve)-increasingInvadedVolume(iPointCurve-1))/volumeTranche;
        else
            saturationProfile(1,2)=increasingInvadedVolume(1)/volumeTranche;
        end
    end


    
end


function [mode,axe,nPointCurve,codeForLiquid,codeForSolid]=ReadCheckInputs( cluster, poreNetwork, options )
    %Checking and reading inputs

    assert(isa(cluster,'ClusterMonophasique'),'First input must be a ClusterMonophasique')
    assert(isa(poreNetwork,'PoreNetwork'),'First input must be a PoreNetwork')

    if isfield(options,'type') 
        if strcmp(options.type,'saturationProfile')
            mode='saturationProfile';

            if isfield(options,'axe')
                assert(ismatrix(options.axe) && length(options.axe)==poreNetwork.GetDimension,'wrong axe input')
                axe=options.axe;
                axe=axe/norm(axe);
            else
                disp('Axe non d�fini pour tracer la courbe de saturation. Choix par d�faut');
                if poreNetwork.Dimension==3;
                    axe=[0 0 1];
                else
                    axe=[0 1];
                end
            end

            if isfield(options,'nPointProfile')
                nPointCurve=options.nPointProfile;
            else
                nPointCurve=20;
            end

        end
    else
        mode='totalSaturation';
    end
    
    if isa(poreNetwork,'PoreNetworkImageBased')
        assert(isfield(options,'codeForLiquid') && isfield(options,'codeForSolid'),'Codes for liquid and solid required when working with images')
        codeForLiquid=options.codeForLiquid;
        codeForSolid=options.codeForSolid;
    
    else
        codeForLiquid='foo';
        codeForSolid='foo';
    end

end
    
function CheckPoreVolume(poreNetwork)
    %V�rification que les volumes des pores sont d�j� calcul�s
    if not(isfield(poreNetwork.GetPoreDataList,'Volume'))
        disp('Calcul des volumes des pores...');
        tic;
        volumePore=poreNetwork.ComputeAllPoreVolume;
        poreNetwork.AddNewPoreData(volumePore,'Volume');
        duree=toc;minutes=floor(duree/60);secondes=duree-60*minutes;
        fprintf('Calcul des volumes des pores termin�. Dur�e : %d minutes %f s.',minutes,secondes);
    end
end
