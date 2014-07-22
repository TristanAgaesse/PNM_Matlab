function conductances = LocalScaleComputeConductancesStokes(poreNetwork,linkInlet,linkOutlet)
    %input : poreNetwork,linkInlet,linkOutlet
    %output : conductances
    viscosite_dyn_water = 1e-3;

    %V�rification si les diametres des liens sont d�j� calcul�s
    if not(isfield(poreNetwork.GetLinkDataList,'Diameter'))
        disp('Calcul du diam�tre des liens...');
        tic;
        nLink = poreNetwork.GetNumberOfLinks;
        diameter = zeros(1,nLink);
       for iLink = 1:nLink
            diameter(iLink) = poreNetwork.ComputeLinkDiameter(iLink);
       end
       poreNetwork.AddNewLinkData(diameter,'Diameter');
       duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
       disp(sprintf('Calcul du diam�tre des liens termin�. Dur�e : %d minutes %f s.',minutes,secondes));
    end

    nLink = poreNetwork.GetNumberOfLinks;
    nPore = poreNetwork.GetNumberOfPores;
    linkDiameter = poreNetwork.GetLinkDataList.Diameter;
    linkSurface = (linkDiameter.^2).*(pi/4);
    dimension = poreNetwork.Dimension;

    poreCenter=poreNetwork.GetPoreCenter(1:nPore);

    linkCenter=poreNetwork.GetLinkCenter(1:nLink);
    

    conductances = zeros(1,nLink);

    internalLinks = setdiff(poreNetwork.GetLinksFrontiere(0),union(linkInlet,linkOutlet));

    
    a=poreCenter(poreNetwork.LinkOwners(internalLinks),:)-linkCenter(internalLinks,:);
    b=poreCenter(poreNetwork.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    surface=linkSurface(internalLinks);
    if dimension==2
        distance1=sqrt(a(:,1).^2+a(:,2).^2) ;
        distance2=sqrt(b(:,1).^2+b(:,2).^2) ;
    elseif dimension==3
        distance1=sqrt(a(:,1).^2+a(:,2).^2+a(:,3).^2) ;
        distance2=sqrt(b(:,1).^2+b(:,2).^2+b(:,3).^2) ;
    end
    conductances(internalLinks)=surface./(8*pi*viscosite_dyn_water*transpose(distance1+distance2));


    for iLink = linkOutlet
        surface = linkSurface(iLink);
        distance = 2*norm(poreCenter(poreNetwork.LinkOwners(iLink),:)-linkCenter(iLink,:));
        conductances(iLink) = surface^2/(8*pi*viscosite_dyn_water*(distance));
    end

    for iLink = linkInlet
        surface = linkSurface(iLink);
        distance = 2*norm(poreCenter(poreNetwork.LinkOwners(iLink),:)-linkCenter(iLink,:));
        conductances(iLink) = surface^2/(8*pi*viscosite_dyn_water*(distance));
    end

end       

        