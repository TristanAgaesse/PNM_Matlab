function conductances = LocalScaleComputeConductancesStokes(network)
    %input : network
    %output : conductances
    
    
    viscosite_dyn_water = 1e-3;

    
    %Get geometric data from the network
    CheckLinkDiameter(network) %Check if link diameters are already computed

    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    linkDiameter = network.GetLinkDataList.Diameter;
    linkSurface = (linkDiameter.^2).*(pi/4);
    dimension = network.Dimension;
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);

    allLinks=1:nLink;
    internalLinks = network.GetLinksFrontiere(0);
    boundaryLinks = GetLinksFrontiere(network,1:network.GetNumberOfBoundaries);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
    
    %Compute conductances
    conductances = zeros(1,nLink);
    conductances(internalLinks)=linkSurface(internalLinks)./(8*pi*viscosite_dyn_water*transpose(distance1(internalLinks)+distance2));
    conductances(boundaryLinks)=diff_O2_dans_N2*linkSurface(boundaryLinks)./(8*pi*viscosite_dyn_water*transpose(2*distance1(boundaryLinks)));
    
end       


function CheckLinkDiameter(network)
    if not(isfield(network.GetLinkDataList,'Diameter'))
        disp('Calcul du diametre des liens...');
        tic;
        nLink = network.GetNumberOfLinks;
        diameter = zeros(1,nLink);
        for iLink = 1:nLink
            diameter(iLink) = network.ComputeLinkDiameter(iLink);
        end
        network.AddNewLinkData(diameter,'Diameter');
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        fprintf('Calcul du diametre des liens termin�. Dur�e : %d minutes %f s.',minutes,secondes);
    end
end


function myNorm=FastNorm(myVect,dimension)
    %Vectorial version of the norm function
    if dimension==2
        myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2) ;
    elseif dimension==3
        myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2+myVect(:,3).^2) ;
    end
end        