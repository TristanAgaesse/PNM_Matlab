 function conductances = LocalScaleComputeConductancesDiffusion(poreNetwork)
    %input : poreNetwork,linkInlet,linkOutlet
    %output : conductances

    CheckLinkDiameter(poreNetwork) %V�rification si les diametres des liens sont d�j� calcul�s

    diff_O2_dans_N2 = 2e-5;

    nLink = poreNetwork.GetNumberOfLinks;
    nPore = poreNetwork.GetNumberOfPores;
    linkDiameter = poreNetwork.GetLinkDataList.Diameter;
    linkSurface = (linkDiameter.^2).*(pi/4);
    dimension = poreNetwork.Dimension;
    
    
    poreCenter=poreNetwork.GetPoreCenter(1:nPore);
    linkCenter=poreNetwork.GetLinkCenter(1:nLink);


    allLinks=1:nLink;
    internalLinks = poreNetwork.GetLinksFrontiere(0);
    boundaryLinks = GetLinksFrontiere(poreNetwork,1:poreNetwork.GetNumberOfBoundaries);
    
    a=poreCenter(poreNetwork.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(poreNetwork.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
    
    conductances = zeros(1,nLink);
    conductances(internalLinks)=diff_O2_dans_N2*linkSurface(internalLinks)./transpose(distance1(internalLinks)+distance2);
    conductances(boundaryLinks)=diff_O2_dans_N2*linkSurface(boundaryLinks)./transpose(2*distance1(boundaryLinks));
    
end

function CheckLinkDiameter(network)
    if not(isfield(network.GetLinkDataList,'Diameter'))
        disp('Calcul du diam�tre des liens...');
        tic;
        nLink = network.GetNumberOfLinks;
        diameter = zeros(1,nLink);
        for iLink = 1:nLink
            diameter(iLink) = network.ComputeLinkDiameter(iLink);
        end
        network.AddNewLinkData(diameter,'Diameter');
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        fprintf('Calcul du diam�tre des liens termin�. Dur�e : %d minutes %f s.',minutes,secondes);
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
