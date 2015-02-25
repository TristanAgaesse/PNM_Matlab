function conductances = LocalScaleComputeConductancesDiffusion(network,diffusivity)
    %input : network, diffusivity
	%output : conductances

    
    %diff_O2_dans_N2 = 2e-5;
    
    %Get geometric data from the network
    CheckLinkDiameter(network) 

    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    linkSurface = network.GetLinkData('Surface');
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
    conductances(internalLinks)=diffusivity*linkSurface(internalLinks)./transpose(distance1(internalLinks)+distance2);
    conductances(boundaryLinks)=diffusivity*linkSurface(boundaryLinks)./transpose(2*distance1(boundaryLinks));
    
    
end


 
function CheckLinkDiameter(network)

    if not(isfield(network.GetLinkDataList,'Diameter'))
        diameter = network.ComputeAllLinkDiameter;
        network.AddNewLinkData(diameter,'Diameter');
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
