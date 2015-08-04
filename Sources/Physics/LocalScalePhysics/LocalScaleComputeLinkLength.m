function linkLength = LocalScaleComputeLinkLength(network)
    %input : network
    %output : conductances
    
    
    %Get geometric data from the network
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
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
    
    %Compute length
    linkLength = zeros(1,nLink);
    linkLength(internalLinks)=transpose(distance1(internalLinks)+distance2);
    linkLength(boundaryLinks)=2*distance1(boundaryLinks);
    
end       


function myNorm=FastNorm(myVect,dimension)
    %Vectorial version of the norm function
    if dimension==2
        myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2) ;
    elseif dimension==3
        myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2+myVect(:,3).^2) ;
    end
end        