function linkCassieBaxterContactAngle = LocalScaleComputeCassieBaxterContactAngle(network,pureContactAngle)
    %input : network
    %output : linkContactAngle
    
    
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
    linkContactAngle = zeros(1,nLink);

    
end       


function myNorm=FastNorm(myVect,dimension)
    %Vectorial version of the norm function
    if dimension==2
        myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2) ;
    elseif dimension==3
        myNorm=sqrt(myVect(:,1).^2+myVect(:,2).^2+myVect(:,3).^2) ;
    end
end        