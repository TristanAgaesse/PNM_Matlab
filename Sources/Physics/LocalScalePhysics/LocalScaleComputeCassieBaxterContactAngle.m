function linkCassieBaxterContactAngle = LocalScaleComputeCassieBaxterContactAngle(network,pureContactAngle)
    %input : network,pureContactAngle
    %output : linkCassieBaxterContactAngle
    
    linkPhases=network.GetLinkData('RawData_NeighborPhases');
    nPhase = size(linkPhases,2);
    
    assert(length(pureContactAngle)==nPhase);
    
    nLink = network.GetNumberOfLinks;
    
    %Compute phase proportional coefficients
    foosum=sum(linkPhases,2);
    validLinks=foosum>0;
    assert(length(validLinks)==sum(validLinks))
    for iPhase=1:nPhase
        linkPhases(validLinks,iPhase)=linkPhases(validLinks,iPhase)./foosum(validLinks);
    end
    
    meanCos = zeros(nLink,1);
    for iPhase=1:nPhase
        meanCos = meanCos+linkPhases(:,iPhase)*pureContactAngle(iPhase);
    end
    
    linkCassieBaxterContactAngle = acos(meanCos);
    
end       

