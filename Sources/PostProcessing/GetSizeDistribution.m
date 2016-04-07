function [ poreSizeDistribution, linkSizeDistribution,pd ] = GetSizeDistribution( network )
%GETSIZEDISTRIBUTION Summary of this function goes here
%   Input : network
%   Output : poreSizeDistribution, linkSizeDistribution
    
    figure
    porePhase=network.GetPoreData('Phase');
    for i=transpose(unique(porePhase))
        hold on
        poreDiameter = network.GetPoreData('Diameter');
        poreDiameter = poreDiameter(porePhase==i);
        poreSizeDistribution=histogram(poreDiameter);%,'Normalization','probability');
        %edges=linspace(min(poreDiameter),max(poreDiameter),20);
        %poreDiameterDistribution = discretize(poreDiameter,edges);
        
%         pd=fitdist(poreDiameter,'Weibull');
%         hold on
%         x=linspace(0,max(poreDiameter));
%         plot(x,pdf(pd,x))
    end
    figure
    
    linkDiameter = network.GetLinkData('Diameter');
    linkSizeDistribution=histogram(linkDiameter,'Normalization','probability');
%     edges=linspace(min(linkDiameter),max(linkDiameter),20);
%     linkDiameterDistribution = discretize(linkDiameter,edges);
    
    
end