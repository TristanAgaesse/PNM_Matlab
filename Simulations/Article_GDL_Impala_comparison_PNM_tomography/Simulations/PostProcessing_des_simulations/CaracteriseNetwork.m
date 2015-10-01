

%Pore Diameter 
figure
hist(network.GetPoreData('Diameter'),100)
title(strcat('Pore Diameter ',titleParameter));

%Constriction Diameter 
figure
hist(2*network.GetLinkData('CapillaryRadius'),100)
title(strcat('Constriction Diameter ',titleParameter));

%Pore number of neighboor
figure
poreVoisins=zeros(1,network.GetNumberOfPores);
for iPore=1:network.GetNumberOfPores
    poreVoisins(iPore)=length(network.GetLinksOfPore(iPore));
end
hist(poreVoisins,max(poreVoisins))
title(strcat('Pore number of neighboor ',titleParameter));

%Correlation between links size and pore size
figure
poreDiameter = network.GetPoreData('Diameter');
linkDiameter = 2*network.GetLinkData('CapillaryRadius');
internalLinks = network.GetLinksFrontiere(0);
sizeCorrelation= zeros(1,2*length(internalLinks));
j=1;
for iInternalLink = internalLinks
    
    p=network.GetPoresOfLink(iInternalLink);
    sizeCorrelation(j)= poreDiameter(p(1))/linkDiameter(iInternalLink);
    sizeCorrelation(j+1)= poreDiameter(p(2))/linkDiameter(iInternalLink);
    j=j+2;
end
hist(sizeCorrelation,100)
title(strcat('PoreDiameter/LinkDiameter ',titleParameter));


%Taille de lien directionnelle
figure
linkDiameter = 2*network.GetLinkData('CapillaryRadius');
internalLinks = network.GetLinksFrontiere(0);
linkDirectionalSize= zeros(length(internalLinks),2);
j=1;
for iInternalLink = internalLinks
    
    p=network.GetPoresOfLink(iInternalLink);
    
    linkDirection=network.GetPoreCenter(p(1))-network.GetPoreCenter(p(2));
    linkDirection=linkDirection./norm(linkDirection);
    
    linkDirectionalSize(j,1)= abs(linkDirection(1));
    linkDirectionalSize(j,2)= linkDiameter(iInternalLink);
    j=j+1;
end


edges = linspace(min(linkDirectionalSize(:,1)),max(linkDirectionalSize(:,1)),11);
[N,bin] = histc(linkDirectionalSize(:,1),edges);
binCenter = (edges(1:end-1)+edges(2:end))/2;

myMean=zeros(length(edges),1);
for j=1:length(linkDirectionalSize(:,1))
    myMean(bin(j))=myMean(bin(j))+linkDirectionalSize(j,2);
end
myMean = myMean(1:end-1)./N(1:end-1);

plot(binCenter,myMean)
title(strcat('Link directional diameter ',titleParameter));


