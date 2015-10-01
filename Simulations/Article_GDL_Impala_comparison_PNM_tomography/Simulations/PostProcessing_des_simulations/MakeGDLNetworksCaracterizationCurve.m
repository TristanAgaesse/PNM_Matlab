function MakeGDLNetworksCaracterizationCurve(nNetwork,legende,networks)
%   Input : nNetwork : number of input networks
%           legende : legende des courbes, legende{i}=name of network_i
%           networks : networks{i}= a network. length(networks) = nNetwork


    function [nPore,nLink,poreDiamData,linkDiamData,poreVoisinsData]=ComputeCaracteristics(network)

        nPore = network.GetNumberOfPores;
        nLink = network.GetNumberOfLinks;

        poreDiameters=network.GetPoreData('Diameter');
        poreDiameters=poreDiameters(poreDiameters<2e-4);
        %[poreDiamBins,poreDiamBinCenters] = hist(poreDiameters,25);
        
        edges = transpose(linspace(0,2e-4,25));
        %edges = linspace(min(poreDiameters),max(poreDiameters),25);
        [poreDiamBinsOccurrences,bin] = histc(poreDiameters,edges);
        poreDiamBinCenters = (edges(1:end-1)+edges(2:end))/2;
        %volume=network.GetPoreData('Volume');
        volume=(4/3)*pi*(poreDiameters/2).^3;
        poreDiamBinsVolume=zeros(length(edges),1);
        for iPore=1:length(poreDiameters);
            poreDiamBinsVolume(bin(iPore)) = poreDiamBinsVolume(bin(iPore))+volume(iPore);
        end
        totalPoreVolume=sum(volume);
        poreDiamBinsVolume=poreDiamBinsVolume./totalPoreVolume;
        poreDiamBinsVolume=poreDiamBinsVolume(1:end-1);
        poreDiamBinsOccurrences=poreDiamBinsOccurrences(1:end-1);
        poreDiamData={poreDiamBinsOccurrences,poreDiamBinsVolume,poreDiamBinCenters};
        
        
        linkDiameters=2*network.GetLinkData('CapillaryRadius');
        linkDiameters=linkDiameters(linkDiameters<2e-4);
        [linkDiamBins,linkDiamBinCenters] = hist(linkDiameters,25);
        linkDiamData={linkDiamBins,linkDiamBinCenters};

        poreVoisins=zeros(1,nPore);
        for iPore=1:nPore
            poreVoisins(iPore)=length(network.GetLinksOfPore(iPore));
        end
        network.AddNewPoreData(poreVoisins,'PoresVoisins');
        poreVoisins=poreVoisins(poreVoisins<=25);
        [poreVoisinsBin,poreVoisinsBinCenters] =hist(poreVoisins,max(poreVoisins));
        poreVoisinsData={poreVoisinsBin,poreVoisinsBinCenters};
    end
    
%     legende={'Hmaxima8','LocalMax4','LocalMax10','LocalMax20','Hmaxima4'};
    
    
    assert(nNetwork==length(networks))
    
    
    for i=1:nNetwork
        [nPore,nLink,poreDiamData,linkDiamData,poreVoisinsData]=ComputeCaracteristics(networks{i});
        nPore_{i}                 = nPore;
        nLink_{i}                 = nLink;
        poreDiamBinsOccurrences_{i} = poreDiamData{1};
        poreDiamBinsVolume_{i} = poreDiamData{2};
        poreDiamBinCenters_{i}    = poreDiamData{3};
        linkDiamBins_{i}          = linkDiamData{1};
        linkDiamBinCenters_{i}    = linkDiamData{2};
        poreVoisinsBin_{i}        = poreVoisinsData{1};
        poreVoisinsBinCenters_{i} = poreVoisinsData{2};
    end
    
    %Plots
    
%     colors={'green',...
%         'blue',...
%         'cyan',...
%         'red',...
%         'magenta',...
%         'black',...
%         'yellow' };
    
    m=colormap(hsv(nNetwork));
    for i=1:nNetwork
        colors{i}=m(i,:);
    end
    
    fig=figure;hold on;
    for i=1:nNetwork
        %plot(poreDiamBinCenters_{i},poreDiamBins_{i}/nPore_{i},'Color',colors{i},'LineWidth',6);
        %plot(poreDiamBinCenters_{i},poreDiamBinsVolume_{i},'Color',colors{i},'LineWidth',6);
        volumeContribution=poreDiamBinsOccurrences_{i}.*(4/3)*pi.*(poreDiamBinCenters_{i}/2).^3;
        plot(poreDiamBinCenters_{i},volumeContribution,'Color',colors{i},'LineWidth',6);
    end
    title('Distribution of pore diameters')
    axe = get(fig, 'Children');
    legend(axe,legende)

    
    fig=figure;hold on;
    for i=1:nNetwork
        %plot(linkDiamBinCenters_{i},linkDiamBins_{i}/nLink_{i},'Color',colors{i},'LineWidth',6);
        plot(linkDiamBinCenters_{i},linkDiamBins_{i},'Color',colors{i},'LineWidth',6);
    end
    title('Distribution of link diameters')
    axe = get(fig, 'Children');
    legend(axe,legende)


    fig=figure;hold on;
    for i=1:nNetwork
        %plot(poreVoisinsBinCenters_{i},poreVoisinsBin_{i}/nPore_{i},'Color',colors{i},'LineWidth',6)
        plot(poreVoisinsBinCenters_{i},poreVoisinsBin_{i},'Color',colors{i},'LineWidth',6)
    end
    title('Distribution of coordination number')
    axe = get(fig, 'Children');
    legend(axe,legende)
    
    
    
    
    
    
%     
%     %Network 1
%     [nPore1,nLink1,poreDiamData1,linkDiamData1,poreVoisinsData1]=ComputeCaracteristics(network1);
%     poreDiamBins1          = poreDiamData1{1};
%     poreDiamBinCenters1    = poreDiamData1{2};
%     linkDiamBins1          = linkDiamData1{1};
%     linkDiamBinCenters1    = linkDiamData1{2};
%     poreVoisinsBin1        = poreVoisinsData1{1};
%     poreVoisinsBinCenters1 = poreVoisinsData1{2};
%     
%     %Network 2 
%     [nPore2,nLink2,poreDiamData2,linkDiamData2,poreVoisinsData2]=ComputeCaracteristics(network2);
%     poreDiamBins2          = poreDiamData2{1};
%     poreDiamBinCenters2    = poreDiamData2{2};
%     linkDiamBins2          = linkDiamData2{1};
%     linkDiamBinCenters2    = linkDiamData2{2};
%     poreVoisinsBin2        = poreVoisinsData2{1};
%     poreVoisinsBinCenters2 = poreVoisinsData2{2};
% 
%     %Network 3 
%     [nPore3,nLink3,poreDiamData3,linkDiamData3,poreVoisinsData3]=ComputeCaracteristics(network3);
%     poreDiamBins3          = poreDiamData3{1};
%     poreDiamBinCenters3    = poreDiamData3{2};
%     linkDiamBins3          = linkDiamData3{1};
%     linkDiamBinCenters3    = linkDiamData3{2};
%     poreVoisinsBin3        = poreVoisinsData3{1};
%     poreVoisinsBinCenters3 = poreVoisinsData3{2};    
%     
%     %Network 4 
%     [nPore4,nLink4,poreDiamData4,linkDiamData4,poreVoisinsData4]=ComputeCaracteristics(network4);
%     poreDiamBins4          = poreDiamData4{1};
%     poreDiamBinCenters4    = poreDiamData4{2};
%     linkDiamBins4          = linkDiamData4{1};
%     linkDiamBinCenters4    = linkDiamData4{2};
%     poreVoisinsBin4        = poreVoisinsData4{1};
%     poreVoisinsBinCenters4 = poreVoisinsData4{2};    
%     
%     %Network 5 
%     [nPore5,nLink5,poreDiamData5,linkDiamData5,poreVoisinsData5]=ComputeCaracteristics(network5);
%     poreDiamBins5          = poreDiamData5{1};
%     poreDiamBinCenters5    = poreDiamData5{2};
%     linkDiamBins5          = linkDiamData5{1};
%     linkDiamBinCenters5    = linkDiamData5{2};
%     poreVoisinsBin5        = poreVoisinsData5{1};
%     poreVoisinsBinCenters5 = poreVoisinsData5{2};  
%     
%     
%     
%     %Plots
%     
%     fig=figure;hold on;
%     plot(poreDiamBinCenters1,poreDiamBins1/nPore1)
%     plot(poreDiamBinCenters2,poreDiamBins2/nPore2,'red')
%     plot(poreDiamBinCenters3,poreDiamBins3/nPore3,'green')
%     plot(poreDiamBinCenters4,poreDiamBins4/nPore4,'black')
%     plot(poreDiamBinCenters5,poreDiamBins5/nPore5,'magenta')
%     title('Pore diameters')
%     axe = get(fig, 'Children');
%     legend(axe,leg)
% 
%     fig=figure;hold on;
%     plot(linkDiamBinCenters1,linkDiamBins1/nLink1)
%     plot(linkDiamBinCenters2,linkDiamBins2/nLink2,'red')
%     plot(linkDiamBinCenters3,linkDiamBins3/nLink3,'green')
%     plot(linkDiamBinCenters4,linkDiamBins4/nLink4,'black')
%     plot(linkDiamBinCenters5,linkDiamBins5/nLink5,'magenta')
%     title('Link diameters')
%     axe = get(fig, 'Children');
%     legend(axe,leg)
% 
% 
%     fig=figure;hold on;
%     plot(poreVoisinsBinCenters1,poreVoisinsBin1/nPore1)
%     plot(poreVoisinsBinCenters2,poreVoisinsBin2/nPore2,'red')
%     plot(poreVoisinsBinCenters3,poreVoisinsBin3/nPore3,'green')
%     plot(poreVoisinsBinCenters4,poreVoisinsBin4/nPore4,'black')
%     plot(poreVoisinsBinCenters5,poreVoisinsBin5/nPore5,'magenta')
%     title('Coordinence')
%     axe = get(fig, 'Children');
%     legend(axe,leg)
    
end