function network = SetUpNetwork()
%SETUPNETWORK Create network with boundary conditions for simulations of inhomogeneous GDL for Yohann Thomas 


    load('GDL_Image/GDLwithoutMembranes_pnmExtract_h8.mat')
    InputScriptImageBasedPNM
    network = CreateImageBasedPoreNetwork(inputContainerMap);

    
    membranePoreLabel =  FindMembranePoreLabel(network);
    network.AddNewPoreData( membranePoreLabel,'MembranePoreLabel'); 
    
    
    [gdlBottomLinks,gdlTopLinks,gdlPores] = FindGDLElements(network);
    
    gdlBottomLinksData=zeros(1,network.GetNumberOfLinks);
    gdlBottomLinksData(gdlBottomLinks)=1;
    network.AddNewLinkData(gdlBottomLinksData,'GdlBottomLinks');
    
    gdlTopLinksData=zeros(1,network.GetNumberOfLinks);
    gdlTopLinksData(gdlTopLinks)=1;
    network.AddNewLinkData(gdlTopLinksData,'GdlTopLinks');
    
    gdlPoreData=zeros(1,network.GetNumberOfPores);
    gdlPoreData(gdlPores)=1;
    network.AddNewPoreData(gdlPoreData,'GdlPores');
        
end


function membranePoreLabel = FindMembranePoreLabel(network)


    PsiFusionImage=ReadTiff('GDL_Image/PSI_FusionImages_2540.tif');
    [fusionlabelEnds,fusionOrderLabels,fusionLabelIndices] = PoreNetworkImageBased.ParseLabeledImage(PsiFusionImage);
    
%     holder=PoreNetworkImageBased.GetVoxelsOfLabel(200,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
%     fibers=PoreNetworkImageBased.GetVoxelsOfLabel(50,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    topMembrane=PoreNetworkImageBased.GetVoxelsOfLabel(180,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    bottomMembrane=PoreNetworkImageBased.GetVoxelsOfLabel(75,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    
    void=PoreNetworkImageBased.GetVoxelsOfLabel(250,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    void1=PoreNetworkImageBased.GetVoxelsOfLabel(114,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    void2=PoreNetworkImageBased.GetVoxelsOfLabel(122,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    void3=PoreNetworkImageBased.GetVoxelsOfLabel(128,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    void4=PoreNetworkImageBased.GetVoxelsOfLabel(139,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    void5=PoreNetworkImageBased.GetVoxelsOfLabel(153,fusionlabelEnds,fusionOrderLabels,fusionLabelIndices);
    
    network.AddNewPoreData( 1:network.GetNumberOfPores,'PoreNumber');
    imagePores = network.GetImagePoreData('PoreNumber');
    
    gdlPores = unique(imagePores([void,void1,void2,void3,void4,void5]));
    gdlPores=gdlPores(gdlPores>0);
    
    topPores = unique(imagePores(topMembrane));
    topPores=topPores(topPores>0);
    bottomPores = unique(imagePores(bottomMembrane));
    bottomPores=bottomPores(bottomPores>0);
    
    membranePoreLabel = zeros(1,network.GetNumberOfPores);
    membranePoreLabel(topPores)=1;
    membranePoreLabel(bottomPores)=2;
    membranePoreLabel(gdlPores)=3;
end



function [gdlBottomLinks,gdlTopLinks,gdlPores] = FindGDLElements(network)

    membraneCodes=network.GetPoreData('MembranePoreLabel'); 
    %membraneCodes=1 si membrane sup, 2 pour membrane inf, 3 si dans GDL
    
    gdlPores = find(membraneCodes==3);
    membraneSupPores = find(membraneCodes==1);
    membraneInfPores = find(membraneCodes==2);
    
    gdlBoundary = network.GetPoreRegionBoundaryLinks(gdlPores);
    membraneSupBoundary = network.GetPoreRegionBoundaryLinks(membraneSupPores);
    membraneInfBoundary = network.GetPoreRegionBoundaryLinks(membraneInfPores);
    
    gdlTopLinks=intersect(gdlBoundary,membraneSupBoundary);
    gdlBottomLinks=intersect(gdlBoundary,membraneInfBoundary);
    
%     
%     internalLinks=network.GetLinksFrontiere(0);
%     linksNeighboors=network.GetPoresOfLink(internalLinks);
%     
%     linksNeighboorsElements=[membraneCodes(linksNeighboors(:,1)),membraneCodes(linksNeighboors(:,2))];
%     linksNeighboorsElements=sort(linksNeighboorsElements,2);
%     
%     gdlBottomLinksIndices= and(linksNeighboorsElements(:,1)==0,linksNeighboorsElements(:,2)==2 );
%     gdlBottomLinks=internalLinks(gdlBottomLinksIndices);
%     
%     gdlTopLinksIndices= and(linksNeighboorsElements(:,1)==0,linksNeighboorsElements(:,2)==1 );
%     gdlTopLinks=internalLinks(gdlTopLinksIndices);
%     
%     gdlPores = find(membraneCodes==0);
end
