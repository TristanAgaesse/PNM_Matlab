function [ diffusionCoefficient,totalSaturation,saturationProfile,cluster,concentrations ] = RunSimu_ImprovedWettabilityHomogeneity( network,contactAngleOption )
%RUNSIMU_IMPROVEDWETTABILITYHOMOGENEITY Summary of this function goes here
%   Detailed explanation goes here
    
    
    linkDiameter=2*network.GetLinkData('CapillaryRadius');
    network.AddNewLinkData(linkDiameter,'Diameter');
    linkSurface = pi*(network.GetLinkData('CapillaryRadius')).^2;
    network.AddNewLinkData(linkSurface,'Surface');
    
    
    gdlBottomLinks = find( network.GetLinkData('GdlBottomLinks') );
    gdlTopLinks = find( network.GetLinkData('GdlTopLinks') );
    gdlPores = find( network.GetPoreData('GdlPores') );
    
    
    %Define contact angle
    %minContactAngle = 100*pi/180;
    %maxContactAngle = 130*pi/180;
    
    
    
    contactAngle = DefineContactAngle(network,gdlPores,contactAngleOption);
    network.AddNewLinkData(contactAngle,'ContactAngle');
    
    
    
    
    % Invasion percolation
    clusterOptions.ThroatPressure='LaplaceCylinder';
    clusterOptions.Coalescence = 'none' ;
    clusterOptions.SurfaceTension=72e-3;
    waterInjectionLinks = network.GetLinksFrontiere(2);
    waterBreakthroughLinks = gdlTopLinks;
    [cluster,~,~]=ComputeInvasionPercolation(network,waterInjectionLinks,waterBreakthroughLinks,'currentWettability',clusterOptions);
    
    % Gas diffusion of the dry GDL
    [dryDiffusionCoefficient,~] = ComputeGasDiffusion(network,network.CreateVoidCluster,gdlBottomLinks,gdlTopLinks,gdlPores);
    
    
    % Gas diffusion in the non-invaded pores
    [wetDiffusionCoefficient,concentrations] = ComputeGasDiffusion(network,cluster,gdlBottomLinks,gdlTopLinks,gdlPores);
    
    saturationOptions.type = 'saturationProfile'; 
    saturationOptions.axe = [1 0 0];
    saturationOptions.codeForLiquid=1;
    saturationOptions.codeForSolid=255;
    [ totalSaturation, saturationProfile ] = ComputeSaturation( cluster, network, saturationOptions );
    
    diffusionCoefficient=wetDiffusionCoefficient/dryDiffusionCoefficient;
    
    %image_IP_Result = uint8(network.GetImagePoreData('ordreInvasionFctPressure'));




end



function contactAngle = DefineContactAngle(network,gdlPores,option)

    disp(option)

    % Trouver la position de chaque lien suivant la profondeur
    [gdlBoundary,gdlInnerLinks]=network.GetPoreRegionBoundaryLinks(gdlPores);
    
    contactAngle = 80*ones(1,network.GetNumberOfLinks);
    
    allGDLLinks=[gdlBoundary,gdlInnerLinks];
    centers = network.GetLinkCenter(allGDLLinks);
    Zlink = centers(:,1);
    
    %Attribuer à chaque lien un angle de contact
    
    Zmin = min(Zlink);
    Zmax = max(Zlink);
    
    
    if strcmp(option.Type,'quadratique')
        
        minContactAngle = option.MinContactAngle;
        maxContactAngle = option.MaxContactAngle;
        
        ZlinkNormalise=(Zlink-Zmin)./(Zmax-Zmin);
        %    loi parabolique entre min et max contact angle
        contactAngle(allGDLLinks)=  minContactAngle + ((ZlinkNormalise-1/2).^2).*4*(maxContactAngle-minContactAngle);
        
    elseif strcmp(option.Type,'PTFE_sur_petite_epaisseur')
        
        epaisseur = option.EpaisseurRelative;
        maxContactAngle = option.MaxContactAngle;
        minContactAngle = option.MinContactAngle;
        
        normalLinks = and(Zlink > Zmin + epaisseur*(Zmax-Zmin), Zlink<Zmax - epaisseur*(Zmax-Zmin));
        
        PTFELinks = setdiff(allGDLLinks,normalLinks);
        
        
        contactAngle(PTFELinks)=maxContactAngle;
        contactAngle(normalLinks)=minContactAngle;
        
    else
        assert(False,'wrong contact angle option')
    end 

end





function [diffusionCoefficient,concentrations]=ComputeGasDiffusion(network,cluster,gdlBottomLinks,gdlTopLinks,gdlPores)

    % Diffusion dans la phase gazeuse restante
    transportPores = intersect(cluster.GetInvadedPoresComplementary,gdlPores) ;

    diffusivity = 2e-5; % 02 in N2 at ambiant conditions
    conductancesDiffusion = LocalScaleComputeConductancesDiffusion(network,diffusivity);


    boundaryConditions.inletLink = gdlTopLinks;
    boundaryConditions.outletLink = gdlBottomLinks;

    boundaryConditions.inletType = 'Dirichlet' ;
    nInletLink=length(boundaryConditions.inletLink);
    boundaryConditions.inletValue = 2 *ones(1,nInletLink);

%     boundaryConditions.inletType = 'Neumann' ;
%     allLinkSurfaces = network.GetLinkData('Surface');   
%     surfacicFlux = 1e-5;
%     boundaryConditions.inletValue = surfacicFlux*allLinkSurfaces(boundaryConditions.inletLink);
%     
    boundaryConditions.outletType = 'Dirichlet' ;
    nOutletLink = length(boundaryConditions.outletLink);
    boundaryConditions.outletValue = 1 *ones(1,nOutletLink);
    


    [ concentrations, ~, diffusionCoefficient ] = ComputeLinearTransport( ...
                                                    network,transportPores, ...
                                                    conductancesDiffusion, ...
                                                    boundaryConditions  ...
                                                    );

end

