function conductances = LocalScaleComputeConductancesConduction(network,parameters)
    % input : - network, 
    %         - parameters : struct with the following optional fields
    %            parameters.GeometricModel.Pore = 'Cylinder' 
    %            parameters.GeometricModel.Link = 'None' or 'SurfaceResistance_RealSurface'
    %            parameters.PoreBulkProp : scalar or array(nPore,1)
    %            parameters.LinkBulkProp : scalar or array(nLink,1)
	% output : conductances
    %
    % Examples : 
    %   1) conductivity = 2e-5;
    %    parameters.PoreBulkProp = conductivity;
    %    conductances = LocalScaleComputeConductancesConduction(network,parameters)
    %     
    %   2) parameters.GeometricModel.Pore = 'Cylinder';
    %    parameters.GeometricModel.Link ='SurfaceResistance_RealSurface';
    %    parameters.PoreBulkProp = 2e-5*rand(1,nPore);
    %    parameters.LinkBulkProp = 1;
    %    conductances = LocalScaleComputeConductancesConduction(network,parameters) 
    %     
    
    
    %Get geometric data from the network
    
    [geometricModel,poreBulkProp,linkBulkProp] = ReadCheckInputs(network,parameters);
    
    
    nLink = network.GetNumberOfLinks;
    
    internalLinks = network.GetLinksFrontiere(0);
    boundaryLinks = network.GetLinksFrontiere(1:network.GetNumberOfBoundaries);
    
    
    %Compute conductances : 
    
    switch geometricModel.Pore
        
        case 'Cylinder'
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_Cylinder(...
                                                    network,poreBulkProp,internalLinks,boundaryLinks);
            
        case 'ConstricivityEquation'
        	[in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_ConstrictivityEquation(...
                                                    network,poreBulkProp,internalLinks,boundaryLinks);
                                                
    end
    
    
    switch geometricModel.Link
        
        case 'None'
            in_resistanceLink = 0;
            bound_resistanceLink = 0;    
        case 'SurfaceResistance_RealSurface'
            [in_resistanceLink,bound_resistanceLink]=ComputeResistanceLink_SurfaceResistance_RealSurface(...
                                                    network,linkBulkProp,internalLinks,boundaryLinks);
            
    end
    
    %pour chaque lien, resistances pores 1, lien et pore 2 en série
    conductances = zeros(nLink,1);
    conductances(internalLinks) = 1./(in_resistanceP1+in_resistanceP2+in_resistanceLink);
    conductances(boundaryLinks) = 1./(bound_resistanceP1+bound_resistanceP2+bound_resistanceLink);
    
end

%conductances(internalLinks)=diffusivity*linkSurface(internalLinks)./(distance1(internalLinks)+distance2);
    %conductances(boundaryLinks)=diffusivity*linkSurface(boundaryLinks)./(2*distance1(boundaryLinks));


%%
function [geometricModel,poreBulkProp,linkBulkProp] = ReadCheckInputs(network,parameter)
    
    assert(isa(network,'PoreNetwork'),'First input must be a pore network')
    nPore = network.GetNumberOfPores;
    nLink = network.GetNumberOfLinks;
    
    % default value
    defaultParameter = struct;
    defaultParameter.GeometricModel.Pore = 'Cylinder';
    defaultParameter.GeometricModel.Link = 'None';
    defaultParameter.PoreBulkProp = 2e-5*ones(1,nPore); %diff_O2_dans_N2
    defaultParameter.LinkBulkProp = zeros(1,nLink);
    
    % authorized values
    valid_Parameter                     = {'GeometricModel','PoreBulkProp','LinkBulkProp'};
    valid_Parameter_GeometricModel      = {'Pore','Link'};
    valid_Parameter_GeometricModel_Pore = {'Cylinder'};
    valid_Parameter_GeometricModel_Link = {'None','SurfaceResistance_RealSurface'};
    
    
    
    %Check parameters
    assert(isa(parameter,'struct'))
    CheckValues(fields(parameter),valid_Parameter)
    
    % Geometric model       
    if isfield(parameter,'GeometricModel')
        assert(isa(parameter.GeometricModel,'struct'))
        CheckValues(fields(parameter.GeometricModel),valid_Parameter_GeometricModel)
        geometricModel=parameter.GeometricModel;
    else
        geometricModel=defaultParameter.GeometricModel;
    end
    
    if isfield(geometricModel,'Pore')
        CheckValues(geometricModel.Pore,valid_Parameter_GeometricModel_Pore)
    else
        geometricModel.Pore=defaultParameter.GeometricModel.Pore;
    end
    
    if isfield(geometricModel,'Link')
        CheckValues(geometricModel.Link,valid_Parameter_GeometricModel_Link)
    else
        geometricModel.Link=defaultParameter.GeometricModel.Link;
    end
    
    
    % PoreBulkProp    
    if isfield(parameter,'PoreBulkProp')
        assert(all(parameter.PoreBulkProp>=0))
        if length(parameter.PoreBulkProp)==1
            %scalar input 
            poreBulkProp = parameter.PoreBulkProp*ones(nPore,1);
        else
            assert(length(parameter.PoreBulkProp)==nPore)
            poreBulkProp = zeros(nPore,1);
            poreBulkProp(:) = parameter.PoreBulkProp(:);
        end
    else
        poreBulkProp = defaultParameter.PoreBulkProp;
    end
    
    
    %LinkBulkProp
    if isfield(parameter,'LinkBulkProp')
        assert(all(parameter.PoreBulkProp>=0))
        if length(parameter.LinkBulkProp)==1
            %scalar input 
            linkBulkProp = parameter.LinkBulkProp*ones(nLink,1);
        else
            assert(length(parameter.LinkBulkProp)==nLink)
            linkBulkProp = zeros(nLink,1);
            linkBulkProp(:) = parameter.LinkBulkProp(:);
        end
    else
        linkBulkProp = defaultParameter.LinkBulkProp;
    end
    
end


%%
function [in_resistance1,in_resistance2,bound_resistance1,bound_resistance2]=ComputeResistancePore_Cylinder(...
                                                            network,poreBulkProp,internalLinks,boundaryLinks)

    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    
    dimension = network.Dimension;
    
    
    allLinks=1:nLink;
    CheckLinkDiameter(network)
    %linkSurface = network.GetLinkData('Surface');
    linkSurface = pi*(network.GetLinkData('Diameter')/2).^2;
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
        %Internal links
    in_resistance1 = distance1(internalLinks)./(poreBulkProp(network.LinkOwners(internalLinks)).*linkSurface(internalLinks));
    in_resistance2 = distance2./(poreBulkProp(network.LinkNeighbours(internalLinks)).*linkSurface(internalLinks));
    
        %Boundary links
    bound_resistance1 = distance1(boundaryLinks)./(poreBulkProp(network.LinkOwners(boundaryLinks)).*linkSurface(boundaryLinks));
    bound_resistance2 = 0;
    
end


function [in_resistance1,in_resistance2,bound_resistance1,bound_resistance2]=ComputeResistancePore_ConstrictivityEquation(...
                                                    network,poreBulkProp,internalLinks,boundaryLinks)
    
    % Effective conductance formula based on the article
    % "Quantitative relationship between microstructure and effective
    % transport properties based on virtual materials testing", 
    % G.Gaiselmann, M.Neumann, V.Schmidt. DOI: 10.1002/aic.14416
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    
    dimension = network.Dimension;
    
    allLinks=1:nLink;
    CheckLinkDiameter(network)
    %linkSurface = network.GetLinkData('Surface');
    linkDiameter = network.GetLinkData('Diameter');
    %linkSurface = pi*(linkDiameter/2).^2;
    poreDiameter = network.GetPoreData('Diameter');
    poreSurface = pi*(poreDiameter/2).^2;
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
    
    [~,porePhase]=unique(poreBulkProp);
    nPhase = max(porePhase);
    phaseVolume=zeros(1,nPhase);
    for i=1:nPhase
        phaseVolume(i)=sum(poreVolume(porePhase==i));
    end
    totalVolume = prod(network.GetImageSize)*(network.GetVoxelEdgeLength)^3;
    globalPorosity = phaseVolume/totalVolume;
    
    
    % Effective geometric parameters of the pore-link-pore
    localPorosity = globalPorosity(porePhase);
    tortuosity = 1;
    S = poreSurface;
    
        %Internal links
    pore1 = network.LinkOwners(internalLinks);
    constrictivity1 = (linkDiameter(internalLinks)/poreDiameter(pore1)).^2;
    conductance1 = poreBulkProp(pore1).*S(pore1)./distance1(internalLinks).*2.03*localPorosity(pore1).^1.57.*constrictivity1.^0.72./tortuosity.^2;    
    pore2 = network.LinkNeighbours(internalLinks);
    constrictivity2 = (linkDiameter(internalLinks)/poreDiameter(pore2)).^2;
    conductance2 = poreBulkProp(pore2).S(pore2)./distance2(internalLinks).*2.03*localPorosity(pore2).^1.57.*constrictivity2.^0.72./tortuosity.^2;
    
    in_resistance1 = 1./conductance1;
    in_resistance2 = 1./conductance2;
    
        %Boundary links
    pore = network.LinkOwners(boundaryLinks);   
    constrictivity = (linkDiameter(boundaryLinks)/poreDiameter(pore)).^2;
    conductance = poreBulkProp(pore).*S(pore)./distance1(boundaryLinks).*2.03*localPorosity(pore).^1.57.*constrictivity.^0.72./tortuosity.^2;
    bound_resistance1 = 1./conductance;
    bound_resistance2 = 0;
    
end 


function [in_resistanceLink,bound_resistanceLink]=ComputeResistanceLink_SurfaceResistance_RealSurface(...
                                                    network,linkBulkProp,internalLinks,boundaryLinks)
    
    
    linkSurface = network.GetLinkData('RealSurface');    

    in_resistanceLink=linkBulkProp(internalLinks)./linkSurface(internalLinks);



    bound_resistanceLink=zeros(1,length(boundaryLinks)); %Here I chose to put no surface resistance on boundary links

end


%%
function CheckValues(myValues,authorizedValues)
    checkValues=ismember(myValues,authorizedValues);
    if ~all(checkValues)
        disp(myValues(~checkValues))
        error('Wrong parameter')
    end
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
