function conductances = LocalScaleComputeConductancesStokes(network,parameters)
% input : - network, 
    %         - parameters : struct with the following optional fields
    %            parameters.GeometricModel.Pore = 'Cylinder' or '2Cubes'
    %            parameters.GeometricModel.Link = 'None' or 'SurfaceResistance_RealSurface'
    %            parameters.PoreBulkProp : bulk dynamic viscosity - scalar or array(nPore,1)
    %            parameters.LinkBulkProp : bulk dynamic viscosity - scalar or array(nLink,1)
	% output : conductances
    %
    % Examples : 
    %   1) dynamicViscosity = 1e-3;
    %    parameters.PoreBulkProp = dynamicViscosity;
    %    conductances = LocalScaleComputeConductancesStokes(network,parameters)
    %     
    %   2) parameters.GeometricModel.Pore = 'Cylinder';
    %    parameters.GeometricModel.Link ='SurfaceResistance_RealSurface';
    %    parameters.PoreBulkProp = 1e-3*rand(1,nPore);
    %    parameters.LinkBulkProp = 1;
    %    conductances = LocalScaleComputeConductancesStokes(network,parameters) 
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
        case '2Cubes'
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_2Cubes(...
                                                    network,poreBulkConductivity,internalLinks,boundaryLinks);
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

    
%%
function [geometricModel,poreBulkProp,linkBulkProp] = ReadCheckInputs(network,parameter)
    
    assert(isa(network,'PoreNetwork'),'First input must be a pore network')
    nPore = network.GetNumberOfPores;
    nLink = network.GetNumberOfLinks;
    
    % default value
    defaultParameter = struct;
    defaultParameter.GeometricModel.Pore = 'Cylinder';
    defaultParameter.GeometricModel.Link = 'None';
    defaultParameter.PoreBulkProp = 1e-3*ones(1,nPore); %viscosite_dyn_water
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
    in_resistance1 = 8*pi*poreBulkProp(network.LinkOwners(internalLinks)).*distance1(internalLinks)./(linkSurface(internalLinks).^2);
    in_resistance2 = 8*pi*poreBulkProp(network.LinkNeighbours(internalLinks)).*distance2./(linkSurface(internalLinks).^2);

        %Boundary links
    bound_resistance1 = 8*pi*poreBulkProp(network.LinkOwners(boundaryLinks)).*distance1(boundaryLinks)./(linkSurface(boundaryLinks).^2);
    bound_resistance2 = 0;
                
end

function [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_2Cubes(...
                                                    network,poreBulkProp,internalLinks,boundaryLinks)
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    
    dimension = network.Dimension;
    
    
    allLinks=1:nLink;
    CheckLinkDiameter(network)
    linkRadius = network.GetLinkData('Diameter')/2;
    poreRadius = network.GetPoreData('Volume').^(1/3);
    %linkSurface = pi*(network.GetLinkData('Diameter')/2).^2;
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
        %Internal links
    lPore = distance1(internalLinks)./(1+linkRadius(internalLinks)/poreRadius(network.LinkOwners(internalLinks)));
    lLink = distance1(internalLinks)-lPore;
    resistanceCubePore = lPore./(poreRadius(network.LinkOwners(internalLinks)).^4);     
    resistanceCubeLink = lLink./(linkRadius(internalLinks).^4); 
    in_resistanceP1 = (resistanceCubePore+resistanceCubeLink)./(8*pi*poreBulkProp(network.LinkOwners(internalLinks)));
    
    lPore = distance2(internalLinks)./(1+linkRadius(internalLinks)/poreRadius(network.LinkNeighbours(internalLinks)));
    lLink = distance2(internalLinks)-lPore;
    resistanceCubePore = lPore./(poreRadius(network.LinkNeighbours(internalLinks)).^4);     
    resistanceCubeLink = lLink./(linkRadius(internalLinks).^4); 
    in_resistanceP2 = (resistanceCubePore+resistanceCubeLink)./(8*pi*poreBulkProp(network.LinkNeighbours(internalLinks)));

        %Boundary links
    lPore = distance1(boundaryLinks)./(1+linkRadius(boundaryLinks)/poreRadius(network.LinkOwners(boundaryLinks)));
    lLink = distance1(boundaryLinks)-lPore;
    resistanceCubePore = lPore./(poreRadius(network.LinkOwners(boundaryLinks)).^4);     
    resistanceCubeLink = lLink./(linkRadius(boundaryLinks).^4); 
    bound_resistanceP1 = (resistanceCubePore+resistanceCubeLink)./(8*pi*poreBulkProp(network.LinkOwners(boundaryLinks)));    
        
    bound_resistanceP2 = 0;                                            
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
