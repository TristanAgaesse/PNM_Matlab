function conductances = LocalScaleComputeConductancesStokes(network,parameters)
% input : - network, 
    %         - parameters : struct with the following optional fields
    %            parameters.GeometricModel.Pore = 'Cylinder' 
    %            parameters.GeometricModel.Link = 'None' or 'SurfaceResistance_RealSurface'
    %            parameters.PoreBulkProp : scalar or array(nPore,1)
    %            parameters.LinkBulkProp : scalar or array(nLink,1)
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
    linkSurface = network.GetLinkData('Surface');
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);

    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);

    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);

        %Internal links
    in_resistance1 = 8*pi*poreBulkProp(network.LinkOwners(internalLinks)).*distance1(internalLinks)./linkSurface(internalLinks);
    in_resistance2 = 8*pi*poreBulkProp(network.LinkNeighbours(internalLinks)).*distance2./linkSurface(internalLinks);

        %Boundary links
    bound_resistance1 = 8*pi*poreBulkProp(network.LinkOwners(boundaryLinks)).*distance1(boundaryLinks)./linkSurface(boundaryLinks);
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
