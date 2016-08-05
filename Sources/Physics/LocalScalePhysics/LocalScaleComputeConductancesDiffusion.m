function conductances = LocalScaleComputeConductancesDiffusion(network,parameters)
    % input : - network, 
    %         - parameters : struct with the following optional fields
    %            parameters.GeometricModel.Pore = 'Cylinder' or '2Cubes' or
    %               'PolynomialProfileDegree1' or 'VolumeConservation_LinearProfile'
    %            parameters.GeometricModel.Link = 'None' or 'SurfaceResistance_RealSurface'
    %            parameters.PoreBulkProp : bulk diffusivity - scalar or array(nPore,1)
    %            parameters.LinkBulkProp : bulk diffusivity - scalar or array(nLink,1)
	% output : conductances
    %
    % Examples : 
    %   1) diffusivity = 2e-5;
    %    parameters.PoreBulkProp = diffusivity;
    %    conductances = LocalScaleComputeConductancesDiffusion(network,parameters)
    %	
    %   2) parameters.GeometricModel.Pore = 'Cylinder';
    %    parameters.GeometricModel.Link ='SurfaceResistance_RealSurface';
    %    parameters.PoreBulkProp = 2e-5*rand(1,nPore);
    %    parameters.LinkBulkProp = 1;
    %    conductances = LocalScaleComputeConductancesDiffusion(network,parameters) 
    %     
    
    
    %Get geometric data from the network
    
    [geometricModel,poreBulkConductivity,linkBulkConductivity] = ReadCheckInputs(network,parameters);
    
    
    nLink = network.GetNumberOfLinks;
    
    internalLinks = network.GetLinksFrontiere(0);
    boundaryLinks = network.GetLinksFrontiere(1:network.GetNumberOfBoundaries);
    
    
    %Compute conductances : 
    
    switch geometricModel.Pore
        
        case 'Cylinder' % constant profile with surface adjusted to link surface
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_Cylinder(...
                                                    network,poreBulkConductivity,internalLinks,boundaryLinks);
                                                
        case 'ConstantProfilePoreRadius' % constant profile with surface adjusted to pore volume
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_ConstantProfilePore(...
                                                            network,poreBulkConductivity,internalLinks,boundaryLinks)    ;                                        
                                                
        case '2Cubes'  % square profil: two cubes
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_2Cubes(...
                                                    network,poreBulkConductivity,internalLinks,boundaryLinks);
                                                
        case 'PolynomialProfileDegree1' %linear profile between link and pore radius
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_PolynomialProfileDegree1(...
                                                    network,poreBulkConductivity,internalLinks,boundaryLinks);
            
        case 'VolumeConservation_LinearProfile'  %linear profile with scaling to have volume conservation
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_VolumeConservation_LinearProfile(...
                                                    network,poreBulkConductivity,internalLinks,boundaryLinks);
            
        
        case 'VolumeConservation_ConstantProfile' %constant profile with scaling to have volume conservation
            [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_VolumeConservation_ConstantProfile(...
                                                            network,poreBulkConductivity,internalLinks,boundaryLinks);
    end
    
    
    switch geometricModel.Link
        
        case 'None'     %no surface resistance at link crossing
            in_resistanceLink = 0;
            bound_resistanceLink = 0;    
        case 'SurfaceResistance_RealSurface' %surface resistance at link crossing proportional to link surface
            [in_resistanceLink,bound_resistanceLink]=ComputeResistanceLink_SurfaceResistance_RealSurface(...
                                                    network,linkBulkConductivity,internalLinks,boundaryLinks);
            
    end
    
    %pour chaque lien, resistances pores 1, lien et pore 2 en série
    conductances = zeros(nLink,1);
    conductances(internalLinks) = 1./(in_resistanceP1+in_resistanceP2+in_resistanceLink);
    conductances(boundaryLinks) = 1./(bound_resistanceP1+bound_resistanceP2+bound_resistanceLink);
    
    
    if sum(isinf(conductances))~=0
        disp('Infinite conductance ! we put them to 0')
        conductances(isinf(conductances))=0;
    end
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
    valid_Parameter_GeometricModel_Pore = {'Cylinder','2Cubes','PolynomialProfileDegree1','VolumeConservation_LinearProfile',...
        'ConstantProfilePoreRadius','VolumeConservation_ConstantProfile'};
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
    %linkSurface = pi*(network.GetLinkData('CapillaryRadius')).^2;
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

function [in_resistance1,in_resistance2,bound_resistance1,bound_resistance2]=ComputeResistancePore_ConstantProfilePore(...
                                                            network,poreBulkProp,internalLinks,boundaryLinks)
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    
    dimension = network.Dimension;
    
    
    allLinks=1:nLink;
    CheckLinkDiameter(network)
    %linkSurface = network.GetLinkData('Surface');
    %poreSurface = pi*(network.GetPoreData('Volume').^(1/3)/2).^2;
    poreVolume=network.GetPoreData('Volume');
    poreDiameter =poreVolume.^(1/3);
    poreSurface = poreDiameter.^2;
    %linkSurface = pi*(network.GetLinkData('Diameter')/2).^2;
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
%     %Compute volume scaling factor for each pore
%     poreVolumeApproximation = zeros(nPore,1);
%     
%     for iPore=1:nPore
%         
%         poreLinks = network.GetLinksOfPore(iPore);
%         nPoreLinks=length(poreLinks);
%         if nPoreLinks==0
%             poreVolumeApproximation(iPore)=poreVolume(iPore);
%         else
%             thisPoreDiameter = poreDiameter(iPore);
%             foo=linkCenter(poreLinks,:)-poreCenter(iPore*ones(1,nPoreLinks),:);
%             L = FastNorm(foo,dimension);
%             poreVolumeApproximation(iPore)=sum(L)*thisPoreDiameter^2-thisPoreDiameter^3*(nPoreLinks-1);%    pi*(thisPoreRadius^3).*(sum(L.*(1+beta+beta.^2)/(thisPoreRadius*3))-nPoreLinks*(2/3)+4/3);
%         end
%     end
%         
%     scaleFactor = (poreVolumeApproximation./poreVolume).^3;
%     scaleFactor(scaleFactor<0)=1;
%     scaleFactor(scaleFactor>1)=1;
%     %assert(all(scaleFactor>0) && all(scaleFactor<=1))
    scaleFactor=ones(nPore,1);
    
        %Internal links
    in_resistance1 = distance1(internalLinks)./(poreBulkProp(network.LinkOwners(internalLinks)).*poreSurface(network.LinkOwners(internalLinks)).*scaleFactor(network.LinkOwners(internalLinks)));%linkSurface(internalLinks));
    in_resistance2 = distance2./(poreBulkProp(network.LinkNeighbours(internalLinks)).*poreSurface(network.LinkNeighbours(internalLinks)).*scaleFactor(network.LinkNeighbours(internalLinks)));%linkSurface(internalLinks));
    
        %Boundary links
    bound_resistance1 = distance1(boundaryLinks)./(poreBulkProp(network.LinkOwners(boundaryLinks)).*poreSurface(network.LinkOwners(boundaryLinks)).*scaleFactor(network.LinkOwners(boundaryLinks)));%linkSurface(boundaryLinks));
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
    %linkRadius = network.GetLinkData('CapillaryRadius')
    %poreRadius = (network.GetPoreData('Volume').^(1/3))/2;
    poreRadius = network.GetPoreData('Diameter')/2;
    %poreRadius = network.GetPoreData('InscribedSphereRadius');
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
    distance1(distance1==0)=min(distance1(distance1>0));
    distance2(distance2==0)=min(distance2(distance2>0));
    
        %Internal links
    lPore = distance1(internalLinks)./(1+linkRadius(internalLinks)./poreRadius(network.LinkOwners(internalLinks)));
    lLink = distance1(internalLinks)-lPore;
    foo=lLink<=0;
    lPore(foo)=distance1(internalLinks(foo))/2;
    lLink(foo)=lPore(foo);
    assert(all(lLink>0))
    
    resistanceCubePore = 0.25*lPore./(poreRadius(network.LinkOwners(internalLinks)).^2);     
    resistanceCubeLink = 0.25*lLink./(linkRadius(internalLinks).^2); 
    in_resistanceP1 = (resistanceCubePore+resistanceCubeLink)./poreBulkProp(network.LinkOwners(internalLinks));
    
    
    lPore = distance2./(1+linkRadius(internalLinks)./poreRadius(network.LinkNeighbours(internalLinks)));
    lLink = distance2-lPore;
    foo=lLink<=0;
    lPore(foo)=distance2(foo)/2;
    lLink(foo)=lPore(foo);
    assert(all(lLink>0))
    
    resistanceCubePore = 0.25*lPore./(poreRadius(network.LinkNeighbours(internalLinks)).^2);     
    resistanceCubeLink = 0.25*lLink./(linkRadius(internalLinks).^2); 
    in_resistanceP2 = (resistanceCubePore+resistanceCubeLink)./poreBulkProp(network.LinkNeighbours(internalLinks));
    
        %Boundary links
    lPore = distance1(boundaryLinks)./(1+linkRadius(boundaryLinks)./poreRadius(network.LinkOwners(boundaryLinks)));
    lLink = distance1(boundaryLinks)-lPore;
    foo=lLink<=0;
    lPore(foo)=distance1(boundaryLinks(foo))/2;
    lLink(foo)=lPore(foo);
    assert(all(lLink>0))
    
    resistanceCubePore = 0.25*lPore./(poreRadius(network.LinkOwners(boundaryLinks)).^2);     
    resistanceCubeLink = 0.25*lLink./(linkRadius(boundaryLinks).^2); 
    bound_resistanceP1 = (resistanceCubePore+resistanceCubeLink)./poreBulkProp(network.LinkOwners(boundaryLinks));    
    
    bound_resistanceP2 = 0;                                            
end

function [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_PolynomialProfileDegree1(...
                                                    network,poreBulkProp,internalLinks,boundaryLinks)
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    
    dimension = network.Dimension;
    
    
    allLinks=1:nLink;
    CheckLinkDiameter(network)
    linkRadius = network.GetLinkData('Diameter')/2;
    poreRadius = network.GetPoreData('Diameter')/2;
    
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
        %Internal links
    effectiveSurface = pi*linkRadius(internalLinks).*poreRadius(network.LinkOwners(internalLinks)) ;
    in_resistanceP1 = distance1(internalLinks)./(poreBulkProp(network.LinkOwners(internalLinks)).*effectiveSurface);
    
    effectiveSurface = pi*linkRadius(internalLinks).*poreRadius(network.LinkNeighbours(internalLinks)) ;
    in_resistanceP2 = distance2./(poreBulkProp(network.LinkNeighbours(internalLinks)).*effectiveSurface);
    
        %Boundary links
    effectiveSurface =  pi*linkRadius(boundaryLinks).*poreRadius(network.LinkOwners(boundaryLinks));
    bound_resistanceP1 = distance1(boundaryLinks)./(poreBulkProp(network.LinkOwners(boundaryLinks)).*effectiveSurface);
    bound_resistanceP2 = 0;
    
end                                                


function [in_resistanceP1,in_resistanceP2,bound_resistanceP1,bound_resistanceP2]=ComputeResistancePore_VolumeConservation_LinearProfile(...
                                                    network,poreBulkProp,internalLinks,boundaryLinks)
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    
    dimension = network.Dimension;
    
    
    allLinks=1:nLink;
    CheckLinkDiameter(network)
    linkRadius = network.GetLinkData('CapillaryRadius');
    poreRadius = network.GetPoreData('InscribedSphereRadius');
    poreVolume = network.GetPoreData('Volume');
    
    %assert(all(4*pi/3*poreRadius.^3<=poreVolume))
    
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
    %Compute volume scaling factor for each pore
    scaleFactor = zeros(nPore,1);
    for iPore=1:nPore
        
        poreLinks = network.GetLinksOfPore(iPore);
        nPoreLinks=length(poreLinks);
        if nPoreLinks==0
            scaleFactor(iPore)=1;
        else
            thisPoreRadius = poreRadius(iPore);
            beta = linkRadius(poreLinks)/thisPoreRadius;
            foo=linkCenter(poreLinks,:)-poreCenter(iPore*ones(1,nPoreLinks),:);
            L = FastNorm(foo,dimension);
            %volume is a2*R^2-a3*R^3
            a2=pi*sum(L.*(1+beta+beta.^2))/3;
            a3=pi*(nPoreLinks*(2/3)-4/3);
            maxAccessibleVolume=(a2^3/a3^2)*(4/9-8/27);
            if not(maxAccessibleVolume>=poreVolume(iPore))
                scaleFactor(iPore)=1;
            else
                %find a conductance typical length which conserves pore volume 
                p=[a3, -a2, 0 , poreVolume(iPore)];
                r = roots(p);
                r = r(isreal(r));

                scaleFactor(iPore)=min(r(r>0))/poreRadius(iPore);
            end
        end
    end
    scaleFactor( scaleFactor>1)=1;
    %scaleFactor = (poreVolumeApproximation./poreVolume).^3;
   % assert(all(scaleFactor>0) && all(scaleFactor<=1))
    
        %Internal links
    effectiveSurface = pi*linkRadius(internalLinks).*poreRadius(network.LinkOwners(internalLinks)) ;
    poreList=network.LinkOwners(internalLinks);
    in_resistanceP1 = distance1(internalLinks)./(poreBulkProp(poreList).*effectiveSurface.*scaleFactor(poreList).^2);
    
    effectiveSurface = pi*linkRadius(internalLinks).*poreRadius(network.LinkNeighbours(internalLinks)) ;
    poreList=network.LinkNeighbours(internalLinks);
    in_resistanceP2 = distance2./(poreBulkProp(poreList).*effectiveSurface.*scaleFactor(poreList).^2);
    
        %Boundary links
    effectiveSurface =  pi*linkRadius(boundaryLinks).*poreRadius(network.LinkOwners(boundaryLinks));
    poreList=network.LinkOwners(boundaryLinks);
    bound_resistanceP1 = distance1(boundaryLinks)./(poreBulkProp(poreList).*effectiveSurface.*scaleFactor(poreList).^2);
    bound_resistanceP2 = 0;                                            
    
end

function [in_resistance1,in_resistance2,bound_resistance1,bound_resistance2]=ComputeResistancePore_VolumeConservation_ConstantProfile(...
                                                            network,poreBulkProp,internalLinks,boundaryLinks)
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    
    dimension = network.Dimension;
    
    
    allLinks=1:nLink;
    CheckLinkDiameter(network)
    %linkSurface = network.GetLinkData('Surface');
    %poreSurface = pi*(network.GetPoreData('Volume').^(1/3)/2).^2;
    poreVolume=network.GetPoreData('Volume');
    poreDiameter =poreVolume.^(1/3);
    poreSurface = poreDiameter.^2;
    %linkSurface = pi*(network.GetLinkData('Diameter')/2).^2;
    poreCenter=network.GetPoreCenter(1:nPore);
    linkCenter=network.GetLinkCenter(1:nLink);
    
    a=poreCenter(network.LinkOwners(allLinks),:)-linkCenter(allLinks,:);
    distance1=FastNorm(a,dimension);
    
    b=poreCenter(network.LinkNeighbours(internalLinks),:)-linkCenter(internalLinks,:);
    distance2=FastNorm(b,dimension);
    
    
    scaleFactor=zeros(nPore,1);
    
    for iPore=1:nPore
        
        poreLinks = network.GetLinksOfPore(iPore);
        nPoreLinks=length(poreLinks);
        if nPoreLinks==0
            scaleFactor(iPore)=1;
        else
            foo=linkCenter(poreLinks,:)-poreCenter(iPore*ones(1,nPoreLinks),:);
            L = FastNorm(foo,dimension);
            %volume is a2*R^2-a3*R^3
            a2=pi*sum(L);
            a3=pi*(nPoreLinks*(2/3)-4/3);
            maxVolume=(a2^3/a3^2)*(4/9-8/27);
            if not(maxVolume>=poreVolume(iPore))
                scaleFactor(iPore)=1;
            else
                p=[a3, -a2, 0 , poreVolume(iPore)];
                r = roots(p);
                r = r(isreal(r));
                foo=min(r(r>0));
                if ~isempty(foo)
                    scaleFactor(iPore)=min(r(r>0))/(poreDiameter(iPore)/2);
                else
                    scaleFactor(iPore)=1;
                end
            end
        end
    end
    scaleFactor( scaleFactor>1)=1;
    
        %Internal links
    in_resistance1 = distance1(internalLinks)./(poreBulkProp(network.LinkOwners(internalLinks)).*poreSurface(network.LinkOwners(internalLinks)).*scaleFactor(network.LinkOwners(internalLinks)).^2);%linkSurface(internalLinks));
    in_resistance2 = distance2./(poreBulkProp(network.LinkNeighbours(internalLinks)).*poreSurface(network.LinkNeighbours(internalLinks)).*scaleFactor(network.LinkNeighbours(internalLinks)).^2);%linkSurface(internalLinks));
    
        %Boundary links
    bound_resistance1 = distance1(boundaryLinks)./(poreBulkProp(network.LinkOwners(boundaryLinks)).*poreSurface(network.LinkOwners(boundaryLinks)).*scaleFactor(network.LinkOwners(boundaryLinks)).^2);%linkSurface(boundaryLinks));
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
