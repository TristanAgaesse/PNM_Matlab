function outputInformation=ComputeHydrophobicityLoss(network,initialCluster,clusterPressure,inletLink,outletLink,options)
%COMPUTEHYDROPHOBICITYLOSS Compute how a cluster evolves when there is a 
%hydrophobicity loss
%
%Algorithm : - Initial step is a water cluster 
%            - then contact angle changes according to a degradation rate
%            (non uniform degradation is possible)
%           - from time to time a pore next to the cluster is invaded if its 
%           capillary pressure becomes smaller than the mean pressure of the 
%           water cluster. This pore invasion can induce several neighbours 
%           pores invasions as a consequence (this is called haines jump or
%           burst). Contact angle and degradation rate are updated after each invasion
%
%Input : network,initialCluster,clusterPressure,inletLink,outletLink,options
%       network : pore network on which the algorithm is run
%       initialCluster : initial water distribution (class=ClusterMonophasique)
%       clusterPressure : mean pressure inside the cluster (given by water drop boundary condition in channel)
%       inletLink,outletLink : specify inlet and outlet
%       options.nIterations= number of iterations
%       options.DeltaContactAngle = maximal loss in contact angle (in radian)
%       options.ClusterGrowth = true, false  : allow degradation to make new invasions or not 
%       options.DegradationMechanism='uniform','uniformInWater','waterSpeed'
%       options.InvasionTimeSearch='exactTimeForLaplaceLaw','localLinearisationOfCapillaryPressure','linearDecreaseOfCapillaryPressure'
%
%
%Output : outputInformation : information to analyse the degradation process 
%         outputInformation.times
%         outputInformation.clusters
%         outputInformation.invasionPressures
%         outputInformation.degradationPercentage
%         outputInformation.fixedPressure
%         outputInformation.inletLink
%         outputInformation.outletLink
%         outputInformation.options
%         outputInformation.fixedPressure

%---------------------------------------------------------------------------------------------    
    
    
    [nIterations,outputInformation,cluster,degradationPercentage,temps] = InitializeAlgorithm(...
                        options,network,initialCluster,clusterPressure,inletLink,outletLink);
    
    
    disp('Begin degradation')
    
    degradationSpeed = ComputeDegradationSpeed(network,cluster,inletLink,outletLink,options);
    
    disp('Iteration : ')
    for iIteration=1:nIterations
        fprintf ('%d ', iIteration)
        
        if options.ClusterGrowth == false
            
            outputInformation.invasionTimeDistribution{iIteration}=[];
            timeStep=1/max(degradationSpeed);
            temps=temps+timeStep;
            
            degradationPercentage = degradationPercentage+timeStep*degradationSpeed;
            UpdateContactAngle(network,degradationPercentage,options.DeltaContactAngle);
            
        elseif options.ClusterGrowth == true
            
            %recherche du prochain pore envahi sous l'effet de la degration
            try
                [indexInvadedLink,timeStep,temps_invasion_potentielle] = FindNextInvadedLink(...
                     cluster,clusterPressure,degradationPercentage,options.DeltaContactAngle,...
                     degradationSpeed,inletLink,outletLink,options);
            catch err
                if (strcmp(err.identifier,'FindNextInvadedLink:EmptyTempsPotentielValable'))
                    disp(err.message)
                    return
                else
                    rethrow(err);
                end
            end

            outputInformation.invasionTimeDistribution{iIteration} = temps_invasion_potentielle+temps;
            temps=temps+timeStep;

            if indexInvadedLink>0
                
                %envahissement de ce pore
                interfaceChangeInformation = cluster.InvadeNewPore(indexInvadedLink);
                cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink);

                outputInformation.invasionPressures{iIteration}=clusterPressure;
                
                %gestion du burst (=haines jump) subsequent
                [indexMinPressureLink,minPressure] = cluster.GetMinimalPressureLink;
                while minPressure<clusterPressure
                    
                    minPressureLink=cluster.InterfaceLinks(indexMinPressureLink);
                    
                    if network.GetFrontiereOfLink(minPressureLink)~=0
                        fprintf('New breakthrough point on boundary %d',network.GetFrontiereOfLink(minPressureLink));
                        cluster.InvadeOutletLink(minPressureLink);
                    else
                        interfaceChangeInformation=cluster.InvadeNewPore(indexMinPressureLink);
                        cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink);
                    end

                    outputInformation.invasionPressures{iIteration} = [outputInformation.invasionPressures{iIteration},minPressure];

                    [indexMinPressureLink,minPressure] = cluster.GetMinimalPressureLink;
                end
            end
        end
        
        outputInformation.nIteration=iIteration;
        outputInformation.times(iIteration)=temps;
        outputInformation.degradationPercentage{iIteration}=degradationPercentage;
        savedCluster=cluster.CopyCluster;
        savedCluster.Network=[];
        outputInformation.clusters{iIteration}=savedCluster;

        %evolution de l'etat de degradation
        degradationPercentage = degradationPercentage+timeStep*degradationSpeed;

        degradationSpeed = ComputeDegradationSpeed(network,cluster,inletLink,outletLink,options);
    end
    
end
    

%---------------------------------------------------------------------------------------------    
function [nIterations,outputInformation,cluster,degradationPercentage,temps] = InitializeAlgorithm(...
                        options,network,initialCluster,clusterPressure,inletLink,outletLink)

    nIterations=options.nIterations;

    outputInformation=struct;
    outputInformation.times=zeros(1,nIterations);
    outputInformation.clusters=cell(1,nIterations);
    outputInformation.invasionPressures=cell(1,nIterations);
    outputInformation.degradationPercentage=cell(1,nIterations);
    outputInformation.fixedPressure=0;

    outputInformation.inletLink=inletLink;
    outputInformation.outletLink=outletLink;
    outputInformation.options=options;

    outputInformation.fixedPressure=clusterPressure;

    cluster = initialCluster.CopyCluster;

    temps=0;

    degradationPercentage=zeros(1,network.GetNumberOfLinks);

end




%---------------------------------------------------------------------------------------------        
function [indexInvadedLink,temps,temps_invasion_potentielle]=FindNextInvadedLink(...
                        cluster,clusterPressure,degradationPercentage,deltaContactAngle,...
                        degradationSpeed,linkInlet,linkOutlet,options)

    %option='localLinearisationOfCapillaryPressure';
    methode=options.InvasionTimeSearch;

    criticalPressures=cluster.GetCriticalPressures;

    if strcmp(methode,'localLinearisationOfCapillaryPressure') || strcmp(methode,'exactTimeForLaplaceLaw')
        %ici on fait varier explicitement les angles de contact. 
        %On recalcule a chaque pas de temps les pressions critiques
        %en fonction des nouveaux angles de contact

        UpdateContactAngle(cluster.Network,degradationPercentage,deltaContactAngle);

        contactAngle = cluster.Network.GetLinkData('ContactAngle');
        
        interfaceUpdateInformation=cluster.GetInterfaceChangeInformation(1:length(cluster.GetInterfaceLinks));  %tous les liens frontiere
        UpdateCriticalPressure(cluster,interfaceUpdateInformation,linkInlet,linkOutlet)

        temps_invasion_potentielle=zeros(1,cluster.Network.GetNumberOfLinks);

        linksToUpdate = cluster.GetInterfaceLinks;
        isDegradating = (degradationSpeed(linksToUpdate)~=0);
        isInternalLink = (cluster.Network.GetFrontiereOfLink(linksToUpdate)==0);
        linksToUpdate = linksToUpdate(and(isDegradating,isInternalLink));

        if strcmp(methode,'localLinearisationOfCapillaryPressure')
            %Pour trouver le prochain pas de temps, on linearise
            %la pression localement en fonction de l'angle de contact

            foo = 1-(clusterPressure*ones(1,cluster.Network.GetNumberOfLinks))./(criticalPressures.*tan(contactAngle));
            temps_invasion_potentielle(linksToUpdate) = foo(linksToUpdate)./degradationSpeed(linksToUpdate);


        elseif strcmp(methode,'exactTimeForLaplaceLaw')
            %Pour trouver le prochain pas de temps, on inverse la loi
            %de Laplace
            R=cluster.Network.GetLinkDataList.Diameter(linksToUpdate)/2;
            gamma=cluster.ClusterOptions.SurfaceTension;
            deltaTheta = contactAngle(linksToUpdate)-acos(-clusterPressure*R/(2*gamma));
            temps_invasion_potentielle(linksToUpdate) = deltaTheta./degradationSpeed(linksToUpdate);
            
        end
        
        assert(isempty(find(isnan(temps_invasion_potentielle),1)),'NaN found !!!')
        
        temps_potentiels_valables=temps_invasion_potentielle(and(temps_invasion_potentielle>0,degradationPercentage<100));
        if isempty(temps_potentiels_valables)
            exception = MException('FindNextInvadedLink:EmptyTempsPotentielValable','Pas de nouvel envahissement possible');
            throw(exception);
        end

        temps=min(temps_potentiels_valables);
        ind=find(abs(temps_invasion_potentielle-temps)<1e-6*abs(temps),1);
        indexInvadedLink=find(cluster.GetInterfaceLinks==ind);
        assert(length(indexInvadedLink)==1)


    elseif strcmp(methode,'linearDecreaseOfCapillaryPressure')
        %ici on suppose que la pression decroit lineairement vis a vis
        %du pourcentage de perte de PTFE. La vitesse de perte de PTFE
        %est recalculee a chaque pas de temps
        temps_invasion_potentielle=zeros(1,cluster.Network.GetNumberOfLinks);
        foo=200*(1-((clusterPressure*ones(1,cluster.Network.GetNumberOfLinks))./criticalPressures))-degradationPercentage;

        for iLink=cluster.GetInterfaceLinks
            if degradationSpeed(iLink)~=0 && cluster.Network.GetFrontiereOfLink(iLink)==0
                temps_invasion_potentielle(iLink)=foo(iLink)/degradationSpeed(iLink);
            end
        end

        temps_potentiels_valables=temps_invasion_potentielle(and(temps_invasion_potentielle>0,degradationPercentage<100));

        temps=min(temps_potentiels_valables);
        ind=find(temps_invasion_potentielle==temps);
        indexInvadedLink=find(cluster.GetInterfaceLinks==ind);

    end
end



%---------------------------------------------------------------------------------------------
function UpdateContactAngle(network,degradationPercentage,deltaContactAngle)

    contactAngle=network.GetLinkDataList.ContactAngle;
    network.RemoveLinkData('ContactAngle');

    if not(isfield(network.GetLinkDataList,'InitialContactAngle'))
        network.AddNewLinkData(contactAngle,'InitialContactAngle');
    end

    initialContactAngle=network.GetLinkDataList.InitialContactAngle;
    
    boundedPercentage=degradationPercentage;
    boundedPercentage(degradationPercentage>100)=100;
    
    contactAngle=initialContactAngle-boundedPercentage/100*deltaContactAngle;
    network.AddNewLinkData(contactAngle,'ContactAngle');
end



%---------------------------------------------------------------------------------------------
function degradationSpeed=ComputeDegradationSpeed(network,clusterLiquide,inletLink,ouletLink,options)

    degradationSpeed=zeros(1,network.GetNumberOfLinks);
    degradationMechanism=options.DegradationMechanism;

    if strcmp(degradationMechanism,'uniform')
        degradationSpeed=ones(1,network.GetNumberOfLinks);

    elseif strcmp(degradationMechanism,'uniformInWater')
        degradationSpeed(clusterLiquide.GetInvadedLinks)=1;
        degradationSpeed(clusterLiquide.GetInterfaceLinks)=1;

    elseif strcmp(degradationMechanism,'waterSpeed')
        %coeff degradation(face)=moyenne des vitesse sur les faces des
        %pores_voisins
        
        dynamicViscosity = 1e-3; %viscosity water at ambiant conditions
        conductancesPermeability = LocalScaleComputeConductancesStokes(network,dynamicViscosity);

        boundaryConditions.inletLink = inletLink;
        boundaryConditions.outletLink = ouletLink;
        boundaryConditions.inletType = 'Dirichlet' ;
        boundaryConditions.outletType = 'Dirichlet' ;
        boundaryConditions.inletValue = 1*ones(1,length(boundaryConditions.inletLink));
        boundaryConditions.outletValue = 0.1*ones(1,length(boundaryConditions.outletLink));

        transportPores = clusterLiquide.GetInvadedPores;
        
        [ ~, waterFlux, ~ ] = ComputeLinearTransport(network,transportPores,conductancesPermeability,boundaryConditions);
        
        assert(isempty(find(isnan(waterFlux),1)),'NaN found !!!')
        
        nPore=network.GetNumberOfPores;
        vitesseMoyennePore=zeros(1,nPore);
        for iPore=1:nPore
            links=network.GetLinksOfPore(iPore);
            vitesseMoyennePore(iPore)=sum(arrayfun(@(x) abs(x),waterFlux(links)))/length(links);
        end

        pores = network.GetPoresOfLink(1:network.GetNumberOfLinks);
        internalLinks = network.GetLinksFrontiere(0);
        boundaryLinks= network.GetLinksFrontiere(1:network.GetNumberOfBoundaries);
        degradationSpeed(internalLinks) = 0.5*(vitesseMoyennePore(pores(internalLinks,1))+vitesseMoyennePore(pores(internalLinks,2)));
        degradationSpeed(boundaryLinks) = vitesseMoyennePore(pores(boundaryLinks,2));

    end
    
    
    assert(isempty(find(isnan(degradationSpeed),1)),'NaN found !!!')
    assert(isempty(find(degradationSpeed<0,1)),'Negative degradation speed found !!!')
end

