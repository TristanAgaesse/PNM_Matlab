function [clusters,invadedPores] = ComputeInvasionPercolationSeveralClusters( network,nCluster,clustersInletLink,clustersOutletLink,wettability,varargin )
%COMPUTEINVASIONPERCOLATIONSEVERALCLUSTERS  Calcule l'invasion de percolation sur un reseau
%de pores, avec plusieurs cluster
%Input : network,nCluster,clustersInletLink,clustersOutletLink,wettability   , ( varargin ) :
%       - network
%       - nCluster : number of cluster at inlet
%       - clustersInletLink : cell(1,nCluster) avec liste des liens d'injection pour chaque cluster   
%       - clustersOutletLink : cell(1,nCluster) avec liste des liens outlet pour chaque cluster
%       - wettability : 'currentWettability', 'hydrophobic', 'hydrophilic' or 'random'.
%       - varargin (optionnel) : clusterOptions 
%               clusterOptions.Coalescence = 'none' or 'numberOfInvadedNeighbours'
%               clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder','PurcellToroid'
%               clusterOptions.SurfaceTension = value of surface tension
%               clusterOptions.StopCondition = 'Breakthrough','OutletPoresReached'
%
%Output : [clusters,invadedPores]


    %Initialisation de l'algorithme.
    disp('Running Invasion Percolation Several Clusters');
    tic;
    
    [clusterOptions,stopCondition]=ReadCheckInputs(network,nCluster,clustersInletLink,clustersOutletLink,wettability,varargin{1} );
    CheckLinkDiameter(network) %Verification si les diametres des liens sont deja calcules
    AssignContactAngle(network,wettability) %Assignation du Contact Angle
 
    
    labelInvadedPores = zeros(1,network.GetNumberOfPores);
    clusters = cell(1,nCluster);
    fusionIndices = 1:nCluster;
    
    %Find connexe component of the network. It's usefull to stop invasion 
    %when the connexe component is fully invaded
    connexComponents = FindSizeOfConnexeComponents(network,cell2mat(clustersInletLink));
    
    
    for iCluster=1:nCluster
    
        cluster = ClusterMonophasique.InitialiseInvasionCluster(clustersInletLink{iCluster},clustersOutletLink{iCluster},network,clusterOptions);

        iteration = 0;
        currentPressure = 0;
        outlet_reached = false;
        outletPores = network.GetPoresFrontiere(clustersOutletLink{iCluster});
        nPore = network.GetNumberOfPores;
        invasionPressureList = zeros(1,nPore);

        %Find accessible pores
        %nPoreAccessible = FindNumberOfAccessiblePores(network,clustersInletLink{iCluster});
        outwardPore = cluster.GetOutwardPore(1:cluster.GetInterfaceLength);
        nPoreAccessible = GetNumberOfAccessiblePores(outwardPore,connexComponents);
        
        poreAllreadyInvaded = 0;
        stop = outlet_reached || iteration>=nPoreAccessible || poreAllreadyInvaded;
        
        %Boucle d'invasion pore par pore
        while not(stop)
            iteration = iteration+1;
            %trouver la face de plus petite pression critique

            [indexInvadedLink,invasionPressure] = cluster.GetMinimalPressureLink;

            if invasionPressure>currentPressure
                currentPressure = invasionPressure;
            end

            invasionPressureList(iteration) = invasionPressure;    

            invadedPore = cluster.GetOutwardPore(indexInvadedLink);

            if invadedPore<=0 && strcmp(stopCondition,'Breakthrough')
                %Stop if breakthrough 
                linkAbsoluteIndex=cluster.GetInterfaceLinkAbsoluteNumber(indexInvadedLink);
                assert(ismember(linkAbsoluteIndex,clustersOutletLink{iCluster}))
                outlet_reached = true;
                cluster.InvadeOutletLink(linkAbsoluteIndex);

            else
                assert(invadedPore >0)

                %envahir le pore associe et update les pressions critiques
                interfaceChangeInformation = cluster.InvadeNewPore(indexInvadedLink);
                cluster.UpdateCriticalPressure(interfaceChangeInformation,clustersInletLink{iCluster},clustersOutletLink{iCluster}); 
                
                poreAllreadyInvaded = labelInvadedPores(invadedPore);
                if not(poreAllreadyInvaded)
                    labelInvadedPores(invadedPore) = iCluster;
                end
                
                
                if strcmp(stopCondition,'OutletPoresReached')
                    %verifier si outlet_reached
                    if ismember(invadedPore,outletPores)
                        outlet_reached = true;
                        breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),clustersOutletLink{iCluster});
                        cluster.InvadeOutletLink(breakthroughLinks);
                    end
                end
            end
            

            stop = outlet_reached || iteration>=nPoreAccessible || poreAllreadyInvaded;
        end

        
        %Fusionner clusters si on a atteint un pore deja envahi
        if poreAllreadyInvaded
            iclusterAllreadyThere = fusionIndices(poreAllreadyInvaded);
            fusionIndices(iCluster)=iclusterAllreadyThere;
            clusters{iclusterAllreadyThere}=clusters{iclusterAllreadyThere}.FuseClusters(cluster);
        
        else
            clusters{iCluster}=cluster.CopyCluster;
        end

    end
        
    
    
    %Mise en forme des output
    newClusters = cell(1,0);
    for iCluster=1:nCluster
        if ~ isempty(clusters{iCluster})
            newClusters{end+1}=clusters{iCluster}.CopyCluster;
        end
    end
    clusters=newClusters;
    
    invadedPores = find(labelInvadedPores);

    
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Invasion percolation finished. Time spent : %d minutes %f s. \n',minutes,secondes);
    
end





%---------------------------------------------------------------------------------------------        
function [clusterOptions,stopCondition]=ReadCheckInputs(network,nCluster,clustersInletLink,clustersOutletLink,wettability,varargin )

    assert(isa(network,'PoreNetwork'),'First argument network must be a Pore Network object')
    assert(isa(clustersInletLink,'cell') &&  length(clustersInletLink)==nCluster,'Third argument clustersInletLink must be a cell of length nCluster (2nd argument)')
    assert(isa(clustersOutletLink,'cell') &&  length(clustersOutletLink)==nCluster,'Fourth argument clustersOutletLink must be a cell of length nCluster (2nd argument)')

    clusterOptions = struct;
    if not(isempty(varargin))
        clusterOptions = varargin{1};        
    end
    
    if isfield(clusterOptions,'StopCondition')
        assert( strcmp(clusterOptions.StopCondition,'Breakthrough') || strcmp(clusterOptions.StopCondition,'OutletPoresReached'),'clusterOptions.StopCondition = Breakthrough or OutletPoresReached')
        
        stopCondition = clusterOptions.StopCondition;
    else
        stopCondition = 'OutletPoresReached';
    end
end



%---------------------------------------------------------------------------------------------        
function connexComponents=FindSizeOfConnexeComponents(network,inletLink)

    fooCluster=network.CreateVoidCluster; 
	totalFloodCluster=fooCluster.GetComplementaryCluster; 
	percoPath=totalFloodCluster.FindPercolationPath(inletLink,1:network.GetNumberOfLinks); 
    
    connexComponents=zeros(1,network.GetNumberOfPores);
    for i=1:length(percoPath)
        thisConnexComponent = percoPath{i}.GetInvadedPores;
        connexComponents(thisConnexComponent) = i;
    end

end

%---------------------------------------------------------------------------------------------        
function nPoreAccessible = GetNumberOfAccessiblePores(outwardPore,connexComponents)

    thisConnexeComponent=unique(connexComponents(outwardPore));
    nPoreAccessible=sum(ismember(connexComponents,thisConnexeComponent));
    
end

%---------------------------------------------------------------------------------------------        
function CheckLinkDiameter(network)
    if not(isfield(network.GetLinkDataList,'Diameter'))
        diameter = network.ComputeAllLinkDiameter;
        network.AddNewLinkData(diameter,'Diameter');
    end
end


%---------------------------------------------------------------------------------------------        
function AssignContactAngle(network,wettability)
    
    if not(strcmp(wettability,'currentWettability'))
        if isfield(network.GetLinkDataList,'ContactAngle')
            network.RemoveLinkData('ContactAngle');
        end
        nLink = network.GetNumberOfLinks;
        contactAngles = zeros(1,nLink);
        for iLink = 1:nLink
            if strcmp(wettability,'hydrophobic')
                contactAngles(iLink) = pi*110/180;
            elseif strcmp(wettability,'hydrophilic')
                contactAngles(iLink) = pi*80/180;
            elseif strcmp(wettability,'random')
                contactAngles(iLink) = pi*(80+30*rand)/110;
            end
        end
        network.AddNewLinkData(contactAngles,'ContactAngle');
    end 

end


