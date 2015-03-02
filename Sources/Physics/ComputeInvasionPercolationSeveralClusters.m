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
%
%Output : [clusters,invadedPores]


    %Initialisation de l'algorithme.
    disp('Running Invasion Percolation Several Clusters');
    tic;
    
    CheckInputs(network,nCluster,clustersInletLink,clustersOutletLink,wettability,varargin )
    CheckLinkDiameter(network) %Verification si les diametres des liens sont deja calcules
    AssignContactAngle(network,wettability) %Assignation du Contact Angle
    

    clusterOptions = struct;
    if not(isempty(varargin))
        clusterOptions = varargin{1};        
    end
    
    
    allInvadedPores=zeros(1,network.GetNumberOfPores);
    clusters = cell(1,nCluster);
    fusionIndices = zeros(1,nCluster);
    
    for iCluster=1:nCluster
    
        cluster = ClusterMonophasique.InitialiseInvasionCluster(clustersInletLink{iCluster},clustersOutletLink{iCluster},network,clusterOptions);

        iteration = 0;
        currentPressure = 0;
        outlet_reached = false;
        outletPores = network.GetPoresFrontiere(clustersOutletLink{iCluster});
        nPore = network.GetNumberOfPores;
        invasionPressureList = zeros(1,nPore);

        %Find accessible pores
        nPoreAccessible=FindNumberOfAccessiblePores(network,clustersInletLink{iCluster});

        poreAllreadyInvaded=0;
        stopCondition = outlet_reached || iteration>=nPoreAccessible || poreAllreadyInvaded;
        
        %Boucle d'invasion pore par pore
        while not(stopCondition)
            iteration = iteration+1;
            %trouver la face de plus petite pression critique

            [indexInvadedLink,invasionPressure] = cluster.GetMinimalPressureLink;

            if invasionPressure>currentPressure
                currentPressure = invasionPressure;
            end

            invasionPressureList(iteration) = invasionPressure;    

            invadedPore = cluster.GetOutwardPore(indexInvadedLink);
           
            poreAllreadyInvaded = allInvadedPores(invadedPore);
            
            %envahir le pore associe et update les pressions critiques

            interfaceChangeInformation = cluster.InvadeNewPore(indexInvadedLink);
            cluster.UpdateCriticalPressure(interfaceChangeInformation,clustersInletLink{iCluster},clustersOutletLink{iCluster}); 

            %verifier si outlet_reached
            if ismember(invadedPore,outletPores)
                outlet_reached = true;
                breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),clustersOutletLink{iCluster});
                cluster.InvadeOutletLink(breakthroughLinks);
            end
            
            if not(poreAllreadyInvaded)
                allInvadedPores(invadedPore)=iCluster;
            end
            
            stopCondition = outlet_reached || iteration>=nPoreAccessible || poreAllreadyInvaded;
        end

        
        %Fusionner clusters si on a atteint un pore deja envahi
        if poreAllreadyInvaded
            iclusterAllreadyThere = fusionIndices(poreAllreadyInvaded);
            fusionIndices(iCluster)=iclusterAllreadyThere;
            clusters{iclusterAllreadyThere}=clusters{iclusterAllreadyThere}.FuseClusters(cluster);
        
        else
            clusters{iCluster}=cluster;
        end

    end
        
    
    
    %Mise en forme des output
    newClusters = cell(1,0);
    for iCluster=1:nCluster
        if ~ isempty(clusters{iCluster})
            newClusters{end+1}=clusters{iCluster};
        end
    end
    clusters=newClusters;
    
    invadedPores = find(allInvadedPores);

    
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Invasion percolation finished. Time spent : %d minutes %f s. \n',minutes,secondes);
    
end





%---------------------------------------------------------------------------------------------        
function CheckInputs(network,nCluster,clustersInletLink,clustersOutletLink,wettability,varargin )

    assert(isa(network,'PoreNetwork'),'First argument network must be a Pore Network object')
    assert(isa(clustersInletLink,'cell') &&  length(clustersInletLink)==nCluster,'Third argument clustersInletLink must be a cell of length nCluster (2nd argument)')
    assert(isa(clustersOutletLink,'cell') &&  length(clustersOutletLink)==nCluster,'Fourth argument clustersOutletLink must be a cell of length nCluster (2nd argument)')


end



%---------------------------------------------------------------------------------------------        
function nPoreAccessible=FindNumberOfAccessiblePores(network,clustersInletLink)

    fooCluster=network.CreateVoidCluster;
    totalFloodCluster=fooCluster.GetComplementaryCluster;
    percoPath=totalFloodCluster.FindPercolationPath(clustersInletLink,1:network.GetNumberOfLinks);
    nPoreAccessible=0;
    for i=1:length(percoPath)
        nPoreAccessible=nPoreAccessible+length(percoPath{i}.GetInvadedPores);
    end

end



%---------------------------------------------------------------------------------------------        
function CheckLinkDiameter(network)
    if not(isfield(network.GetLinkDataList,'Diameter'))
        disp('Calcul du diam�tre des liens...');
        tic;
        nLink = network.GetNumberOfLinks;
        diameter = zeros(1,nLink);
        for iLink = 1:nLink
            diameter(iLink) = network.ComputeLinkDiameter(iLink);
        end
        network.AddNewLinkData(diameter,'Diameter');
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        fprintf('Calcul du diam�tre des liens termin�. Dur�e : %d minutes %f s. \n',minutes,secondes);
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

