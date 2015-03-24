function [cluster,breakthroughPressure,invasionPressureList]  =  ComputeInvasionPercolation(network,inletLink,outletLink,wettability,varargin)
%ComputeInvasionPercolation Calcule l'invasion de percolation sur un r�seau
%de pores, avec un unique cluster
%Input : network,inletLink,outletLink,wettability   , ( varargin ) :
%       - network
%       - inletLink : liste des liens d'injection   
%       - outletLink : liste des liens de percée possible
%       - wettability : 'currentWettability', 'hydrophobic', 'hydrophilic' or 'random'.
%       - varargin (optionnel) : clusterOptions 
%               clusterOptions.Coalescence = 'none' or 'numberOfInvadedNeighbours'
%               clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder','PurcellToroid'
%               clusterOptions.SurfaceTension = value of surface tension
%
%Output : [cluster,breakthroughPressure,invasionPressureList]


    %Initialisation de l'algorithme.
    disp('Running Invasion Percolation');
    tic;
    CheckLinkDiameter(network) %Verification si les diametres des liens sont deja calcules
    AssignContactAngle(network,wettability) %Assignation du Contact Angle


    clusterOptions = struct;
    if not(isempty(varargin))
        clusterOptions = varargin{1};        
    end
    cluster = ClusterMonophasique.InitialiseInvasionCluster(inletLink,outletLink,network,clusterOptions);
    
    iteration = 0;
    currentPressure = 0;
    outlet_reached = false;
    outletPores = network.GetPoresFrontiere(outletLink);
    nPore = network.GetNumberOfPores;
    invasionPressureList = zeros(1,nPore);
    
    %Find accessible pores
    nPoreAccessible=FindNumberOfAccessiblePores(network,inletLink);
    
    stopCondition = outlet_reached || iteration>=nPoreAccessible ;
    
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
        
        %envahir le pore associe et update les pressions critiques
        
        interfaceChangeInformation = cluster.InvadeNewPore(indexInvadedLink);
        cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink); 
        
        %verifier si outlet_reached
        if ismember(invadedPore,outletPores)
            outlet_reached = true;
            breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),outletLink);
            cluster.InvadeOutletLink(breakthroughLinks);
        end
        
        stopCondition = outlet_reached || iteration>=nPoreAccessible ;
    end
    
    %Mise en forme des output
    breakthroughPressure = currentPressure;
    
    invasionPressureList = invasionPressureList(1:iteration);
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Invasion percolation finished. Time spent : %d minutes %f s. \n',minutes,secondes);
    
end



%---------------------------------------------------------------------------------------------        
function nPoreAccessible=FindNumberOfAccessiblePores(network,inletLink)

    fooCluster=network.CreateVoidCluster;
    totalFloodCluster=fooCluster.GetComplementaryCluster;
    percoPath=totalFloodCluster.FindPercolationPath(inletLink,1:network.GetNumberOfLinks);
    nPoreAccessible=0;
    for i=1:length(percoPath)
        nPoreAccessible=nPoreAccessible+length(percoPath{i}.GetInvadedPores);
    end

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

