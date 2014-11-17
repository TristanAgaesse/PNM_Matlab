function [cluster,breakthroughPressure,invasionPressureList]  =  ComputeInvasionPercolation(network,inletLink,outletLink,wettability,varargin)
%ComputeInvasionPercolation Calcule l'invasion de percolation sur un r�seau
%de pores, avec un unique cluster
%Input : network,inletLink,outletLink,wettability   , ( varargin ) :
%       - network
%       - inletLink : liste des liens d'injection   
%       - outletLink : liste des liens de percée possible
%       - wettability : 'currentWettability', 'hydrophobic', 'hydrophilic' or 'random'.
%       - varargin (optionnel) : clusterOptions (voir ClusterMonophasique)
%
%Output : [cluster,breakthroughPressure,invasionPressureList]


    %V�rification si les diametres des liens sont d�j� calcul�s
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
        disp(sprintf('Calcul du diam�tre des liens termin�. Dur�e : %d minutes %f s.',minutes,secondes));
    end
    
    %Assignation du Contact Angle
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
    
    
    
    
    %Pr�paration de l'invasion : initialisation du cluster.
    disp('Calcul d''invasion percolation...');
    tic;
    
    clusterOptions = struct;
    if not(isempty(varargin))
        clusterOptions = varargin{1};        
    end
    cluster = ClusterMonophasique.InitialiseInvasionCluster(inletLink,outletLink,network,clusterOptions);
    
    %Boucle d'invasion pore par pore
    time = 0;
    currentPressure = 0;
    outlet_reached = false;
    outletPores = network.GetPoresFrontiere(outletLink);
    nPore = network.GetNumberOfPores;
    invasionPressureList = zeros(1,nPore);
    
    
    while not(outlet_reached) && time<nPore
        time = time+1;
        %trouver la face de plus petite pression critique
        
        [indexInvadedLink,invasionPressure] = cluster.GetMinimalPressureLink;
        
        if invasionPressure>currentPressure
            currentPressure = invasionPressure;
        end
        
        invasionPressureList(time) = invasionPressure;    
        
        invadedPore = cluster.GetOutwardPore(indexInvadedLink);
        
        %envahir le pore associ� et mettre � jour les pressions critiques
        
        interfaceChangeInformation = cluster.InvadeNewPore(indexInvadedLink);
        cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink); 
        
        %v�rifier si outlet_reached
        if ismember(invadedPore,outletPores)
            outlet_reached = true;
            breakthroughLinks = intersect(cluster.Network.GetLinksOfPore(invadedPore),outletLink);
            cluster.InvadeOutletLink(breakthroughLinks);
        end
    end
    breakthroughPressure = currentPressure;
    
    invasionPressureList = invasionPressureList(1:time);
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    disp(sprintf('Calcul d''invasion percolation termin�. Dur�e : %d minutes %f s.',minutes,secondes));
    
end