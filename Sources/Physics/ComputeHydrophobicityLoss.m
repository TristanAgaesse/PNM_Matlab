function floodingStepInformation=ComputeHydrophobicityLoss(network,inletLink,outletLink,options,varargin)
%COMPUTEHYDROPHOBICITYLOSS High level function which compute how a
%cluster is changed when there is a hydrophobicity loss which rate is 
%correlated to the presence of water
%
%Input : network,inletLink,ouletLink,options,varargin (optionnel)
%       options.nIterations= nombre d'itérations
%       options.MechanismeDegradation='uniforme','uniformeDansEau','sommeVitesses'
%       options.RechercheNextInvadedLink='localLinearisationOfCapillaryPressure','linearDecreaseOfCapillaryPressure'
%       varargin (optionnel) :clusterOptions (voir ClusterMonophasique)
%Output : floodingStepInformation
    
    
    nIterations=options.nIterations;
    
    floodingStepInformation=struct;
    floodingStepInformation.times=zeros(1,nIterations+1);
    floodingStepInformation.clusters=cell(1,nIterations+1);
    floodingStepInformation.invasionPressures=cell(1,nIterations+1);
    floodingStepInformation.pourcentageDegradation=cell(1,nIterations+1);
    floodingStepInformation.fixedPressure=0;
    
    floodingStepInformation.inletLink=inletLink;
    floodingStepInformation.outletLink=outletLink;
    floodingStepInformation.options=options;
    
    %Invasion initiale du r�seau
    temps=0;
    clusterOptions=struct;
    if not(isempty(varargin))
        clusterOptions=varargin{1};        
    end
    [cluster,breakthroughPressure,invasionPressureList] = ComputeInvasionPercolation(network,inletLink,outletLink,'hydrophobic',clusterOptions);
    
    floodingStepInformation.times(1)=0;
    floodingStepInformation.clusters{1}=cluster.CopyCluster;
    floodingStepInformation.invasionPressures{1}=invasionPressureList;
    
    pression_reference=0.9*breakthroughPressure;
    floodingStepInformation.fixedPressure=pression_reference;
    
    %Initialisation des param�tres de d�gradation
    pourcentageDegradation=zeros(1,network.GetNumberOfLinks);
    floodingStepInformation.pourcentageDegradation{1}=pourcentageDegradation;
    
    vitesseDegradation=ComputeVitesseDegradation(network,cluster,inletLink,outletLink,options);
    
    disp('Envahissement secondaire sous l''effet de la perte d''hydrophobie')

    for iIteration=1:nIterations
        disp(iIteration);
        floodingStepInformation.nIteration=iIteration;
        try
            [indexInvadedLink,timeStep]=FindNextInvadedLink(cluster,pression_reference,pourcentageDegradation,vitesseDegradation,inletLink,outletLink,options);
        catch err
            if (strcmp(err.identifier,'FindNextInvadedLink:EmptyTempsPotentielValable'))
                %error('MATLAB:myCode:dimensions', err.message);
                disp(err.message)
                return
                % Display any other errors as usual.
            else
                rethrow(err);
            end

        end
        temps=temps+timeStep;
        
        %en cas d'envahissement :
        if indexInvadedLink>0
            %mise � jour localisation de l'eau :
            interfaceChangeInformation=cluster.InvadeNewPore(indexInvadedLink);
            cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink);
            
            floodingStepInformation.invasionPressures{iIteration+1}=pression_reference;
            %envahissement �ventuel des pores voisins
            
            [indexMinPressureLink,minPressure]=cluster.GetMinimalPressureLink;
            while minPressure<pression_reference
                
                minPressureLink=cluster.InterfaceLinks(indexMinPressureLink);
                if network.GetFrontiereOfLink(minPressureLink)~=0
                    disp(sprintf('New breakthrough point on boundary %d',network.GetFrontiereOfLink(minPressureLink)));
                    cluster.InvadeOutletLink(minPressureLink);
                else
                    interfaceChangeInformation=cluster.InvadeNewPore(indexMinPressureLink);
                    cluster.UpdateCriticalPressure(interfaceChangeInformation,inletLink,outletLink);
                end
                
                floodingStepInformation.invasionPressures{iIteration+1}=[floodingStepInformation.invasionPressures{iIteration+1},minPressure];
                
                [indexMinPressureLink,minPressure]=cluster.GetMinimalPressureLink;
            end
            
            floodingStepInformation.times(iIteration+1)=temps;
            floodingStepInformation.clusters{iIteration+1}=cluster.CopyCluster;
            floodingStepInformation.pourcentageDegradation{iIteration+1}=pourcentageDegradation;
            
            %mise � jour contactAngle et parametres degradation
            pourcentageDegradation=pourcentageDegradation+timeStep*vitesseDegradation;
            pourcentageDegradation=min(pourcentageDegradation,100*ones(1,length(pourcentageDegradation)));
            
            deltaContactAngle=30/180*pi;
            contactAngle=cluster.Network.GetLinkDataList.ContactAngle;
            cluster.Network.RemoveLinkData('ContactAngle');
            
            if not(isfield(cluster.Network.GetLinkDataList,'InitialContactAngle'))
                cluster.Network.AddNewLinkData(contactAngle,'InitialContactAngle');
            end
            
            initialContactAngle=cluster.Network.GetLinkDataList.InitialContactAngle;
            
            contactAngle=initialContactAngle-pourcentageDegradation/100*deltaContactAngle;
            cluster.Network.AddNewLinkData(contactAngle,'ContactAngle');
            
            
            
            vitesseDegradation=ComputeVitesseDegradation(network,cluster,inletLink,outletLink,options);
        end
        
    end
    
    
    
    function [indexInvadedLink,temps]=FindNextInvadedLink(cluster,pression_reference,pourcentageDegradation,vitesseDegradation,linkInlet,linkOutlet,options)
        
        %option='localLinearisationOfCapillaryPressure';
        methode=options.RechercheNextInvadedLink;
        
        criticalPressures=cluster.GetCriticalPressures;
        
        if strcmp(methode,'localLinearisationOfCapillaryPressure')
            %ici on fait varier explicitement les angles de contact. 
            %On recalcule a chaque pas de temps les
            %pressions critiques en fonction des nouveaux angles de
            %contact. Pour trouver le prochain pas de temps, on linearise
            %la pression localement en fonction de l'angle de contact
            
            
            
            interfaceUpdateInformation=cluster.GetInterfaceChangeInformation(1:length(cluster.GetInterfaceLinks));  %tous les liens frontiere
            UpdateCriticalPressure(cluster,interfaceUpdateInformation,linkInlet,linkOutlet)
            
            temps_invasion_potentielle=zeros(1,cluster.Network.GetNumberOfLinks);
            foo=1-(pression_reference*ones(1,cluster.Network.GetNumberOfLinks))./(criticalPressures.*tan(cluster.Network.GetLinkDataList.ContactAngle));
            
            for iLink=cluster.GetInterfaceLinks
                if vitesseDegradation(iLink)~=0 && cluster.Network.GetFrontiereOfLink(iLink)==0
                    temps_invasion_potentielle(iLink)=foo(iLink)/vitesseDegradation(iLink);
                end
            end
            
            temps_potentiels_valables=temps_invasion_potentielle(and(temps_invasion_potentielle>0,pourcentageDegradation<100));
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
            foo=200*(1-((pression_reference*ones(1,cluster.Network.GetNumberOfLinks))./criticalPressures))-pourcentageDegradation;

            for iLink=cluster.GetInterfaceLinks
                if vitesseDegradation(iLink)~=0 && cluster.Network.GetFrontiereOfLink(iLink)==0
                    temps_invasion_potentielle(iLink)=foo(iLink)/vitesseDegradation(iLink);
                end
            end

            temps_potentiels_valables=temps_invasion_potentielle(and(temps_invasion_potentielle>0,pourcentageDegradation<100));

            temps=min(temps_potentiels_valables);
            ind=find(temps_invasion_potentielle==temps);
            indexInvadedLink=find(cluster.GetInterfaceLinks==ind);
        end
    end



    function vitesseDegradation=ComputeVitesseDegradation(network,clusterLiquide,inletLink,ouletLink,options)
        
        vitesseDegradation=zeros(1,network.GetNumberOfLinks);
        mechanismeDegradation=options.MechanismeDegradation;
        
        if strcmp(mechanismeDegradation,'uniforme')
            vitesseDegradation=ones(1,network.GetNumberOfLinks);
        
        elseif strcmp(mechanismeDegradation,'uniformeDansEau')
            vitesseDegradation(clusterLiquide.GetInvadedLinks)=1;
            vitesseDegradation(clusterLiquide.GetInterfaceLinks)=1;
        
        elseif strcmp(mechanismeDegradation,'sommeVitesses')
            %coeff degradation(face)=moyenne des vitesse sur les faces des
            %pores_voisins
            [ ~ , ~, fluidVelocity, ~ ] = ComputePermeability(network,clusterLiquide,inletLink,ouletLink);
            
            nPore=network.GetNumberOfPores;
            vitesseMoyennePore=zeros(1,nPore);
            for iPore=1:nPore
                links=network.GetLinksOfPore(iPore);
                vitesseMoyennePore(iPore)=sum(arrayfun(@(x) abs(x),fluidVelocity(links)))/length(links);
            end
            
            for iLink=1:network.GetNumberOfLinks
                pores=network.GetPoresOfLink(iLink);
                if pores(1)==-1
                    vitesseDegradation(iLink)=vitesseMoyennePore(pores(2))/mean(vitesseMoyennePore(not(isnan(vitesseMoyennePore))));
                else
                    vitesseDegradation(iLink)=0.5*(vitesseMoyennePore(pores(1))+vitesseMoyennePore(pores(2)))/mean(vitesseMoyennePore(not(isnan(vitesseMoyennePore))));
                end
            end
            vitesseDegradation(isnan(vitesseDegradation))=0;
            
        elseif strcmp(mechanismeDegradation,'sommeDebits')
%             for num_face=cluster.GetInterfaceLinks
%                faces_of_owner=network.GetLinksOfPore(network.LinkOwners(num_face));
%                 moyenne_vitesses_owner=sum(arrayfun(@(x) abs(x),fluidVelocity(faces_of_owner)))/length(faces_of_owner);
% 
%                 if network.LinkNeighbours(num_face)~=-1
%                     faces_of_neighbour=network.GetLinksOfPore(network.LinkNeighbours(num_face));
%                     moyenne_vitesses_owner=sum(arrayfun(@(x) abs(x),fluidVelocity(faces_of_neighbour)))/length(faces_of_neighbour);
%                 else
%                     moyenne_vitesse_neighbour=0;
%                 end
%                 vitesseDegradation(num_face)=0.5*(moyenne_vitesses_owner+moyenne_vitesse_neighbour);
%             end
        end
    end

end