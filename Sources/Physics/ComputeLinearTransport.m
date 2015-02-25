function [ fieldValue, flux, fluxSurfaciques, effectiveTransportProperty ]  =  ComputeLinearTransport(network,transportPores,conductances,boundaryConditions)
%ComputeLinearTransport Resout un probleme de transport lineaire dans un
%ensemble de pores. La physique est specifiee par des conductances et des
%conditions limites.
%
%Input : network, transportPores, conductances,boundaryConditions
%       - network : pore network object
%       - transportPores : list of pores where there is transport
%       - conductances : conductances adaptees au transport lineaire
%       considere et a la geometrie ; length(conductances)=network.GetNumberOfLinks
%       - boundaryConditions=struct with the following fields
%       	boundaryConditions.inletLink = list of inletLinks
%       	boundaryConditions.outletLink = list of outletLinks
%       	boundaryConditions.inletType = 'Neumann' or 'Dirichlet'
%       	boundaryConditions.outletType = 'Neumann' or 'Dirichlet'
%       	boundaryConditions.inletValue = value of boundary condition 
%       	boundaryConditions.outletValue = value of boundary condition 
%
%Output : [ fieldValue, flux, fluxSurfaciques, effectiveTransportProperty ]

    disp('Computing linear transport ')
    tic;

    %Checking inputs and initial state
    [inletLink,outletLink]=ReadCheckInputs(network,transportPores,conductances,boundaryConditions);
    %CheckDiffusionConductances(network) ;

    
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    fieldValue = zeros(1,nPore);
    fluxSurfaciques = zeros(1,nLink);    %flux oriente de owner a neighbour
    flux = zeros(1,nLink); %flux oriente de owner a neighbour   
    
    %Decompose transportPores into percolation paths
    cluster = network.CreateVoidCluster;
    cluster.InvadedPores = transportPores;
    cluster.SetClusterOptions = struct;
    composantesConnexesPercolation = cluster.FindPercolationPath(inletLink,outletLink);
    
    %Resolution composante connexe par composante connexe
    for iComposanteConnexe = 1:length(composantesConnexesPercolation)
        
        
        %Extraction des liens envahis (liens envahis internes, inlet et
        %outlet)
        percolatingCluster = composantesConnexesPercolation{iComposanteConnexe};
        poresPercolants = percolatingCluster.GetInvadedPores;
        nPorePercolant = length(poresPercolants);

        liens_envahis = percolatingCluster.GetInvadedLinks;
        interfaceLinks = percolatingCluster.GetInterfaceLinks;
        liens_internes_envahis = setdiff(liens_envahis,interfaceLinks);
        liens_inlet_envahis = intersect(interfaceLinks,inletLink);
        liens_outlet_envahis = intersect(interfaceLinks,outletLink);
    
        poresPercolantsIndices=zeros(1,nPore);
        poresPercolantsIndices(poresPercolants)=1:nPorePercolant;
                
        
        %Remplissage matrice
        stiffnessMatrix = FillMatrix(network,conductances,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant);
        
        %Remplissage du terme de droite 
        rigthHandSide = FillRigthHandSide(network,conductances,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant);
        

        %Resolution du systeme lineaire

        if network.GetDimension==2
            concentr=mldivide(stiffnessMatrix,rigthHandSide);
        else
            L = ichol(stiffnessMatrix); %preconditionnement
          concentr = minres(stiffnessMatrix,rigthHandSide,1e-4,100,L,L');
        end
        
        
        
        clear matrice;
        
        %Reassemblage de l'output concentration
        [fieldValue,flux,fluxSurfaciques]=ReassembleComputedField(fieldValue,flux,fluxSurfaciques,network,poresPercolants,concentr,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions);
        
    end
    

    %Check computation with a mass balance
    
    totalInletDebit = sum(flux(inletLink));
    totalOutletDebit = sum(flux(outletLink));
    
    if(abs(totalOutletDebit)>0)
        assert(abs(totalInletDebit-totalOutletDebit)/abs(totalOutletDebit)<1e-1,'Non conservation de la matière !');
    end
       
    
    effectiveTransportProperty = ComputeDiffusionCoefficient(fieldValue,flux,inletLink,outletLink);
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Computing linear transport finished. Time spent : %d minutes %f s. ',minutes,secondes);
    
end



%---------------------------------------------------------------------------------------------
% function CheckDiffusionConductances(network)
%     
%     if not(isfield(network.GetLinkDataList,'ConductancesDiffusion'))
%         disp('Calcul des conductances de Diffusion...');
%         tic;
%         conductances = LocalScaleComputeConductancesDiffusion(network);
%         network.AddNewLinkData(conductances,'ConductancesDiffusion');
%         duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
%         fprintf('Calcul des conductances de diffusion termine. Duree : %d minutes %f s.',minutes,secondes);
%    
%     end
% end



%---------------------------------------------------------------------------------------------
function matrice = FillMatrix(network,conductances,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant)

    
    value_diag=zeros(1,nPorePercolant);

    indiceI_nonDiag=zeros(2*length(liens_internes_envahis),1);
    indiceJ_nonDiag=zeros(2*length(liens_internes_envahis),1);
    value_nonDiag=zeros(2*length(liens_internes_envahis),1);

    
    
%     %Contribution des liens internes
    numOwner = network.LinkOwners(liens_internes_envahis);
    numNeighbour = network.LinkNeighbours(liens_internes_envahis);
    indiceOwner = poresPercolantsIndices(numOwner);
    indiceNeighbour = poresPercolantsIndices(numNeighbour);
    %termes diagonaux
    for i=1:length(liens_internes_envahis)
         value_diag(indiceOwner(i))=value_diag(indiceOwner(i))+conductances(liens_internes_envahis(i));
         value_diag(indiceNeighbour(i))=value_diag(indiceNeighbour(i))+conductances(liens_internes_envahis(i));
    end
    %termes non diagonaux
    matrixIndex=1:2:2*length(liens_internes_envahis);
    indiceI_nonDiag(matrixIndex)=indiceOwner;
    indiceJ_nonDiag(matrixIndex)=indiceNeighbour;
    value_nonDiag(matrixIndex)=-conductances(liens_internes_envahis);

    matrixIndex=2:2:2*length(liens_internes_envahis);
    indiceI_nonDiag(matrixIndex)=indiceNeighbour;
    indiceJ_nonDiag(matrixIndex)=indiceOwner;
    value_nonDiag(matrixIndex)=-conductances(liens_internes_envahis);
        

    %Contribution des liens inlet et outlet
    if strcmp(boundaryConditions.inletType,'Dirichlet')

      for numLien = liens_inlet_envahis
         numOwner = network.LinkOwners(numLien);
         indiceOwner = poresPercolantsIndices(numOwner);
         %compl�ments aux termes diagonaux 
        value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
      end
    end

    if strcmp(boundaryConditions.outletType,'Dirichlet')

      for numLien = liens_outlet_envahis
         numOwner = network.LinkOwners(numLien);
         indiceOwner = poresPercolantsIndices(numOwner);
         %compl�ments aux termes diagonaux 
        value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
      end
    end

    value_diag=transpose(value_diag);
    indice_diag=transpose(1:nPorePercolant);

    I=vertcat(indiceI_nonDiag,indice_diag);
    J=vertcat(indiceJ_nonDiag,indice_diag);
    Value=vertcat(value_nonDiag,value_diag);
    matrice = sparse(I,J,Value,nPorePercolant,nPorePercolant);

end



%---------------------------------------------------------------------------------------------
function terme_droite = FillRigthHandSide(network,conductances,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant)
    terme_droite = zeros(1,nPorePercolant);
    
    linkDiameters = network.GetLinkDataList.Diameter;
    conductances = network.GetLinkDataList.('ConductancesDiffusion');
    
    if strcmp(boundaryConditions.inletType,'Dirichlet')
        
    	indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_inlet_envahis));
        for i=1:length(liens_inlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))+boundaryConditions.inletValue*conductances(liens_inlet_envahis(i));
        end
        
        
    elseif  strcmp(boundaryConditions.inletType,'Neumann')
        
        surface_face = pi/4*(linkDiameters(liens_inlet_envahis)).^2;
        indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_inlet_envahis));
        for i =1:length(liens_inlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))+boundaryConditions.inletValue*surface_face(i);
        end
    end

    if strcmp(boundaryConditions.outletType,'Dirichlet')
        
        indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_outlet_envahis));
        for i =1:length(liens_outlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))+boundaryConditions.outletValue*conductances(liens_outlet_envahis(i));
        end
        
    elseif  strcmp(boundaryConditions.outletType,'Neumann')
        
        surface_face = pi/4*(linkDiameters(liens_outlet_envahis)).^2;
        indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_outlet_envahis));
        for i =1:length(liens_outlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))-boundaryConditions.outletValue*surface_face(i);        
        end
    end
    
    terme_droite=transpose(terme_droite);
end



%---------------------------------------------------------------------------------------------
function [fieldValue,flux,fluxSurfaciques]=ReassembleComputedField(fieldValue,flux,fluxSurfaciques,network,poresPercolants,concentr,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions)

    for i = 1:length(poresPercolants)
        fieldValue(poresPercolants(i)) = concentr(i);
    end

    linkDiameters = network.GetLinkDataList.Diameter;
    conductances = network.GetLinkDataList.('ConductancesDiffusion');
    
    %Calcul des vitesses internes en fonction des fieldValue
    numOwner = network.LinkOwners(liens_internes_envahis);
    numNeighbour = network.LinkNeighbours(liens_internes_envahis);
    debit(liens_internes_envahis) = conductances(liens_internes_envahis).*(fieldValue(numOwner)-fieldValue(numNeighbour));
    surface_lien_equivalente = pi/4*(linkDiameters(liens_internes_envahis)).^2;
    fluxSurfaciques(liens_internes_envahis) = debit(liens_internes_envahis)./surface_lien_equivalente;
    
    %ajout des flux inlet
    if strcmp(boundaryConditions.inletType,'Dirichlet')
        for numLien = liens_inlet_envahis
            numOwner = network.LinkOwners(numLien);
            debit = conductances(numLien)*(boundaryConditions.inletValue-fieldValue(numOwner));
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
            flux(numLien) = debit;
            fluxSurfaciques(numLien) = debit/surface_lien_equivalente;    
        end

    elseif  strcmp(boundaryConditions.inletType,'Neumann')
      for numLien = liens_inlet_envahis 

        surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
        debit = boundaryConditions.inletValue*surface_lien_equivalente;
        flux(numLien) = debit;
        fluxSurfaciques(numLien) = debit/surface_lien_equivalente;   
      end
    end

    %ajout des flux outlet
    if strcmp(boundaryConditions.outletType,'Dirichlet')
        for numLien = liens_outlet_envahis

            numOwner = network.LinkOwners(numLien);
            debit = conductances(numLien)*(fieldValue(numOwner)-boundaryConditions.outletValue);
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
            flux(numLien) = debit;
            fluxSurfaciques(numLien) = debit/surface_lien_equivalente;    
        end

    elseif  strcmp(boundaryConditions.outletType,'Neumann')
      for numLien = liens_outlet_envahis 

        surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
        debit = boundaryConditions.outletValue*surface_lien_equivalente;
        flux(numLien) = debit;
        fluxSurfaciques(numLien) = debit/surface_lien_equivalente;   
      end
    end

end



%---------------------------------------------------------------------------------------------
function effectiveTransportProperty = ComputeDiffusionCoefficient(fieldValue,flux,inletLink,outletLink)

    totalOutletDebit = sum(flux(outletLink));
    deltaConcentration = mean(fieldValue(inletLink))-mean(fieldValue(outletLink));
    
    effectiveTransportProperty = totalOutletDebit/deltaConcentration;   
end



%---------------------------------------------------------------------------------------------
function [inletLink,outletLink]=ReadCheckInputs(network,transportPores,conductances,boundaryConditions)

    assert( isa(network,'PoreNetwork'),'LinearTransport : first input must be a PoreNetwork object')
    assert( length(transportPores)<=network.GetNumberOfPore,'LinearTransport : second input transportPores must be of length <=network.GetNumberOfPore')
    assert( length(conductances)==network.GetNumberOfLinks,'LinearTransport : third input conductance must be of length network.GetNumberOfLinks')

    %bcCheck1 = isa(boundaryConditions,'struct') ;
    %inletLink 
    
    
%       	boundaryConditions.outletLink = list of outletLinks
%       	boundaryConditions.inletType = 'Neumann' or 'Dirichlet'
%       	boundaryConditions.outletType = 'Neumann' or 'Dirichlet'
%       	boundaryConditions.inletValue = value of boundary condition
%       	boundaryConditions.outletValue 
     
         
     inletLink = boundaryConditions.inletLink;
     outletLink = boundaryConditions.outletLink;

end



