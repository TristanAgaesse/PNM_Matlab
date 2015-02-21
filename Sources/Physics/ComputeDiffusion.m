function [ concentrations, debits, fluxSurfaciques, diffusionCoefficient ]  =  ComputeDiffusion(poreNetwork,cluster,boundaryConditions)
%COMPUTEDIFFUSION Resout un probleme de diffusion en statique dans un
%cluster
%Input : - network, cluster, boundaryConditions
%       - boundaryConditions=struct with the following (facultative) fields
%       - boundaryConditions.inletLink = list of inletLinks
%       - boundaryConditions.outletLink = list of outletLinks
%       - boundaryConditions.inletType = 'Neumann' or 'Dirichlet'
%       - boundaryConditions.outletType = 'Neumann' or 'Dirichlet'
%       - boundaryConditions.inletValue = value of boundary condition
%       - boundaryConditions.outletValue = value of boundary condition
%Output : [ concentrations, debits, fluxSurfaciques, diffusionCoefficient ]


    %Checking inputs and initial state
    [inletLink,outletLink]=ReadCheckInputs(poreNetwork,cluster,boundaryConditions);
    CheckDiffusionConductances(poreNetwork) ;

    
    
    nLink = poreNetwork.GetNumberOfLinks;
    nPore = poreNetwork.GetNumberOfPores;
    concentrations = zeros(1,nPore);
    fluxSurfaciques = zeros(1,nLink);    %flux oriente de owner a neighbour
    debits = zeros(1,nLink); %flux oriente de owner a neighbour   
    
    
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
        stiffnessMatrix = FillMatrix(poreNetwork,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant);
        
        %Remplissage du terme de droite 
        rigthHandSide = FillRigthHandSide(poreNetwork,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant);
        

        %Resolution du systeme lineaire
        disp('Solving linear system')
        tic;
        
        if poreNetwork.GetDimension==2
            concentr=mldivide(stiffnessMatrix,rigthHandSide);
        else
            L = ichol(stiffnessMatrix); %preconditionnement
          concentr = minres(stiffnessMatrix,rigthHandSide,1e-4,100,L,L');
        end
        
        
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        fprintf('Solving linear system finished. Time spent : %d minutes %f s.',minutes,secondes);
        clear matrice;
        
        %Reassemblage de l'output concentration
        [concentrations,debits,fluxSurfaciques]=ReassembleComputedField(concentrations,debits,fluxSurfaciques,poreNetwork,poresPercolants,concentr,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions);
        
    end
    

    %Check computation with a mass balance
    
    totalInletDebit = sum(debits(inletLink));
    totalOutletDebit = sum(debits(outletLink));
    
    if(abs(totalOutletDebit)>0)
        %assert(abs(totalInletDebit-totalOutletDebit)/abs(totalOutletDebit)<1e-1,'Non conservation de la matière !');
    end
       
    
    diffusionCoefficient = ComputeDiffusionCoefficient(concentrations,debits,inletLink,outletLink);
    
end



%---------------------------------------------------------------------------------------------
function CheckDiffusionConductances(poreNetwork)
    
    if not(isfield(poreNetwork.GetLinkDataList,'ConductancesDiffusion'))
        disp('Calcul des conductances de Diffusion...');
        tic;
        conductances = LocalScaleComputeConductancesDiffusion(poreNetwork);
        poreNetwork.AddNewLinkData(conductances,'ConductancesDiffusion');
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        fprintf('Calcul des conductances de diffusion termine. Duree : %d minutes %f s.',minutes,secondes);
   
    end
end



%---------------------------------------------------------------------------------------------
function matrice = FillMatrix(poreNetwork,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant)

    conductances = poreNetwork.GetLinkDataList.('ConductancesDiffusion');
    
    value_diag=zeros(1,nPorePercolant);

    indiceI_nonDiag=zeros(2*length(liens_internes_envahis),1);
    indiceJ_nonDiag=zeros(2*length(liens_internes_envahis),1);
    value_nonDiag=zeros(2*length(liens_internes_envahis),1);

    
    
    %Contribution des liens internes
    numOwner = poreNetwork.LinkOwners(liens_internes_envahis);
    numNeighbour = poreNetwork.LinkNeighbours(liens_internes_envahis);
    indiceOwner = poresPercolantsIndices(numOwner);
    indiceNeighbour = poresPercolantsIndices(numNeighbour);
    %termes diagonaux
    value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(liens_internes_envahis);
    value_diag(indiceNeighbour)=value_diag(indiceNeighbour)+conductances(liens_internes_envahis);
    %termes non diagonaux
    position=1:2:2*length(liens_internes_envahis);
    indiceI_nonDiag(position)=indiceOwner;
    indiceJ_nonDiag(position)=indiceNeighbour;
    value_nonDiag(position)=-conductances(liens_internes_envahis);

    position=2:2:2*length(liens_internes_envahis);
    indiceI_nonDiag(position)=indiceNeighbour;
    indiceJ_nonDiag(position)=indiceOwner;
    value_nonDiag(position)=-conductances(liens_internes_envahis);
        
    %position=0;    
%     for numLien = liens_internes_envahis
%         numOwner = poreNetwork.LinkOwners(numLien);
%         numNeighbour = poreNetwork.LinkNeighbours(numLien);
%         indiceOwner = poresPercolantsIndices(numOwner);
%         assert(numNeighbour ~= -1);
%         indiceNeighbour = poresPercolantsIndices(numNeighbour);
% 
%         %termes diagonaux
%         value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
%         value_diag(indiceNeighbour)=value_diag(indiceNeighbour)+conductances(numLien);
% 
%         %termes non diagonaux
% 
%         indiceI_nonDiag(position+1)=indiceOwner;
%         indiceJ_nonDiag(position+1)=indiceNeighbour;
%         value_nonDiag(position+1)=-conductances(numLien);
% 
%         indiceI_nonDiag(position+2)=indiceNeighbour;
%         indiceJ_nonDiag(position+2)=indiceOwner;
%         value_nonDiag(position+2)=-conductances(numLien);
% 
%         position=position+2;
%     end

    %Contribution des liens inlet et outlet
    if strcmp(boundaryConditions.inletType,'Dirichlet')

      for numLien = liens_inlet_envahis
         numOwner = poreNetwork.LinkOwners(numLien);
         indiceOwner = poresPercolantsIndices(numOwner);
         %compl�ments aux termes diagonaux 
        value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
      end
    end

    if strcmp(boundaryConditions.outletType,'Dirichlet')

      for numLien = liens_outlet_envahis
         numOwner = poreNetwork.LinkOwners(numLien);
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
function terme_droite = FillRigthHandSide(poreNetwork,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions,poresPercolantsIndices,nPorePercolant)
    terme_droite = zeros(1,nPorePercolant);
    
    linkDiameters = poreNetwork.GetLinkDataList.Diameter;
    conductances = poreNetwork.GetLinkDataList.('ConductancesDiffusion');
    
    if strcmp(boundaryConditions.inletType,'Dirichlet')
        
      indiceOwner = poresPercolantsIndices(poreNetwork.LinkOwners(liens_inlet_envahis));
      terme_droite(indiceOwner) = terme_droite(indiceOwner)+boundaryConditions.inletValue*conductances(liens_inlet_envahis);

    elseif  strcmp(boundaryConditions.inletType,'Neumann')
        
        surface_face = pi/4*(linkDiameters(liens_inlet_envahis)).^2;
        indiceOwner = poresPercolantsIndices(poreNetwork.LinkOwners(liens_inlet_envahis));
        terme_droite(indiceOwner) = terme_droite(indiceOwner)+boundaryConditions.inletValue*surface_face;
        
    end

    if strcmp(boundaryConditions.outletType,'Dirichlet')
        
        indiceOwner = poresPercolantsIndices(poreNetwork.LinkOwners(liens_outlet_envahis));
        terme_droite(indiceOwner) = terme_droite(indiceOwner)+boundaryConditions.outletValue*conductances(liens_outlet_envahis);
        
    elseif  strcmp(boundaryConditions.outletType,'Neumann')
        
        surface_face = pi/4*(linkDiameters(liens_outlet_envahis)).^2;
        indiceOwner = poresPercolantsIndices(poreNetwork.LinkOwners(liens_outlet_envahis));
        terme_droite(indiceOwner) = terme_droite(indiceOwner)-boundaryConditions.outletValue*surface_face;        
    
    end
    
    terme_droite=transpose(terme_droite);
end



%---------------------------------------------------------------------------------------------
function [concentrations,debits,fluxSurfaciques]=ReassembleComputedField(concentrations,debits,fluxSurfaciques,poreNetwork,poresPercolants,concentr,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,boundaryConditions)

    for i = 1:length(poresPercolants)
        concentrations(poresPercolants(i)) = concentr(i);
    end

    linkDiameters = poreNetwork.GetLinkDataList.Diameter;
    conductances = poreNetwork.GetLinkDataList.('ConductancesDiffusion');
    
    %Calcul des vitesses internes en fonction des concentrations
    numOwner = poreNetwork.LinkOwners(liens_internes_envahis);
    numNeighbour = poreNetwork.LinkNeighbours(liens_internes_envahis);
    debit(liens_internes_envahis) = conductances(liens_internes_envahis).*(concentrations(numOwner)-concentrations(numNeighbour));
    surface_lien_equivalente = pi/4*(linkDiameters(liens_internes_envahis)).^2;
    fluxSurfaciques(liens_internes_envahis) = debit(liens_internes_envahis)./surface_lien_equivalente;
    
    %ajout des debits inlet
    if strcmp(boundaryConditions.inletType,'Dirichlet')
        for numLien = liens_inlet_envahis
            numOwner = poreNetwork.LinkOwners(numLien);
            debit = conductances(numLien)*(boundaryConditions.inletValue-concentrations(numOwner));
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
            debits(numLien) = debit;
            fluxSurfaciques(numLien) = debit/surface_lien_equivalente;    
        end

    elseif  strcmp(boundaryConditions.inletType,'Neumann')
      for numLien = liens_inlet_envahis 

        surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
        debit = boundaryConditions.inletValue*surface_lien_equivalente;
        debits(numLien) = debit;
        fluxSurfaciques(numLien) = debit/surface_lien_equivalente;   
      end
    end

    %ajout des debits outlet
    if strcmp(boundaryConditions.outletType,'Dirichlet')
        for numLien = liens_outlet_envahis

            numOwner = poreNetwork.LinkOwners(numLien);
            debit = conductances(numLien)*(concentrations(numOwner)-boundaryConditions.outletValue);
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
            debits(numLien) = debit;
            fluxSurfaciques(numLien) = debit/surface_lien_equivalente;    
        end

    elseif  strcmp(boundaryConditions.outletType,'Neumann')
      for numLien = liens_outlet_envahis 

        surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
        debit = boundaryConditions.outletValue*surface_lien_equivalente;
        debits(numLien) = debit;
        fluxSurfaciques(numLien) = debit/surface_lien_equivalente;   
      end
    end

end



%---------------------------------------------------------------------------------------------
function diffusionCoefficient = ComputeDiffusionCoefficient(concentrations,debits,inletLink,outletLink)

    totalOutletDebit = sum(debits(outletLink));
    deltaConcentration = mean(concentrations(inletLink))-mean(concentrations(outletLink));
    
    diffusionCoefficient = totalOutletDebit/deltaConcentration;   
end



%---------------------------------------------------------------------------------------------
function [inletLink,outletLink]=ReadCheckInputs(poreNetwork,cluster,boundaryConditions)


     inletLink = boundaryConditions.inletLink;
     outletLink = boundaryConditions.outletLink;
%     
%     boundaryConditions.inletType
%     boundaryConditions.outletType
%     
%     boundaryConditions.inletValue
%     boundaryConditions.outletValue
%     
%     
%     % Boundary Condition : 'Neumann' or 'Dirichlet'
%     outletBoundaryCondition = 'Dirichlet';
%     inletBoundaryCondition = 'Dirichlet';
% 
%     %Value of Boundary Condition
%     concentrationInlet = 1;
%     concentrationOutlet = 0.1;
%     %debitSurfaciqueInlet = 0.11;

end



