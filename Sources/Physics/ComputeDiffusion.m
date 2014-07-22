function [ concentrations, debits, fluxSurfaciques, diffusionCoefficient ]  =  ComputeDiffusion(poreNetwork,cluster,inletLink,outletLink)
%COMPUTEDIFFUSION Resout un probleme de diffusion en statique dans un
%cluster
%Input : - network, cluster, inletLink, outletLink
%Output : [ concentrations, debits, fluxSurfaciques, diffusionCoefficient ]



    % Boundary Condition : 'Neumann' or 'Dirichlet'
    outletBoundaryCondition = 'Dirichlet';
    inletBoundaryCondition = 'Dirichlet';

    %Value of Boundary Condition
    concentrationInlet = 1;
    concentrationOutlet = 0.1;
    %debitSurfaciqueInlet = 0.11;
    
    
    %Verification si les conductances de diffusion sont deja calculees
    if not(isfield(poreNetwork.GetLinkDataList,'ConductancesDiffusion'))
        disp('Calcul des conductances de Diffusion...');
        tic;
        conductances = LocalScaleComputeConductancesDiffusion(poreNetwork,inletLink,outletLink);
        poreNetwork.AddNewLinkData(conductances,'ConductancesDiffusion');
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        disp(sprintf('Calcul des conductances de diffusion termine. Duree : %d minutes %f s.',minutes,secondes));
    else
        conductances = poreNetwork.GetLinkDataList.('ConductancesDiffusion');
    end
    
    nLink = poreNetwork.GetNumberOfLinks;
    nPore = poreNetwork.GetNumberOfPores;
    concentrations = zeros(1,nPore);
    fluxSurfaciques = zeros(1,nLink);    %flux oriente de owner a neighbour
    debits = zeros(1,nLink); %flux oriente de owner a neighbour   
    
    
    linkDiameters = poreNetwork.GetLinkDataList.Diameter;
    
    composantesConnexesPercolation = cluster.FindPercolationPath(inletLink,outletLink);
    
    for iComposanteConnexe = 1:length(composantesConnexesPercolation)
        percolatingCluster = composantesConnexesPercolation{iComposanteConnexe};
        
        %Extraction des liens envahis (liens envahis internes, inlet et
        %outlet)
        poresPercolants = percolatingCluster.GetInvadedPores;
        nPorePercolant = length(poresPercolants);

        liens_envahis = percolatingCluster.GetInvadedLinks;
        interfaceLinks = percolatingCluster.GetInterfaceLinks;
        liens_internes_envahis = setdiff(liens_envahis,interfaceLinks);
        liens_inlet_envahis = intersect(interfaceLinks,inletLink);
        liens_outlet_envahis = intersect(interfaceLinks,outletLink);
    
        poresPercolantsIndices=zeros(1,nPore);
        for i = 1:nPorePercolant
            poresPercolantsIndices( poresPercolants(i) ) = i;
        end
        
        
        %Remplissage matrice

        value_diag=zeros(1,nPorePercolant);
        
        indiceI_nonDiag=zeros(2*length(liens_internes_envahis),1);
        indiceJ_nonDiag=zeros(2*length(liens_internes_envahis),1);
        value_nonDiag=zeros(2*length(liens_internes_envahis),1);
        
        position=0;
        for numLien = liens_internes_envahis
            numOwner = poreNetwork.LinkOwners(numLien);
            numNeighbour = poreNetwork.LinkNeighbours(numLien);
            indiceOwner = poresPercolantsIndices(numOwner);
            assert(numNeighbour ~= -1);
            indiceNeighbour = poresPercolantsIndices(numNeighbour);
            
            %termes diagonaux
            %matrice(indiceOwner,indiceOwner) = matrice(indiceOwner,indiceOwner)+conductances(numLien);
            value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
            %matrice(indiceNeighbour,indiceNeighbour) = matrice(indiceNeighbour,indiceNeighbour)+conductances(numLien);
            value_diag(indiceNeighbour)=value_diag(indiceNeighbour)+conductances(numLien);
            
            %termes non diagonaux
            %matrice(indiceOwner,indiceNeighbour) = -conductances(numLien);
            indiceI_nonDiag(position+1)=indiceOwner;
            indiceJ_nonDiag(position+1)=indiceNeighbour;
            value_nonDiag(position+1)=-conductances(numLien);
            %matrice(indiceNeighbour,indiceOwner) = -conductances(numLien);
            indiceI_nonDiag(position+2)=indiceNeighbour;
            indiceJ_nonDiag(position+2)=indiceOwner;
            value_nonDiag(position+2)=-conductances(numLien);
            
            position=position+2;
        end

        
        if strcmp(inletBoundaryCondition,'Dirichlet')
            
          for numLien = liens_inlet_envahis
             numOwner = poreNetwork.LinkOwners(numLien);
             indiceOwner = poresPercolantsIndices(numOwner);
             %compl�ments aux termes diagonaux 
            %matrice(indiceOwner,indiceOwner) = matrice(indiceOwner,indiceOwner)+conductances(numLien);
            value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
          end
        end
        
        if strcmp(outletBoundaryCondition,'Dirichlet')
            
          for numLien = liens_outlet_envahis
             numOwner = poreNetwork.LinkOwners(numLien);
             indiceOwner = poresPercolantsIndices(numOwner);
             %compl�ments aux termes diagonaux 
            %matrice(indiceOwner,indiceOwner) = matrice(indiceOwner,indiceOwner)+conductances(numLien);
            value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
          end
        end
        
        value_diag=transpose(value_diag);
        indice_diag=transpose(1:nPorePercolant);
        
        I=vertcat(indiceI_nonDiag,indice_diag);
        J=vertcat(indiceJ_nonDiag,indice_diag);
        Value=vertcat(value_nonDiag,value_diag);
        matrice = sparse(I,J,Value,nPorePercolant,nPorePercolant);
        
        %Remplissage du terme de droite connu
        terme_droite = zeros(length(poresPercolants),1);

        if strcmp(inletBoundaryCondition,'Dirichlet')
          for numLien = liens_inlet_envahis
            
            numOwner = poreNetwork.LinkOwners(numLien);
            indiceOwner = poresPercolantsIndices(numOwner);
            assert(poreNetwork.LinkNeighbours(numLien) == -1);

            terme_droite(indiceOwner) = terme_droite(indiceOwner)+concentrationInlet*conductances(numLien);
          end

        elseif  strcmp(inletBoundaryCondition,'Neumann')
          for numLien = liens_inlet_envahis
              
              surface_face = pi/4*(linkDiameters(numLien))^2;
              numOwner = poreNetwork.LinkOwners(numLien);
              indiceOwner = poresPercolantsIndices(numOwner);
              assert(poreNetwork.LinkNeighbours(numLien) == -1);

              terme_droite(indiceOwner) = terme_droite(indiceOwner)+debitSurfaciqueInlet*surface_face;
          end
        end
        
        if strcmp(outletBoundaryCondition,'Dirichlet')
          for numLien = liens_outlet_envahis
            
            numOwner = poreNetwork.LinkOwners(numLien);
            indiceOwner = poresPercolantsIndices(numOwner);
            assert(poreNetwork.LinkNeighbours(numLien) == -1);

            terme_droite(indiceOwner) = terme_droite(indiceOwner)+concentrationOutlet*conductances(numLien);
          end

        elseif  strcmp(outletBoundaryCondition,'Neumann')
          for numLien = liens_outlet_envahis
              
              surface_face = pi/4*(linkDiameters(numLien))^2;
              numOwner = poreNetwork.LinkOwners(numLien);
              indiceOwner = poresPercolantsIndices(numOwner);
              assert(poreNetwork.LinkNeighbours(numLien) == -1);

              terme_droite(indiceOwner) = terme_droite(indiceOwner)-debitSurfaciqueOutlet*surface_face;
          end
        end       

        %Resolution du systeme lineaire
        disp('Solving linear system')
        tic;
        if poreNetwork.GetDimension==2
            concentr=mldivide(matrice,terme_droite);
        else
            L = ichol(matrice); %preconditionnement
          concentr = minres(matrice,terme_droite,1e-4,100,L,L');
        end
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        disp(sprintf('Solving linear system finished. Time spent : %d minutes %f s.',minutes,secondes));


        clear matrice;
        
        %Mise en forme de l'output concentration
        
        for i = 1:length(poresPercolants)
            concentrations(poresPercolants(i)) = concentr(i);
        end


        for numLien = liens_internes_envahis
            %Calcul des vitesses internes en fonction des concentrations
            numOwner = poreNetwork.LinkOwners(numLien);
            numNeighbour = poreNetwork.LinkNeighbours(numLien);
            debit = conductances(numLien)*(concentrations(numOwner)-concentrations(numNeighbour));
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;  
            debits(numLien) = debit;
            fluxSurfaciques(numLien) = debit/surface_lien_equivalente;   
        end
        
        %ajout des debits inlet
        if strcmp(inletBoundaryCondition,'Dirichlet')
            for numLien = liens_inlet_envahis
                numOwner = poreNetwork.LinkOwners(numLien);
                debit = conductances(numLien)*(concentrationInlet-concentrations(numOwner));
                surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
                debits(numLien) = debit;
                fluxSurfaciques(numLien) = debit/surface_lien_equivalente;    
            end
        
        elseif  strcmp(inletBoundaryCondition,'Neumann')
          for numLien = liens_inlet_envahis 
            
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
            debit = debitSurfaciqueInlet*surface_lien_equivalente;
            debits(numLien) = debit;
            fluxSurfaciques(numLien) = debit/surface_lien_equivalente;   
          end
        end
        
        %ajout des debits outlet
        if strcmp(outletBoundaryCondition,'Dirichlet')
            for numLien = liens_outlet_envahis
             
                numOwner = poreNetwork.LinkOwners(numLien);
                debit = conductances(numLien)*(concentrations(numOwner)-concentrationOutlet);
                surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
                debits(numLien) = debit;
                fluxSurfaciques(numLien) = debit/surface_lien_equivalente;    
            end
        
        elseif  strcmp(outletBoundaryCondition,'Neumann')
          for numLien = liens_outlet_envahis 
            
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
            debit = debitSurfaciqueInlet*surface_lien_equivalente;
            debits(numLien) = debit;
            fluxSurfaciques(numLien) = debit/surface_lien_equivalente;   
          end
        end


    end
    
    

    %Calcul de la diffusivit� effective

    liens_envahis = cluster.GetInvadedLinks;
    liens_inlet_envahis = intersect(liens_envahis,inletLink);
    liens_outlet_envahis = intersect(liens_envahis,outletLink);
    
    totalOutletDebit =  0;
    for iLink = liens_outlet_envahis
        totalOutletDebit = totalOutletDebit+debits(iLink);
    end
    totalInletDebit =  0;
    for iLink = liens_inlet_envahis
        totalInletDebit = totalInletDebit+debits(iLink);
    end

    
    if(abs(totalOutletDebit)>0)
        assert(abs(totalInletDebit-totalOutletDebit)/abs(totalOutletDebit)<1e-1,'Non conservation de la matière !');
    end
    
    deltaConcentration = concentrationInlet-concentrationOutlet;
    diffusionCoefficient = totalOutletDebit/deltaConcentration;   
    
end

