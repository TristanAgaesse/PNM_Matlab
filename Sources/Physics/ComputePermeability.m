function [ pressions, debits, vitessesMoyennes, permeabilityCoefficient ]  =  ComputePermeability(poreNetwork,cluster,inletLink,outletLink)
%ComputePermeability R�sout les �quations de Stokes dans un cluster
%Input :   network, cluster, inletLink, outletLink
%Output : [ pressions, debits, vitessesMoyennes, permeabilityCoefficient ]


    % Boundary Condition : 'Neumann' or 'Dirichlet'
    %inletBoundaryCondition = 'Dirichlet';
    inletBoundaryCondition = 'Neumann';
    outletBoundaryCondition = 'Dirichlet';
    %outletBoundaryCondition = 'Neumann';
    
    %Value of Boundary Condition
    %pressionInlet = 10e5;
    debitSurfaciqueInlet = 1e-3;
    pressionOutlet = 1e5;
    %debitSurfaciqueOutlet = 1e-11;
    
    
    
    
    %V�rification si les conductances de Stokes sont d�j� calcul�s
    if not(isfield(poreNetwork.GetLinkDataList,'ConductancesStokes'))
        disp('Calcul des conductances Stokes...');
        tic;
        conductances = LocalScaleComputeConductancesStokes(poreNetwork,inletLink,outletLink);
        poreNetwork.AddNewLinkData(conductances,'ConductancesStokes');
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        disp(sprintf('Calcul des conductances Stokes termine. Duree : %d minutes %f s.',minutes,secondes));
    else
        conductances = poreNetwork.GetLinkDataList.('ConductancesStokes');
    end
    
    nLink = poreNetwork.GetNumberOfLinks;
    nPore = poreNetwork.GetNumberOfPores;
    pressions = zeros(1,nPore);
    vitessesMoyennes = zeros(1,nLink);    %vitesse orient�e de owner � neighbour
    debits = zeros(1,nLink); %flux orient� de owner � neighbour
    
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
            indiceOwner = find(poresPercolants == numOwner);
            assert(poreNetwork.LinkNeighbours(numLien) == -1);

            terme_droite(indiceOwner) = terme_droite(indiceOwner)+pressionInlet*conductances(numLien);
          end

        elseif  strcmp(inletBoundaryCondition,'Neumann')
          for numLien = liens_inlet_envahis
              
              surface_face = pi/4*(linkDiameters(numLien))^2;
              numOwner = poreNetwork.LinkOwners(numLien);
              indiceOwner = find(poresPercolants == numOwner);
              assert(poreNetwork.LinkNeighbours(numLien) == -1);

              terme_droite(indiceOwner) = terme_droite(indiceOwner)+debitSurfaciqueInlet*surface_face;
          end
        end
    
        if strcmp(outletBoundaryCondition,'Dirichlet')
          for numLien = liens_outlet_envahis
            
            numOwner = poreNetwork.LinkOwners(numLien);
            indiceOwner = find(poresPercolants == numOwner);
            assert(poreNetwork.LinkNeighbours(numLien) == -1);

            terme_droite(indiceOwner) = terme_droite(indiceOwner)+pressionOutlet*conductances(numLien);
          end

        elseif  strcmp(outletBoundaryCondition,'Neumann')
          for numLien = liens_outlet_envahis
              
              surface_face = pi/4*(linkDiameters(numLien))^2;
              numOwner = poreNetwork.LinkOwners(numLien);
              indiceOwner = find(poresPercolants == numOwner);
              assert(poreNetwork.LinkNeighbours(numLien) == -1);

              terme_droite(indiceOwner) = terme_droite(indiceOwner)-debitSurfaciqueOutlet*surface_face;
          end
        end  


        %R�solution du syst�me lin�aire
        disp('Solving linear system')
        tic;
        
        if poreNetwork.GetDimension==2
            press=mldivide(matrice,terme_droite);
        else
            L = ichol(matrice); %preconditionnement
          press = minres(matrice,terme_droite,1e-4,100,L,L');
        end
        
        duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
        disp(sprintf('Solving linear system finished. Time spent : %d minutes %f s.',minutes,secondes));
        
        
%         L = ichol(matrice); %pr�conditionnement
%         [press,flag] = pcg(matrice,terme_droite,1e-4,100,L,L');
%         if flag ~= 0
%             disp('Probleme de convergence du systeme lineaire. Flag : ');
%             disp(flag);
%         end


        %Mise en forme de l'output pression
        for i = 1:length(poresPercolants)
            pressions(poresPercolants(i)) = press(i);
        end



        for numLien = liens_internes_envahis
            %Calcul des vitesses internes en fonction des pressions
            numOwner = poreNetwork.LinkOwners(numLien);
            numNeighbour = poreNetwork.LinkNeighbours(numLien);
            debit = conductances(numLien)*(pressions(numOwner)-pressions(numNeighbour));
            debits(numLien) = debit;
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;   
            vitessesMoyennes(numLien) = debit/surface_lien_equivalente;      
        end
        

        %ajout des vitesses inlet
        if strcmp(inletBoundaryCondition,'Dirichlet')
            for numLien = liens_inlet_envahis
                numOwner = poreNetwork.LinkOwners(numLien);
                debit = conductances(numLien)*(pressionInlet-pressions(numOwner));
                surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
                debits(numLien) = debit;
                vitessesMoyennes(numLien) = debit/surface_lien_equivalente;    
            end
        
        elseif  strcmp(inletBoundaryCondition,'Neumann')
          for numLien = liens_inlet_envahis 
            
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
            debit = debitSurfaciqueInlet*surface_lien_equivalente;
            debits(numLien) = debit;
            vitessesMoyennes(numLien) = debit/surface_lien_equivalente;   
          end
        end
        
        %ajout des debits outlet
        if strcmp(outletBoundaryCondition,'Dirichlet')
            for numLien = liens_outlet_envahis
             
                numOwner = poreNetwork.LinkOwners(numLien);
                debit = conductances(numLien)*(pressions(numOwner)-pressionOutlet);
                surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2; 
                debits(numLien) = debit;
                vitessesMoyennes(numLien) = debit/surface_lien_equivalente;    
            end
        
        elseif  strcmp(outletBoundaryCondition,'Neumann')
          for numLien = liens_outlet_envahis 
            
            surface_lien_equivalente = pi/4*(linkDiameters(numLien))^2;
            debit = debitSurfaciqueOutlet*surface_lien_equivalente;
            debits(numLien) = debit;
            vitessesMoyennes(numLien) = debit/surface_lien_equivalente;   
          end
        end
        
 
    end

    
    %calcul de la permeability

    liens_envahis = cluster.GetInvadedLinks;
    liens_inlet_envahis = intersect(liens_envahis,inletLink);
    liens_outlet_envahis = intersect(liens_envahis,inletLink);
    
    totalOutletDebit =  0;
    for iLink = liens_outlet_envahis
        totalOutletDebit = totalOutletDebit+debits(iLink);
    end
    totalInletDebit =  0;
    for iLink = liens_inlet_envahis
        totalInletDebit = totalInletDebit+debits(iLink);
    end
    

    if(abs(totalOutletDebit)>0)
        assert(abs(totalInletDebit-totalOutletDebit)/abs(totalOutletDebit)<5e-2,'Non conservation de la matière !');
    end
    pressionInlet = max(pressions) ; % TO DO
    
    deltaConcentration = pressionInlet-pressionOutlet;
    permeabilityCoefficient = totalOutletDebit/deltaConcentration;   
    
    
end

