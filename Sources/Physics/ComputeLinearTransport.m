function [ fieldValue, flux,  effectiveConductance ]  =  ComputeLinearTransport(network,transportPores,conductances,boundaryConditions)
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
%       	boundaryConditions.inletValue = one value for each inlet links (if Neuman, flux oriente de owner a neighbour)
%       	boundaryConditions.outletValue = one value for each inlet links (if Neuman, flux oriente de owner a neighbour)
%           boundaryConditions.solver (optional) = 'mldivide','minres','pcg' or 'gmres'
%
%Output : [ fieldValue, flux, effectiveConductance ]
    
    
    disp('Computing linear transport ')
    tic;
    
    
    %Checking inputs
    [inletLink,outletLink,inletValue,outletValue,inletType,outletType,fieldValue,flux,conductances,transportPores,boundaryLinkInnerPore,solver,effectiveConductanceFormula]=...
                    CheckInputs(network,transportPores,conductances,boundaryConditions);
    
    
    %Decompose transportPores into percolation paths = connexes components
    cluster = network.CreateVoidCluster;
    cluster.InvadedPores = transportPores;
    cluster.SetClusterOptions(struct);
    composantesConnexesPercolation = cluster.FindPercolationPath(inletLink,outletLink);
    
    
    %Resolution percolation path by percolation path
    
    for iComposanteConnexe = 1:length(composantesConnexesPercolation)
        
        %Extraction des liens envahis (liens envahis internes, inlet et outlet)
        
        percolatingCluster = composantesConnexesPercolation{iComposanteConnexe};
        poresPercolants = percolatingCluster.GetInvadedPores;
        nPorePercolant = length(poresPercolants);
        
        liens_envahis = percolatingCluster.GetInvadedLinks;
        interfaceLinks = percolatingCluster.GetInterfaceLinks;
        liens_internes_envahis = setdiff(liens_envahis,interfaceLinks);
        liens_inlet_envahis = intersect(interfaceLinks,inletLink);
        liens_outlet_envahis = intersect(interfaceLinks,outletLink);
        
        nPore = network.GetNumberOfPores;
        poresPercolantsIndices=zeros(1,nPore);
        poresPercolantsIndices(poresPercolants)=1:nPorePercolant;
        
        
        %Remplissage matrice
        
        stiffnessMatrix = FillMatrix(network,conductances,liens_internes_envahis,...
                            liens_inlet_envahis,liens_outlet_envahis,...
                            inletType,outletType,poresPercolantsIndices,nPorePercolant,boundaryLinkInnerPore);
        
        
        %Remplissage du terme de droite 
        
        rigthHandSide = FillRigthHandSide(network,conductances,...
                        liens_inlet_envahis,liens_outlet_envahis,inletType,outletType,...
                        inletValue,outletValue,poresPercolantsIndices,nPorePercolant,boundaryLinkInnerPore);
        
        
        %Resolution du systeme lineaire
        
        if strcmp(solver,'mldivide')
            fieldVal=mldivide(stiffnessMatrix,rigthHandSide);
            
        elseif strcmp(solver,'minres')
            L = ichol(stiffnessMatrix); %preconditionnement
        	fieldVal = minres(stiffnessMatrix,rigthHandSide,1e-4,100,L,L');
            
        elseif strcmp(solver,'gmres')
        	fieldVal = gmres(stiffnessMatrix,rigthHandSide);       
            
        elseif strcmp(solver,'pcg')
            L = ichol(stiffnessMatrix); %preconditionnement
        	fieldVal = pcg(stiffnessMatrix,rigthHandSide,1e-4,100,L,L');            
        end
        clear stiffnessMatrix;
        
        
        %Reassemblage de l'output concentration
        [fieldValue,flux]=ReassembleComputedField(fieldValue,flux,network,...
                            poresPercolants,conductances,fieldVal,...
                            liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,...
                            inletType,outletType,inletValue,outletValue,boundaryLinkInnerPore);
        
    end
    
    %Check computation
    CheckComputation(fieldValue,flux,inletLink,outletLink)
          
    
    %Compute effective transport property
    effectiveConductance = ComputeEffectiveConductance(network,inletLink,outletLink,inletValue,outletValue,inletType,outletType,fieldValue,flux,conductances,transportPores,boundaryLinkInnerPore,effectiveConductanceFormula);
    
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Computing linear transport finished. Time spent : %d minutes %f s. \n',minutes,secondes);
    
end





%---------------------------------------------------------------------------------------------
function [inletLink,outletLink,inletValue,outletValue,inletType,outletType,fieldValue,flux,conductances,transportPores,boundaryLinkInnerPore,solver,effectiveConductanceFormula]=...
                        CheckInputs(network,transportPores,conductances,boundaryConditions)
                    
    %Check inputs and initialize algorithm
    conductances=transpose(conductances);
    
    nPore=network.GetNumberOfPores;
    nLink=network.GetNumberOfLinks;
    
    assert( isa(network,'PoreNetwork'),'LinearTransport : first input must be a PoreNetwork object')
    assert( length(transportPores)<=nPore,...
            'LinearTransport : second input transportPores must be of length <=network.GetNumberOfPore')
    assert( ~isempty(transportPores),'LinearTransport : second input transportPores must not be empty')
    assert( size(conductances,2)==nLink,...
            'LinearTransport : third input conductance must be of length network.GetNumberOfLinks')

    if size(transportPores,1)~=1
        assert(size(transportPores,2)==1);
        transportPores=transpose(transportPores);
    end    
    transportPores=unique(transportPores);
    assert(min(transportPores)>0 && max(transportPores)<=nPore,'Wrong pore number in transport pores')
        
    assert( isa(boundaryConditions,'struct') );
    assert( max(boundaryConditions.inletLink)<=network.GetNumberOfLinks )
    assert( max(boundaryConditions.outletLink)<=network.GetNumberOfLinks )
    assert( strcmp(boundaryConditions.inletType,'Neumann') || strcmp(boundaryConditions.inletType,'Dirichlet'))
    assert( strcmp(boundaryConditions.outletType,'Neumann') || strcmp(boundaryConditions.outletType,'Dirichlet'))
	assert( length(boundaryConditions.inletValue)==length(boundaryConditions.inletLink))
	assert( length(boundaryConditions.outletValue)==length(boundaryConditions.outletLink))     
         
    
    if isfield(boundaryConditions,'solver')
        solver=boundaryConditions.solver;
        assert(strcmp(solver,'mldivide') ||strcmp(solver,'minres') ||strcmp(solver,'pcg') ||strcmp(solver,'gmres'),'Solver type not supported'); 
    else 
        if network.GetDimension==2
            solver = 'mldivide';
        else
            solver = 'minres';
        end
    end
    
    
    [boundaryLinkInnerPore,boundaryConditions]=FindBoundaryLinkInnerPore(transportPores,boundaryConditions,network); 
    
    
    nLink = network.GetNumberOfLinks;
    nPore = network.GetNumberOfPores;
    fieldValue = zeros(1,nPore);
    flux = zeros(1,nLink);               %flux oriente de owner a neighbour   

    inletValue=zeros(1,network.GetNumberOfLinks);
    inletValue(boundaryConditions.inletLink)=boundaryConditions.inletValue;
    outletValue=zeros(1,network.GetNumberOfLinks);
    outletValue(boundaryConditions.outletLink)=boundaryConditions.outletValue;

    inletType = boundaryConditions.inletType;
    outletType = boundaryConditions.outletType;
    
    inletLink = boundaryConditions.inletLink;
    outletLink = boundaryConditions.outletLink;
      
    
    defaultFormulaType='fluxOverPotential';
    defaultFormulaDirection='None';
    
    if isfield(boundaryConditions,'effectiveConductanceFormula')
        formula=boundaryConditions.effectiveConductanceFormula;
        assert(isfield(formula,'Type'))
        type=formula.Type;
        assert(strcmp(type,'fluxOverPotential') || strcmp(type,'integralOverVolume'),'effectiveConductanceFormula type not supported');
        if strcmp(type,'integralOverVolume')
            assert(isfield(formula,'Direction'))
            direction=formula.Direction;
            assert(strcmp(direction,'x') || strcmp(direction,'y') || strcmp(direction,'z'),'effectiveConductanceFormula direction not supported');
        else
            effectiveConductanceFormula.Direction=defaultFormulaDirection;
        end
    else 
        effectiveConductanceFormula.Type=defaultFormulaType;
        effectiveConductanceFormula.Direction=defaultFormulaDirection;
    end

end


%---------------------------------------------------------------------------------------------
function [boundaryLinkInnerPore,boundaryConditions]=FindBoundaryLinkInnerPore(transportPores,boundaryConditions,network)
    % output : boundaryLinkInnerPore : donne pour chaque lien inlet et 
    %       outlet le numero du pore qui donne vers l'interieur du domaine.    
    
    inletLink = boundaryConditions.inletLink;
    outletLink = boundaryConditions.outletLink;
    inletValue = boundaryConditions.inletValue;
    outletValue = boundaryConditions.outletValue;
    
    if size(inletLink,1)~=1
        inletLink=transpose(inletLink);
    end
    if size(outletLink,1)~=1
        outletLink=transpose(outletLink);
    end
    if size(inletValue,1)~=1
        inletValue=transpose(inletValue);
    end
    if size(outletValue,1)~=1
        outletValue=transpose(outletValue);
    end
    
    assert(isempty(intersect(inletLink,outletLink)),'inletLink and outletLink have some link in common !');
    
    transportRegionBoundary = network.GetPoreRegionBoundaryLinks(transportPores);
    
    inletLinkToKeep=ismember(inletLink,transportRegionBoundary);
    outletLinkToKeep=ismember(outletLink,transportRegionBoundary);
    
    inletLink=inletLink(inletLinkToKeep);
    inletValue=inletValue(inletLinkToKeep);
    outletLink=outletLink(outletLinkToKeep);
    outletValue=outletValue(outletLinkToKeep);
    
    assert(~isempty(inletLink) && ~isempty(outletLink),'No inlet link or outlet link is on the transport region boundary !')
    
    boundaryLinkInnerPore=zeros(1,network.GetNumberOfLinks);
    
    boundaryLink=[inletLink,outletLink];
    linksNeighboors=network.GetPoresOfLink(boundaryLink);
    
    % Check that all boundary link has only one neighbor pore inside the 
    % transport domain and note the number of this pore
    
    
    % 1 : boundary link qui sont internes au reseau 
    
    isInternal = ismember(boundaryLink,network.GetLinksFrontiere(0));
    internalLinksNeighboors= linksNeighboors(isInternal,:);
    
    booleanTransportPores=zeros(network.GetNumberOfPores,1);
    booleanTransportPores(transportPores)=1;
    
    linksNeighboorsElements=horzcat(booleanTransportPores(internalLinksNeighboors(:,1)),booleanTransportPores(internalLinksNeighboors(:,2)));
    [linksNeighboorsElements,order]=sort(linksNeighboorsElements,2);
    
    assert(all(linksNeighboorsElements(:,2)),'un link inlet ou outlet non situe sur la frontiere de transport pores !')
    assert(all(not(linksNeighboorsElements(:,1))),'un link inlet ou outlet non situe sur la frontiere de transport pores !')
    
    boundaryLinkInnerPore( boundaryLink(isInternal) )=order(:,end);

    % 2 : boundary link qui sont sur la frontiere du reseau
    
    boudaryLink_surface = boundaryLink(not(isInternal));
    innerPore = network.LinkOwners(boudaryLink_surface);
    assert(all(booleanTransportPores(innerPore)),'un link inlet ou outlet non situe sur la frontiere de transport pores !');
    boundaryLinkInnerPore(boudaryLink_surface) = innerPore;
    
    assert( all( booleanTransportPores(boundaryLinkInnerPore(boundaryLinkInnerPore>0))))
    
    boundaryConditions.inletLink = inletLink ;
    boundaryConditions.outletLink = outletLink ;
    
    boundaryConditions.inletValue = inletValue;
    boundaryConditions.outletValue = outletValue;
end


%---------------------------------------------------------------------------------------------
function CheckComputation(fieldValue,flux,inletLink,outletLink)
    %Check computation with a mass balance
    
    totalInletDebit = sum(flux(inletLink));
    totalOutletDebit = sum(flux(outletLink));
    
    if(abs(totalOutletDebit)>0)
        fluxRelativeError = abs(totalInletDebit-totalOutletDebit)/abs(totalOutletDebit);
        assert(fluxRelativeError<1e-1,'Non conservation de la matière !');
    end
    
end 


%---------------------------------------------------------------------------------------------
function effectiveConductance = ComputeEffectiveConductance(network,inletLink,outletLink,inletValue,outletValue,inletType,outletType,...
                                                            fieldValue,flux,conductances,transportPores,boundaryLinkInnerPore,effectiveConductanceFormula)
    %Compute effective transport property
    
    formulaType = effectiveConductanceFormula.Type;
    direction = effectiveConductanceFormula.Direction;
    
    
    switch formulaType
    
        case 'fluxOverPotential'
            % a very simple formula is used here
            switch inletType
                case 'Dirichlet'
                    inletMeanField =  mean(inletValue(inletLink));
                case 'Neumann'
                    inletFlux = inletValue;
                    inletField = fieldValue(boundaryLinkInnerPore(inletLink))-conductances(inletLink).*inletFlux(inletLink);
                    inletMeanField =  mean(inletField);
            end
            
            switch outletType
                case 'Dirichlet'
                    outletMeanField =  mean(outletValue(outletLink));
                case 'Neumann'
                    outletFlux = outletValue;
                    outletField = fieldValue(boundaryLinkInnerPore(outletLink))-conductances(outletLink).*outletFlux(outletLink);
                    outletMeanField =  mean(outletField);
            end
            
            deltaConcentration = inletMeanField-outletMeanField;
            
            totalOutletDebit = sum(flux(outletLink));
            
            effectiveConductance = abs( totalOutletDebit/deltaConcentration );   
            
        case 'integralOverVolume'
            % integral over volume of beta*gradxT
            
            % sum over internal link
            
            % V_half_link=V_pore/nLinkOfPore
            
            
            % order pore & link centers in x direction
            
            % Tout = fieldValue(elementXmax)
            % Tin = fieldValue(elementXmin)
            % Lx = 
            % betaPore = 
            % betaPore*((Tout-Tin)/Lx)*V_half_link
            
            
            
            
            % add inlet and outlet links
            
            
    end
end




%---------------------------------------------------------------------------------------------
function matrice = FillMatrix(network,conductances,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,...
                              inletType,outletType,poresPercolantsIndices,nPorePercolant,boundaryLinkInnerPore)
    %Fill sparse rigidity matrix
    
    value_diag=zeros(1,nPorePercolant);

    indiceI_nonDiag=zeros(2*length(liens_internes_envahis),1);
    indiceJ_nonDiag=zeros(2*length(liens_internes_envahis),1);
    value_nonDiag=zeros(2*length(liens_internes_envahis),1);

    
    %Contribution des liens internes
    
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
    if strcmp(inletType,'Dirichlet')
        for numLien = liens_inlet_envahis
            %numOwner = network.LinkOwners(numLien);
            numOwner = boundaryLinkInnerPore(numLien);
            indiceOwner = poresPercolantsIndices(numOwner);
            %complements aux termes diagonaux 
            value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
        end
    end

    if strcmp(outletType,'Dirichlet')
    	for numLien = liens_outlet_envahis
        	%numOwner = network.LinkOwners(numLien);
            numOwner = boundaryLinkInnerPore(numLien);
        	indiceOwner = poresPercolantsIndices(numOwner);
        	%complements aux termes diagonaux 
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
function terme_droite = FillRigthHandSide(network,conductances,liens_inlet_envahis,liens_outlet_envahis,...
                            inletType,outletType,inletValue,outletValue,poresPercolantsIndices,nPorePercolant,boundaryLinkInnerPore)
                        
    %Fill rigth hand side

    terme_droite = zeros(1,nPorePercolant);
        
    %Contribution of inlet links
            
    if strcmp(inletType,'Dirichlet')
        
        indiceOwner = poresPercolantsIndices(boundaryLinkInnerPore(liens_inlet_envahis));
    	%indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_inlet_envahis));
        for i=1:length(liens_inlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))...
                            +inletValue(liens_inlet_envahis(i))*conductances(liens_inlet_envahis(i));
        end
        
    elseif  strcmp(inletType,'Neumann')
        
        %indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_inlet_envahis));
        indiceOwner = poresPercolantsIndices(boundaryLinkInnerPore(liens_inlet_envahis));
        for i =1:length(liens_inlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))+inletValue(liens_inlet_envahis(i));
        end
    end
    
    %Contribution of outlet links
    
    if strcmp(outletType,'Dirichlet')
        
        %indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_outlet_envahis));
        indiceOwner = poresPercolantsIndices(boundaryLinkInnerPore(liens_outlet_envahis));
        for i =1:length(liens_outlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))...
                    +outletValue(liens_outlet_envahis(i))*conductances(liens_outlet_envahis(i));
        end
        
    elseif  strcmp(outletType,'Neumann')
        
        %indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_outlet_envahis));
        indiceOwner = poresPercolantsIndices(boundaryLinkInnerPore(liens_outlet_envahis));
        for i =1:length(liens_outlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))-outletValue(liens_outlet_envahis(i));       
        end
    end
    
    terme_droite=transpose(terme_droite);
end



%---------------------------------------------------------------------------------------------
function [fieldValue,flux]=ReassembleComputedField(fieldValue,flux,network,poresPercolants,conductances,fieldVal,...
                        liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,...
                        inletType,outletType,inletValue,outletValue,boundaryLinkInnerPore)
                    
    %Reassemble the solution of the differents connexes components

    fieldValue(poresPercolants) = fieldVal;
    
    %Calcul des flux internes
    numOwner = network.LinkOwners(liens_internes_envahis);
    numNeighbour = network.LinkNeighbours(liens_internes_envahis);
    flux(liens_internes_envahis) = conductances(liens_internes_envahis).*(fieldValue(numOwner)-fieldValue(numNeighbour));

    %calcul des flux inlet
        
    if strcmp(inletType,'Dirichlet')
        %numOwner = network.LinkOwners(liens_inlet_envahis);
        numOwner = boundaryLinkInnerPore(liens_inlet_envahis);
        flux(liens_inlet_envahis) = conductances(liens_inlet_envahis).*(inletValue(liens_inlet_envahis)-fieldValue(numOwner));

    elseif  strcmp(inletType,'Neumann')
    	flux(liens_inlet_envahis) = inletValue(liens_inlet_envahis);
    end

    %calcul des flux outlet
    
    if strcmp(outletType,'Dirichlet')
        %numOwner = network.LinkOwners(liens_outlet_envahis);
        numOwner = boundaryLinkInnerPore(liens_outlet_envahis);
        flux(liens_outlet_envahis) = conductances(liens_outlet_envahis).*(fieldValue(numOwner)-outletValue(liens_outlet_envahis));

    elseif  strcmp(outletType,'Neumann')
        flux(liens_outlet_envahis) = outletValue(liens_outlet_envahis);
    end

end



