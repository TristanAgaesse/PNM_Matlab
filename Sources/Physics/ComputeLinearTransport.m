function [ fieldValue, flux,  effectiveTransportProperty ]  =  ComputeLinearTransport(network,transportPores,conductances,boundaryConditions)
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
%
%Output : [ fieldValue, flux, effectiveTransportProperty ]


    disp('Computing linear transport ')
    tic;

    
    %Checking inputs
    [inletLink,outletLink,inletValue,outletValue,inletType,outletType,fieldValue,flux]=...
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
                            inletType,outletType,poresPercolantsIndices,nPorePercolant);
        
        
        %Remplissage du terme de droite 
        
        rigthHandSide = FillRigthHandSide(network,conductances,...
                        liens_inlet_envahis,liens_outlet_envahis,inletType,outletType,...
                        inletValue,outletValue,poresPercolantsIndices,nPorePercolant);
        
        
        %Resolution du systeme lineaire

        if network.GetDimension==2
            fieldVal=mldivide(stiffnessMatrix,rigthHandSide);
        else
            L = ichol(stiffnessMatrix); %preconditionnement
        	fieldVal = minres(stiffnessMatrix,rigthHandSide,1e-4,100,L,L');
        end
        clear stiffnessMatrix;
        
        
        %Reassemblage de l'output concentration
        [fieldValue,flux]=ReassembleComputedField(fieldValue,flux,network,...
                            poresPercolants,conductances,fieldVal,...
                            liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,...
                            inletType,outletType,inletValue,outletValue);
        
    end
    
    %Check computation
    CheckComputation(fieldValue,flux,inletLink,outletLink)
          
    
    %Compute effective transport property
    effectiveTransportProperty = ComputeEffectiveTransportProperty(network,fieldValue,flux,transportPores,inletLink,outletLink);
    
    
    duree = toc;minutes = floor(duree/60);secondes = duree-60*minutes;
    fprintf('Computing linear transport finished. Time spent : %d minutes %f s. \n',minutes,secondes);
    
end





%---------------------------------------------------------------------------------------------
function [inletLink,outletLink,inletValue,outletValue,inletType,outletType,fieldValue,flux]=...
                        CheckInputs(network,transportPores,conductances,boundaryConditions)
                    
    %Check inputs and initialize algorithm

    assert( isa(network,'PoreNetwork'),'LinearTransport : first input must be a PoreNetwork object')
    assert( length(transportPores)<=network.GetNumberOfPores,...
            'LinearTransport : second input transportPores must be of length <=network.GetNumberOfPore')
    assert( ~isempty(transportPores),'LinearTransport : second input transportPores must not be empty')
    assert( length(conductances)==network.GetNumberOfLinks,...
            'LinearTransport : third input conductance must be of length network.GetNumberOfLinks')

    assert( isa(boundaryConditions,'struct') );
    assert( max(boundaryConditions.inletLink)<=network.GetNumberOfLinks )
    assert( max(boundaryConditions.outletLink)<=network.GetNumberOfLinks )
    assert( strcmp(boundaryConditions.inletType,'Neumann') || strcmp(boundaryConditions.inletType,'Dirichlet'))
    assert( strcmp(boundaryConditions.outletType,'Neumann') || strcmp(boundaryConditions.outletType,'Dirichlet'))
	assert( length(boundaryConditions.inletValue)==length(boundaryConditions.inletLink))
	assert( length(boundaryConditions.outletValue)==length(boundaryConditions.outletLink))     
         
    inletLink = boundaryConditions.inletLink;
    outletLink = boundaryConditions.outletLink;

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
end


%---------------------------------------------------------------------------------------------
function CheckComputation(fieldValue,flux,inletLink,outletLink)
    %Check computation with a mass balance
    
    totalInletDebit = sum(flux(inletLink));
    totalOutletDebit = sum(flux(outletLink));
    
    if(abs(totalOutletDebit)>0)
        assert(abs(totalInletDebit-totalOutletDebit)/abs(totalOutletDebit)<1e-1,'Non conservation de la matière !');
    end

end 


%---------------------------------------------------------------------------------------------
function effectiveTransportProperty = ComputeEffectiveTransportProperty(network,fieldValue,flux,transportPores,inletLink,outletLink)
    %Compute effective transport property
    
    %very simplistic formula here, should be checked
    totalOutletDebit = sum(flux(outletLink));
    
    
    inletPores=setdiff(network.GetPoresOfLink(inletLink),-1);
    inletPores=inletPores(ismember(inletPores,transportPores));
    
    outletPores=setdiff(network.GetPoresOfLink(outletLink),-1);
    outletPores=outletPores(ismember(outletPores,transportPores));
    
    deltaConcentration = mean(fieldValue(inletPores))-mean(fieldValue(outletPores));
    
    
    effectiveTransportProperty = abs( totalOutletDebit/deltaConcentration );   
end




%---------------------------------------------------------------------------------------------
function matrice = FillMatrix(network,conductances,liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,...
                              inletType,outletType,poresPercolantsIndices,nPorePercolant)
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
            numOwner = network.LinkOwners(numLien);
            indiceOwner = poresPercolantsIndices(numOwner);
            %complements aux termes diagonaux 
            value_diag(indiceOwner)=value_diag(indiceOwner)+conductances(numLien);
        end
    end

    if strcmp(outletType,'Dirichlet')
    	for numLien = liens_outlet_envahis
        	numOwner = network.LinkOwners(numLien);
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
                            inletType,outletType,inletValue,outletValue,poresPercolantsIndices,nPorePercolant)
                        
    %Fill rigth hand side

    terme_droite = zeros(1,nPorePercolant);
        
    %Contribution of inlet links
            
    if strcmp(inletType,'Dirichlet')
        
    	indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_inlet_envahis));
        for i=1:length(liens_inlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))...
                            +inletValue(liens_inlet_envahis(i))*conductances(liens_inlet_envahis(i));
        end
        
    elseif  strcmp(inletType,'Neumann')
        
        indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_inlet_envahis));
        for i =1:length(liens_inlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))+inletValue(liens_inlet_envahis(i));
        end
    end
    
    %Contribution of outlet links
    
    if strcmp(outletType,'Dirichlet')
        
        indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_outlet_envahis));
        for i =1:length(liens_outlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))...
                    +outletValue(liens_outlet_envahis(i))*conductances(liens_outlet_envahis(i));
        end
        
    elseif  strcmp(outletType,'Neumann')
        
        indiceOwner = poresPercolantsIndices(network.LinkOwners(liens_outlet_envahis));
        for i =1:length(liens_outlet_envahis)
            terme_droite(indiceOwner(i)) = terme_droite(indiceOwner(i))-outletValue(liens_outlet_envahis(i));       
        end
    end
    
    terme_droite=transpose(terme_droite);
end



%---------------------------------------------------------------------------------------------
function [fieldValue,flux]=ReassembleComputedField(fieldValue,flux,network,poresPercolants,conductances,fieldVal,...
                        liens_internes_envahis,liens_inlet_envahis,liens_outlet_envahis,...
                        inletType,outletType,inletValue,outletValue)
                    
    %Reassemble the solution of the differents connexes components

    fieldValue(poresPercolants) = fieldVal;
    
    %Calcul des flux internes
    numOwner = network.LinkOwners(liens_internes_envahis);
    numNeighbour = network.LinkNeighbours(liens_internes_envahis);
    flux(liens_internes_envahis) = transpose(conductances(liens_internes_envahis)).*(fieldValue(numOwner)-fieldValue(numNeighbour));

    %calcul des flux inlet
        
    if strcmp(inletType,'Dirichlet')
        numOwner = network.LinkOwners(liens_inlet_envahis);
        flux(liens_inlet_envahis) = transpose(conductances(liens_inlet_envahis)).*(inletValue(liens_inlet_envahis)-fieldValue(numOwner));

    elseif  strcmp(inletType,'Neumann')
    	flux(liens_inlet_envahis) = inletValue(liens_inlet_envahis);
    end

    %calcul des flux outlet
    
    if strcmp(outletType,'Dirichlet')
        numOwner = network.LinkOwners(liens_outlet_envahis);
        flux(liens_outlet_envahis) = transpose(conductances(liens_outlet_envahis)).*(fieldValue(numOwner)-outletValue(liens_outlet_envahis));

    elseif  strcmp(outletType,'Neumann')
        flux(liens_outlet_envahis) = outletValue(liens_outlet_envahis);
    end

end



