%%
%Dans ce tutorial nous allons apprendre a utiliser les physiques de base.
%Les physiques sont codees dans des fonctions dont le nom commence par 
%Compute. Pour voir la liste des physiques disponibles, taper Compute+autocomplétion.

%Commencons par contruire un reseau de pores.
[ network,viewer ]=CreateNetwork('3block2D');

viewer.View('Network')

figure;colorbar;
viewer.View('Boundaries')

%%
%Calculons maintenant une invasion percolation.
%La fonction ComputeInvasionPercolation prend en arguement le reseau, une
%liste de liens d'entree pour l'eau (inletLink), une liste de liens de 
%sortie (outletLink), une option concernant les angles de contact et 
%facultativement une option detaillant les regles d'invasion. 
%Ici nous choisissons de faire rentrer l'eau par la frontiere du bas et la 
%faire sortir par la frontiere du haut. Les numeros de ces frontieres sont
%accessibles soit dans le fichier de geometrie macroscopique, soit avec la
%fonction viewer.View('Boundaries'). Les liens correspondants à ces
%frontieres sont donnes par la methode GetLinksFrontiere du reseau.

help ComputeInvasionPercolation

inletLink=network.GetLinksFrontiere(2);
outletLink=network.GetLinksFrontiere([4,5,6]);

clusterOptions=struct;

cluster = ComputeInvasionPercolation(network,inletLink,outletLink,'hydrophobic',clusterOptions);

class(cluster)

%%
%ComputeInvasionPercolation renvoit un objet de type ClusterMonophasique.
%Ce type d'objet decrit la repartition d'un phase dans le reseau. La classe
%ClusterMonophasique est celle qui gere les deplacements des frontieres
%entre phases. Pour avoir la liste des pores envahis, utiliser la methode
%GetInvadedPores du cluster. 

viewer.View('PoreList',cluster.GetInvadedPores)


%%
%Calculons maintenant la permeabilite du reseau de pores. 
%
%En utilisant l'autocompletion a partir du mot Compute, on trouve la fonction
%ComputeLinearTransport. help ou doc ComputeLinearTransport nous donne la liste 
%des parametres d'entree et de sortie de cette fonction. Nous allons 
%utiliser la  fonction ComputeLinearTransport avec des conductances de 
%permeabilite, dans une region correspondant au reseau tout entier. 

transportPores = 1:network.GetNumberOfPores ;

dynamicViscosity = 1e-3; %viscosity water at ambiant conditions
conductancesPermeability = LocalScaleComputeConductancesStokes(network,dynamicViscosity);


boundaryConditions.inletLink = network.GetLinksFrontiere([4,5,6]);
boundaryConditions.outletLink = network.GetLinksFrontiere([1,2,3]);
boundaryConditions.inletType = 'Dirichlet' ;
boundaryConditions.outletType = 'Dirichlet' ;

nInletLink = length(boundaryConditions.inletLink);
nOutletLink = length(boundaryConditions.outletLink);

boundaryConditions.inletValue = 1*ones(1,nInletLink);
boundaryConditions.outletValue = 0.1*ones(1,nOutletLink);

    
[ pressure, ~, permeabilityCoefficient ] = ComputeLinearTransport( ...
                                            network,transportPores, ...
                                           	conductancesPermeability, ...
                                           	boundaryConditions  ...
                                           	);


disp('Permeability coefficient :')
disp(permeabilityCoefficient)

figure;
viewer.View('PoreField',pressure)


%%
%Calculons maintenant la diffusion relative. 
%
%Nous allons calculer la diffusion dans les pores qui n'ont pas ete envahis
%lors de l'invasion percolation. Il est possible d'indiquer la region dans
%laquelle a lieu la diffusion avec une liste de pores envahis. Ici ce sera
%le complementaire des pores envahis par le cluster liquide.

transportPores = cluster.GetInvadedPoresComplementary ;

diffusivity = 2e-5; % 02 in N2 at ambiant conditions
conductancesDiffusion = LocalScaleComputeConductancesDiffusion(network,diffusivity);


boundaryConditions.inletLink = network.GetLinksFrontiere([4,5,6]);
boundaryConditions.outletLink = network.GetLinksFrontiere([1,2,3]);

boundaryConditions.inletType = 'Neumann' ;
allLinkSurfaces = network.GetLinkData('Surface');   
surfacicFlux = 1;
boundaryConditions.inletValue = surfacicFlux*allLinkSurfaces(boundaryConditions.inletLink);

nOutletLink = length(boundaryConditions.outletLink);
boundaryConditions.outletType = 'Dirichlet' ;
boundaryConditions.outletValue = -1 *ones(1,nOutletLink);


[ concentrations, ~, diffusionCoefficient ] = ComputeLinearTransport( ...
                                                network,transportPores, ...
                                                conductancesDiffusion, ...
                                                boundaryConditions  ...
                                                );
                                            
disp('Diffusion coefficient :')
disp(diffusionCoefficient)

figure
viewer.View('PoreField',concentrations)

