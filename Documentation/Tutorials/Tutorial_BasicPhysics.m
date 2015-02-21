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
%sortie (outletLink) et une option concernant les angles de contact. 
%Ici nous choisissons de faire rentrer l'eau par la frontiere du bas et la 
%faire sortir par la frontiere du haut. Les numeros de ces frontieres sont
%accessibles soit dans le fichier de geometrie macroscopique, soit avec la
%fonction viewer.View('Boundaries'). Les liens correspondants à ces
%frontieres sont donnes par la methode GetLinksFrontiere du reseau.

help ComputeInvasionPercolation

inletLink=network.GetLinksFrontiere(2);
outletLink=network.GetLinksFrontiere([4,5,6]);
[cluster,breakthroughPressure,invasionPressureList]=ComputeInvasionPercolation(network,inletLink,outletLink,'hydrophobic');

class(cluster)

%%
%ComputeInvasionPercolation renvoit un objet de type ClusterMonophasique.
%Ce type d'objet decrit la repartition d'un phase dans le reseau. La classe
%ClusterMonophasique est celle qui gere les deplacements des frontieres
%entre phases. Pour avoir la liste des pores envahis, utiliser la methode
%GetInvadedPores du cluster. 

viewer.View('PoreList',cluster.GetInvadedPores)


%%
%Nous allons maintenant faire un calcul de diffusion. En utilisant 
%l'autocompletion a partir du mot Compute, on trouve la fonction
%ComputeDiffusion. help ou doc ComputeDiffusion nous donne la liste des parametres
%d'entree et de sortie de cette fonction.
%
%Nous allons calculer la diffusion dans les pores qui n'ont pas ete envahis
%lors de l'invasion percolation. Il est possible d'indiquer la region dans
%laquelle a lieu la diffusion avec un objet cluster, qui contient une
%liste de pores envahis. Ici le cluster a donner a ComputeDiffusion est le
%complementaire du cluster obtenu par invasion percolation.


complementaryCluster=cluster.GetComplementaryCluster;

boundaryConditions=struct;
boundaryConditions.inletLink = network.GetLinksFrontiere([4,5,6]);
boundaryConditions.outletLink = network.GetLinksFrontiere([1,2,3]);
boundaryConditions.inletType = 'Dirichlet' ;
boundaryConditions.outletType = 'Dirichlet' ;
boundaryConditions.inletValue = 1;
boundaryConditions.outletValue = 0.1;


[ concentrations, ~, ~, diffusionCoefficient ]=ComputeDiffusion(network,complementaryCluster, boundaryConditions);
disp(diffusionCoefficient)

figure
viewer.View('PoreField',concentrations)

%%
%Calculons maintenant la permeabilite du reseau de pores. Nous allons 
%utiliser la  fonction ComputePermeabilite dans une region correspondant au
%reseau tout entier. Pour cela il nous faut un cluster qui recouvre tout le
%reseau. Le moyen le plus simple de l'obtenir utilise la fonction 
%CreateFullCluster du reseau.

fullCluster=network.CreateFullCluster;

[ concentrations, ~, ~, diffusionCoefficient ]=ComputeDiffusion(network,fullCluster, boundaryConditions);
disp(diffusionCoefficient)

figure;
viewer.View('PoreField',concentrations)