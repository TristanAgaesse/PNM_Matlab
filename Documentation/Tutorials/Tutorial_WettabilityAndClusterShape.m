%% Tutoriel Wettability and ClusterShape
%Dans ce tutoriel nous allons voir comment faire une invasion percolation
%en définissant nous même les angles de contact. Nous allons apprendre à
%gérer les poreData et les linkData qui contiennent les propriétés des
%pores et des liens. Enfin nous allons exporter le réseau et les résultats 
%de simulation vers Paraview pour visualiser les résultats. 

%%
%Commençons par créer un réseau 3D. 
[ network,viewer ]=CreateNetwork('1block3D');

inletLink=network.GetLinksFrontiere(1);outletLink=network.GetLinksFrontiere(2);


%%
%L'objet network contient des data associées aux pores et liens. Ce peut 
%etre les propriétées géométriques des pores et liens, des résultats de 
%simulations... 
%On peut accéder aux poresData et linkData directement dans l'interface
%graphique de Matlab en double cliquant sur l'objet network. Il est aussi
%possible de les manipuler en ligne de commande. 
%
%Ces données peuvent servir dans les algorithmes physiques 
%(tailles...). C'est le cas pour l'invasion percolation qui
%utilise les diametres des liens et les angles de contact dans les liens.
%Pour faire une invasion percolation avec des angles de contact
%personnalises, il faut ajouter a la liste des link data un tableau
%contenant les angles de contact. On utilise la fonction AddNewLinkData
%avec le nom 'ContactAngle'. 

contactAngle=80;
theta=contactAngle*pi/180*ones(1,network.GetNumberOfLinks);
network.AddNewLinkData(theta,'ContactAngle');

network.GetLinkDataList 

network.GetLinkData('ContactAngle');


%%
%Calculons donc une invasion percolation avec l'angle de contact de 80°
%defini a l'etape precedente. Pour visualiser les pores envahis dans
%Paraview ultérieurement, il faut rajouter une liste de pores envahis aux 
%poreData du réseau. On voit aussi que la fonction
%ComputeInvasionPercolation a rajoute la liste des diametres de liens qui
%lui manquait.

[cluster,breakthroughPressure,invasionPressureList]=ComputeInvasionPercolation(network,inletLink,outletLink,'currentWettability');

network.AddNewPoreData(cluster.GetInvadedPoresBooleans,'InvadedPores_80');

network.GetPoreDataList



%%
%Faisons maintenant une invasion percolation avec un angle de contact
%uniforme de 110°. On retire d'abord l'ancien linkData 'ContactAngle', on
%ajoute le nouveau puis . 

contactAngle=110;
network.RemoveLinkData('ContactAngle');
theta=contactAngle*pi/180*ones(1,network.GetNumberOfLinks);
network.AddNewLinkData(theta,'ContactAngle');

[cluster,breakthroughPressure,invasionPressureList]=ComputeInvasionPercolation(network,inletLink,outletLink,'currentWettability');

network.AddNewPoreData(cluster.GetInvadedPoresBooleans,'InvadedPores_110');


network.GetPoreDataList


%%
%Calculons maintenant la diffusion dans le cluster. Pour visualiser les concentrations dans
%Paraview ultérieurement, il faut rajouter les concentrations aux data du
%réseau. 
    
transportPores = cluster.GetInvadedPores ;

conductancesDiffusion = LocalScaleComputeConductancesDiffusion(network);

boundaryConditions.inletLink = inletLink;
boundaryConditions.outletLink = outletLink;
boundaryConditions.inletType = 'Dirichlet';
boundaryConditions.outletType = 'Dirichlet';
boundaryConditions.inletValue = 1*ones(1,length(boundaryConditions.inletLink));
boundaryConditions.outletValue = 0.1*ones(1,length(boundaryConditions.outletLink));

concentrations = ComputeLinearTransport(network,transportPores,conductancesDiffusion,boundaryConditions);
network.AddNewPoreData(concentrations,'DiffusionConcentrations');


%%
%La fonction ComputeDiffusion utilise les conductances de diffusion des
%liens. Si 'ConductancesDiffusion' n'est pas encore defini dans les
%linkData du reseau, la fonction la calcule et l'ajoute au reseau pour une
%future utilisation. 

network.GetLinkDataList


%%
%On peut visualiser le reseau dans Matlab avec le viewer. Neanmoins les 
%fonctionnalites du viewer sont moins etendues en 3D qu'en 2D et le rendu est moins bon.
%Cela est du aux limitations des outils de visualisation de Matlab. 

viewer.View('Network')


%%
%Il est possible d'exporter le reseau dans Paraview pour une meilleure visualisation. Pour cela on ecrit un
%fichier .vtk. Dans ce fichier il y aura des informations sur la geometrie
%du reseau et les data associes aux pores, liens etc... Il y a deux types
%de rendu d'un reseau, Ball and Stick et maillage. Le rendu maillage est
%disponible seulement pour les reseaux des classes PoreNetworkMesh et
%classes derivees.

network.ExportToParaview('NetworkTutorial_Wettability')

