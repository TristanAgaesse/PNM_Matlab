%%
%Dans ce tutorial nous allons apprendre a utiliser calculer une invasion
%percolation avec de multiples points d'injection independants.


%Commencons par contruire un reseau de pores.
[ network,viewer ]=CreateNetwork('GDL_2D');


%Definissons les lois d'invasion et la mouillabilite
clusterOptions.Coalescence = 'none';
clusterOptions.CapillaryPressureLaw = 'LaplaceCylinder';
clusterOptions.SurfaceTension = 72e-3 ; 

contactAngle=110;
theta=contactAngle*pi/180*ones(1,network.GetNumberOfLinks);
network.AddNewLinkData(theta,'ContactAngle');

%Define the injection points : here we are choosing that 100% of inlet links  
%are independant injection points
inletLink = network.GetLinksFrontiere([1,2,3]);
nCluster = length(inletLink);

clustersInletLink = cell(1,nCluster);
clustersOutletLink = cell(1,nCluster);
for iCluster=1:nCluster
    clustersInletLink{iCluster} = inletLink(iCluster);
    clustersOutletLink{iCluster} = network.GetLinksFrontiere([4,5,6]);
end

[clusters,invadedPores] = ComputeInvasionPercolationSeveralClusters( ...
                          	network,nCluster,...
                          	clustersInletLink,clustersOutletLink,...
                          	'currentWettability',clusterOptions );



%Compute number of breakthrough points

nCluster = length(clusters);

figure
clusterPores=zeros(1,network.GetNumberOfPores);
for iCluster=1:nCluster
    clusterPores(clusters{iCluster}.GetInvadedPores)=iCluster;
end
colorScale=[0,nCluster];
viewer.ViewPoreData(clusterPores,colorScale)
