%Nous allons appprendre a utiliser la classe MacroscopicGeometry pour creer des
%reseaux de pores de forme variees.

%%
%Une geometrie peut etre definie dans un fichier ou par un script qui cree un 
%objet MacroscopicGeometry.

myGeometry=MacroscopicGeometry();

%%
%Nous allons utiliser la librairie math geom (utilisee par PNM_Matlab) pour creer
%geometries exotiques.
%L'element de base d'une geometrie est le block. Celui-ci peut-etre n'importe 
%quel polynome convexe. C'est inspire de l'outil blockMesh de openFoam.
%Pour le moment il n'y a que les reseaux Voronoi qui peuvent avoir des formes
%exotiques, les reseaux structures doivent etre dans un pave aligne avec les axes

nPore = 1000;


[V, E, F] = createOctahedron;

myGeometry.AddVertices(V);

for iBoundary=1:size(F,1)
    boundary=struct;
    boundary.Face=F(iBoundary,:);
    boundary.Rugosity='flat';
    boundary.Type='surface';
    boundary.Name='';
    myGeometry.AddBoundary(boundary);
end

block=struct; 
block.VertexNumbers=1:size(V,1);
block.Filling.PoreNumber=nPore;
block.Filling.FiberThickness=0.0001;
myGeometry.AddBlock(block)


myGeometry.SetDimension(3);
myGeometry.SetConvertToMeters(1e-3);  %scale unit =  1e-3 m
myGeometry.SetNetworkType('PoreNetworkMeshFibrous') ;
myGeometry.SetPavageType('RandomVoronoi');

myGeometry.WriteGeometryFile('Octahedron3D.txt');

networkBuilder=NetworkBuilder(myGeometry);

network=networkBuilder.BuildNetwork();


viewer=Viewer(network.PrivateInternalOutputStruct);
viewer.View('Network')
