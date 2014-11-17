%Nous allons appprendre � utiliser la classe MacroscopicGeometry pour creer des
%geometries.

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

[V E F] = createSoccerBall;

myGeometry.SetDimension(3);

myGeometry.AddVertices(V);

block=struct; 
block.VertexNumbers=1:60;
block.Filling.PoreNumber=10;
block.Filling.FiberThickness=0.0001;
myGeometry.AddBlock(block)

%Need to add the faces too...
nBoundary=length(F);
for iBoundary=1:nBoundary
    boundary=struct;
    boundary.Face=F{iBoundary};
    boundary.Rugosity='flat';
    boundary.Type='surface';
    
    myGeometry.AddBoundary(boundary);
end




myGeometry.SetNetworkType('PoreNetworkMeshFibrous') ;
myGeometry.SetPavageType('RandomVoronoi');





networkBuilder=NetworkBuilder(myGeometry);

network=networkBuilder.BuildNetwork();

myGeometry.WriteGeometryFile('SoccerBall3D');
