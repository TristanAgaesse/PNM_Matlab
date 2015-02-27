%Nous allons appprendre a utiliser la classe MacroscopicGeometry pour creer des
%reseaux de pores de forme variees.

%%
%Une geometrie peut etre definie dans un fichier ou par un script qui cree un 
%objet MacroscopicGeometry.
%Dans cet exemple nous allons utiliser la librairie math geom pour creer
%une geometrie exotique, l'octahedron .

[V, E, F] = createOctahedron;

%Creons maintenant un objet MacroscopicGeometry correspondant a l'octahedron.
%L'element de base d'une geometrie est le block. Celui-ci peut-etre n'importe 
%quel polynome convexe. Cette syntaxe est inspiree de l'outil blockMesh de openFoam.

myGeometry=MacroscopicGeometry();

myGeometry.SetDimension(3);
myGeometry.SetConvertToMeters(1e-3);  %scale unit =  1e-3 m
myGeometry.SetNetworkType('PoreNetworkMeshFibrous') ;
myGeometry.SetPavageType('RandomVoronoi');

myGeometry.AddVertices(V);

nPore = 1000;
block=struct; 
block.VertexNumbers=1:size(V,1);
block.Filling.PoreNumber=nPore;
block.Filling.FiberThickness=0.0001;
myGeometry.AddBlock(block)

for iBoundary=1:size(F,1)
    boundary=struct;
    boundary.Face=F(iBoundary,:);
    boundary.Rugosity='flat';
    boundary.Type='surface';
    boundary.Name='';
    myGeometry.AddBoundary(boundary);
end

%Une fois la geometrie definie, le reseau est construit avec un NetworkBuilder
networkBuilder=NetworkBuilder(myGeometry);

network=networkBuilder.BuildNetwork();

%Visualiser le reseau dans Matlab avec un viewer et exportons vers paraview.
network.ExportToParaview('tutorialMacroxcopicGeometry.vtk')

viewer=Viewer(network.PrivateInternalOutputStruct);
viewer.View('Network')


%Il est possible de creer un fichier contenant les parametres de la
%geometrie macroscopique. Il suffira d'utiliser la fonction CreateNetwork
%avec ce fichier pour creer le reseau, ou d'utiliser
%MacroscopicGeometry.LoadGeometry pour recupere l'objet MacroscopicGeometry.

myGeometry.WriteGeometryFile('Octahedron3D.txt');

myNewGeometry=MacroscopicGeometry();
myNewGeometry.LoadGeometry('Octahedron3D.txt');


%Pour le moment il n'y a que les reseaux Voronoi qui peuvent avoir des formes
%exotiques, les reseaux structures doivent etre dans un pave aligne avec les axes
