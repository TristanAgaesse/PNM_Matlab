%Nous allons appprendre a utiliser la classe MacroscopicGeometry pour creer des
%reseaux de pores de forme variees.


%La fonction basique pour creer un reseau est CreateNetwork qui lit un 
%fichier de geometrie. En fait cette fonction cree en memoire un objet 
%MacroscopicGeometry qui contient les informations du fichier puis
%construit le reseau avec un objet NetworkBuilder qui interprete l'objet 
%MacroscopicGeometry. 
%Exemple avec un reseau regulier : creons un reseau a partir du fichier
%1block3D_Structured et visualisons le avec paraview.

folderName=fileparts(which('Tutorial_GeometryDefinition.m'));

myGeometry=MacroscopicGeometry();
filename=strcat(folderName,'/1block3D_Structured');
myGeometry.LoadGeometry(filename);

disp(myGeometry)

networkBuilder=NetworkBuilder(myGeometry);

network=networkBuilder.BuildNetwork();
filename=strcat(folderName,'/tutorialMacroscopicGeometry_StructuredNetwork');
network.ExportToParaview(filename)

clear myGeometry

%%
%Une geometrie peut etre definie dans un fichier ou par un script qui cree un 
%objet MacroscopicGeometry.
%Dans cet exemple nous allons utiliser un script pour creer une geometrie 
%complexe, l'octahedron.

[V, E, F] = createOctahedron; %utilise la librairie math_geom

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
    boundary.Rugosity='rough';
    boundary.Type='surface';
    boundary.Name='';
    myGeometry.AddBoundary(boundary);
end

%Une fois la geometrie definie, le reseau est construit avec un NetworkBuilder
networkBuilder=NetworkBuilder(myGeometry);

network=networkBuilder.BuildNetwork();

%Visualiser le reseau dans Matlab avec un viewer et exportons vers paraview.
filename=strcat(folderName,'/tutorialMacroscopicGeometry_octahedron');
network.ExportToParaview(filename)

viewer=Viewer(network.PrivateInternalOutputStruct);
figure
viewer.View('Network')


%Il est possible de creer un fichier contenant les parametres de la
%geometrie macroscopique. 
filename=strcat(folderName,'/Octahedron3D.txt');
myGeometry.WriteGeometryFile(filename);



%Pour le moment il n'y a que les reseaux Voronoi qui peuvent avoir des formes
%exotiques, les reseaux reguliers doivent etre dans un pave aligne avec les axes



