%Script construisant l'input de CreateImageBasedPoreNetwork a partir des
%informations fournies par le script python PNMWatershedAnGeometryAnalysis


inputContainerMap=containers.Map;

%Proprietes de l'image 3D du materiau
inputContainerMap('MaterialImage') = myImage;               %image 3D du materiau (tableau de voxels)
inputContainerMap('VoxelEdgeLength') = 1e-6;                 %taille d'un voxel en mètres

%Proprietes des pores (pores = espaces vides séparés par des watershed lines)
inputContainerMap('PoresImage') =  imagePores;             %image des pores labelises
inputContainerMap('PorePropertyVolume') = poreVolumes;              %tableau nPore contenant le volume de chaque pore
inputContainerMap('PorePropertyCenter') = poreCenters;             %tableau (nPore,3) contenant le barycentre de chaque pore

%Proprietes des internal links (internal link = interface entre deux pores : watershed lines)
%inputContainerMap('InternalLinkImage') = imageLiensDilates ;               %image des internal links labelises puis dilates d'un pixel
inputContainerMap('InterfaceToPore') = interfaceToPore;
inputContainerMap('InternalLinkPropertyDiameter') = internalLinkDiameters;                %tableau contenant le diametre de chaque internal links
inputContainerMap('InternalLinkPropertyCenter') = internalLinkBarycenters;              %tableau (nInternalLink,3) contenant le barycentre de chaque internal links

% Proprietes des liens frontieres (liens frontieres = slices des pores sur les bords de l'image)

boundaryLinkPropertyCenter = cell(1,6);
boundaryLinkPropertyCenter{1} = boundaryCenters0;
boundaryLinkPropertyCenter{2} = boundaryCenters1;
boundaryLinkPropertyCenter{3} = boundaryCenters2;
boundaryLinkPropertyCenter{4} = boundaryCenters3;
boundaryLinkPropertyCenter{5} = boundaryCenters4;
boundaryLinkPropertyCenter{6} = boundaryCenters5;
inputContainerMap('BoundaryLinkPropertyCenter') = boundaryLinkPropertyCenter;    %cell (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) contenant les centres des liens frontieres (tableaux (nPore,3) avec NaN si le pore i n'intersecte pas la slice)

boundaryLinkPropertyDiameter = cell(1,6);
boundaryLinkPropertyDiameter{1} = boundaryDiameters0;
boundaryLinkPropertyDiameter{2} = boundaryDiameters1;
boundaryLinkPropertyDiameter{3} = boundaryDiameters2;
boundaryLinkPropertyDiameter{4} = boundaryDiameters3;
boundaryLinkPropertyDiameter{5} = boundaryDiameters4;
boundaryLinkPropertyDiameter{6} = boundaryDiameters5;
inputContainerMap('BoundaryLinkPropertyDiameter') = boundaryLinkPropertyDiameter;  %cell (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) contenant les diametres des liens frontieres (tableaux (1,nPore) avec NaN si le pore i n'intersecte pas la slice)


 
 
 





 