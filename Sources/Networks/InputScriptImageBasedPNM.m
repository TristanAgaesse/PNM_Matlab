%Script construisant l'input de CreateImageBasedPoreNetwork a partir des
%informations fournies par le script python PNMWatershedAnGeometryAnalysis


inputContainerMap=containers.Map;

%Proprietes de l'image 3D du materiau
inputContainerMap('MaterialImage') = myImage;               %image 3D du materiau (tableau de voxels)
inputContainerMap('VoxelEdgeLength') = 1e-6;                 %taille d'un voxel en mètres


%Proprietes des pores (pores = espaces vides séparés par des watershed lines)
inputContainerMap('PoresImage') =  imagePores;             %image des pores labelises
inputContainerMap('PorePropertyVolume') = poreVolumes;              %tableau nPore contenant le volume de chaque pore (en nombre de voxels)
inputContainerMap('PorePropertyCenter') = poreCenters;             %tableau (nPore,3) contenant le barycentre de chaque pore
inputContainerMap('PorePhase') = porePhase;         %tableau nPore contenant la phase à laquelle appartient chaque pore
%Other pore properties
inputContainerMap('OtherPoreProperties') = struct;


%Proprietes des internal links (internal link = interface entre deux pores : watershed lines)
inputContainerMap('InterfaceToPore') = interfaceToPore;
inputContainerMap('InternalLinkCapillaryRadius') = internalLinkCapillaryRadius; %tableau contenant le rayon déduit de la distance map de chaque internal links
inputContainerMap('InternalLinkPropertyCenter') = internalLinkBarycenters;              %tableau (nInternalLink,3) contenant le barycentre de chaque internal links
%Other internal links properties
myStruct=struct;
myStruct.('GeometricSurface') = internalLinkGeometricSurface; %tableau contenant la surface chaque internal links (en nombre de voxels)
inputContainerMap('OtherInternalLinkProperties') = myStruct;


%Proprietes des liens frontieres (liens frontieres = slices des pores sur les bords de l'image)
boundaryLinkPropertyCenter = cell(1,6);
for i = 0:5
    boundaryLinkPropertyCenter{i+1} = eval(sprintf('boundaryCenters%d',i));
end
inputContainerMap('BoundaryLinkPropertyCenter') = boundaryLinkPropertyCenter;    %cell (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) contenant les centres des liens frontieres (tableaux (nPore,3) avec NaN si le pore i n'intersecte pas la slice)

boundaryLinkPropertyCapillaryRadius = cell(1,6);
for i = 0:5
    boundaryLinkPropertyCapillaryRadius{i+1} = eval(sprintf('boundaryCapillaryRadius%d',i));
end
inputContainerMap('BoundaryLinkPropertyCapillaryRadius') = boundaryLinkPropertyCapillaryRadius;  %cell (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) contenant les rayons capillaires des liens frontieres (tableaux (1,nPore) avec NaN si le pore i n'intersecte pas la slice)
%Other internal links properties
boundaryLinkPropertyGeometricSurface = cell(1,6);
for i = 0:5
    boundaryLinkPropertyGeometricSurface{i+1} = eval(sprintf('boundaryGeometricSurface%d',i));
end
myStruct=struct;
myStruct.('GeometricSurface') = boundaryLinkPropertyGeometricSurface;  %cell (1,6) (slices dans l'ordre Xmin Xmax Ymin Ymax Zmin Zmax) contenant la surface des liens frontieres en nombre de voxels (tableaux (1,nPore) avec NaN si le pore i n'intersecte pas la slice)
inputContainerMap('OtherBoundaryLinkProperties') = myStruct;


