Scripts utilises pour l'article : 
Comparison between pore network simulations and X-ray tomography images of ex-situ 
water distributions in a gas diffusion layer of proton exchange membrane fuel cell 

Tristan Agaesse
01/10/2015


Pretraitement des images experimentales:
  - Detection Boundary conditions : identifier GDL, sample holder, membrane 
  superieure et inferieure
  
  - Rassembler (fusion) les images expérimentales  donnant repartitions d'eau a 
  differentes pressions pour avoir tout sur une seule image
  
  - detection des fibres : separer binder et fibres. N'a pas ete utilise dans l'article
  
  
Lancement des simulations :
  - Extraction du reseau de pores
  
  - Simulations Matlab sur le reseau de pores
  
  - Simulations Full Morphology


Post-processing des simulations :
  - Courbes de pression capillaire
  - Voxel Match
  - Rendu 3D des repartitions d'eau (paraview calculator contour et VirtualMaterial.Resize)
 
  