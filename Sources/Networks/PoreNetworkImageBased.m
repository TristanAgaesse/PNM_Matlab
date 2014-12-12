classdef PoreNetworkImageBased < PoreNetworkEuclidien
    %PoreNetworkImageBased : Subclass of PoreNetworkEuclidien
    % L'element en plus par rapport à PoreNetworkEuclidien est qu'il y a une 
    % image 3D (voxels) associé au réseau de pores. Le réseau de pores est 
    % construit à partir de cette image par watershed segmentation des pores.
    
    properties 
        VoxelEdgeLength
        MaterialImage
        ParsedPore
        OrderPore
    end
    
    methods

        function pore_network_image_based=PoreNetworkImageBased(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry,voxelEdgeLength,materialImage,parsedPores,orderPore)
            %constructeur
            %input : dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter
            %output : pore_network_image_based
            pore_network_image_based=pore_network_image_based@PoreNetworkEuclidien(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry);
            
            pore_network_image_based.VoxelEdgeLength=voxelEdgeLength;
            pore_network_image_based.MaterialImage=materialImage;
            pore_network_image_based.ParsedPore=parsedPores;
            pore_network_image_based.OrderPore=orderPore;
            
        end
        
        function image=GetImage(network)
            image=network.MaterialImage;
        end
        
        function imageSize=GetImageSize(network)
            imageSize=size(network.MaterialImage);
        end
        
        function linearIndices=GetVoxelOfPore(network,iPore)
            linearIndices=network.OrderPore(network.ParsedPore{iPore}(2):network.ParsedPore{iPore}(3));
        end
        
        function image=GetImagePoreData(network,poreDataName)
            %input :network,poreDataName
            %output: image
            assert(isfield(network.GetPoreDataList,poreDataName),'wrong pore data name')
            
            image=zeros(network.GetImageSize);
            
            for iPore=1:network.GetNumberOfPores
                image(network.GetVoxelOfPore(iPore))=network.GetPoreDataList.(poreDataName)(iPore);
            end
        end
        
        
    end

end