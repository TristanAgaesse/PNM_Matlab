classdef PoreNetworkImageBased < PoreNetworkEuclidien
    %PoreNetworkImageBased : Subclass of PoreNetworkEuclidien
    % L'element en plus par rapport à PoreNetworkEuclidien est qu'il y a une 
    % image 3D (voxels) associé au réseau de pores. Le réseau de pores est 
    % construit à partir de cette image par watershed segmentation des pores.
    
    properties 
        VoxelEdgeLength
        MaterialImage
        PoreVoxelEnds
        PoreVoxelOrder
    end
    
    methods

        function pore_network_image_based=PoreNetworkImageBased(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry,voxelEdgeLength,materialImage,poresImage)
            %constructeur
            %input : dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter
            %output : pore_network_image_based
            pore_network_image_based=pore_network_image_based@PoreNetworkEuclidien(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry);
            
            pore_network_image_based.VoxelEdgeLength=voxelEdgeLength;
            pore_network_image_based.MaterialImage=materialImage;
            [labelEnds,orderLabels]=PoreNetworkImageBased.ParseLabeledImage(poresImage);
            pore_network_image_based.PoreVoxelEnds=labelEnds;
            pore_network_image_based.PoreVoxelOrder=orderLabels;
        end
        
        function image=GetImage(network)
            image=network.MaterialImage;
        end
        
        function imageSize=GetImageSize(network)
            imageSize=size(network.MaterialImage);
        end
        
        function linearIndices=GetVoxelOfPore(network,iPore)
            vRange=network.PoreVoxelEnds(iPore)+1:network.PoreVoxelEnds(iPore+1);
            linearIndices=network.PoreVoxelOrder(vRange);
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
    
    methods (Static=true )
        function [labelEnds,orderLabels] = ParseLabeledImage(labelImage)
            
            labelImage = reshape(labelImage,[1,numel(labelImage)]);
        	[sortedLabels,orderLabels] = sort(labelImage);
            labelEnds=find(sortedLabels([2:end,1])-sortedLabels) ;   
        end
        
    end
end