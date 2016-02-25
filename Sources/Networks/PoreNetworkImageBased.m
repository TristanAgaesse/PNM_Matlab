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
        PoreVoxelLabelIndices
    end
    
    methods

        function pore_network_image_based=PoreNetworkImageBased(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry,voxelEdgeLength,materialImage,poresImage)
            %constructeur
            %input : dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter
            %output : pore_network_image_based
            pore_network_image_based=pore_network_image_based@PoreNetworkEuclidien(dimension,pores,owners,neighbours,boundaries,poreCenter,linkCenter,myGeometry);
            
            pore_network_image_based.VoxelEdgeLength=voxelEdgeLength;
            pore_network_image_based.MaterialImage=materialImage;
            [labelEnds,orderLabels,labelIndices]=PoreNetworkImageBased.ParseLabeledImage(poresImage);
            pore_network_image_based.PoreVoxelEnds=labelEnds;
            pore_network_image_based.PoreVoxelOrder=orderLabels;
            pore_network_image_based.PoreVoxelLabelIndices=labelIndices;
        end
        
        function image=GetImage(network)
            image=network.MaterialImage;
        end
        
        function imageSize=GetImageSize(network)
            imageSize=size(network.MaterialImage);
        end
        
        function linearIndices=GetVoxelOfPore(network,iPore)
%             numLabel = network.PoreVoxelLabelIndices(iPore);
%             vRange=network.PoreVoxelEnds(numLabel)+1:network.PoreVoxelEnds(numLabel+1);
%             linearIndices=network.PoreVoxelOrder(vRange);
            linearIndices=PoreNetworkImageBased.GetVoxelsOfLabel(iPore,network.PoreVoxelEnds,network.PoreVoxelOrder,network.PoreVoxelLabelIndices);
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
        
        
        function diameter = ComputeAllLinkDiameter(network)
            diameter = 2*network.GetLinkData('CapillaryRadius');
        end
        
    end
    
    methods (Static=true )
        function [labelEnds,orderLabels,labelIndices] = ParseLabeledImage(labelImage)
            
            labelImage = reshape(labelImage,[1,numel(labelImage)]);
        	[sortedLabels,orderLabels] = sort(labelImage);
            labelEnds=find(sortedLabels([2:end,1])-sortedLabels) ; 
%             assert(labelEnds(end)<numel(labelImage))
            labelEnds=[0,labelEnds];
            
            labelIndices = zeros(1,max(labelImage(:))+1);
            labelIndices(sortedLabels(labelEnds(2:end))+1)=1:length(labelEnds)-1;
        end
        
        function linearIndices=GetVoxelsOfLabel(numLabel,labelEnds,orderLabels,labelIndices)
            
            if numLabel>=0 && numLabel<=length(labelIndices)
                
                iLabel = labelIndices(numLabel+1);
                if iLabel>0 || numLabel==0
                	vRange=labelEnds(iLabel)+1:labelEnds(iLabel+1);
                	linearIndices=orderLabels(vRange);
            
                else
                    linearIndices=[];
                end
            
            else
                linearIndices=[];
            end
            
        end
        
        
        function Test_ParseLabeledImage()

            image=[0, 0, 2, 1,2,1, 4];
            [labelEnds,orderLabels,labelIndices] = PoreNetworkImageBased.ParseLabeledImage(image);
            voxels0=PoreNetworkImageBased.GetVoxelsOfLabel(0,labelEnds,orderLabels,labelIndices);
            voxels1=PoreNetworkImageBased.GetVoxelsOfLabel(1,labelEnds,orderLabels,labelIndices);
            voxels2=PoreNetworkImageBased.GetVoxelsOfLabel(2,labelEnds,orderLabels,labelIndices);
            voxels3=PoreNetworkImageBased.GetVoxelsOfLabel(3,labelEnds,orderLabels,labelIndices);
            voxels4=PoreNetworkImageBased.GetVoxelsOfLabel(4,labelEnds,orderLabels,labelIndices);

            assert( isequal(voxels0,[1, 2]))
            assert( isequal(voxels1,[4, 6]))  
            assert( isequal(voxels2,[3,5])) 
            assert( isequal(voxels3,[]))  
            assert( isequal(voxels4,7))  
            
            image=[[0, 0, 2]; [1,2,1];[4 5 5] ];
            [labelEnds,orderLabels,labelIndices] = PoreNetworkImageBased.ParseLabeledImage(image);
            voxels0=PoreNetworkImageBased.GetVoxelsOfLabel(0,labelEnds,orderLabels,labelIndices);
            voxels1=PoreNetworkImageBased.GetVoxelsOfLabel(1,labelEnds,orderLabels,labelIndices);
            voxels2=PoreNetworkImageBased.GetVoxelsOfLabel(2,labelEnds,orderLabels,labelIndices);
            voxels3=PoreNetworkImageBased.GetVoxelsOfLabel(3,labelEnds,orderLabels,labelIndices);
            voxels4=PoreNetworkImageBased.GetVoxelsOfLabel(4,labelEnds,orderLabels,labelIndices);
            voxels5=PoreNetworkImageBased.GetVoxelsOfLabel(5,labelEnds,orderLabels,labelIndices);
            
            assert( isequal(voxels0,[1, 4]))
            assert( isequal(voxels1,[2, 8]))  
            assert( isequal(voxels2,[5,7])) 
            assert( isequal(voxels3,[]))  
            assert( isequal(voxels4,3)) 
            assert( isequal(voxels5,[6,9]))
            
            
            
            disp('Test OK')
        end

    end

    
end