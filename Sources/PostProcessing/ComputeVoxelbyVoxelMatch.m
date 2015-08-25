function voxelByVoxelError=ComputeVoxelbyVoxelMatch(image1,image2,labels)

    nLabel = length(labels);
    voxelByVoxelError = zeros(nLabel,1);

    [expLabelEnds,expOrderLabels,expLabelIndices] = PoreNetworkImageBased.ParseLabeledImage(image1);
    [simuLabelEnds,simuOrderLabels,simuLabelIndices] = PoreNetworkImageBased.ParseLabeledImage(image2);
    
    for i = 1:nLabel
        
        iLabel = labels(i);
        expVoxels = PoreNetworkImageBased.GetVoxelsOfLabel(iLabel,expLabelEnds,expOrderLabels,expLabelIndices);
        simuVoxels = PoreNetworkImageBased.GetVoxelsOfLabel(iLabel,simuLabelEnds,simuOrderLabels,simuLabelIndices);
        
        voxelByVoxelError(i) = numel(intersect(expVoxels,simuVoxels))/numel(expVoxels);
        
    end
    
    
    
end

