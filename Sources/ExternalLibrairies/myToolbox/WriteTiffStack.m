function WriteTiffStack( matlab3dImage, filename )
%WriteTiffStack Write a tiff stack sile from a Matlab 3D array
%Input : - matlab3dImage
%Output : - filename 

    imwrite(matlab3dImage(:,:,1), filename)
    for k = 2:size(matlab3dImage,3)
        imwrite(matlab3dImage(:,:,k), filename, 'writemode', 'append');
    end
end

