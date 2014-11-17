function WriteTiffStack(filename,image)
%input : filename,image

%Writes a 3D voxels matrix to a tiff stack file.

    imwrite(image(:,:,1), filename)
    for k = 2:size(image,3)
        imwrite(image(:,:,k), filename, 'writemode', 'append');
    end

end