function WriteTiffStackNew( matlab3dImage, filename )
%WriteTiffStack Write a tiff stack sile from a Matlab 3D array
%Input : - matlab3dImage, filename
    
    matlab3dImage=uint32(matlab3dImage);

    t = Tiff(filename,'w');
    t.close()
    
    t = Tiff(filename,'a');
        
    for k = 1:size(matlab3dImage,3)
        
        tagstruct.ImageLength = size(matlab3dImage,1);
        tagstruct.ImageWidth = size(matlab3dImage,2);
        tagstruct.Photometric = Tiff.Photometric.MinIsWhite;
        tagstruct.BitsPerSample = 32;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Software = 'MATLAB';
        t.setTag(tagstruct)
        
        t.write(matlab3dImage(:,:,k));
        
        t.writeDirectory();
    end
    
    t.close();
end

