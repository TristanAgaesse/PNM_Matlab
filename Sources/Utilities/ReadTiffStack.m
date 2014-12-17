function imageStack=ReadTiffStack(filename)
%ReadTiffStack Reads a tiff stack which represents a 3D image
%Input : - filename
%Output : - imageStack : 3D matrix

    info = imfinfo(filename);
    imageStack = [];
    numberOfImages = length(info);
    for k = 1:numberOfImages
    currentImage = imread(filename, k, 'Info', info);
    imageStack(:,:,k) = currentImage;
    end

end