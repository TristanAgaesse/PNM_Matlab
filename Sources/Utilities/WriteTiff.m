function WriteTiff(filename,image)
%WRITETIFF Write a 3D image to a tiff file. Uses the python library tifffile.py  
%   input : filename,image

    disp('Writing tiff file')

    image=int16(image);
    
    %Get directory to WriteTiff.py
    scriptDirectory = fileparts(mfilename('fullpath'));
    
    pythonScript=fullfile(scriptDirectory,'WriteTiff.py');
    
    %Get absolute name of filename
    curDir = pwd;
    [pathstr, name] = fileparts(filename);
    if ~isempty(pathstr)
        cd(pathstr)
    end
    outputDir=pwd; % full (absolute) path
    cd(curDir) % get back to where you were
    
    
    %write the image in a .mat file
    matFileName=fullfile(outputDir,strcat('temp_',name,'.mat'));
    save(matFileName,'image','-v7.3')
    
    %call python script which reads the image in the .mat file and write it
    %to tiff format
    system(['python',' ',pythonScript,' ',scriptDirectory,' ',outputDir,' ',name])
    
    delete(matFileName)
    
end

