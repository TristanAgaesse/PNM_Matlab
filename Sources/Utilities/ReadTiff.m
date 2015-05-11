function image=ReadTiff(filename)
%READTIFF Read a 3D image from a tiff file. Uses the python library tifffile.py  
%   input : filename
%   output : image

    disp('Reading tiff file')
    
    %Get directory to ReadTiff.py
    scriptDirectory = fileparts(mfilename('fullpath'));
    
    pythonScript=fullfile(scriptDirectory,'ReadTiff.py');
    
    %Get absolute name of filename
    curDir = pwd;
    [pathstr, tiffFileName] = fileparts(filename);
    if ~isempty(pathstr)
        cd(pathstr)
    end
    tiffFileDir=pwd; % full (absolute) path
    cd(curDir) % get back to where you were
    
    
    
    %call python script which reads the image in the .tiff file and writes
    %a temporary .mat file
    system(['python',' ',pythonScript,' ',scriptDirectory,' ',tiffFileDir,' ',tiffFileName])
    
    %load the image in the temporary .mat file
    matFileName=fullfile(tiffFileDir,strcat('temp_',tiffFileName,'.mat'));
    load(matFileName)
    
    delete(matFileName)
end
