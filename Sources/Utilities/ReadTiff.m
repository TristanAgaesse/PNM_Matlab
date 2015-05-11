function image=ReadTiff(filename)
%READTIFF Read a 3D image from a tiff file. Uses the python library tifffile.py  
%   input : filename
%   output : image
    
    assert(ischar(filename),'Input argument of ReadTiff must be the name of a tiff file')
    disp('Reading tiff file')
    
    %Get directory to ReadTiff.py
    scriptDirectory = fileparts(mfilename('fullpath'));
    
    pythonScript=fullfile(scriptDirectory,'ReadTiff.py');
    
    %Get absolute name of filename
    curDir = pwd;
    [pathstr, tiffFileName,ext] = fileparts(filename);
    if ~isempty(pathstr)
        cd(pathstr)
    end
    tiffFileDir=pwd; % full (absolute) path
    cd(curDir) % get back to where you were
    
    
    
    %call python script which reads the image in the .tiff file and writes
    %a temporary .mat file
    
    system(['python',' ',pythonScript,' ',scriptDirectory,' ',tiffFileDir,' ',tiffFileName,ext])
    
    %load the image in the temporary .mat file
    matFileName=fullfile(tiffFileDir,strcat('temp_',tiffFileName,'.mat'));
    load(matFileName)
    
    image=temp_image;    
    
    delete(matFileName)
end
