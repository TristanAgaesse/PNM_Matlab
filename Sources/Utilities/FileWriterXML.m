classdef FileWriterXML <handle
    %FILEWRITERXML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = private)
        FileID
        FileName
    end
    
    methods
        % Construct an object and 
        % save the file ID  
        function file_writer = FileWriterXML(filename) 
            file_writer.FileName=filename;
            file_writer.FileID = fopen(filename,'a');
        end
        
        function Write(fileWriter,contenu)
            wPref.StructItem = false;
            xml_write(fileWriter.FileName,contenu,'PNM_Matlab',wPref); %utilise xml_write
        end
        
        
        % Delete methods are always called before a object 
        % of the class is destroyed 
        function delete(obj)
            fclose(obj.FileID);
        end 
        
    end
    
end

