classdef FileReaderXML < handle
    %XMLFILEREADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = private)
        FileID
        FileName
    end
    
    methods
        function file_reader = FileReaderXML(filename) 
            file_reader.FileID = fopen(filename,'r');
            file_reader.FileName=filename;
        end
        
        
        function file_content=Read(file_reader)
            file_content = xml_read(file_reader.FileName);
        end
        
        function delete(obj)
         fclose(obj.FileID);
        end 
        
    end
    
end

