classdef FileWriterFreecad <handle
    %FileWriterFreecad Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = private)
        FolderName
        ParentFolderName
    end
    
    methods
        % Construct an object and 
        % save the file ID  
        function file_writer = FileWriterFreecad(folderName) 
            file_writer.FolderName=folderName;
            file_writer.ParentFolderName=pwd;
            mkdir(folderName)
            cd(folderName)
%             a=dir(file_writer.FileName);
%             if ~isempty(a)
%                 if a.bytes~=0
%                     user_input=input('Fichier d�j� existant. Taper 1 pour effacer le contenu, 2 pour �crire � la suite. \n');
%                     if user_input==1
%                         delete(file_writer.FileName);
%                     end
%                 end
%             end
%             fid=fopen(file_writer.FileName,'a');
%             
%             if fid == -1
%                 error('Cannot open file for writing.');
%             else
%                 file_writer.FileID = fid;
%             end
        end
        
        
        
        function Write(file_writer,vtk_struct)
                        
            %points.txt
            fid=fopen('points.txt','w');
            data=vtk_struct.Points;
            nPoints=length(data(:,1));
            dlmwrite('points.txt', data,'-append', 'delimiter', ' ')
            fclose(fid);

            %vertices.txt
            fid=fopen('vertices.txt','w');
            data=vtk_struct.Vertices;
            nVertices=length(data(:,1));
            nChiffresSignificatifs=ceil(log(nPoints)/log(10))+1;
            dlmwrite('vertices.txt', data,'-append', 'delimiter', ' ','precision', nChiffresSignificatifs)
            fclose(fid);

            %fibres.txt
            fid=fopen('fibres.txt','w');
            data=vtk_struct.Lines;
            nLines=length(data(:,1));
            nChiffresSignificatifs=ceil(log(nPoints)/log(10))+1;
            dlmwrite('fibres.txt', data,'-append','delimiter', ' ','precision', nChiffresSignificatifs)
            fclose(fid);            
            
            %radius.txt
            fid=fopen('radius.txt','w');
            list_data=vtk_struct.PointData;
            data=list_data.Edge_FiberDiameter;
            dlmwrite('radius.txt', data,'-append','delimiter', ' ')
            fclose(fid);   
            
            cd(file_writer.ParentFolderName)
        end
        
        

    end
    
end

