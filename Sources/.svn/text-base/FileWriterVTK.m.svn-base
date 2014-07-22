classdef FileWriterVTK <handle
    %FILEWRITERXML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = private)
        FileID
        FileName
    end
    
    methods
        % Construct an object and 
        % save the file ID  
        function file_writer = FileWriterVTK(filename) 
            file_writer.FileName=strcat(filename,'.vtk');
            a=dir(file_writer.FileName);
            if ~isempty(a)
                if a.bytes~=0
                    user_input=input('Fichier déjà existant. Taper 1 pour effacer le contenu, 2 pour écrire à la suite. \n');
                    if user_input==1
                        delete(file_writer.FileName);
                    end
                end
            end
            fid=fopen(file_writer.FileName,'a');
            
            if fid == -1
                error('Cannot open file for writing.');
            else
                file_writer.FileID = fid;
            end
            
        end
        
        function Write(file_writer,vtk_struct)
            
            nl = sprintf('\n'); %code pour new line
            fid=file_writer.FileID;
            %EN TETE
            fwrite(fid, ['# vtk DataFile Version 3.0' ...
                nl file_writer.FileName ...
                nl 'ASCII' ...
                nl 'DATASET POLYDATA' nl ]);
         
            %POINTS
            data=vtk_struct.Points;
            nPoints=length(data(:,1));
            fwrite(fid, ['POINTS ' num2str(nPoints) ' double'  ]);
            dlmwrite(file_writer.FileName, data,'-append','roffset', 1, 'delimiter', ' ')
%             for i=1:nPoints
%                 fwrite(fid, num2str(data(i,:),'%G '));
%                 fwrite(fid, nl);
%             end
            
            %VERTICES
            data=vtk_struct.Vertices;
            nVertices=length(data(:,1));
            fwrite(fid, [nl 'VERTICES ' num2str(nVertices) ' ' num2str(2*nVertices)  ]);
            nChiffresSignificatifs=ceil(log(nPoints)/log(10))+1;
            dlmwrite(file_writer.FileName, data,'-append','roffset', 1, 'delimiter', ' ','precision', nChiffresSignificatifs)


            %LINES
            data=vtk_struct.Lines;
            nLines=length(data(:,1));
            fwrite(fid, [nl 'LINES ' num2str(nLines) ' ' num2str(3*nLines)  ]);
            nChiffresSignificatifs=ceil(log(nPoints)/log(10))+1;
            dlmwrite(file_writer.FileName, data,'-append','roffset', 1, 'delimiter', ' ','precision', nChiffresSignificatifs)
%             for i=1:nLines
%                 fwrite(fid, num2str(data(i,:)));
%                 fwrite(fid, nl);
%             end
            
            %POLYGONS
            data=vtk_struct.Polygons;
            nPolygons=length(data(:,1));
            n=sum(cellfun('length',data));
            fwrite(fid, ['POLYGONS ' num2str(nPolygons) ' ' num2str(n) nl  ]);
            for i=1:nPolygons
                fwrite(fid, num2str(data{i}));
                fwrite(fid, nl);
            end
            
            %POINT DATA
            fwrite(fid, ['POINT_DATA ' num2str(nPoints) nl  ]);
            list_data=vtk_struct.PointData;
            names=fieldnames(list_data);
            for i=1:length(names)
                fwrite(fid, ['SCALARS ' names{i} ' double 1' nl  ]);
                fwrite(fid, ['LOOKUP_TABLE default' nl  ]);
                data=list_data.(names{i});
                dlmwrite(file_writer.FileName, data,'-append','roffset', 1, 'delimiter', ' ')
%                 for j=1:length(data(:,1))
%                     fwrite(fid, num2str(data(j,:),'%G '));
%                     fwrite(fid, nl);
%                 end
            end
            
            %CELL DATA
            fwrite(fid, ['CELL_DATA ' num2str(nLines+nPolygons+nVertices) nl  ]);
            list_data=vtk_struct.CellData;
            names=fieldnames(list_data);
            for i=1:length(names)
                fwrite(fid, ['SCALARS ' names{i} ' double 1' nl  ]);
                fwrite(fid, ['LOOKUP_TABLE default' nl  ]);
                data=list_data.(names{i});
                dlmwrite(file_writer.FileName, data,'-append','roffset', 1, 'delimiter', ' ')
%                 for j=1:length(data(:,1))
%                     fwrite(fid, num2str(data(j,:),'%G '));
%                     fwrite(fid, nl);
%                 end
            end
            
           
        end
        
        
        % Delete methods are always called before a object 
        % of the class is destroyed 
        function delete(obj)
            fclose(obj.FileID);
        end 
        
    end
    
end

