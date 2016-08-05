classdef FileWriterPNM_CEAcpp <handle
    %FILEWRITERXML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = private)
        FileID
        FileName
    end
    
    methods
        % Construct an object and 
        % save the file ID  
        function file_writer = FileWriterPNM_CEAcpp(filename) 
            file_writer.FileName=strcat(filename,'.txt');
            a=dir(file_writer.FileName);
            if ~isempty(a)
                if a.bytes~=0
                    user_input=input('Fichier deja existant. Taper 1 pour effacer le contenu, 2 pour ecrire a la suite. \n');
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
        
        function Write(file_writer,babe_struct)
            
            fid=file_writer.FileID;
            nl = sprintf('\n'); %code pour new line
            
            nPore=babe_struct.NumberOfPores;
            nThroat=babe_struct.NumberOfThroats;
            nAgglomerate=babe_struct.NumberOfAgglomerates;
            nConnexion=babe_struct.NumberOfConnexions;
            
            %Partie du fichier EN TETE
            fwrite(fid, ['Geometry File' nl]);
            
            
            %Partie du fichier NetworkSize
            data=babe_struct.NetworkPosition;
            fwrite(fid, [nl 'Network Position' nl]);
            %dlmwrite(file_writer.FileName, data,'-append','roffset', 1, 'delimiter', ' ')
            fwrite(fid, [num2str(data) nl]);
            
            data=babe_struct.NetworkLength;
            fwrite(fid, [nl 'Network Length' nl]);
            %dlmwrite(file_writer.FileName, data,'-append','roffset', 1, 'delimiter', ' ')
            fwrite(fid, [num2str(data) nl]);
            
            %Partie du fichier Pores
            fwrite(fid, [nl 'Number Of Pores :' nl]);
            fwrite(fid, [num2str(nPore) nl ]);
            fwrite(fid, [nl '(Id Diam X Y Z Domain)'  nl ]);
            
            for iPore=1:nPore
                
                data=babe_struct.Pores.NumericData{iPore};
                fwrite(fid, [num2str(data(1)) ' ' num2str(data(2:end),' %e')]);
                string1=babe_struct.Pores.StringData{iPore}{1};
                string2=babe_struct.Pores.StringData{iPore}{2};
                fwrite(fid, [' ' string1 ' ' string2 nl]);
            end
            
            
            %Partie du fichier Throats
            fwrite(fid, [nl 'Number Of Throats :' nl]);
            fwrite(fid, [num2str(nThroat) nl ]);
            fwrite(fid, [nl '(Id Diam X Y Z Domain Label)'  nl ]);
            
            for iThroat=1:nThroat
                
                data=babe_struct.Throats.NumericData{iThroat};
                fwrite(fid, [num2str(data(1)) ' ' num2str(data(2:end),' %e')]);
                string1=babe_struct.Throats.StringData{iThroat}{1};
                string2=babe_struct.Throats.StringData{iThroat}{2};
                string3=babe_struct.Throats.StringData{iThroat}{3};
                fwrite(fid, [' ' string1 ' ' string2 ' ' string3 nl]);
            end
            
            
            %Partie du fichier Agglomerates
            fwrite(fid, [nl 'Number Of Agglomerates :' nl]);
            fwrite(fid, [num2str(babe_struct.NumberOfAgglomerates) nl ]);
            fwrite(fid, [nl '(Id Diam X Y Z Domain)'  nl ]);
            
            for iAgglomerate=1:nAgglomerate
                
                data=babe_struct.Pores.NumericData{iAgglomerate};
                fwrite(fid, [num2str(data(1)) ' ' num2str(data(2:end),' %e')]);
                string1=babe_struct.Agglomerates.StringData{iAgglomerate};
                fwrite(fid, [string1 nl]);
            end
            
            
            
            %Partie du fichier Connexions
            fwrite(fid, [nl 'Number Of Connexions :' nl]);
            fwrite(fid, [num2str(babe_struct.NumberOfConnexions) nl ]);
            fwrite(fid, [nl '(Id Diam X Y Z Domain)'  nl ]);

            for iConnexion=1:nConnexion
                
                data=babe_struct.Pores.NumericData{iConnexion};
                fwrite(fid, [num2str(data(1)) ' ' num2str(data(2:end),' %e')]);
                string1=babe_struct.Connexions.StringData{iConnexion};
                fwrite(fid, [string1 nl]);
            end

            
            %Partie du fichier PoreLinking
            
            fwrite(fid, [nl 'Pore Linking' nl]);
            for iPore=1:nPore
                data=babe_struct.PoreLinking{iPore};
                dlmwrite(file_writer.FileName, data,'-append','delimiter', ' ');
            end
            
            
            %Partie du fichier ThroatLinking
            
            fwrite(fid, [nl 'Throat Linking' nl]);
            for iThroat=1:nThroat
                data=babe_struct.ThroatLinking{iThroat};
                dlmwrite(file_writer.FileName, data,'-append','delimiter', ' ');
            end
            
            
            %Partie du fichier AgglomerateToPoreLinking
            
            fwrite(fid, [nl 'Agglomerate to Pore Linking' nl]);
            for iAgglomerate=1:nAgglomerate
                data=babe_struct.AgglomeratesToPoreLinking{iAgglomerate};
                dlmwrite(file_writer.FileName, data,'-append','delimiter', ' ');
            end
            
            
            %Partie du fichier ConnexionLinking 
            
            fwrite(fid, [nl 'Connexion Linking' nl]);
            for iConnexion=1:nConnexion
                data=babe_struct.ConnexionLinking{iConnexion};
                dlmwrite(file_writer.FileName, data,'-append','delimiter', ' ');
            end

        end
        
        
        % Delete methods are always called before a object 
        % of the class is destroyed 
        function delete(obj)
            fclose(obj.FileID);
        end 
        
    end
    
end

