classdef DataPoreList <  Data
    %DATAPORELIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        PoreDatas
        NombrePore %nombre de pores
    end
    
    methods
        function data_pore_list=DataPoreList(nPore)
            data_pore_list.NombrePore=nPore;
            data_pore_list.PoreDatas=struct;
        end
        
        function AddData(data_pore_list,data,name)
            
            if size(data,1)~=data_pore_list.NombrePore
                data=transpose(data);
                assert(size(data,1)==data_pore_list.NombrePore,'Un pore data doit etre un tableau de taille NombrePores')
            end            
            data_pore_list.PoreDatas.(name)=data;
        end

        function RemoveData(data_list,name)
            data_list.PoreDatas=rmfield(data_list.PoreDatas,name);
        end
        
        
    end
    
end

