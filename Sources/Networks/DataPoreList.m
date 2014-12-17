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
            assert(length(data)==data_pore_list.NombrePore,'Un pore data doit �tre un tableau de taille NombrePores');
            data_pore_list.PoreDatas.(name)=data;
        end

        function RemoveData(data_list,name)
            data_list.PoreDatas=rmfield(data_list.PoreDatas,name);
        end
        
        
    end
    
end

