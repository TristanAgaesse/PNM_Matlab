classdef DataLinkList <  Data
    %DATALINKLIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LinkDatas
        NombreLink %nombre de liens
    end
    
    methods
        function data_link_list=DataLinkList(nLink)
            data_link_list.NombreLink=nLink;
            data_link_list.LinkDatas=struct;
        end
        
        function AddData(data_link_list,data,name)
            
            if size(data,1)~=data_link_list.NombreLink
                data=transpose(data);
                assert(size(data,1)==data_link_list.NombreLink,'n link data doit etre un tableau de taille NombreLinks')
            end            
            data_link_list.LinkDatas.(name)=data;
        end
       
        function RemoveData(data_link_list,name)
            data_link_list.LinkDatas=rmfield(data_link_list.LinkDatas,name);
        end
        
    end
    
end

