function fileContent= read_XLSX( filename )
    %READ_XLSX Renvoie une structure contenant le contenu d'une feuille de calcul au
    %format .xlsx 

    [ndata,headertext,raw]=xlsread(filename);

    numericData=cellfun(@str2double,raw);
    numericData=numericData(2:end,:);

    %gestion d'un bug qui renvoie des nan à certains endroits
    nanIndices=~isnan(ndata);
    numericData(nanIndices)=ndata(nanIndices);

    %creation de fileContent
    nColumn=length(numericData(1,:));
    fileContent=struct;
    for iColumn=1:nColumn
        columName=headertext{1,iColumn};
        columName(columName==' ') = '';
        fileContent.(columName)=numericData(:,iColumn);
    end

end


