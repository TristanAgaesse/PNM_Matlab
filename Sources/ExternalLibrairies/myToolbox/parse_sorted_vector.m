function C = parse_sorted_vector( sortedVector )
%parse_sorted_vector Donne les valeurs uniques d'un tableau ordonné avec
%les indices min et max de ces valeurs
% input : sortedVector = un tableau ordonné
% output : - C : C{i}(1)= valeur unique numero i,  
%               [C{i}(2):C{i}(3)]=indices des valeurs C{i}(1) dans sortedVector

    
    currentValue=sortedVector(1);
    numCurrentValue=1;
    
    C=cell(1,length(sortedVector));
    
    %une recherche binaire rapide pour la première valeur, utile si la première valeur
    %est présente dans la majorité du tableau
    [first,last]= find_sorted_vector(sortedVector,currentValue);
    C{numCurrentValue}(1)=currentValue;
    C{numCurrentValue}(2)=first;
    C{numCurrentValue}(3)=last;
    
    if last<length(sortedVector)
    
        for i=(last+1):length(sortedVector)
            %puis un parcours linéaire du tableau pour les autres valeurs
            
            if sortedVector(i)~=currentValue
                
                numCurrentValue=numCurrentValue+1;
                currentValue=sortedVector(i);
                
                C{numCurrentValue}(1)=currentValue;
                C{numCurrentValue}(2)=i;
                
                C{numCurrentValue-1}(3)=i-1;

            end
            
        end
        
        C{numCurrentValue}(3)=length(sortedVector);
        
    end
 
    C=C(1:(numCurrentValue));
end

