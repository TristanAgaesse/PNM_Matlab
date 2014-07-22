function  [first,last]= find_sorted_vector(sortedVector,searchfor)
%find_sorted_vector Searches for all occurrences of searchfor in a sorted vector
%input: sortedVector,searchfor
%output : [first,last] : intervale des indices où est searchfor
    a=1;
    first=numel(sortedVector);
    last=1;
    d=numel(sortedVector);
    while (a+1<first||last+1<d)
        lw=(floor((a+first)/2));
        if (sortedVector(lw)<searchfor)
            a=lw;
        else
            first=lw;
        end
        lw=(floor((last+d)/2));
        if (sortedVector(lw)<=searchfor)
            last=lw;
        else
            d=lw;
        end
    end

end

