function SEM_value = SEM(values,dim)

if exist('dim','var') && dim==2
    values = values';
    SEM_value = SEM(values);
    SEM_value = SEM_value';
    return
end

if all(size(values)>1)
    SEM_value = NaN(1,size(values,2));
    for i=1:size(values,2)
        idx = ~isnan(values(:,i));
        if sum(idx)<2, continue, end
        SEM_value(i) = SEM(values(idx,i));

    end
    return
end

values(isnan(values)) = [];

if length(values)>2
    SEM_value = std(ToColumn(values))/sqrt(length(values));
elseif length(values)==2
    SEM_value = abs(diff(values)/2);
else
    SEM_value = NaN;
end
