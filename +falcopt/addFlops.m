function c = addFlops(a, b)

    fieldsA = fieldnames(a);
    fieldsB = fieldnames(b);
    fields = unique({fieldsA{:}, fieldsB{:}});
    c = struct();
    for f=1:length(fields)
        if any(strcmp(fields{f}, fieldsA)) && any(strcmp(fields{f}, fieldsB))
            c.(fields{f}) = a.(fields{f}) + b.(fields{f});
        elseif any(strcmp(fields{f}, fieldsA))
            c.(fields{f}) = a.(fields{f});
        else
            c.(fields{f}) = b.(fields{f});
        end
    end

end