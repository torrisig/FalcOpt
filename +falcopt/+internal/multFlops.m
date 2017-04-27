%% multFlops multiply each element of struct a with scalar b.
%
% c = multFlops( a, b)
%
% Inputs:
% a    - struct
% b    - scalar
% Outputs:
% c    - struct
function c = multFlops(a, b)
    if(~isscalar(b))
        error('second term must be scalar')
    end
    fields = fieldnames(a);
    for f=1:length(fields)
        c.(fields{f}) = b*a.(fields{f});
    end
end