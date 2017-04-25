function [data, info] = fcn2struct(varargin)
% FCN2STRUCT given a function handle f generate c-code equivalent
%
% Input:
% - matlab function handle
% - options structure containing: nx,nu,nw,precision
%
% Options:
% - 'name': name of variable returned by the function
% - 'structure': can be 'unique' or 'sparse'. Default 'unique'
%
% Outputs:
% - data: the c-code data for the function
% - info: .structure : the structure of the matrix
%         .flops: the flops of the function
%
p = inputParser;
p.addRequired('function');
p.addRequired('opt');
p.addParameter('structure', 'unique');
p.addParameter('name','x');
p.parse(varargin{:});
f = p.Results.function;
o = p.Results.opt;

data = [];
info.flops.add = 0;
info.flops.mul = 0;
info.flops.div= 0;

x = sym('x', [o.nx,1]);
u = sym('u', [o.nu,1]);
w = sym('w', [o.nw,1]);
sl = sym('sl', [50,1]); % ToDo Tommaso change size of slack

static = strcmp(class(f),'double');
switch o.precision
    case 'float'
        digits(9);
    case 'double'
        digits(17);
    otherwise
        
end

if( ~static)
    M = vpa(simplify(f));
else
    M = vpa(f);
end
[n,m] = size(M);
values = [];
indeces = [];
rows = [];
cols = [];
ind = 1;
S = zeros(n,m);

% extract the function/matrix structure
for i= 1:n
    for j= 1:m
        if( M(i,j) ~= 0)
            if ( isempty(values))
                values = [ values; M(i,j)];
                indeces = [indeces, ind];
                S(i,j) = ind;
            elseif( ~any(has(values,M(i,j))) )||(~isequal(p.Results.structure,'unique'))
                values = [ values; M(i,j)];
                ind = ind +1;
                indeces = [indeces, ind];
                S(i,j) = ind;
            else
                S(i,j) = indeces(find(has(values,M(i,j))== 1));
            end
            rows = [ rows; i];
            cols = [ cols; j];
        end
    end
end

info.structure.mat = S;
info.structure.values = values;
info.structure.rows = rows;
info.structure.cols = cols;
info.structure.num = size(values,1);
info.static = static;

fun = @replaceEsp;
for i= 1:size(values)
    
    tmp_c = ccode(values(i));
    tmp = strrep(tmp_c, 't0 = ','');
    
    % substitute the pow(.,n) c function with (.)*(.)*(.)... n-times
    tmp = regexprep(tmp,'pow\((.*?,\d\.\d)\)','${fun($1)}');
   
    % w1 -> w[0] for all variable x,u,w
    for j = o.nw:-1:1
        tmp = strrep(tmp, char(w(j)), sprintf('w[%i]',j-1));
    end
    for j = o.nx:-1:1
        tmp = strrep(tmp, char(x(j)), sprintf('x[%i]',j-1));
    end
    for j = o.nu:-1:1
        tmp = strrep(tmp, char(u(j)), sprintf('u[%i]',j-1));
    end
    %ToDo Tommaso: remove after no slack needed
    for j = 50:-1:1
        tmp = strrep(tmp, char(sl(j)), sprintf('sl[%i]',j-1));
    end
    
    if( static)
        tmp = regexprep( tmp,';','');
        data = [ data sprintf([tmp ', '])];
    else
        data = [data sprintf(['\t' p.Results.name '[%i] = ' tmp '\n'], i-1)];
    end
    
    
end

% change c-function from math.h in case the precision is float
switch o.precision
    case 'float'
        data = regexprep( data, {'sin[^h]' 'sinh' 'asin[^h]' 'asinh[^f]' 'cos[^h]' 'cosh' 'acos[^h]' 'acosh[^f]' 'tan[^h]' ...
                                'tanh' 'atan[^h]' 'atanh[^f]' 'atan2' 'sqrt' 'log[^21]' 'log10' 'exp[^m]' 'expm1' ...
                                'log1p' 'log2' 'hypot' 'erf[^c]' 'erfc' 'gamma[^f]' 'ceil' 'floor' 'mod' ...
                                'round'},...
                                {'sinf(' 'sinhf' 'asinf' 'asinhf' 'cosf(' 'coshf' 'acosf(' 'acoshf' 'tanf(' ...
                                'tanhf' 'atanf'  'atanhf' 'atan2f(' 'sqrtf' 'logf(' 'log10f' 'expf(' 'expm1f' ...
                                'log1pf' 'log2f' 'hypotf' 'erff(' 'erfcf' 'tgammaf(' 'ceilf' 'floorf' 'fmodf' ...
                                'roundf'});      
end

%flops
info.flops.add = info.flops.add + size(strfind(data,'+'),2);
info.flops.mul = info.flops.mul + size(strfind(data,'*'),2);
info.flops.div = info.flops.div + size(strfind(data,'/'),2);

minus = regexp( data,'([^\E\(\s\\\*][\-])?', 'tokens'); % count minus operations avoiding counting sign
info.flops.add = info.flops.add + length(minus);
    
end

function res = replaceEsp(c)
ex = strsplit(char(c),',');
esponente = str2num(ex{2});
num = sprintf(['(' char(ex{1}) ')']);
res = num;
if ~isempty(regexp(num,'pow'))
    res = sprintf(['pow(' c ')']);
elseif ~rem(esponente,1)&& esponente < 9 && esponente > 0
    for j = 1:esponente-1
        res = sprintf([res '*' num]);
    end
else
    res = sprintf(['pow(' num ', %i)'],esponente);
end
end
