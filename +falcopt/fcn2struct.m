%% fcn2struct given a function handle f generate c-code equivalent
%
% Input:
% - sym (or numeric) matrix. (see matlab symbolic toolbox 'sym' variables)
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
function [data, info] = fcn2struct(varargin)
    p = inputParser;
    p.addRequired('function');
    p.addRequired('opt');
    p.addParameter('structure', 'unique');
    p.addParameter('name','x');
    p.parse(varargin{:});
    f = p.Results.function;
    o = p.Results.opt;

    data = [];
    subr = {};
    info.flops.add = 0;
    info.flops.mul = 0;
    info.flops.div= 0;

    x = sym('x', [o.nx,1],'real');
    u = sym('u', [o.nu,1],'real');
    w = sym('w', [o.nw,1],'real');
    sl = sym('sl', [50,1],'real'); % ToDo Tommaso change size of slack
    
    % check parameter .function
    if( isa(f,'function_handle'))
        error_s = sprintf(['Parameter "function" must be a sym-matrix (see matlab symbolic toolbox) or numeric matrix.\n'...
            'use (for example) x = sym(''x'',[5,1]); and call fcn2struct(f(x),opt);']);
        error(error_s);
    end
    
    % check if sym matrix or numeric matrix. If error assume sym matrix
    try
        static = isempty(symvar(f));
    catch
        static = 0;
    end


    if( ~static)
        % check existence of subexpression in f entries
        ind = 1;
        string = 's1';
        for i=1:10                  % limit to a max of 10 subexpression
            [f,subr{ind}] = subexpr(f,string);
            if isempty(subr{ind})
                subr(ind) = [];
                break;
            else
                subr{ind} = vpa(simplify(subr{ind}));
                ind = ind +1;
                string = sprintf('s%i',ind);
            end
            
        end
        
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
                %elseif( ~any(has(values,M(i,j))) )||(~isequal(p.Results.structure,'unique'))
                elseif( ~any(find(values==M(i,j))) )||(~isequal(p.Results.structure,'unique'))
                    values = [ values; M(i,j)];
                    ind = ind +1;
                    indeces = [indeces, ind];
                    S(i,j) = ind;
                else
                    %S(i,j) = indeces(find(has(values,M(i,j))== 1));
                    S(i,j) = indeces(find(values==M(i,j)));
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

    fun = @replaceEsp; % initialize function called by regexprep (to substitute 'pow')
    
    % substitute 'pow' in subexpression where possible and save  subexpressions to data
    for i = 1:length(subr)
        tmp_c = ccode(subr{i});
        tmp = strrep(tmp_c, 't0 = ','');
        
        % substitute the pow(.,n) c function with (.)*(.)*(.)... n-times
        tmp = regexprep(tmp,'pow\((.*?,\d\.\d)\)','${fun($1)}');
        
        % convert numbers into single/double format
        fun_2 = @(x)falcopt.internal.num2str(str2double(x),'precision',o.precision);
        tmp = regexprep(tmp,'(\d+\.\d+E[\-\+]*\d+|\d+\.\d+\d|\d\.\d)','${fun_2($1)}');
       
        % save into data
        data = [ data, sprintf(['\t' o.real ' s%i = %s\n'],i,tmp)];
    end
    
    % subs. 'pow' in values where possible and save to data
    for i= 1:size(values)
        
        tmp_c = ccode(values(i));
        tmp = strrep(tmp_c, 't0 = ','');
        
        % substitute the pow(.,n) c function with (.)*(.)*(.)... n-times
        tmp = regexprep(tmp,'pow\((.*?,\d\.\d)\)','${fun($1)}');
        
        % convert numbers into single/double format
        fun_2 = @(x)falcopt.internal.num2str(str2double(x),'precision',o.precision);
        tmp = regexprep(tmp,'(\d+\.\d+E[\-\+]*\d+|\d+\.\d+\d|\d\.\d)','${fun_2($1)}');
        
        
        if( static)
            tmp = regexprep( tmp,';','');
            data = [ data sprintf([tmp ', '])];
        else
            data = [data sprintf(['\t' p.Results.name '[%i] = ' tmp '\n'], i-1)];
        end 
    end
    
    % replace x1 with x[0] (the same for all x,u,w,sl vectors)
    for j = o.nw:-1:1
        data = strrep(data, char(w(j)), sprintf('w[%i]',j-1));
    end
    for j = o.nx:-1:1
        data = strrep(data, char(x(j)), sprintf('x[%i]',j-1));
    end
    for j = o.nu:-1:1
        data = strrep(data, char(u(j)), sprintf('u[%i]',j-1));
    end
    %ToDo Tommaso
    for j = 50:-1:1
        data = strrep(data, char(sl(j)), sprintf('sl[%i]',j-1));
    end
    
    % change c-functions (math.h) in case of single precision
    switch o.precision
        case 'single'
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
    
    % flops for principal functions
    info.flops.add = info.flops.add + size(regexp(data,'sqrt|sqrtf'),2)*4; %approx 4x an addition
    info.flops.add = info.flops.add + size(regexp(data,'sin|sinf|cos|cosf'),2)*15;  %approx 15x an addition
    info.flops.add = info.flops.add + size(regexp(data,'tan|tanf'),2)*20;  %approx 20x an addition
    info.flops.add = info.flops.add + size(regexp(data,'atan|atanf'),2)*23;  %approx 23x an addition
    info.flops.add = info.flops.add + size(regexp(data,'exp|expf'),2)*10;  %approx 10x an addition
    
    minus = regexp( data,'([^\E\(\s\\\*][\-])?', 'tokens'); % count minus operations avoiding counting sign
    info.flops.add = info.flops.add + length(minus);
    
end

function res = replaceEsp(c)
ex = strsplit(char(c),',');
esponente = str2num(ex{2});
num = sprintf(['(' char(ex{1}) ')']);
res = num;
if ~isempty(regexp(num,'pow')) % avoid errors in case of expression like pow(pow(...))
    res = sprintf(['pow(' c ')']);
elseif ~rem(esponente,1)&& esponente < 40 && esponente > 0  % only in case of integer exponent \in [0,40]
    for j = 1:esponente-1
        res = sprintf([res '*' num]);
    end
else
    res = sprintf(['pow(' num ', %i)'],esponente);
end
end
