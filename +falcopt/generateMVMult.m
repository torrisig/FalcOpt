%% generateMVMult Generate data and code for a C-function for 
% matrix-vector multiplication and addition of the following form:
%  r = M1*v1 + ... + Ml*vl + m1 + ... + mk
% 
% [data, code, info] = generateMVMult(M)
% or
% [data, code, info] = generateMVMult(M, ...) [with options]
% or
% [data, code, info] = generateMVMult(M, m, ...) [with options]
% 
% Returns two strings [data, code] containing the static data and the code, respectively, and
%  returns an info struct containing information about the number of static
%  data points in 'size', the number of FLOPS in 'flops' and the used names in 'names' 
%
% the following options are available:
%  .names     - A struct with function and variable names to be used. Default is
%                struct('fun', 'mvmul', 'M', 'M', 'm', 'm', 'v', 'v', 'r', 'r', 'prefix', '')
%  .structure - A string or struct that determines the structure of M and m. Default is 'sparse'.
%                if a struct, then of the form struct('M', {{'dense', 'sparse'}}, 'm', {{'unique'}}).
%                only affects how the structure of non-static data is interpreted and utilized.
%                the possibilities are:
%                'dense' for no structure,
%                'sparse' for ommitted zero elements (only non-zero elements of non-static matrices are considered)
%                'unique' for ommitted zero and reccurring elements (recurring elements of non-static matrices are assumed to be the same)
%                'ordered' for ommitted zero and reccurring elements including permutation of elements (elements of non-static matrices are stored in ascending order of their values)
%                'indexed' for ommitted zero and recurring elements including permutation of elements (elements of non-static matrices are stored according to their values as indices starting from 1)
%  .add       - A boolean, if true r = r + M1*v1 + ... + Ml*vl + m1 + ... + mk is computed
%  .sub       - A boolean, if true r = r - M1*v1 - ... - Ml*vl + m1 + ... + mk is computed
%  .stripped  - A boolean, if true then code only contains the computation, no function and data is a struct with elements .M and .m which are cells containing only the static data.
%  .precision - A string that determines the precision of computation and stored data in the generated code. 
%                Needs to be either 'single' or 'double'. Default: '' (unspecified, will try to extract from type);
%  .types     - A string or struct. If a string then it indicates the data type,
%                can be 'float', 'double' or a custom type. Needs to be compatible with the floating-point type of .precision.
%                If a struct it also indicates the function type. Default is struct('fun', 'static inline void', 'data', 'double').
%  .symmetric - A boolean, if true then all matrices are considered symmetric and code is generate to exploit this.
%                Can also be a boolean vector. Default: false.
%                NOTE: Static matrices defined as symmetric are always symmetrized via (M+M')/2
%  .transpose - A boolean, if true then all matrices are transposed.
%                Can also be a boolean vector. Default: false.
%  .static    - A boolean, if true all matrices M are considered static and are defined in 'data'
%                If a boolean vector this is applied to the matrices M individually.
%                Can also be a struct, default: struct('M', true(1,l), 'm', true(1,k)).
%  .scale     - A static scalar that scales the result. Default is 1.
%                Can also be a struct struct('M', [..], 'm', [..]) to scale individually.
%  .indent    - Indentation to be used in code generation
%  .inline    - Inline keyword to be uised. Default is inline
%  .verbose   - Level of procedural output of this function
%  .test      - Level of tests performed
%  .epsFactor - Factor of machine-precision that is considered. Default is 100
%  .seed      - Seed used for random number generation in testing
% If fields of options structs are ommitted, the default is used.
% 

% Copyright (c) 2017, ETH Zurich, Automatic Control Laboratory 
%                    Damian Frick <falcopt@damianfrick.com>
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
function [data, code, info] = generateMVMult(varargin)
    functionId = 'generateMVMult';

    defaultNames = struct('fun', 'mvmul', ...
                          'M', 'M', ...
                          'm', 'm', ...
                          'v', 'v', ...
                          'r', 'r', ...
                          'prefix', '');
    structures = {'dense', 'sparse', 'unique', 'ordered', 'indexed'};
    indentTypes = {'generic', 'data', 'code'};
    defaultStructure = 'sparse';
    p = inputParser;
    p.CaseSensitive = true;
    p.addRequired('M', @(x)(isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x)))));
    p.addOptional('m', [], @(x)(isempty(x) || isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x)))));
    p.addParameter('names', defaultNames, @isstruct);
    p.addParameter('structure', defaultStructure, @(x)((ischar(x) && any(strcmp(x, structures))) || (isstruct(x) && any(isfield(x, {'M', 'm'})))));
    p.addParameter('add', false, @islogical);
    p.addParameter('sub', false, @islogical);
    p.addParameter('stripped', false, @islogical);
    p.addParameter('precision', '', @(s)(ischar(s) && any(strcmp(s, {'single', 'double'}))));
    p.addParameter('types', 'double', @(x)(ischar(x) || (isstruct(x) && any(isfield(x, {'fun', 'data'}) && (~isfield(x, 'fun') || ischar(x.fun)) && (~isfield(x, 'data') || ischar(x.data))))));
    p.addParameter('symmetric', false, @islogical);
    p.addParameter('transpose', false, @islogical);
    p.addParameter('static', true, @(x)(islogical(x) || (isstruct(x) && any(isfield(x, {'M', 'm'})))));
    p.addParameter('scale', true, @(x)((isnumeric(x) && numel(x)==1) || (isstruct(x) && any(isfield(x, {'M', 'm'})) && (~isfield(x, 'M') || isnumeric(x.M) ) && (~isfield(x, 'm') || isnumeric(x.m) ))));
    p.addParameter('indent', '    ', @(x)(ischar(x) || (isstruct(x) && isfield(x, 'generic') && all(cellfun(@(y)(~isfield(x, y) || ischar(x.(y))), indentTypes)))));
    p.addParameter('inline', 'inline', @ischar);
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && x >= 0));
    p.addParameter('test', 0, @(x)(isnumeric(x) && x >= 0));
    p.addParameter('epsFactor', 100, @(x)(isnumeric(x) && x > 0));
    p.addParameter('seed', 1, @isnumeric);
    p.parse(varargin{:});
    options = p.Results;
    M = options.M;
    m = options.m;
    names = options.names;
    
    if options.verbose >= 2
        fprintf('Code generation started...\n');
        fprintf('. processing options\n');
    end
    
    % Additional parameters
    options.lineSplitThreshold = 5;
    
    dims = struct();
    
    %% Bring into canonical form
    if ~iscell(M) % M
        M = {M};
    end
    dims.l = length(M);
    if ~iscell(m) % m
        if isempty(m)
            m = {};
        else
            m = {m};
        end
    end
    dims.k = length(m);
    
    %% Check dimension of M and m
    % first check transpose option
    if length(options.transpose) == 1
        options.transpose = logical(options.transpose.*true(1,dims.l));
    elseif length(options.transpose) ~= dims.l
        throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(options.transpose)) ' of transpose flags does not match the number ' num2str(dims.l) ' of matrices M.']));
    end
    % transpose all transposed matrices
    for i=1:dims.l
        if options.transpose(i)
            M{i} = M{i}';
        end
    end
    % Store dimensions
    dims.m = size(M{1},1);
    % Check dimensions
    for i=1:dims.l % M
        dims.n(i) = size(M{i},2);
        if dims.m ~= size(M{i},1)
            if ~options.transpose(i)
                throw(MException([functionId ':InvalidDimension'], ['The ' num2str(i) '-th input matrix M has ' num2str(size(M{i},1)) ' rows, which does not match the number ' num2str(dims.m) ' of required rows.']));
            else
                throw(MException([functionId ':InvalidDimension'], ['The ' num2str(i) '-th (transposed) input matrix M'' has ' num2str(size(M{i},1)) ' columns, which does not match the number ' num2str(dims.m) ' of required colums.']));
            end
        end
    end
    for i=1:length(m) % m
        if any([dims.m, 1] ~= size(m{i}))
            throw(MException([functionId ':InvalidDimension'], ['The ' num2str(i) '-th input vector m has dimension ' mat2str(size(m{i})) ', which does not match the required dimension ' mat2str([dims.m, 1]) '.']));
        end
    end
    
    %% Check names and bring into canonical form
    fields = fieldnames(defaultNames);
    for f=1:length(fields)
        % Set defaults
        if ~isfield(names, fields{f})
            names.(fields{f}) = defaultNames.(fields{f});
        end
        if ~any(strcmp(fields{f}, {'M', 'm', 'v'}))
            continue;
        end
        if ~iscell(names.(fields{f}))
            names.(fields{f}) = {names.(fields{f})};
        end
        if any(strcmp(fields{f}, {'M', 'v'}))
            if length(names.(fields{f})) == 1 && dims.l > 1
                if iscell(names.(fields{f}))
                	base = names.(fields{f}){1};
                else
                    base = names.(fields{f});
                end
                for i=1:dims.l
                    names.(fields{f}){i} = sprintf('%s%i', base, i);
                end
            end
            if length(names.(fields{f})) ~= dims.l
                throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(names.(fields{f}))) ' of names for ' fields{f} ' does not match the number ' num2str(dims.l) ' of required names. Cannot extrapolate, since there is more than one name supplied']));
            end
        elseif strcmp(fields{f}, 'm')
            if dims.k == 0
                names.(fields{f}) = {};
            end
            if length(names.(fields{f})) == 1 && dims.k > 1
                if iscell(names.(fields{f}))
                	base = names.(fields{f}){1};
                else
                    base = names.(fields{f});
                end
                for i=1:dims.k
                    names.(fields{f}){i} = sprintf('%s%i', base, i);
                end
            end
            if length(names.(fields{f})) ~= dims.k
                throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(names.(fields{f}))) ' of names for ' fields{f} ' does not match the number ' num2str(dims.k) ' of required names. Cannot extrapolate, since there is more than one name supplied']));
            end
        end
    end
    
    %% Check optional properties
    % Scale
    if ~isstruct(options.scale)
        options.scale = struct('M', repmat(options.scale, 1, dims.l), 'm', repmat(options.scale, 1, dims.k));
    end
    if ~isfield(options.scale, 'M')
        options.scale.M = ones(1, dims.l); % Set scaling to 1
    end
    if ~isfield(options.scale, 'm')
        options.scale.m = ones(1, dims.k); % Set scaling to 1
    end
    if length(options.scale.M) ~= dims.l
        throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(options.scale.M)) ' of scaling factors for the matrices M does not match the number ' num2str(dims.l) ' of matrices M.']));
    end
    if length(options.scale.m) ~= dims.k
        throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(options.scale.m)) ' of scaling factors for the vectors m does not match the number ' num2str(dims.k) ' of vectors m.']));
    end
    
    % Static vs. non-static data
    if ~isstruct(options.static)
        options.static = struct('M', options.static, 'm', true(1,dims.k));
    else
        if ~isfield(options.static, 'M')
            options.static.M = true(1,dims.l);
        end
        if ~isfield(options.static, 'm')
            options.static.m = true(1,dims.k);
        end
    end
    if length(options.static.M) == 1
        options.static.M = logical(options.static.M.*true(1,dims.l));
    elseif length(options.static.M) ~= dims.l
        throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(options.static.M)) ' of static flags for the matrices M does not match the number ' num2str(dims.l) ' of matrices M.']));
    end
    if length(options.static.m) == 1
        options.static.m = logical(options.static.m.*true(1,dims.k));
    elseif length(options.static.m) ~= dims.k
        throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(options.static.m)) ' of static flags of vectors m does not match the number ' num2str(dims.k) ' of vectors m.']));
    end
    
    % Ensure add and sub options are not enabled simultaneuously
    if options.add && options.sub
        throw(MException([functionId ':InvalidParameters'], '"Sub" and "add" cannot be enabled at the same time.'));
    end
    % Transform sub into add by changing the sign of the scaling
    if options.sub
        options.add = true;
        options.sub = false;
        options.scale.M = -options.scale.M;
        options.scale.m = -options.scale.m;
    end
    % Directly add scaling to all static matrices and vectors
    for i=1:dims.l
        if options.static.M(i)
            M{i} = options.scale.M(i)*M{i};
            options.scale.M(i) = 1;
        end
    end
    for i=1:dims.k
        if options.static.m(i)
            m{i} = options.scale.m(i)*m{i};
            options.scale.m(i) = 1;
        end
    end
    % Prefix (add prefix to every externally defined element)
    names.fun = [names.prefix names.fun];
    for i=1:dims.l
        if options.static.M(i)
            names.M{i} = [names.prefix names.M{i}];
        end
    end
    for i=1:dims.k
        if options.static.m(i)
            names.m{i} = [names.prefix names.m{i}];
        end
    end
    names.prefix = ''; % Remove prefix, since it has been incorporated
    info.names = names;
    
    % Structure
    if ~isstruct(options.structure)
        structure = options.structure;
        options.structure = struct('M', {cell(1,dims.l)}, 'm', {cell(1,dims.k)});
        options.structure.M(:) = {structure};
        options.structure.m(:) = {structure};
    end
    if ~isfield(options.structure, 'M')
        options.structure.M = cell(1,dims.l);
        options.structure.M(:) = {defaultStructure}; % Setting the default
    else
        if ~iscell(options.structure.M)
            if ~any(strcmp(options.structure.M, structures))
                throw(MException([functionId ':InvalidOption'], ['The structure flag ''' options.structure.M ''' for M does not exist, needs to be one of {' strjoin(structures, ', ') '}.']));
            end
            structure = options.structure.M;
            options.structure.M = cell(1,dims.l);
            options.structure.M(:) = {structure};
        else
            if length(options.structure.M) == 1
                if ~any(strcmp(options.structure.M, structures))
                    throw(MException([functionId ':InvalidOption'], ['The structure flag ''' options.structure.M ''' for M does not exist, needs to be one of {' strjoin(structures, ', ') '}.']));
                end
                structure = options.structure.M{1};
                options.structure.M = cell(1,dims.l);
                options.structure.M(:) = {structure};
            elseif length(options.structure.M) ~= dims.l
                throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(length(options.structure.M) )) ' of structure flags of matrices M does not match the number ' num2str(dims.l) ' of matrices M.']));
            elseif ~all(cellfun(@(s)(any(strcmp(s, structures))), options.structure.M))
                indices = find(cellfun(@(s)(any(strcmp(s, structures))), options.structure.M));
                throw(MException([functionId ':InvalidOption'], ['The structure flag ''' options.structure.M{indices(1)} ''' for ' names.M{indices(1)} '(#' num2str(indices(1)) ') does not exist, needs to be one of {' strjoin(structures, ', ') '}.']));
            end
        end
    end
    if ~isfield(options.structure, 'm')
        options.structure.m = cell(1,dims.k);
        options.structure.m(:) = {defaultStructure}; % Setting the default
        if ~iscell(options.structure.m)
            if ~any(strcmp(options.structure.m, structures))
                throw(MException([functionId ':InvalidOption'], ['The structure flag ''' options.structure.m ''' for m does not exist, needs to be one of {' strjoin(structures, ', ') '}.']));
            end
            structure = options.structure.m;
            options.structure.m = cell(1,dims.k);
            options.structure.m(:) = {structure};
        else
            if length(options.structure.m) == 1
                if ~any(strcmp(options.structure.m, structures))
                    throw(MException([functionId ':InvalidOption'], ['The structure flag ''' options.structure.m ''' for m does not exist, needs to be one of {' strjoin(structures, ', ') '}.']));
                end
                structure = options.structure.m{1};
                options.structure.m = cell(1,dims.k);
                options.structure.m(:) = {structure};
            elseif length(options.structure.m) ~= dims.k
                throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(length(options.structure.m) )) ' of structure flags of vectors m does not match the number ' num2str(dims.k) ' of vectors m.']));
            elseif ~all(cellfun(@(s)(any(strcmp(s, structures))), options.structure.m))
                indices = find(cellfun(@(s)(any(strcmp(s, structures))), options.structure.m));
                throw(MException([functionId ':InvalidOption'], ['The structure flag ''' options.structure.m{indices(1)} ''' for ' names.m{indices(1)} '(#' num2str(indices(1)) ') does not exist, needs to be one of {' strjoin(structures, ', ') '}.']));
            end
        end
    end
    
    % Ensure and exploit symmetry
    if any(options.symmetric)
        warning([functionId ':MissingImplemementation'], 'The special treatment of symmetric matrices has not yet been implemented. Ignoring...');
    end
    if length(options.symmetric) == 1
        options.symmetric = logical(options.symmetric.*true(1,dims.l));
    elseif length(options.symmetric) ~= dims.l
        throw(MException([functionId ':InvalidDimension'], ['The number ' num2str(length(options.symmetric)) ' of symmetry flags does not match the number ' num2str(dims.l) ' of matrices M.']));
    end
    for i=1:length(options.symmetric)
        if options.symmetric(i) && norm(M{i} - (M{i}+M{i}')/2, inf) > eps
            if options.static(i)
                warning([functionId ':InvalidData'], ['The matrix ' names.M{i} ' (#'  num2str(i) ') is not symmetric (max. error: ' num2str(norm(M{i} - (M{i}+M{i}')/2, inf)) '). Symmetrizing...']);
            else
                throw(MException([functionId ':InvalidData'], ['The matrix ' names.M{i} ' (#'  num2str(i) ') is not symmetric (max. error: ' num2str(norm(M{i} - (M{i}+M{i}')/2, inf)) '). Cannot be symmetrized, since not a static matrix.']));
            end
        end
        % Set transpose flag to zero if also symmetric
        if options.symmetric(i) && options.transpose(i)
            options.transpose(i) = false;
        end
        if options.symmetric(i) && options.static.M(i)
            M{i} = (M{i}+M{i}')/2;
        end
    end
    
    % Types
    if ischar(options.types)
        options.types = struct('data', options.types);
    end
    if ~isfield(options.types, 'data')
        options.types.data = 'double';
    end
    if ~isfield(options.types, 'fun')
        options.types.fun = 'static ';
        if ~isempty(options.inline)
            options.types.fun = [options.types.fun options.inline ' '];
        end
        options.types.fun = [options.types.fun 'void'];
    end
    
    % Precision (due to legacy usage)
    if isempty(options.precision)
        warning([functionId ':MissingPrecision'], 'Precision was not specified. Will try to determine based on supplied type. NOTE: this feature is depricated and will be removed in future versions.');
        switch options.types.data
            case 'float'
                options.precision = 'single';
            case 'single'
                options.precision = 'single';
            case 'double'
                options.precision = 'double';
            otherwise
                warning([functionId ':InvalidPrecision'], ['Precision was not specified, could not infer proper precision from given type "' options.types.data '". Setting to "double".']);
                options.precision = 'double';
        end
    end
    
    % Indentation
    if ~isstruct(options.indent)
        indent = options.indent;
        options.indent = struct();
        for i=1:length(indentTypes)
            options.indent.(indentTypes{i}) = indent;
        end
    end
    for i=1:length(indentTypes)
        if ~isfield(options.indent, indentTypes{i})
            options.indent.(indentTypes{i}) = options.indent.generic;
        end
        options.indent.(indentTypes{i}) = options.indent.(indentTypes{i})(:)'; % Make sure is row vector
    end
    
    if options.verbose == 1
        fprintf(['Generating code for ' names.fun '()\n']);
    end
    
    %% Generate data
    if options.verbose >= 2
        fprintf('. processing data\n');
    end
    info.size.data = 0;
    data = '';
    mexTransformCode = '';
    for i=1:dims.l
        if options.static.M(i)
            [d, ~, in] = falcopt.generateData(M{i}, 'names', struct('M', names.M{i}, 'fromDense', ['transform_' names.M{i}]), ...
                                                  'type', options.types.data, 'precision', options.precision, 'structure', 'unique', 'noones', true, 'symmetric', options.symmetric(i), ...
                                                  'stripped', options.stripped, 'static', true, 'const', true, ...
                                                  'indent', options.indent, 'verbose', max(0,options.verbose-1));
            elements.M.ones.sign{i} = in.ones.sign;
            elements.M.ones.rows{i} = in.ones.rows;
            elements.M.ones.cols{i} = in.ones.cols;
            elements.M.ones.num(i) = in.ones.num;
            elements.M.ones.mat{i} = in.ones.mat;
            if options.stripped
                data.M{i} = d;
            elseif in.structure.stored.num > 0
                if info.size.data > 0
                    data = [data, sprintf('\n')]; %#ok
                end
                data = [data, d]; %#ok
            end
            % Bookkeeping of number of static variables
            info.size.data = info.size.data + in.size;
        else
            if ~options.transpose(i)
                [~, c, in] = falcopt.generateData(M{i}, 'names', struct('M', names.M{i}, 'fromDense', ['transform_' names.M{i}]), ...
                                                      'type', options.types.data, 'precision', options.precision, 'structure', options.structure.M{i}, 'noones', false, 'symmetric', options.symmetric(i), ...
                                                      'indent', options.indent, 'verbose', max(0,options.verbose-1));
            else
                [~, c, in] = falcopt.generateData(M{i}', 'names', struct('M', names.M{i}, 'fromDense', ['transform_' names.M{i}]), 'transpose', true, ...
                                                       'type', options.types.data, 'precision', options.precision, 'structure', options.structure.M{i}, 'noones', false, 'symmetric', options.symmetric(i), ...
                                                       'indent', options.indent, 'verbose', max(0,options.verbose-1));
            end
            elements.M.ones.sign{i} = [];
            elements.M.ones.rows{i} = [];
            elements.M.ones.cols{i} = [];
            elements.M.ones.num(i) = 0;
            elements.M.ones.mat{i} = zeros(dims.m, dims.n(i));
            % Add code for MEX transformation
            mexTransformCode = [mexTransformCode c.fromDense sprintf(['\n'])]; %#ok
        end
        elements.M.stored.values{i} = in.structure.stored.values;
        elements.M.stored.rows{i} = in.structure.stored.rows;
        elements.M.stored.cols{i} = in.structure.stored.cols;
        elements.M.stored.num(i) = in.structure.stored.num;
        elements.M.stored.mat{i} = in.structure.stored.mat;
        elements.M.access.rows{i} = in.structure.access.rows;
        elements.M.access.cols{i} = in.structure.access.cols;
        elements.M.access.indices{i} = in.structure.access.indices;
        elements.M.access.num(i) = in.structure.access.num;
        if ~options.stripped
            elements.M.utils.rand{i} = in.utils.rand;
        end
    end
    for i=1:dims.k
        if options.static.m(i)
            [d, ~, in] = falcopt.generateData(m{i}, 'names', struct('M', names.m{i}, 'fromDense', ['transform_' names.m{i}]), ...
                                                  'type', options.types.data, 'precision', options.precision, 'structure', 'unique', ...
                                                  'stripped', options.stripped, 'indent', options.indent, 'verbose', max(0,options.verbose-1));
            % Add data
            if options.stripped
                data.m{i} = d;
            elseif in.structure.stored.num > 0
                if info.size.data > 0
                    data = [data, sprintf('\n')]; %#ok
                end
                data = [data, d]; %#ok
            end
            % Bookkeeping of number of static variables
            info.size.data = info.size.data + in.size;
        else
            [~, c, in] = falcopt.generateData(m{i}, 'names', struct('M', names.m{i}, 'fromDense', ['transform_' names.m{i}]), ...
                                                  'type', options.types.data, 'precision', options.precision, 'structure', options.structure.m{i}, ...
                                                  'indent', options.indent, 'verbose', max(0,options.verbose-1));
            % Add code for MEX transformation
            mexTransformCode = [mexTransformCode c.fromDense sprintf(['\n'])]; %#ok
        end
        elements.m.stored.values{i} = in.structure.stored.values;
        %elements.m.stored.rows{i} = in.structure.stored.rows;
        elements.m.stored.rows{i} = in.structure.stored.rows;
        elements.m.stored.num(i) = in.structure.stored.num;
        elements.m.stored.mat{i} = in.structure.stored.mat;
        %elements.m.structure.access.rows{i} = in.structure.access.rows;
        elements.m.access.rows{i} = in.structure.access.rows;
        elements.m.access.indices{i} = in.structure.access.indices;
        elements.m.access.num(i) = in.structure.access.num;
        if ~options.stripped
            elements.m.utils.rand{i} = in.utils.rand;
        end
    end
    if info.size.data > 0 && ~options.stripped
        data = [sprintf([options.indent.data '/* Static data for ' names.fun '() */' '\n']), data];
    end
    info.elements = elements;
    
    %% Generate code
    if options.verbose >= 2
        fprintf('. generating code\n');
    end
    info.flops.add = 0;
    info.flops.mul = 0;
    % Function header
    code = '';
    if ~options.stripped
        code = sprintf([options.indent.code options.types.fun ' ' names.fun '(' options.types.data '* ' names.r]);
        funIndent = length([options.types.fun ' ' names.fun '(']);
        if dims.l >= options.lineSplitThreshold
            code = [code, sprintf([', \n' options.indent.code repmat(' ', 1 , funIndent)])];
        else
            code = [code, sprintf(', ')];
        end
        code = [code, sprintf(['const ' options.types.data '* ' strjoin(names.v, [', const ' options.types.data '* '])])];
        if dims.l > 0 && any(~options.static.M & (elements.M.stored.num > 0)) % Whether there are non-static matrices that have to be included in the function interface
            nonStaticNonEmptyM = ~options.static.M & (elements.M.stored.num > 0);
            if dims.l >= options.lineSplitThreshold || sum(nonStaticNonEmptyM) >= options.lineSplitThreshold
                code = [code, sprintf([', \n' options.indent.code repmat(' ', 1, funIndent)])];
            else
                code = [code, sprintf(', ')];
            end
            code = [code, sprintf(['const ' options.types.data '* ' strjoin(names.M(nonStaticNonEmptyM), [', const ' options.types.data '* '])])];
        end
        if dims.k > 0 && any(~options.static.m & (elements.m.stored.num > 0))
            nonStaticNonEmptym = ~options.static.m & (elements.m.stored.num > 0);
            if sum(nonStaticNonEmptym) >= options.lineSplitThreshold || (dims.l > 0 && sum(nonStaticNonEmptyM) >= options.lineSplitThreshold)
                code = [code, sprintf([', \n' options.indent.code repmat(' ', 1, funIndent)])];
            else
                code = [code, sprintf(', ')];
            end
            code = [code, sprintf(['const ' options.types.data '* ' strjoin(names.m(nonStaticNonEmptym), [', const ' options.types.data '* '])])];
        end
        code = [code, sprintf([') {' '\n'])];
    end
    % Code
    % TODO deal with scale = 0!
    for i=1:dims.m
        % Iterate over elements of r[i]
        if (dims.l > 0 && any(cellfun(@(x)(any(x==i)), elements.M.access.rows))) || (dims.l > 0 && any(cellfun(@(x)(any(x==i)), elements.M.ones.rows))) || (dims.k > 0 && any(cellfun(@(x)(any(x==i)), elements.m.access.rows))) % Whether r[i] is affected by M's or m's
            if options.add
                code = [code, sprintf([options.indent.code options.indent.generic names.r '[' num2str(i-1) '] = ' names.r '[' num2str(i-1) ']'])]; %#ok
                info.flops.add = info.flops.add+1;
            else
                code = [code, sprintf([options.indent.code options.indent.generic names.r '[' num2str(i-1) '] = '])]; %#ok
            end
            bFirst = true;
            % M
            if dims.l > 0
                % Find all M's that contribute to r[i]
                J = find(cellfun(@(x)(any(x==i)), elements.M.access.rows) | cellfun(@(x)(any(x==i)), elements.M.ones.rows));
                for j=J
                    K = find(elements.M.access.rows{j} == i)'; % Find all elements in row i for the jth M
                    Kone = find(elements.M.ones.rows{j} == i)';
                    cols = unique([elements.M.access.cols{j}(K) elements.M.ones.cols{j}(Kone)], 'sorted');
                    for k=cols
                        if any(elements.M.access.cols{j}(K) == k) % If the factor is a regular factor (not -1, 1)
                            if ~bFirst
                                if options.scale.M(j) >= 0
                                    code = [code, sprintf(' + ')]; %#ok
                                else
                                    code = [code, sprintf(' - ')]; %#ok
                                end
                                info.flops.add = info.flops.add+1;
                            else
                                if options.add
                                    if options.scale.M(j) >= 0
                                        code = [code, sprintf(' + ')]; %#ok
                                    else
                                        code = [code, sprintf(' - ')]; %#ok
                                    end
                                elseif options.scale.M(j) < 0
                                    code = [code, sprintf('-')]; %#ok
                                end
                            end
                            if abs(options.scale.M(j)) == 1
                                code = [code, sprintf([names.M{j} '[' num2str(elements.M.access.indices{j}(elements.M.access.rows{j} == i & elements.M.access.cols{j} == k)-1) ']*' names.v{j} '[' num2str(k-1) ']'])]; %#ok
                            else
                                code = [code, sprintf([falcopt.internal.num2str(abs(options.scale.M(j)), options.precision) '*' names.M{j} '[' num2str(elements.M.access.indices{j}(elements.M.access.rows{j} == i & elements.M.access.cols{j} == k)-1) ']*' names.v{j} '[' num2str(k-1) ']'])]; %#ok
                                info.flops.mul = info.flops.mul+1;
                            end
                            info.flops.mul = info.flops.mul+1;
                        else
                            if bFirst % Add minus-sign if first element is a -1
                                if elements.M.ones.mat{j}(i,k) < 0
                                    code = [code, sprintf('-')]; %#ok
                                end
                            else % Add minus or plus, depending on whether the next operation is a multiplication with -1
                                if elements.M.ones.mat{j}(i,k) < 0
                                     code = [code, sprintf(' - ')]; %#ok
                                else
                                    code = [code, sprintf(' + ')]; %#ok
                                end
                                info.flops.add = info.flops.add+1;
                            end
                            code = [code, sprintf([names.v{j} '[%i]'], k-1)]; %#ok
                        end
                        bFirst = false;
                    end
                end
            end
            % m
            if dims.k > 0
                J = find(cellfun(@(x)(any(x==i)), elements.m.access.rows));
                for j=J
                    if bFirst
                        bFirst = false;
                        if options.scale.m(j) < 0
                            code = [code, sprintf('-')];
                        end
                    else
                        if options.scale.m(j) >= 0
                            code = [code, sprintf(' + ')]; %#ok
                        else
                            code = [code, sprintf(' - ')]; %#ok
                        end
                        info.flops.add = info.flops.add+1;
                    end
                    if abs(options.scale.M(j)) == 1
                        code = [code, sprintf([names.m{j} '[%i]'], elements.m.stored.mat{j}(i)-1)]; %#ok
                    else
                        code = [code, sprintf([falcopt.internal.num2str(abs(options.scale.m(j), 'precision', options.precision), options.types.data) '*' names.m{j} '[' num2str(elements.m.stored.mat{j}(i)-1) ']'])]; %#ok
                        info.flops.mul = info.flops.mul+1;
                    end
                end
            end
            code = [code, sprintf(';\n')]; %#ok
        elseif ~options.add
            code = [code, sprintf([options.indent.code options.indent.generic names.r '[%i] = 0.0'], i-1)]; %#ok
            code = [code, sprintf(';\n')]; %#ok
        end
    end
    if ~options.stripped
        code = [code, sprintf([options.indent.code '}'])];
    else
        code = code(1:end-1);
    end
    
    if options.test > 0 && ~options.stripped
        if options.verbose >= 2
            fprintf('. generating MEX file\n');
        end
        %% Generate MEX file
        filename = [tempname '.c'];
        f = fopen(filename, 'w+');
        fprintf(f, ['/* Temporary MEX file for testing ' names.fun '() */' '\n' ...
                    '\n' ...
                    '#include "mex.h"' '\n' '\n']);
        % Static data
        fprintf(f, ['\n' '/***************' '\n' ...
                         ' * Static Data *' '\n' ...
                         ' ***************/' '\n']);
        fprintf(f, data);
        % Code generated function
        fprintf(f, ['\n' '/***************************' '\n' ...
                         ' * Code-generated Function *' '\n' ...
                         ' ***************************/' '\n']);
        fprintf(f, [code '\n']);
        % Auxiliary functions
        fprintf(f, ['\n' '/***********************' '\n' ...
                         ' * Auxiliary Functions *' '\n' ...
                         ' ***********************/' '\n']);
        [~, c] = falcopt.generateData(ones(dims.m,1), 'structure', 'dense', ...
                                                    'names', struct('M', 'r', 'toDense', 'transform_R_toOutput', 'fromDense', 'copy_R'), ...
                                                    'indent', struct('generic', options.indent.generic, 'code', ''), ...
                                                    'verbose', 0);
        if options.add % Copy and transform r into appropriate data type
            fprintf(f, [c.fromDense '\n']);
        end
        if ~strcmp(options.types.data, 'double') % Transform r into appropriate data type for output
            fprintf(f, [c.toDense '\n']);
        end
        if ~strcmp(options.types.data, 'double') % Transform v's into appropriate data type
            for i=1:dims.l
                [~, c] = falcopt.generateData(ones(dims.n(i),1), 'structure', 'dense', ...
                                                               'names', struct('M', 'v', 'fromDense', ['transform_' names.v{i}]), ...
                                                               'indent', struct('generic', options.indent.generic, 'code', ''), ...
                                                               'verbose', 0);
                fprintf(f, [c.fromDense '\n']);
            end
        end
        if any(options.transpose & ~options.static.M)
            fprintf(f, [options.indent.code 'static ' options.inline ' void transpose(const double* M, unsigned int n, unsigned int m, double* Mt) {' '\n']);
            if ~isempty(options.inline)
                fprintf(f, [options.inline ' ']);
            end
            fprintf(f, ['void transpose(const double* M, unsigned int n, unsigned int m, double* Mt) {' '\n']);
            fprintf(f, [options.indent.code options.indent.generic 'unsigned int i,j;' '\n']);
            fprintf(f, [options.indent.code options.indent.generic 'for(i=0; i<m; i++) { /* iterate over columns */' '\n']);
            fprintf(f, [options.indent.code options.indent.generic options.indent.generic 'for(j=0; j<n; j++) { Mt[j*m+i] = M[i*n+j]; }' '\n']);
            fprintf(f, [options.indent.code options.indent.generic '}' '\n']);
            fprintf(f, [options.indent.code '}' '\n']);
        end
        % Transform code for M's and m's
        fprintf(f, [mexTransformCode '\n']);
        
        fprintf(f, ['\n' '/****************' '\n' ...
                         ' * MEX Function *' '\n' ...
                         ' ****************/' '\n']);
        fprintf(f, ['\n' 'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {' '\n']);
        fprintf(f, [options.indent.generic '/* Processing in- and outputs */' '\n']);
        if dims.k>0
            nInputs = dims.l+sum(~options.static.M & (elements.M.stored.num > 0))+sum(~options.static.m & (elements.m.stored.num > 0));
        else
            nInputs = dims.l+sum(~options.static.M & (elements.M.stored.num > 0));
        end
        if options.add
            nInputs = nInputs+1;
        end
        if nInputs == 1
            fprintf(f, [options.indent.generic 'if(nrhs != 1) { mexErrMsgIdAndTxt("' names.fun ':InvalidInput", "Requires 1 input."); }' '\n']);
        else
            fprintf(f, [options.indent.generic 'if(nrhs != %i) { mexErrMsgIdAndTxt("' names.fun ':InvalidInput", "Requires %i inputs."); }' '\n'], nInputs, nInputs);
        end
        i=1;
        % Return vector r
        fprintf(f, [options.indent.generic 'if(nlhs != 1) { mexErrMsgIdAndTxt("' names.fun ':InvalidInput", "Requires 1 output."); }' '\n']);
        fprintf(f, [options.indent.generic 'plhs[0] = mxCreateDoubleMatrix(%i,1,mxREAL); /* Output vector r */' '\n'], dims.m);
        if options.add
            fprintf(f, [options.indent.generic '/* Result ' names.r ' (initial value) */' '\n']);
            fprintf(f, [options.indent.generic 'if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != %i) { mexErrMsgIdAndTxt("' names.fun ':InvalidDimension", "Return vector ' names.r ' has invalid dimension."); }' '\n'], dims.m);
            fprintf(f, [options.indent.generic options.types.data ' ' names.r '[%i]; double* ' names.r '_in = mxGetPr(prhs[0]);' '\n'], dims.m);
            fprintf(f, [options.indent.generic 'copy_R(' names.r '_in, ' names.r ');' '\n']);
            i=i+1;
        else
            if strcmp(options.types.data, 'double')
                fprintf(f, [options.indent.generic 'double* ' names.r ' = mxGetPr(plhs[0]);' '\n']);
            else
                fprintf(f, [options.indent.generic options.types.data ' ' names.r '[%i];' '\n']);
            end
        end
        % Vectors v
        fprintf(f, [options.indent.generic '/* Vectors ' strjoin(names.v, ', ') ' */' '\n']);
        for j=1:dims.l
            fprintf(f, [options.indent.generic 'if(!mxIsDouble(prhs[%i]) || mxIsComplex(prhs[%i])) { mexErrMsgIdAndTxt("' names.fun ':InvalidType", "Vector ' names.v{j} ' has invalid type."); }' '\n'], i-1, i-1);
            fprintf(f, [options.indent.generic 'if(mxGetM(prhs[%i]) != %i || mxGetN(prhs[%i]) != 1) { mexErrMsgIdAndTxt("' names.fun ':InvalidDimension", "Vector ' names.v{j} ' has invalid dimension."); }' '\n'], i-1, dims.n(j), i-1);
            if strcmp(options.types.data, 'double')
                fprintf(f, [options.indent.generic 'double* ' names.v{j} ' = mxGetPr(prhs[%i]);' '\n'], i-1);
            else
                fprintf(f, [options.indent.generic options.types.data ' ' names.v{j} '[%i]; double* ' names.v{j} '_in = mxGetPr(prhs[%i]);' '\n'], dims.n(j), i-1);
                fprintf(f, [options.indent.generic 'transform_' names.v{i} '(' names.v{j} '_in, ' names.v{j} ');' '\n']);
            end
            i=i+1;
        end
        % Non-static matrices M
        if any(~options.static.M & (elements.M.stored.num > 0))
            fprintf(f, [options.indent.generic '/* Non-static matrices ' strjoin(names.M(~options.static.M), ', ') ' */' '\n']);
            if any(options.transpose & ~options.static.M)
                fprintf(f, [options.indent.generic 'double ttemp[%i]; /* Temporary array for transposing */' '\n'], dims.m*max(dims.n(options.transpose & ~options.static.M)));
            end
            for j=1:dims.l
                if ~options.static.M(j) && (elements.M.stored.num(j) > 0)
                    if ~options.transpose(j)
                        fprintf(f, [options.indent.generic 'if(!mxIsDouble(prhs[%i]) || mxIsComplex(prhs[%i]) || mxGetM(prhs[%i]) != %i || mxGetN(prhs[%i]) != %i) { mexErrMsgIdAndTxt("' names.fun ':InvalidDimension", "Matrix ' names.M{j} ' has invalid dimension."); }' '\n'], i-1, i-1, i-1, dims.m, i-1, dims.n(j));
                    else
                        fprintf(f, [options.indent.generic 'if(!mxIsDouble(prhs[%i]) || mxIsComplex(prhs[%i]) || mxGetM(prhs[%i]) != %i || mxGetN(prhs[%i]) != %i) { mexErrMsgIdAndTxt("' names.fun ':InvalidDimension", "Matrix ' names.M{j} ' has invalid dimension."); }' '\n'], i-1, i-1, i-1, dims.n(j), i-1, dims.m);
                    end
                    fprintf(f, [options.indent.generic options.types.data ' ' names.M{j} '[%i]; '], elements.M.stored.num(j));
                    fprintf(f, ['double* ' names.M{j} '_in = mxGetPr(prhs[%i]);' '\n'], i-1);
                    if options.transpose(j)
                        fprintf(f, [options.indent.generic 'transpose(' names.M{j} '_in, %i, %i, ttemp);' '\n'], dims.n(j), dims.m);
                        fprintf(f, [options.indent.generic 'transform_' names.M{j} '(ttemp, ' names.M{j} ');' '\n']);
                    else
                        fprintf(f, [options.indent.generic 'transform_' names.M{j} '(' names.M{j} '_in, ' names.M{j} ');' '\n']);
                    end
                    i = i+1;
                end
            end
        end
        % Non-static vectors m
        if dims.k > 0 && any(~options.static.m & (elements.m.stored.num > 0))
            fprintf(f, [options.indent.generic '/* Non-static vectors ' strjoin(names.m(~options.static.m), ', ') ' */' '\n']);
            for j=1:dims.k
                if ~options.static.m(j) && (elements.m.stored.num(j) > 0)
                    fprintf(f, [options.indent.generic 'if(!mxIsDouble(prhs[%i]) || mxIsComplex(prhs[%i]) || mxGetM(prhs[%i]) != %i || mxGetN(prhs[%i]) != 1) { mexErrMsgIdAndTxt("' names.fun ':InvalidDimension", "Vector ' names.m{j} ' has invalid dimension."); }' '\n'], i-1, i-1, i-1, dims.m, i-1);
                    if strcmp(options.structure.M{j}, 'dense')
                        fprintf(f, [options.indent.generic options.types.data ' ' names.m{j} '[%i]; '], dims.n(j));
                    else
                        fprintf(f, [options.indent.generic options.types.data ' ' names.m{j} '[%i]; '], elements.m.stored.num(j));
                    end
                    fprintf(f, ['double* ' names.m{j} '_in = mxGetPr(prhs[%i]);' '\n'], i-1);
                    fprintf(f, [options.indent.generic 'transform_' names.m{j} '(' names.m{j} '_in, ' names.m{j} ');' '\n']);
                    i = i+1;
                end
            end
        end
        % Call function
        fprintf(f, [options.indent.generic '/* Call function */' '\n']);
        fprintf(f, [options.indent.generic names.fun '(' names.r]);
        funIndent = length([names.fun '(']);
        if dims.l >= options.lineSplitThreshold
            fprintf(f, [', \n' options.indent.generic repmat(' ', 1 , funIndent)]);
        else
            fprintf(f, ', ');
        end
        fprintf(f, strjoin(names.v, ', '));
        if dims.l>0 && any(~options.static.M & (elements.M.stored.num > 0)) % Whether there are non-static matrices that have to be included in the function interface
            if dims.l >= options.lineSplitThreshold || sum(~options.static.M & (elements.M.stored.num > 0)) >= options.lineSplitThreshold
                fprintf(f, [', \n' options.indent.generic repmat(' ', 1, funIndent)]);
            else
                fprintf(f, ', ');
            end
            fprintf(f, strjoin(names.M(~options.static.M & (elements.M.stored.num > 0)), ', '));
        end
        if dims.k>0 && any(~options.static.m & (elements.m.stored.num > 0))
            if sum(~options.static.m & (elements.m.stored.num > 0)) >= options.lineSplitThreshold || (dims.l>0 && sum(~options.static.M & (elements.M.stored.num > 0)) >= options.lineSplitThreshold)
                fprintf(f, [', \n' options.indent.generic repmat(' ', 1, funIndent)]);
            else
                fprintf(f, ', ');
            end
            fprintf(f, strjoin(names.m(~options.static.m & (elements.m.stored.num > 0)), ', '));
        end
        fprintf(f, ');\n');
        % Outputs
        if ~strcmp(options.types.data, 'double')
            fprintf(f, [options.indent.generic '/* Prepare outputs */' '\n']);
            fprintf(f, [options.indent.generic 'transform_R_toOutput(' names.r ', mxGetPr(plhs[0]));' '\n']);
        end
        fprintf(f, ['}' '\n']);
        fclose(f);
        if options.verbose >= 3
            open(filename);
        end
        % Compile
        if options.verbose <= 1
            compile = 'mex -silent';
        elseif options.verbose <= 2
            compile = 'mex';
        else
            compile = 'mex -v';
        end
        if options.verbose >= 2
            fprintf('. compiling MEX file\n');
        end
        compile = [compile ' CFLAGS="$CFLAGS -Wall" ' filename ' -outdir ./ -output ' names.fun '_mex'];
        if options.verbose >= 3
            fprintf(compile);
        end
        eval(compile);
        % Test
        if options.verbose == 1
            fprintf('Testing code 0%%\n');
        elseif options.verbose >= 2
            fprintf('. testing code 0%%\n');
        end
        % Run tests
        errors = 0;
        maxError = 0;
        prevlen = length(sprintf('%i%%\n', 0));
        stream = RandStream('mt19937ar', 'Seed', options.seed);
        for i=1:floor(options.test)
            % Generate data
            vs = {};
            for j=1:dims.l % v's
                vs = [vs, {randn(dims.n(j),1)}]; %#ok
            end
            Ms = {};
            for j=1:dims.l
                if ~options.static.M(j)
                    if ~options.transpose(j)
                        Ms = [Ms, {elements.M.utils.rand{j}(stream, -2,2)}]; %#ok
                    else
                        Ms = [Ms, {elements.M.utils.rand{j}(stream, -2,2)'}]; %#ok
                    end
                end
            end
            ms = {};
            for j=1:dims.k % m's
                if ~options.static.m(j)
                    ms = [ms, {elements.m.utils.rand{i}(stream, -2,2)}]; %#ok
                end
            end
            % Compute r using code-generated function
            eval(['rs = ' names.fun '_mex(vs{:}, Ms{:}, ms{:});']);
            % Compute correct value
            r = zeros(dims.m,1);
            idx = 1;
            for j=1:dims.l
                if options.static.M(j)
                    r = r + M{j}*vs{j};
                else
                    if ~options.transpose(j)
                        r = r + Ms{idx}*vs{j};
                    else
                        r = r + Ms{idx}'*vs{j};
                    end
                    idx = idx+1;
                end
            end
            idx = 1;
            for j=1:dims.k
                if options.static.m(j)
                    r = r + m{j};
                else
                    r = r + ms{idx};
                    idx = idx+1;
                end
            end
            if isnan(norm(r-rs,inf)) || norm(r-rs,inf) > options.epsFactor*eps(options.types.data)
                errors = errors + 1;
            end
            if options.verbose >= 1
                if errors > 0
                    fprintf([repmat('\b', 1, prevlen) '%i%% (%i%% failed)\n'], floor(100*i/floor(options.test)), floor(100*errors/floor(options.test)));
                    prevlen = length(sprintf('%i%% (%i%% failed)\n', floor(100*i/floor(options.test)), floor(100*errors/floor(options.test))));
                else
                    fprintf([repmat('\b', 1, prevlen) '%i%%\n'], floor(100*i/floor(options.test)));
                    prevlen = length(sprintf('%i%%\n', floor(100*i/floor(options.test))));
                end
            end
        end
        if options.verbose >= 1 && errors == 0
            fprintf('\b success\n');
        end
        % Cleanup
        if options.verbose >= 2
            fprintf('. cleaning up\n');
        end
        delete([names.fun '_mex*']);
    end
    if options.verbose >= 2
        fprintf('...done\n');
    end
    if options.verbose >= 1
        if options.test == 0 || options.stripped || errors == 0
            fprintf('Code successfully generated\n');
        else
            fprintf('Code generated, testing revealed %i%% errors, the maximum error is %1.5g\n', floor(100*errors/floor(options.test)), maxError);
        end
    end
    if options.test > 0 && ~options.stripped && errors == options.test
        throw(MException([functionId ':FaultyImplementation'], 'All tests failed. Likely bug in implementation.'));
    end

end

