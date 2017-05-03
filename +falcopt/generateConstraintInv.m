%% generateConstraintInv Generate code for a C-function for 
% solving the system of linear constraints induced by the jacobian:
%  r = (dg'*dg+diag(s)^2)^-1*v
% 
% [code, info] = generateConstraintInv(dp, dt)
% or
% [code, info] = generateConstraintInv(dp, dt, ...) [with options]
% 
% Returns a string 'code' containing the code, and
%  returns an info struct containing information about the number of static
%  data points in 'size', the number of FLOPS in 'flops' and the used names in 'names' 
%
% the following options are available:
%  .N         - The prediction horizon.
%  .names     - A struct with function and variable names to be used. Default is
%                struct('fun', 'solveConstraintSystem', 'dg', 'dg', 'dt', 'dt', 'prefix', '')
%  .bounds    - A struct with two fields .lb and .ub each containing a logical matrix.
%                The logical matrices must be of dimension (nu x N)
%                 with the logicals indicating that there is a lower/upper bound (if true)
%                 for the particular stage k in 1,...,N and component i in 1,...,nu of u.
%                Alternatively if the lower and upper bound are the same for all k in 1,...,N
%                 then only a logical vector of dimension nu can be supplied.
%                Default: struct('lb', true(nu,1), 'ub', true(nu,1))
%  .structure - A string that determines the structure of dp and dt. Default is 'sparse'.
%                only affects how the structure of dp, dt is interpreted and utilized.
%                the possibilities are:
%                'dense' for no structure,
%                'sparse' for ommitted zero elements (only non-zero elements are considered)
%                'unique' for ommitted zero and reccurring elements (recurring elements are assumed to be the same)
%                'ordered' for ommitted zero and reccurring elements including permutation of elements (elements are stored in ascending order of their values)
%                'indexed' for ommitted zero and recurring elements including permutation of element (elements are stored according to their values as indices starting from 1)
%               dp and dt can be treated differently if it is a struct e.g. struct('dp', 'dense', 'dt' 'sparse').
%  .precision - A string that determines the precision of computation and stored data in the generated code. 
%                Needs to be either 'single' or 'double'. Default: '' (unspecified, will try to extract from type);
%  .types     - A string or struct. If a string then it indicates the data type, either 'single' or 'double',
%                if a struct it also indicates the function type. Default is struct('fun', 'static inline void', 'data', 'double')
%  .indent    - Indentation to be used in code generation
%  .inline    - Inline keyword to be uised. Default is inline
%  .verbose   - Level of procedural output of this function
%  .test      - Level of tests performed
%  .epsFactor - Factor of machine-precision that is considered. Default is 100
% If fields of options structs are ommitted, the default is used.
% 

% Copyright (c) 2017 Damian Frick <falcopt@damianfrick.com>
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
function [code, info] = generateConstraintInv(varargin)
    defaultNames = struct('fun', 'solveConstraintSystem', ...
                          'r', 'r', ...
                          'dp', 'dp', ...
                          'dt', 'dt', ...
                          'v', 'v', ...
                          'slacks', 'sq', ...
                          'M', 'M', ...
                          'Mi', 'Mi', ...
                          'tmp', 'tmp', ...
                          'm1', 'm1', ...
                          'm2', 'm2', ...
                          'D1', 'D1', ...
                          'D2', 'D2', ...
                          'r1', 'r1', ...
                          'r2', 'r2', ...
                          'c', 'c', ...
                          'c1', 'c1', ...
                          'c2', 'c2', ...
                          'tau', 'tau', ...
                          'gamma', 'gamma', ...
                          'invert', 'invert', ...
                          'mtMult', 'mtMult', ...
                          'mtMultdt', 'mtMultdt', ...
                          'DMult', 'DMult', ...
                          'DMultdt', 'DMultdt', ...
                          'MiMult', 'MiMult', ...
                          'mMult', 'mMult', ...
                          'dptMultdt', 'dptMultdt', ...
                          'mMultAdd', 'mMultAdd', ...
                          'mMultSub', 'mMultSub', ...
                          'dttMultAdd', 'dttMultAdd', ...
                          'dttMultSub', 'dttMultSub', ...
                          'processSlacks', 'processSlacks', ...
                          'prefix', '');
    structures = {'dense', 'sparse', 'unique', 'ordered', 'indexed'};
    defaultStructure = struct('dp', 'sparse', 'dt', 'sparse');
    p = inputParser;
    p.CaseSensitive = true;
    p.addRequired('dp', @(x)(isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x)))));
    p.addRequired('dt', @(x)(isnumeric(x) || (iscell(x) && all(cellfun(@isnumeric, x)))));
    p.addParameter('N', 1, @isnumeric);
    p.addParameter('names', defaultNames, @isstruct);
    p.addParameter('bounds', struct(), @isstruct);
    p.addParameter('structure', defaultStructure, @(x)((ischar(x) && any(strcmp(x, structures))) || (isstruct(x) && all(isfield(x, {'dp', 'dt'})) && ischar(x.dp) && any(strcmp(x.dp, structures)) && ischar(x.dt) && any(strcmp(x.dt, strctures)))));
    p.addParameter('precision', '', @(s)(ischar(s) && any(strcmp(s, {'single', 'double'}))));
    p.addParameter('types', 'double', @(x)((ischar(x) && any(strcmp(x, {'float', 'single', 'double'}))) || (isstruct(x) && any(isfield(x, {'fun', 'data'}) && (~isfield(x, 'data') || any(strcmp(x.data, {'float', 'single', 'double'})))))));
    p.addParameter('indent', '    ', @(x)(ischar(x) || (isstruct(x) && all(isfield(x, {'generic', 'code'})) && ischar(x.generic) && ischar(x.data) && ischar(x.code))));
    p.addParameter('inline', 'inline', @ischar);
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && x >= 0));
    p.addParameter('test', 0, @(x)(isnumeric(x) && x >= 0));
    p.addParameter('epsFactor', 100, @(x)(isnumeric(x) && x > 0));
    p.parse(varargin{:});
    options = p.Results;
    dp = options.dp;
    dt = options.dt;
    names = options.names;
    
    if options.verbose >= 2
        fprintf('Code generation started...\n');
        fprintf('. processing options\n');
    end
    
    %% Bring dp and dt into canonical form
    if ~iscell(dp)
        dpAll = dp; clear dp;
        for k=1:options.N
            dp{k} = dpAll;
        end
    end
    if ~iscell(dt)
        if isempty(dt)
            clear dt;
            for k=1:options.N
                dt{k} = zeros(size(dp{k},1),0);
            end
        else
            dtAll = dt; clear dt;
            for k=1:options.N
                dt{k} = dtAll;
            end
        end
    end
    
    %% Extract and check dimensions
    if length(dp) ~= length(dt)
        throw(MException('falcopt:generateConstraintInv:InvalidDimension', 'Number of stages N and the number of cell-elements of dp and dt need to match.'));
    end
    options.N = length(dp);
    % Extract input and constraint dimension
    for k=1:options.N
        if ~isempty(dp{k})
            if ~isempty(dp{k})
                dims.nu(k) = size(dp{k},1); % Number of inputs of stage k
                dims.np(k) = size(dp{k},2); % Number of non-linear constraint (constraints of p) of stage k
            else
                dims.nu(k) = size(dt{k},1);
                dims.np(k) = 0;
            end
        else
            dims.nu(k) = size(dt{k},1);
            dims.np(k) = 0;
        end
        dims.nt(k) = size(dt{k},2);
        if isempty(dt{k})
            dt{k} = zeros(dims.nu(k),0);
        end
        if isempty(dp{k})
            dp{k} = zeros(dims.nu(k),0);
        end
        if dims.nu(k) ~= size(dt{k},1)
            throw(MException('falcopt:generateConstraintInv:InvalidDimension', ['Number of rows of dp and dt for k=' num2str(k) ' do not match.']));
        end
    end
    if any(dims.nt > 1)
        throw(MException('falcopt:generateConstraintInv:InvalidDimension', 'Jacobian dt can be at most a vector (at most one column).'));
    end
    % Only one u for all stages k
    if any(dims.nu ~= dims.nu(1))
        throw(MException('falcopt:generateConstraintInv:MissingImplementation', 'The dimension of the input needs to be the same for all stages.'));
    end
    
    %% Check and process bounds
    if ~isfield(options.bounds, 'lb')
        options.bounds.lb = true(dims.nu(1),options.N);
    elseif size(options.bounds.lb,2) == 1 && options.N > 1
        options.bounds.lb = repmat(options.bounds.lb, 1, options.N);
    end
    if ~isfield(options.bounds, 'ub')
        options.bounds.ub = true(dims.nu(1),options.N);
    elseif size(options.bounds.ub,2) == 1 && options.N > 1
        options.bounds.ub = repmat(options.bounds.ub, 1, options.N);
    end
    % Check
    if any(size(options.bounds.lb) ~= [dims.nu(1), options.N])
        throw(MException('falcopt:generateConstraintInv:InvalidDimension', 'Bound matrix .lb has invalid dimensions.'));
    elseif any(any(~islogical(options.bounds.lb)))
        throw(MException('falcopt:generateConstraintInv:InvalidInput', 'Bound matrix .lb needs to be a logical matrix.'));
    end
    if any(size(options.bounds.ub) ~= [dims.nu(1), options.N])
        throw(MException('falcopt:generateConstraintInv:InvalidDimension', 'Bound matrix .ub has invalid dimensions.'));
    elseif any(any(~islogical(options.bounds.ub)))
        throw(MException('falcopt:generateConstraintInv:InvalidInput', 'Bound matrix .ub needs to be a logical matrix.'));
    end
    % Determine number of bounds
    for k=1:options.N
        dims.lb(k) = sum(options.bounds.lb(:,k));
        dims.ub(k) = sum(options.bounds.ub(:,k));
        dims.bd(k) = sum(options.bounds.lb(:,k) | options.bounds.ub(:,k)); % Number of components of u that have at least one bound
    end
    
    %% Compute overall constraint dimension
    for k=1:options.N
        dims.nc(k) = dims.lb(k) + dims.ub(k) + dims.np(k);
    end
    
    %% Check names and bring into canonical form
    fields = fieldnames(defaultNames);
    for f=1:length(fields)
        % Set defaults
        if ~isfield(names, fields{f})
            names.(fields{f}) = defaultNames.(fields{f});
        end
    end
    info.names = names;
    
    %% Check optional properties
    % Structure
    if ~strcmp(options.structure.dp, 'sparse') || ~strcmp(options.structure.dt, 'sparse')
        throw(MException('falcopt:generateConstraintInv:InvalidParameter', 'For now generateConstraintInv can only handle ''sparse'' structures.'));
    end
    % Types
    if ischar(options.types)
        options.types = struct('data', options.types);
    end
    if ~isfield(options.types, 'data')
        options.types.data = 'double';
    end
    if ~isfield(options.types, 'fun')
        %options.types.fun = ['static ' options.inline ' void'];
        options.types.fun = 'static void';
    end
    if strcmp(options.types.data, 'single')
        warning('falcopt:generateConstraintInv:InvalidType', 'The type "single" is depricated and will not be allowed in a future version.');
        options.types.data = 'float';
    end
    % Precision (due to legacy usage)
    if isempty(options.precision)
        warning('falcopt:generateConstraintInv:MissingPrecision', 'Precision was not specified. Will try to determine based on supplied type. NOTE: this feature is depricated and will be removed in future versions.');
        switch options.types.data
            case 'float'
                options.precision = 'single';
            case 'single'
                options.precision = 'single';
            case 'double'
                options.precision = 'double';
            otherwise
                warning('falcopt:generateConstraintInv:InvalidPrecision', ['Precision was not specified, could not infer proper precision from given type "' options.types.data '". Setting to "double".']);
                options.precision = 'double';
        end
    end
    
    % Indentation
    if ~isstruct(options.indent)
    	options.indent = struct('generic', options.indent, 'data', options.indent, 'code', options.indent);
    end
    fields = fieldnames(options.indent);
    for i=1:length(fields)
        options.indent.(fields{i}) = options.indent.(fields{i})(:)'; % Make sure is row vector
    end
    
    % Prefix (add prefix to every externally defined element)
    names.fun = [names.prefix names.fun];
    names.prefix = ''; % Remove prefix, since it has been incorporated
    
    if options.verbose == 1
        fprintf(['Generating code for ' names.fun '()\n']);
    end
    
    %% Pre-processing, determining unique combinations of bounds, dp and dt
    [indices, structure] = extractUnique(dp, dt, dims, options);
    
    %% Generating code
    if options.verbose >= 2
        fprintf('. generating code\n');
    end
    info.flops.add = 0;
    info.flops.mul = 0;
    info.flops.inv = 0;
    code = '';
    
    %% Generating inverse/inverses
    inverseIndices = nan(1,options.N);
    idx = 1;
    while any(isnan(inverseIndices))
        k = find(isnan(inverseIndices), 1, 'first'); % Get index of first dp for which no inverse has been generated, yet
        I = cellfun(@(c)(size(dp{k},2) == size(c,2) && all(all(((abs(dp{k})'*abs(dp{k}) + eye(dims.np(k)))~=0) == ((abs(c)'*abs(c) + eye(dims.np(k)))~=0)))), dp) & isnan(inverseIndices);
        inverseIndices(I) = idx;
        if dims.np(k) > 0
            [invCode, invInfo] = falcopt.generateInverse(double(((abs(dp{k})'*abs(dp{k})+ eye(dims.np(k)))~=0)), 'names', struct('fun', [names.fun '_' names.invert num2str(idx)]), 'symmetric', true, ...
                                                       'indent', options.indent, 'inline', options.inline, 'types', options.types, 'precision', options.precision, ...
                                                       'verbose', max(0,options.verbose-1), 'test', options.test);
            structure.inv{idx} = invInfo;
            code = [code, invCode, sprintf('\n')]; %#ok
        else
            structure.inv{idx} = [];
        end
        idx = idx+1;
    end
    indices.inv.n = length(unique(inverseIndices));
    indices.inv.map = cell(1,indices.inv.n);
    for i=1:indices.inv.n
        indices.inv.map{i} = find(inverseIndices == i); % Indices of stages belonging to the i-th unique dp
    end
    indices.inv.k = @(i)(indices.inv.map{i}(1));
    indices.inv.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.inv.map)));
    clear inverseIndices;
    
    %% Compute number of elements of dp and dt for each stage
    for k=1:options.N
        % dt
        if ~isempty(structure.dt{indices.dt.i(k)})
            dims.mt(k) = structure.dt{indices.dt.i(k)}.structure.M.num;
        else
            dims.mt(k) = 0;
        end
        if (~isempty(structure.dt{indices.dt.i(k)}) && dims.nt(k) == 0) || (isempty(structure.dt{indices.dt.i(k)}) && dims.nt(k) > 0)
            throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
        end
        % dp
        if ~isempty(structure.dp{indices.dp.i(k)})
            dims.mp(k) = structure.dp{indices.dp.i(k)}.structure.M.num;
        else
            dims.mp(k) = 0;
        end
        if (~isempty(structure.dp{indices.dp.i(k)}) && dims.np(k) == 0) || (isempty(structure.dp{indices.dp.i(k)}) && dims.np(k) > 0)
            throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
        end
        % M, Mi
        if ~isempty(structure.inv{indices.inv.i(k)})
            dims.mM(k) = structure.inv{indices.inv.i(k)}.structure.M.num;
            dims.mMi(k) = structure.inv{indices.inv.i(k)}.structure.Mi.num;
        else
            dims.mM(k) = 0;
            dims.mMi(k) = 0;
        end
        % m1
        if ~isempty(structure.m1{indices.lbAndDp.i(k)})
            dims.mm1(k) = structure.m1{indices.lbAndDp.i(k)}.structure.M.num;
        else
            dims.mm1(k) = 0;
        end
        % m2
        if ~isempty(structure.m2{indices.ubAndDp.i(k)})
            dims.mm2(k) = structure.m2{indices.ubAndDp.i(k)}.structure.M.num;
        else
            dims.mm2(k) = 0;
        end
        % D1
        if ~isempty(structure.D1{indices.bounds.i(k)})
            dims.mD1(k) = structure.D1{indices.bounds.i(k)}.structure.M.num;
        else
            dims.mD1(k) = 0;
        end
        % D2
        if ~isempty(structure.D2{indices.bounds.i(k)})
            dims.mD2(k) = structure.D2{indices.bounds.i(k)}.structure.M.num;
        else
            dims.mD2(k) = 0;
        end
    end
    % Checking some assumptions
    nt = 0;
    np = 0;
    for k=1:options.N
        if nt ~= sum(dims.mt(1:k-1))
            throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
        elseif ~isempty(structure.dt{indices.dt.i(k)})
            nt = nt + structure.dt{indices.dt.i(k)}.structure.M.num;
        end
        if np ~= sum(dims.mp(1:k-1))
            throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
        elseif ~isempty(structure.dp{indices.dp.i(k)})
            np = np + structure.dp{indices.dp.i(k)}.structure.M.num;
        end
    end
    if max(dims.nt) > 0 && nt ~= sum(dims.mt)
        throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
    end
    if max(dims.np) > 0 && np ~= sum(dims.mp)
        throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
    end
    
    %% Generating matrix-vector multiplication functions
    % Generate mtMult for computing [m1' m2']*v
    % TODO generate documentation
    for i=1:indices.bdAndDp.n
        m1Structure = structure.m1{indices.lbAndDp.i(indices.bdAndDp.k(i))};
        m2Structure = structure.m2{indices.ubAndDp.i(indices.bdAndDp.k(i))};
        if ~isempty(m1Structure) && ~isempty(m2Structure)
            [~, multCode] = falcopt.generateMVMult({m1Structure.structure.M.mat, m2Structure.structure.M.mat}, 'structure', m1Structure.structure.M.type, ...
                                                                                                                     'names', struct('fun', [names.fun '_' names.mtMult num2str(i)]), ...
                                                                                                                     'types', options.types.data, 'precision', options.precision, ...
                                                                                                                     'transpose', true, 'static', false, 'indent', options.indent, ...
                                                                                                                     'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(m1Structure)
            [~, multCode] = falcopt.generateMVMult(m1Structure.structure.M.mat, 'structure', m1Structure.structure.M.type, ...
                                                                                           'names', struct('fun', [names.fun '_' names.mtMult num2str(i)]), ...
                                                                                           'types', options.types.data, 'precision', options.precision, ...
                                                                                           'transpose', true, 'static', false, 'indent', options.indent, ...
                                                                                           'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(m2Structure)
            [~, multCode] = falcopt.generateMVMult(m2Structure.structure.M.mat, 'structure', m2Structure.structure.M.type, ...
                                                                                           'names', struct('fun', [names.fun '_' names.mtMult num2str(i)]), ...
                                                                                           'types', options.types.data, 'precision', options.precision, ...
                                                                                           'transpose', true, 'static', false, 'indent', options.indent, ...
                                                                                           'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            % If matrices are all zeros
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    
    % Generate mMultAdd for computing r += [m1; m2]*x
    % TODO generate documentation
    for i=1:indices.bdAndDp.n
        m1Structure = structure.m1{indices.lbAndDp.i(indices.bdAndDp.k(i))};
        m2Structure = structure.m2{indices.ubAndDp.i(indices.bdAndDp.k(i))};
        if ~isempty(m1Structure) && ~isempty(m2Structure)
            [~, multCode] = falcopt.generateMVMult({[m1Structure.structure.M.mat; zeros(size(m2Structure.structure.M.mat))], [zeros(size(m1Structure.structure.M.mat)); m2Structure.structure.M.mat]}, ...
                                                 'structure', m1Structure.structure.M.type, 'add', true, ...
                                                 'names', struct('fun', [names.fun '_' names.mMultAdd num2str(i)]), ...
                                                 'types', options.types.data, 'precision', options.precision, ...
                                                 'static', false, 'indent', options.indent, ...
                                                 'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(m1Structure)
            [~, multCode] = falcopt.generateMVMult(m1Structure.structure.M.mat, 'structure', m1Structure.structure.M.type, 'add', true, ...
                                                                              'names', struct('fun', [names.fun '_' names.mMultAdd num2str(i)]), ...
                                                                              'types', options.types.data, 'precision', options.precision, ...
                                                                              'static', false, 'indent', options.indent, ...
                                                                              'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(m2Structure)
            [~, multCode] = falcopt.generateMVMult(m2Structure.structure.M.mat, 'structure', m2Structure.structure.M.type, 'add', true, ...
                                                                              'names', struct('fun', [names.fun '_' names.mMultAdd num2str(i)]), ...
                                                                              'types', options.types.data, 'precision', options.precision, ...
                                                                              'static', false, 'indent', options.indent, ...
                                                                              'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            % If matrices are all zeros
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    % Generate mMultSub for computing r -= [m1; m2]*x
    % TODO generate documentation
    for i=1:indices.bdAndDp.n
        m1Structure = structure.m1{indices.lbAndDp.i(indices.bdAndDp.k(i))};
        m2Structure = structure.m2{indices.ubAndDp.i(indices.bdAndDp.k(i))};
        if ~isempty(m1Structure) && ~isempty(m2Structure)
            [~, multCode] = falcopt.generateMVMult({[m1Structure.structure.M.mat; zeros(size(m2Structure.structure.M.mat))], [zeros(size(m1Structure.structure.M.mat)); m2Structure.structure.M.mat]}, ...
                                                 'structure', m1Structure.structure.M.type, 'sub', true, ...
                                                 'names', struct('fun', [names.fun '_' names.mMultSub num2str(i)]), ...
                                                 'types', options.types.data, 'precision', options.precision, ...
                                                 'static', false, 'indent', options.indent, ...
                                                 'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(m1Structure)
            [~, multCode] = falcopt.generateMVMult(m1Structure.structure.M.mat, 'structure', m1Structure.structure.M.type, 'sub', true, ...
                                                                              'names', struct('fun', [names.fun '_' names.mMultSub num2str(i)]), ...
                                                                              'types', options.types.data, 'precision', options.precision, ...
                                                                              'static', false, 'indent', options.indent, ...
                                                                              'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(m2Structure)
            [~, multCode] = falcopt.generateMVMult(m2Structure.structure.M.mat, 'structure', m2Structure.structure.M.type, 'sub', true, ...
                                                                              'names', struct('fun', [names.fun '_' names.mMultSub num2str(i)]), ...
                                                                              'types', options.types.data, 'precision', options.precision, ...
                                                                              'static', false, 'indent', options.indent, ...
                                                                              'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            % If matrices are all zeros
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    
    % Generate mtMultdt for computing [m1' m2']*[-dt(L) dt(U)]
    % TODO generate documentation
    for i=1:indices.bdAndDpAndDt.n
        if ~isempty(structure.bdAndDpAndDt.lb{i}) && ~isempty(structure.bdAndDpAndDt.ub{i})
            [~, multCode] = falcopt.generateMVMult({structure.bdAndDpAndDt.lb{i}.structure.M.mat, structure.bdAndDpAndDt.ub{i}.structure.M.mat}, 'structure', structure.bdAndDpAndDt.lb{i}.structure.M.type, ...
                                                                                                                     'scale', struct('M', [-1 1]), ...
                                                                                                                     'names', struct('fun', [names.fun '_' names.mtMultdt num2str(i)]), ...
                                                                                                                     'types', options.types.data, 'precision', options.precision, ...
                                                                                                                     'transpose', true, 'static', false, 'indent', options.indent, ...
                                                                                                                     'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDpAndDt.lb{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDpAndDt.lb{i}.structure.M.mat, 'structure', structure.bdAndDpAndDt.lb{i}.structure.M.type, ...
                                                                                           'scale', struct('M', -1), ...
                                                                                           'names', struct('fun', [names.fun '_' names.mtMultdt num2str(i)]), ...
                                                                                           'types', options.types.data, 'precision', options.precision, ...
                                                                                           'transpose', true, 'static', false, 'indent', options.indent, ...
                                                                                           'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDpAndDt.ub{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDpAndDt.ub{i}.structure.M.mat, 'structure', structure.bdAndDpAndDt.ub{i}.structure.M.type, ...
                                                                                           'names', struct('fun', [names.fun '_' names.mtMultdt num2str(i)]), ...
                                                                                           'types', options.types.data, 'precision', options.precision, ...
                                                                                           'transpose', true, 'static', false, 'indent', options.indent, ...
                                                                                           'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            % If matrices are all zeros
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    
    % Generate DMult for computing D*v
    % TODO generate documentation
    for i=1:indices.bounds.n
        if ~isempty(structure.D1{i}) && ~isempty(structure.D2{i})
            [~, multCode] = falcopt.generateMVMult({structure.D1{i}.structure.M.mat, structure.D2{i}.structure.M.mat}, 'structure', structure.D1{i}.structure.M.type, ...
                                                                                        'names', struct('fun', [names.fun '_' names.DMult num2str(i)]), ...
                                                                                        'types', options.types.data, 'precision', options.precision, ...
                                                                                        'static', false, 'indent', options.indent, ...
                                                                                        'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.D1{i})
            [~, multCode] = falcopt.generateMVMult(structure.D1{i}.structure.M.mat, 'structure', structure.D1{i}.structure.M.type, ...
                                                                                        'names', struct('fun', [names.fun '_' names.DMult num2str(i)]), ...
                                                                                        'types', options.types.data, 'precision', options.precision, ...
                                                                                        'static', false, 'indent', options.indent, ...
                                                                                        'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.D2{i})
            [~, multCode] = falcopt.generateMVMult(structure.D2{i}.structure.M.mat, 'structure', structure.D2{i}.structure.M.type, ...
                                                                                        'names', struct('fun', [names.fun '_' names.DMult num2str(i)]), ...
                                                                                        'types', options.types.data, 'precision', options.precision, ...
                                                                                        'static', false, 'indent', options.indent, ...
                                                                                        'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    
    % Generate DMultdt for computing D*[-dt(L); dt(U)]
    % TODO generate documentation
    for i=1:indices.bdAndDt.n
        if ~isempty(structure.bdAndDt.D1{i}) && ~isempty(structure.bdAndDt.D2{i})
            [~, multCode] = falcopt.generateMVMult({structure.bdAndDt.D1{i}.structure.M.mat, structure.bdAndDt.D2{i}.structure.M.mat}, 'structure', structure.bdAndDt.D1{i}.structure.M.type, ...
                                                                                        'names', struct('fun', [names.fun '_' names.DMultdt num2str(i)]), ...
                                                                                        'scale', struct('M', [-1 1]), 'types', options.types.data, 'precision', options.precision, ...
                                                                                        'static', false, 'indent', options.indent, ...
                                                                                        'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDt.D1{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDt.D1{i}.structure.M.mat, 'structure', structure.bdAndDt.D1{i}.structure.M.type, ...
                                                                                        'names', struct('fun', [names.fun '_' names.DMultdt num2str(i)]), ...
                                                                                        'scale', -1, 'types', options.types.data, 'precision', options.precision, ...
                                                                                        'static', false, 'indent', options.indent, ...
                                                                                        'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDt.D2{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDt.D2{i}.structure.M.mat, 'structure', structure.bdAndDt.D2{i}.structure.M.type, ...
                                                                                        'names', struct('fun', [names.fun '_' names.DMultdt num2str(i)]), ...
                                                                                        'types', options.types.data, 'precision', options.precision, ...
                                                                                        'static', false, 'indent', options.indent, ...
                                                                                        'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    
    % Generate MiMult for computing Mi*x
    % TODO generate documentation
    for i=1:indices.inv.n
        if ~isempty(structure.inv{i})
            [~, multCode] = falcopt.generateMVMult(structure.inv{i}.structure.Mi.mat, 'structure', structure.inv{i}.structure.Mi.type, ...
                                                                                    'names', struct('fun', [names.fun '_' names.MiMult num2str(i)]), ...
                                                                                    'types', options.types.data, 'precision', options.precision, ...
                                                                                    'transpose', true, 'static', false, 'indent', options.indent, ...
                                                                                    'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
            code = [code, multCode, sprintf('\n')]; %#ok
        end
    end
    
    % Generate dp'*dt
    % TODO generate documentation
    for i=1:indices.dpAndDt.n
        if ~isempty(structure.dpAndDt{i})
            [~, multCode] = falcopt.generateMVMult(structure.dpAndDt{i}.structure.M.mat, 'structure', structure.dpAndDt{i}.structure.M.type, ...
                                                                                       'names', struct('fun', [names.fun '_' names.dptMultdt num2str(i)]), ...
                                                                                       'types', options.types.data, 'precision', options.precision, ...
                                                                                       'static', false, 'indent', options.indent, ...
                                                                                       'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
            code = [code, multCode, sprintf('\n')]; %#ok
        end
    end
    
    % Generate dttMultAdd for r = r+[-dt(L)' dt(U)']
    % TODO generate documentation
    for i=1:indices.bdAndDt.n
        if ~isempty(structure.bdAndDt.dtL{i}) && ~isempty(structure.bdAndDt.dtU{i})
            [~, multCode] = falcopt.generateMVMult({structure.bdAndDt.dtL{i}.structure.M.mat, structure.bdAndDt.dtU{i}.structure.M.mat}, 'structure', structure.bdAndDt.dtL{i}.structure.M.type, ...
                                                                                          'names', struct('fun', [names.fun '_' names.dttMultAdd num2str(i)]), ...
                                                                                          'scale', struct('M', [-1 1]), 'add', true, 'types', options.types.data, 'precision', options.precision, ...
                                                                                          'static', false, 'indent', options.indent, ...
                                                                                          'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDt.dtL{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDt.dtL{i}.structure.M.mat, 'structure', structure.bdAndDt.dtL{i}.structure.M.type, ...
                                                                                          'names', struct('fun', [names.fun '_' names.dttMultAdd num2str(i)]), ...
                                                                                          'scale', struct('M', -1), 'add', true, 'types', options.types.data, 'precision', options.precision, ...
                                                                                          'static', false, 'indent', options.indent, ...
                                                                                          'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDt.dtU{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDt.dtU{i}.structure.M.mat, 'structure', structure.bdAndDt.dtU{i}.structure.M.type, ...
                                                                                          'names', struct('fun', [names.fun '_' names.dttMultAdd num2str(i)]), ...
                                                                                          'add', true, 'types', options.types.data, 'precision', options.precision, ...
                                                                                          'static', false, 'indent', options.indent, ...
                                                                                          'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    % Generate dttMultSub for r = r-[-dt(L)' dt(U)']
    % TODO generate documentation
    for i=1:indices.bdAndDt.n
        if ~isempty(structure.bdAndDt.dtL{i}) && ~isempty(structure.bdAndDt.dtU{i})
            [~, multCode] = falcopt.generateMVMult({structure.bdAndDt.dtL{i}.structure.M.mat, structure.bdAndDt.dtU{i}.structure.M.mat}, 'structure', structure.bdAndDt.dtL{i}.structure.M.type, ...
                                                                                          'names', struct('fun', [names.fun '_' names.dttMultSub num2str(i)]), ...
                                                                                          'scale', struct('M', [-1 1]), 'sub', true, 'types', options.types.data, 'precision', options.precision, ...
                                                                                          'static', false, 'indent', options.indent, ...
                                                                                          'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDt.dtL{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDt.dtL{i}.structure.M.mat, 'structure', structure.bdAndDt.dtL{i}.structure.M.type, ...
                                                                                          'names', struct('fun', [names.fun '_' names.dttMultSub num2str(i)]), ...
                                                                                          'scale', struct('M', -1), 'sub', true, 'types', options.types.data, 'precision', options.precision, ...
                                                                                          'static', false, 'indent', options.indent, ...
                                                                                          'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        elseif ~isempty(structure.bdAndDt.dtU{i})
            [~, multCode] = falcopt.generateMVMult(structure.bdAndDt.dtU{i}.structure.M.mat, 'structure', structure.bdAndDt.dtU{i}.structure.M.type, ...
                                                                                          'names', struct('fun', [names.fun '_' names.dttMultSub num2str(i)]), ...
                                                                                          'sub', true, 'types', options.types.data, 'precision', options.precision, ...
                                                                                          'static', false, 'indent', options.indent, ...
                                                                                          'inline', options.inline, 'verbose', max(0,options.verbose-1), 'test', options.test);
        else
            continue;
        end
        code = [code, multCode, sprintf('\n')]; %#ok
    end
    
    %% Generate code to construct various matrices from slacks
    [c, i] = generateSlacksCode(structure, indices, names, dims, options);
    code = [code, sprintf('\n') c sprintf('\n')];
    info.flops = falcopt.internal.addFlops(info.flops, i.flops);
    
    %% Generate code
    code = [code, sprintf([options.indent.code '\n'])];
    % TODO: Add documentation
    code = [code, sprintf([options.indent.code '/** ' '\n'])];
    code = [code, sprintf([options.indent.code ' */' '\n'])];
    % Function header
    code = [code, sprintf([options.indent.code options.types.fun ' ' names.fun '(' options.types.data '* ' names.r])];
    if max(dims.np) > 0
        code = [code, sprintf([', const ' options.types.data '* ' names.dp])];
    end
    if max(dims.nt) > 0
        code = [code, sprintf([', const ' options.types.data '* ' names.dt])];
    end
    code = [code, sprintf([', const ' options.types.data '* ' names.v ', const ' options.types.data '* ' names.slacks ') {' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
    code = [code, sprintf([options.indent.code options.indent.generic '/* Variables to store intermediate results */' '\n'])];
    % TODO make statically allocated?
    code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.D1 '[' num2str(max(dims.mD1)) '];' ' '])];
    code = [code, sprintf([options.types.data ' ' names.D2 '[' num2str(max(dims.mD2)) '];' '\n'])];
    if max(dims.np) > 0
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.r1 '[' num2str(max(dims.np)) '];' ' '])];
        code = [code, sprintf([options.types.data ' ' names.r2 '[' num2str(max(dims.np)) '];' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.M '[' num2str(max(dims.mM)) ']; /* Temporary variable for storing matrix */' '\n'])]; 
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.Mi '[' num2str(max(max(dims.mMi), max(dims.mp))) ']; /* Temporary variable for storing inverse of ' names.M ' */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.m1 '[' num2str(max(dims.mm1)) '];' ' '])];
        code = [code, sprintf([options.types.data ' ' names.m2 '[' num2str(max(dims.mm2)) '];' '\n'])];
    end
    code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.tmp '[' num2str(max(sum(options.bounds.lb | options.bounds.ub,1))) '];' '\n'])];
    if max(dims.nt) > 0
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.gamma ' = -' names.v '[' num2str(sum(dims.lb)+sum(dims.ub)+sum(dims.np)) ']; /* Initialize ' names.gamma ' with last element of ' names.v ' */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.tau ' = ' names.slacks '[' num2str(sum(dims.lb)+sum(dims.ub)+sum(dims.np)) ']; /* Initialize ' names.tau ' with the last element of ' names.slacks ' plus the squared 2-norm of ' names.dt ' (later) */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' ' names.c '[' num2str(sum(dims.lb)+sum(dims.ub)+sum(dims.np)) '];'])];
        if max(dims.np) > 0
            code = [code, sprintf([' ' options.types.data ' ' names.c1 '[%i];'], max(dims.np))];
            code = [code, sprintf([' ' options.types.data ' ' names.c2 '[%i];'], max(dims.np))];
        end
        code = [code, sprintf('\n')];
    end
    if max(dims.nt) > 0
        code = [code, sprintf([options.indent.code options.indent.generic '/* Update ' names.tau ' with squared 2-norm of ' names.dt ' */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<' num2str(sum(dims.mt)) '; i++) {' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic names.tau ' += ' names.dt '[i]*' names.dt '[i]; /* Adding squared 2-norm of ' names.dt ' */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic '}' '\n'])];
        info.flops.add = info.flops.add + sum(dims.mt);
        info.flops.mul = info.flops.mul + sum(dims.mt);
    end
    % TODO add option to not unroll in special cases
    % TODO finalize FLOPS count
    for k=1:options.N
        code = [code, sprintf(['\n' options.indent.code options.indent.generic '/** k = ' num2str(k-1) ' **/' '\n'])]; %#ok
        % Build auxiliary matrices
        code = [code, sprintf([options.indent.code options.indent.generic '/* Build auxiliary matrices */' '\n'])]; %#ok
        % Slacks
        code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.processSlacks '(&' names.slacks '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) ']'])]; %#ok
        % dp
        if max(dims.np) > 0
            code = [code, sprintf([', &' names.dp '[' num2str(sum(dims.mp(1:k-1))) ']'])]; %#ok
        end
        % stage index k
        if max([indices.bounds.n, indices.bdAndDp.n, indices.lbAndDp.n, indices.ubAndDp.n]) > 1
            code = [code, sprintf([', ' num2str(k-1)])]; %#ok
        end
        if max(dims.np) > 0
            code = [code, sprintf([', ' names.M])]; %#ok
            code = [code, sprintf([', ' names.Mi])]; %#ok
            code = [code, sprintf([', ' names.m1])]; %#ok
            code = [code, sprintf([', ' names.m2])]; %#ok
        end
        code = [code, sprintf([', ' names.D1])]; %#ok
        code = [code, sprintf([', ' names.D2])]; %#ok
        code = [code, sprintf([', ' names.tmp])]; %#ok
        if max(dims.np) > 0
            code = [code, sprintf([', ' names.Mi])]; %#ok Is not a bug, we use Mi also as temporary matrix
        end
        code = [code, sprintf([');' '\n'])]; %#ok
        
        % Compute initial value of rk
        if dims.nt(k) > 0
            code = [code, sprintf([options.indent.code options.indent.generic '/* Initialize vectors ' names.r '_k^{a,b} and ' names.c '_k^{a,b} */' '\n'])]; %#ok
        else
            code = [code, sprintf([options.indent.code options.indent.generic '/* Initialize vector ' names.r '_k^{a,b} */' '\n'])]; %#ok
        end
        if ~isempty(structure.D1{indices.bounds.i(k)}) && ~isempty(structure.D2{indices.bounds.i(k)})
            if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
            code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.DMult num2str(indices.bounds.i(k)) '(' ... 
                                   '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                   '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                   '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                   names.D1 ', ' names.D2 ');' '\n'])]; %#ok
        elseif ~isempty(structure.D1{indices.bounds.i(k)})
            if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
            code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.DMult num2str(indices.bounds.i(k)) '(' ... 
                                   '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                   '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                   names.D1 ');' '\n'])]; %#ok
        elseif ~isempty(structure.D2{indices.bounds.i(k)})
            if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
            code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.DMult num2str(indices.bounds.i(k)) '(' ... 
                                   '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ... 
                                   '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                   names.D2 ');' '\n'])]; %#ok
        end
        % Compute initial value of ck
        if dims.nt(k) > 0
            if ~isempty(structure.bdAndDt.D1{indices.bdAndDt.i(k)}) && ~isempty(structure.bdAndDt.D2{indices.bdAndDt.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.DMultdt num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                       names.D1 ', ' names.D2 ');' '\n'])]; %#ok
            elseif ~isempty(structure.bdAndDt.D1{indices.bdAndDt.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.DMultdt num2str(indices.bdAndDt.i(k)) '(' ... 
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                       names.D1 ');' '\n'])]; %#ok
            elseif ~isempty(structure.bdAndDt.D2{indices.bdAndDt.i(k)})
                if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.DMultdt num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                       names.D2 ');' '\n'])]; %#ok
            end
        end
        
        % Compute auxiliary vectors r1, r2 and c1, c2
        if dims.np(k) > 0
            if dims.nt(k) > 0
                code = [code, sprintf([options.indent.code options.indent.generic '/* Compute auxiliary vectors ' names.r1 ', ' names.r2 ' and ' names.c1 ', ' names.c2 ' */' '\n'])]; %#ok
            else
                code = [code, sprintf([options.indent.code options.indent.generic '/* Compute auxiliary vectors ' names.r1 ' and ' names.r2 ' */' '\n'])]; %#ok
            end
            % First compute r1
            if ~isempty(structure.m1{indices.lbAndDp.i(k)}) && ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mtMult num2str(indices.bdAndDp.i(k)) '(' ... 
                                       names.r1 ', ' ... 
                                       '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                       '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                       names.m1 ', ' names.m2 '); /* ' names.r1 ' */' '\n'])]; %#ok
            elseif ~isempty(structure.m1{indices.lbAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mtMult num2str(indices.bdAndDp.i(k)) '(' ...
                                       names.r1 ', ' ... 
                                       '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                       names.m1 '); /* ' names.r1 ' */' '\n'])]; %#ok
            elseif ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mtMult num2str(indices.bdAndDp.i(k)) '(' ...
                                       names.r1 ', ' ... 
                                       '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                       names.m2 '); /* ' names.r1 ' */' '\n'])]; %#ok
            end
            % and c1
            if dims.nt(k) > 0
                if ~isempty(structure.bdAndDpAndDt.lb{indices.bdAndDpAndDt.i(k)}) && ~isempty(structure.bdAndDpAndDt.ub{indices.bdAndDpAndDt.i(k)})
                    if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                        throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                    end
                    code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mtMultdt num2str(indices.bdAndDpAndDt.i(k)) '(' ...
                                           names.c1 ', ' ... 
                                           '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                           '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                           names.m1 ', ' names.m2 '); /* ' names.c1 ' */' '\n'])]; %#ok
                elseif ~isempty(structure.bdAndDpAndDt.lb{indices.bdAndDpAndDt.i(k)})
                    if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                        throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                    end
                    code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mtMultdt num2str(i) '(' ...
                                           names.c1 ', ' ... 
                                           '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                           names.m1 '); /* ' names.c1 ' */' '\n'])]; %#ok
                elseif ~isempty(structure.bdAndDpAndDt.ub{indices.bdAndDpAndDt.i(k)})
                    if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                        throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                    end
                    code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mtMultdt num2str(i) '(' ...
                                           names.c1 ', ' ... 
                                           '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                           names.m2 '); /* ' names.c1 ' */' '\n'])]; %#ok
                end
            end
            % Then compute r2
            code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.MiMult num2str(indices.inv.i(k)) '(' names.r2 ', ' names.r1 ', ' names.Mi '); /* ' names.r2 ' = ' names.Mi '*' names.r1 ' */ ' '\n'])]; %#ok
            % and c2
            if dims.nt(k) > 0
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.MiMult num2str(indices.inv.i(k)) '(' names.c2 ', ' names.c1 ', ' names.Mi '); /* ' names.c2 ' = ' names.Mi '*' names.c1 ' */ ' '\n'])]; %#ok
            end
        end
        
        % Finalize rk
        if dims.np(k) > 0
            code = [code, sprintf([options.indent.code options.indent.generic '/* Complete computation of ' names.r '_k^{a,b,c} */' '\n'])]; %#ok
            % Update rk^{a,b} using r2
            if ~isempty(structure.m1{indices.lbAndDp.i(k)}) && ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultAdd num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       names.r2 ', ' names.r2 ', ' ...
                                       names.m1 ', ' names.m2 '); /* Update ' names.r '_k^{a,b} using ' names.r2 ' */' '\n'])]; %#ok
            elseif ~isempty(structure.m1{indices.lbAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultAdd num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       names.r2 ', ' ...
                                       names.m1 ') /* Update ' names.r '_k^{a,b} using ' names.r2 ' */;' '\n'])]; %#ok
            elseif ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultAdd num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       names.r2 ', ' ...
                                       names.m2 ') /* Update ' names.r '_k^{a,b} using ' names.r2 ' */;' '\n'])]; %#ok
            end
            % Compute initial value of rk^c
            code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.MiMult num2str(indices.inv.i(k)) '(' ... 
                                   '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                   '&' names.v '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                   names.Mi '); /* r_k^c = ' names.Mi '*v_k^c */ ' '\n'])]; %#ok
            % Update rk^{a,b} using rk^c
            if ~isempty(structure.m1{indices.lbAndDp.i(k)}) && ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultSub num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       names.m1 ', ' names.m2 '); /* Update ' names.r '_k^{a,b} using ' names.r '_k^c */' '\n'])]; %#ok
            elseif ~isempty(structure.m1{indices.lbAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultSub num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       names.m1 ') /* Update ' names.r '_k^{a,b} using ' names.r '_k^c */;' '\n'])]; %#ok
            elseif ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultSub num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       names.m2 ') /* Update ' names.r '_k^{a,b} using ' names.r '_k^c */;' '\n'])]; %#ok
            end
            % Finalize rk^c
            code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<' num2str(dims.np(k)) '; i++) { ' ...
                                   names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '+i] -= ' names.r2 '[i]; }' ... 
                                   ' /* Finalize ' names.r '_k^c by subtracting ' names.r2 ' */' '\n'])]; %#ok
        end
        
        % Finalize ck
        if dims.np(k) > 0 && dims.nt(k) > 0
            code = [code, sprintf([options.indent.code options.indent.generic '/* Complete computation of ' names.c '_k^{a,b,c} */' '\n'])]; %#ok
            % Update ck^{a,b} using c2
            if ~isempty(structure.m1{indices.lbAndDp.i(k)}) && ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultAdd num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       names.c2 ', ' names.c2 ', ' ...
                                       names.m1 ', ' names.m2 '); /* Update ' names.c '_k^{a,b} using ' names.c2 ' */' '\n'])]; %#ok
            elseif ~isempty(structure.m1{indices.lbAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultAdd num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       names.c2 ', ' ...
                                       names.m1 ') /* Update ' names.c '_k^{a,b} using ' names.c2 ' */;' '\n'])]; %#ok
            elseif ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultAdd num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       names.c2 ', ' ...
                                       names.m2 ') /* Update ' names.c '_k^{a,b} using ' names.c2 ' */;' '\n'])]; %#ok
            end
            % Compute dp'*dt and store in c1
            code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.dptMultdt num2str(indices.dpAndDt.i(k)) '(' ...
                                   names.c1 ', ' ...
                                   '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                   '&' names.dp '[' num2str(sum(dims.mp(1:k-1))) ']); /* ' names.c1 ' = dp_k''*dt_k */' '\n'])]; %#ok
            % Compute initial value of ck^c using c1
            code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.MiMult num2str(indices.inv.i(k)) '(' ... 
                                   '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                   names.c1 ', ' ...
                                   names.Mi '); /* r_k^c = ' names.Mi '*' names.c1 ' */ ' '\n'])]; %#ok
            % Update ck^{a,b} using ck^c
            if ~isempty(structure.m1{indices.lbAndDp.i(k)}) && ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) == 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultSub num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       names.m1 ', ' names.m2 '); /* Update ' names.c '_k^{a,b} using ' names.c '_k^c */' '\n'])]; %#ok
            elseif ~isempty(structure.m1{indices.lbAndDp.i(k)})
                if dims.lb(k) == 0 || dims.ub(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultSub num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       names.m1 ') /* Update ' names.c '_k^{a,b} using ' names.c '_k^c */;' '\n'])]; %#ok
            elseif ~isempty(structure.m2{indices.ubAndDp.i(k)})
                if dims.ub(k) == 0 || dims.lb(k) ~= 0 % Sanity check
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.mMultSub num2str(indices.bdAndDp.i(k)) '(' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ... 
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '], ' ...
                                       names.m2 ') /* Update ' names.c '_k^{a,b} using ' names.c '_k^c */;' '\n'])]; %#ok
            end
            % Finalize ck^c
            code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<' num2str(dims.np(k)) '; i++) { ' ...
                                   names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '+i] -= ' names.c2 '[i]; }' ... 
                                   ' /* Finalize ' names.c '_k^c by subtracting ' names.c2 ' */' '\n'])]; %#ok
        end
        
        % Update tau and gamma
        if dims.nt(k) > 0
            code = [code, sprintf([options.indent.code options.indent.generic '/* Update ' names.gamma ' and ' names.tau ' */' '\n'])]; %#ok
            if ~isempty(structure.bdAndDt.dtL{indices.bdAndDt.i(k)}) && ~isempty(structure.bdAndDt.dtU{indices.bdAndDt.i(k)})
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.dttMultSub num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.gamma ', ' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) ']); /* First part of update of ' names.gamma ' */' '\n'])]; %#ok
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.dttMultAdd num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.tau ', ' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) ']); /* First part of update of ' names.tau ' */' '\n'])]; %#ok
            elseif ~isempty(structure.bdAndDt.dtL{indices.bdAndDt.i(k)})
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.dttMultSub num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.gamma ', ' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) ']); /* First part of update of ' names.gamma ' */' '\n'])]; %#ok
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.dttMultAdd num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.tau ', ' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) ']); /* First part of update of ' names.tau ' */' '\n'])]; %#ok
                
            elseif ~isempty(structure.bdAndDt.dtU{indices.bdAndDt.i(k)})
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.dttMultSub num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.gamma ', ' ...
                                       '&' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) ']); /* First part of update of ' names.gamma ' */' '\n'])]; %#ok
                code = [code, sprintf([options.indent.code options.indent.generic names.fun '_' names.dttMultAdd num2str(indices.bdAndDt.i(k)) '(' ...
                                       '&' names.tau ', ' ...
                                       '&' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)) '], ' ...
                                       '&' names.dt '[' num2str(sum(dims.mt(1:k-1))) ']); /* First part of update of ' names.tau ' */' '\n'])]; %#ok
            end
            if dims.np(k) > 0
                code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<' num2str(dims.np(k)) '; i++) {' '\n'])]; %#ok
                code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic names.gamma ' += ' names.c1 '[i]*' names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '+i]; /* Finalize update of ' names.gamma ' */' '\n'])]; %#ok
                code = [code, sprintf([options.indent.code options.indent.generic options.indent.generic names.tau ' += ' names.c1 '[i]*' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+' num2str(dims.lb(k)+dims.ub(k)) '+i]; /* Finalize update of ' names.tau ' */' '\n'])]; %#ok
                code = [code, sprintf([options.indent.code options.indent.generic '}' '\n'])]; %#ok
            end
        end
    end
    
    % Use ck to finalize rk
    if max(dims.nt) > 0
        code = [code, sprintf(['\n' options.indent.code options.indent.generic '/* Finalize ' names.r '_k using ' names.c '_k */' '\n'])];
        code = [code, sprintf([options.indent.code options.indent.generic names.tau ' = ' names.gamma '/' names.tau '; /* Only ' names.gamma ' divided by ' names.tau ' is needed */' '\n'])];
        for k=1:options.N
            if dims.nt(k) > 0
                code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<' num2str(dims.lb(k)+dims.ub(k)+dims.np(k)) '; i++) { ' ...
                                       names.r '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+i] += ' names.tau '*' names.c '[' num2str(sum(dims.lb(1:k-1)+dims.ub(1:k-1)+dims.np(1:k-1))) '+i];' ...                                     
                                       ' } /* k = ' num2str(k-1) ' */' '\n'])]; %#ok
            end
        end
        code = [code, sprintf([options.indent.code options.indent.generic names.r '[' num2str(sum(dims.lb+dims.ub+dims.np)) '] = -' names.tau '; /* Compute ' names.r '_N */' '\n'])];
    end
    
    % Close function
    code = [code, sprintf([options.indent.code '}' '\n'])];
end

% Function to generate code that processes slacks
%  builds matrices M, m1, m2, C and D 
function [code, info] = generateSlacksCode(structure, indices, names, dims, options)
    localNames.T = 'T';
    localNames.tmp = 'tmp';
    localNames.s = 's';

    % Initialize flop counts
    info.flops.mul = 0;
    info.flops.add = 0;
    info.flops.div = 0;

    %% Function interface
    docu = sprintf([options.indent.code '/** ' '\n' ...
                    options.indent.code ' * @brief Processe slacks and generate matrices ']);
    if max(dims.np) > 0
        docu = [docu sprintf([names.M ', ' names.m1 ', ' names.m2 ' as well as '])];
    end
    docu = [docu sprintf([names.D1 ' and ' names.D2 '.' '\n'])];
    code = sprintf([options.indent.code 'static ' options.inline ' void ' names.fun '_' names.processSlacks '(']);
    bFirst = true;
    % Slacks
	docu = [docu sprintf([options.indent.code ' * @param ' localNames.s ' The (squared) slacks for stage k.'])];
    if ~all((dims.lb+dims.ub+dims.np) == (dims.lb(1)+dims.ub(1)+dims.np(1)))
        docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' localNames.s ' '])) ' A vector of length '])];
        for k=1:options.N-1
            docu = [docu sprintf([num2str(dims.lb(k)+dims.ub(k)+dims.np(k)) '(k=' num2str(k-1) '), '])]; %#ok
        end
        docu = [docu sprintf(['or ' num2str(dims.lb(end)+dims.ub(end)+dims.np(end)) '(k=' num2str(options.N-1) ').' '\n'])];
    else
        docu = [docu sprintf([' A vector of dimension ' num2str(dims.lb(1)+dims.ub(1)+dims.np(1)) '.' '\n'])];
    end
    if ~bFirst
        code = [code sprintf(', ')];
    else
        bFirst = false;
    end
    code = [code sprintf(['const ' options.types.data '* ' localNames.s])];
    % Jacobian dp
    if max(dims.np) > 0
        docu = [docu sprintf([options.indent.code ' * @param ' names.dp ' Jacobian of the constraints p_k.'])];
        if ~all(cellfun(@(c)(isfield(c,'elements')),structure.dp)) || ~all(cellfun(@(c)(c.elements.M.num),structure.dp) == structure.dp{1}.elements.M.num)
            docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' names.dp ' '])) ' A matrix with '])];
            for k=1:options.N-1
                if ~isempty(structure.dp{indices.dp.i(k)})
                    docu = [docu sprintf([num2str(structure.dp{indices.dp.i(k)}.elements.M.num) '(k=' num2str(k-1) '), '])]; %#ok
                else
                    docu = [docu sprintf(['0(k=' num2str(k-1) '), '])]; %#ok
                end
            end
            if ~isempty(structure.dp{indices.dp.i(options.N)})
                docu = [docu sprintf(['or ' num2str(structure.dp{indices.dp.i(options.N)}.elements.M.num) '(k=' num2str(options.N-1) ') elements.' '\n'])];
            else
                docu = [docu sprintf(['or 0(k=' num2str(options.N-1) ') elements.' '\n'])];
            end
        else
            docu = [docu sprintf([' A matrix with ' num2str(structure.dp{1}.elements.M.num) ' elements.' '\n'])];
        end
        if ~bFirst
            code = [code sprintf(', ')];
        else
            bFirst = false;
        end
        code = [code sprintf(['const ' options.types.data '* ' names.dp])];
    end
    % Stage index k (if needed)
    if max([indices.bounds.n, indices.bdAndDp.n, indices.lbAndDp.n, indices.ubAndDp.n]) > 1
        docu = [docu sprintf([options.indent.code ' * @param k The stage for which the matrices are computed. Needs to be an unsigned integer value from 0 to ' num2str(options.N-1) '.' '\n'])];
        if ~bFirst
        code = [code sprintf(', ')];
        else
            bFirst = false;
        end
        code = [code sprintf('unsigned int k')];
    end
    % M
    if max(dims.np) > 0
        docu = [docu sprintf([options.indent.code ' * @param ' names.M ' Return value.'])];
        if ~all(cellfun(@(c)(isfield(c,'elements')),structure.inv)) || ~all(cellfun(@(c)(c.elements.M.num),structure.inv) == structure.inv{1}.elements.M.num)
            docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' names.M ' '])) ' A matrix with '])];
            for k=1:options.N-1
                if ~isempty(structure.inv{indices.inv.i(k)})
                    docu = [docu sprintf([num2str(structure.inv{indices.inv.i(k)}.elements.M.num) '(k=' num2str(k-1) '), '])]; %#ok
                else
                    docu = [docu sprintf(['0(k=' num2str(k-1) '), '])]; %#ok
                end
            end
            if ~isempty(structure.inv{indices.inv.i(options.N)})
                docu = [docu sprintf(['or ' num2str(structure.inv{indices.inv.i(options.N)}.elements.M.num) '(k=' num2str(options.N-1) ') elements.' '\n'])];
            else
                docu = [docu sprintf(['or 0(k=' num2str(options.N-1) ') elements.' '\n'])];
            end
        else
            docu = [docu sprintf([' A matrix with ' num2str(structure.inv{1}.elements.M.num) ' elements.' '\n'])];
        end
        if ~bFirst
            code = [code sprintf(', ')];
        else
            bFirst = false;
        end
        code = [code sprintf([options.types.data '* ' names.M])];
    end
    % Mi
    if max(dims.np) > 0
        docu = [docu sprintf([options.indent.code ' * @param ' names.Mi ' Return value. The inverse of ' names.M '.'])];
        if ~all(cellfun(@(c)(isfield(c,'elements')),structure.inv)) || ~all(cellfun(@(c)(c.elements.Mi.num),structure.inv) == structure.inv{1}.elements.Mi.num)
            docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' names.Mi ' '])) ' A matrix with '])];
            for k=1:options.N-1
                if ~isempty(structure.inv{indices.inv.i(k)})
                    docu = [docu sprintf([num2str(structure.inv{indices.inv.i(k)}.elements.Mi.num) '(k=' num2str(k-1) '), '])]; %#ok
                else
                    docu = [docu sprintf(['0(k=' num2str(k-1) '), '])]; %#ok
                end
            end
            if ~isempty(structure.inv{indices.inv.i(options.N)})
                docu = [docu sprintf(['or ' num2str(structure.inv{indices.inv.i(options.N)}.elements.Mi.num) '(k=' num2str(options.N-1) ') elements.' '\n'])];
            else
                docu = [docu sprintf(['or 0(k=' num2str(options.N-1) ') elements.' '\n'])];
            end
        else
            docu = [docu sprintf([' A matrix with ' num2str(structure.inv{1}.elements.Mi.num) ' elements.' '\n'])];
        end
        if ~bFirst
            code = [code sprintf(', ')];
        else
            bFirst = false;
        end
        code = [code sprintf([options.types.data '* ' names.Mi])];
    end
    % m1
    if max(dims.np) > 0
        docu = [docu sprintf([options.indent.code ' * @param ' names.m1 ' Return value.'])];
        if ~all(cellfun(@(c)(isfield(c,'elements')),structure.m1)) || ~all(cellfun(@(c)(c.elements.M.num),structure.m1) == structure.m1{1}.elements.M.num)
            docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' names.m1 ' '])) ' A matrix with '])];
            for k=1:options.N-1
                if ~isempty(structure.m1{indices.lbAndDp.i(k)})
                    docu = [docu sprintf([num2str(structure.m1{indices.lbAndDp.i(k)}.elements.M.num) '(k=' num2str(k) '), '])]; %#ok
                else
                    docu = [docu sprintf(['0(k=' num2str(k) '), '])]; %#ok
                end
            end
            if ~isempty(structure.m1{indices.lbAndDp.i(options.N)})
                docu = [docu sprintf(['or ' num2str(structure.m1{indices.lbAndDp.i(options.N)}.elements.M.num) '(k=' num2str(options.N) ') elements.' '\n'])];
            else
                docu = [docu sprintf(['or 0(k=' num2str(options.N) ') elements.' '\n'])];
            end
        else
            docu = [docu sprintf([' A matrix with ' num2str(structure.m1{1}.elements.M.num) ' elements.' '\n'])];
        end
        if ~bFirst
            code = [code sprintf(', ')];
        else
            bFirst = false;
        end
        code = [code sprintf([options.types.data '* ' names.m1])];
    end
    % m2
    if max(dims.np) > 0
        docu = [docu sprintf([options.indent.code ' * @param ' names.m2 ' Return value.'])];
        if ~all(cellfun(@(c)(isfield(c,'elements')),structure.m2)) || ~all(cellfun(@(c)(c.elements.M.num),structure.m2) == structure.m2{1}.elements.M.num)
            docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' names.m2 ' '])) ' A matrix with '])];
            for k=1:options.N-1
                if ~isempty(structure.m2{indices.ubAndDp.i(k)})
                    docu = [docu sprintf([num2str(structure.m2{indices.ubAndDp.i(k)}.elements.M.num) '(k=' num2str(k-1) '), '])]; %#ok
                else
                    docu = [docu sprintf(['0(k=' num2str(k-1) '), '])]; %#ok
                end
            end
            if ~isempty(structure.m2{indices.ubAndDp.i(options.N)})
                docu = [docu sprintf(['or ' num2str(structure.m2{indices.ubAndDp.i(options.N)}.elements.M.num) '(k=' num2str(options.N-1) ') elements.' '\n'])];
            else
                docu = [docu sprintf(['or 0(k=' num2str(options.N-1) ') elements.' '\n'])];
            end
        else
            docu = [docu sprintf([' A matrix with ' num2str(structure.m2{1}.elements.M.num) ' elements.' '\n'])];
        end
        if ~bFirst
            code = [code sprintf(', ')];
        else
            bFirst = false;
        end
        code = [code sprintf([options.types.data '* ' names.m2])];
    end
    % D1
    docu = [docu sprintf([options.indent.code ' * @param ' names.D1 ' Return value.'])];
    if ~all(cellfun(@(c)(isfield(c,'elements')),structure.D1)) || ~all(cellfun(@(c)(c.elements.M.num),structure.D1) == structure.D1{1}.elements.M.num)
        docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' names.D1 ' '])) ' A matrix with '])];
        for k=1:options.N-1
            if ~isempty(structure.D1{indices.bounds.i(k)})
                docu = [docu sprintf([num2str(structure.D1{indices.bounds.i(k)}.elements.M.num) '(k=' num2str(k-1) '), '])]; %#ok
            else
                docu = [docu sprintf(['0(k=' num2str(k-1) '), '])]; %#ok
            end
        end
        if ~isempty(structure.D1{indices.bounds.i(options.N)})
            docu = [docu sprintf(['or ' num2str(structure.D1{indices.bounds.i(options.N)}.elements.M.num) '(k=' num2str(options.N-1) ') elements.' '\n'])];
        else
            docu = [docu sprintf(['or 0(k=' num2str(options.N-1) ') elements.' '\n'])];
        end
    else
        docu = [docu sprintf([' A matrix with ' num2str(structure.D1{1}.elements.M.num) ' elements.' '\n'])];
    end
    if ~bFirst
        code = [code sprintf(', ')];
    else
        bFirst = false;
    end
    code = [code sprintf([options.types.data '* ' names.D1])];
    % D2
    docu = [docu sprintf([options.indent.code ' * @param ' names.D2 ' Return value.'])];
    if ~all(cellfun(@(c)(isfield(c,'elements')),structure.D2)) || ~all(cellfun(@(c)(c.elements.M.num),structure.D2) == structure.D2{1}.elements.M.num)
        docu = [docu sprintf(['\n' options.indent.code ' * ' repmat(' ',1,length(['@param ' names.D2 ' '])) ' A matrix with '])];
        for k=1:options.N-1
            if ~isempty(structure.D2{indices.bounds.i(k)})
                docu = [docu sprintf([num2str(structure.D2{indices.bounds.i(k)}.elements.M.num) '(k=' num2str(k-1) '), '])]; %#ok
            else
                docu = [docu sprintf(['0(k=' num2str(k-1) '), '])]; %#ok
            end
        end
        if ~isempty(structure.D2{indices.bounds.i(options.N)})
            docu = [docu sprintf(['or ' num2str(structure.D2{indices.bounds.i(options.N)}.elements.M.num) '(k=' num2str(options.N-1) ') elements.' '\n'])];
        else
            docu = [docu sprintf(['or 0(k=' num2str(options.N-1) ') elements.' '\n'])];
        end
    else
        docu = [docu sprintf([' A matrix with ' num2str(structure.D2{1}.elements.M.num) ' elements.' '\n'])];
    end
    if ~bFirst
        code = [code sprintf(', ')];
    else
        bFirst = false;
    end
    code = [code sprintf([options.types.data '* ' names.D2])];
    % tmp vector
    docu = [docu sprintf([options.indent.code ' * @param ' localNames.tmp ' A vector of dimension ' num2str(max(dims.bd)) ' used for storing intermediate results.' '\n'])];
    if ~bFirst
        code = [code sprintf(', ')];
    else
        bFirst = false;
    end
    code = [code sprintf([options.types.data '* ' localNames.tmp])];
    if max(dims.np) > 0
        % T matrix
        docu = [docu sprintf([options.indent.code ' * @param ' localNames.T ' A vector (matrix) of dimension ' num2str(max(dims.mp)) ' used for storing intermediate results.' '\n'])];
        if ~bFirst
            code = [code sprintf(', ')];
        else
            bFirst = false; %#ok
        end
        code = [code sprintf([options.types.data '* ' localNames.T])];
    end
    % Combine docu and function definition
    code = [docu sprintf([options.indent.code ' */' '\n']) code sprintf([') {' '\n'])];
    
    %% Generate code to compute auxilary vector tmp
    code = [code sprintf([options.indent.code options.indent.generic '/** Compute intermediate results involving divisions **/' '\n'])];
    % Iterate over all unique combinations of bounds with the same structure
    for i=1:indices.bounds.n
        % Construct check
        if indices.bounds.n > 1
            check = generateCheck(indices.bounds.map{i});
            code = [code, sprintf([options.indent.code options.indent.generic]) check sprintf([' /* Unique bounds combination #' num2str(i) ' */ ' '\n'])]; %#ok
            indent = [options.indent.code options.indent.generic options.indent.generic];
        else
            indent = [options.indent.code options.indent.generic];
        end
        k = indices.bounds.k(i); % Index of first bounds corresponding to this unique structure
        idx = 1;
        % For all components with lower and/or upper bounds
        J = find(options.bounds.lb(:,k) | options.bounds.ub(:,k))';
        for j=J
            li = sum(options.bounds.lb(1:j,k)); % Index of lower bound corresponding to j-th component (if there is a bound)
            lu = sum(options.bounds.ub(1:j,k)); % Index of upper bound corresponding to j-th component (if there is a bound)
            % TODO make 1.0 to 1.0f in float case using improved num2str, once available. Also apply to other pieces of code.
            if options.bounds.lb(j,k) && options.bounds.ub(j,k) % Both lower and upper bound
                code = [code sprintf([indent localNames.tmp '[' num2str(idx-1) '] = 1.0/(' localNames.s '[' num2str(li-1) ']*' localNames.s '[' num2str(dims.lb(k)) '+' num2str(lu-1) ']' ...
                                             ' + ' localNames.s '[' num2str(li-1) '] + ' localNames.s '[' num2str(dims.lb(k)) '+' num2str(lu-1) ']);' '\n'])]; %#ok
                info.flops.mul = info.flops.mul+1*length(indices.bounds.map{i});
                info.flops.add = info.flops.add+3*length(indices.bounds.map{i});
                info.flops.div = info.flops.div+1*length(indices.bounds.map{i});
            elseif options.bounds.lb(j,k)
                code = [code sprintf([indent localNames.tmp '[' num2str(idx-1) '] = 1.0/(1.0 + ' localNames.s '[' num2str(li-1) ']);' '\n'])]; %#ok
                info.flops.add = info.flops.add+1*length(indices.bounds.map{i});
                info.flops.div = info.flops.div+1*length(indices.bounds.map{i});
            elseif options.bounds.ub(j,k)
                code = [code sprintf([indent localNames.tmp '[' num2str(idx-1) '] = 1.0/(1.0 + ' localNames.s '[' num2str(dims.lb(k)) '+' num2str(lu-1) ']);' '\n'])]; %#ok
                info.flops.add = info.flops.add+1*length(indices.bounds.map{i});
                info.flops.div = info.flops.div+1*length(indices.bounds.map{i});
            else
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
            idx = idx+1;
        end
        if indices.bounds.n > 1
            code = [code sprintf([options.indent.code options.indent.generic '}' '\n'])]; %#ok
        end
    end
    
    %% Generate matrix M
    code = [code sprintf([options.indent.code options.indent.generic '/** Construct ' names.M '_k := [' names.dp '_k]''*diag({S_i}_{i=1}^nu)*[' names.dp '_k]+diag({' localNames.s '_{k,i+L+U}}_{i=1}^np)' '\n'...
                          options.indent.code options.indent.generic '     where S_i = 1 if i-th component of u has neither lower nor upper bounds' '\n' ...
                          options.indent.code options.indent.generic '               = ' localNames.s '_{k,j}/(1+' localNames.s '_{k,j}) if i-th component of u has only lower bounds (j is index of lower bound)' '\n' ...
                          options.indent.code options.indent.generic '               = ' localNames.s '_{k,L+j}/(1+' localNames.s '_{k,L+j}) if i-th component of u has only upper bounds (j is index of upper bound)' '\n' ...
                          options.indent.code options.indent.generic '               = ' localNames.s '_{k,j}*' localNames.s '_{k,L+l}/(' localNames.s '_{k,j}*' localNames.s '_{k,L+l} + ' localNames.s '_{k,j} + ' localNames.s '_{k,L+l})' ... 
                                                                                       'if i-th component of u has lower and upper bounds (j is index of lower bound, l is index of upper bound)' ...
                          options.indent.code options.indent.generic '     and nu is the number of inputs, np is the number of constraints, L is the number of lower bounds and U is the number of upper bounds. **/' '\n'])];
    % Iterate over combinations of unique upper/lower bound combinations and structures of dp
    for i=1:indices.bdAndDp.n
        % Construct check
        if indices.bdAndDp.n > 1
            check = generateCheck(indices.bdAndDp.map{i});
            code = [code sprintf([options.indent.code options.indent.generic]) check sprintf([' /* Unique bounds and dp combination #' num2str(i) ' */ ' '\n'])]; %#ok
            indent = [options.indent.code options.indent.generic options.indent.generic];
        else
            indent = [options.indent.code options.indent.generic];
        end
        k = indices.bdAndDp.k(i); % Index of first bounds corresponding to this unique structure
        dpIdx = indices.dp.i(k);
        if ~isempty(structure.dp{dpIdx})
            % Construct M
            % Iterate over (non-zero) elements of dp
            code = [code, sprintf([indent '/* Construct temporary matrix ' localNames.T ' */' '\n'])]; %#ok
            for j=1:structure.dp{dpIdx}.elements.M.num
                li = sum(options.bounds.lb(1:structure.dp{dpIdx}.elements.M.row(j),k)); % Index of lower bound corresponding to j-th component (if there is a bound)
                lu = sum(options.bounds.ub(1:structure.dp{dpIdx}.elements.M.row(j),k)); % Index of upper bound corresponding to j-th component (if there is a bound)
                lb = sum(options.bounds.lb(1:structure.dp{dpIdx}.elements.M.row(j),k) | options.bounds.ub(1:structure.dp{dpIdx}.elements.M.row(j),k)); % Index of lower or upper bound corresponding to j-th component (if there are bounds)
                % If it has lower and upper bound
                if options.bounds.lb(structure.dp{dpIdx}.elements.M.row(j),k) && options.bounds.ub(structure.dp{dpIdx}.elements.M.row(j),k)
                    code = [code, sprintf([indent localNames.T '[' num2str(j-1) '] = ' names.dp '[' num2str(j-1) ']*' localNames.s '[' num2str(li-1) ']*' localNames.s '[' num2str(dims.lb(k)) '+' num2str(lu-1) ']*' localNames.tmp '[' num2str(lb-1) '];'])]; %#ok
                    info.flops.mul = info.flops.mul+3*length(indices.bdAndDp.map{i});
                elseif options.bounds.lb(structure.dp{dpIdx}.elements.M.row(j),k)
                    code = [code, sprintf([indent localNames.T '[' num2str(j-1) '] = ' names.dp '[' num2str(j-1) ']*' localNames.s '[' num2str(li-1) ']*' localNames.tmp '[' num2str(lb-1) '];'])]; %#ok
                    info.flops.mul = info.flops.mul+2*length(indices.bdAndDp.map{i});
                elseif options.bounds.ub(structure.dp{dpIdx}.elements.M.row(j),k)
                    code = [code, sprintf([indent localNames.T '[' num2str(j-1) '] = ' names.dp '[' num2str(j-1) ']*' localNames.s '[' num2str(dims.lb(k)) '+' num2str(lu-1) ']*' localNames.tmp '[' num2str(lb-1) '];'])]; %#ok
                    info.flops.mul = info.flops.mul+2*length(indices.bdAndDp.map{i});
                else
                    code = [code, sprintf([indent localNames.T '[' num2str(j-1) '] = ' names.dp '[' num2str(j-1) '];'])]; %#ok
                end
                code = [code, sprintf(['/* Element (' num2str(structure.dp{dpIdx}.elements.M.row(j)) ',' num2str(structure.dp{dpIdx}.elements.M.col(j)) ') */' '\n'])]; %#ok
            end
            % Iterate over (non-zero) elements of M
            invIdx = indices.inv.i(k); % Index of unique inverse
            code = [code, sprintf([indent '/* Construct matrix ' names.M '_k */' '\n'])]; %#ok
            for j=1:structure.inv{invIdx}.structure.M.num
                code = [code sprintf([indent names.M '[' num2str(j-1) '] ='])]; %#ok
                % Extract row indices for dp and T to compute j-th element of M
                I = find(structure.dp{dpIdx}.structure.M.mat(:,structure.inv{invIdx}.elements.M.row(j))' & structure.dp{dpIdx}.structure.M.mat(:,structure.inv{invIdx}.elements.M.col(j))');
                bFirst = true;
                for l=I
                    if bFirst
                        bFirst = false;
                    else
                        code = [code sprintf(' +')];  %#ok
                        info.flops.add = info.flops.add+1*length(indices.bdAndDp.map{i});
                    end
                    idx1 = find(structure.dp{dpIdx}.elements.M.row == l & structure.dp{dpIdx}.elements.M.col == structure.inv{invIdx}.elements.M.row(j));
                    idx2 = find(structure.dp{dpIdx}.elements.M.row == l & structure.dp{dpIdx}.elements.M.col == structure.inv{invIdx}.elements.M.col(j));
                    code = [code sprintf([' ' names.dp '[' num2str(idx1-1) ']*' localNames.T '[' num2str(idx2-1) ']'])]; %#ok
                    info.flops.mul = info.flops.mul+1*length(indices.bdAndDp.map{i});
                end
                if structure.inv{invIdx}.elements.M.row(j) == structure.inv{invIdx}.elements.M.col(j)
                    code = [code sprintf([' + ' localNames.s '[' num2str(dims.lb(k)+dims.ub(k)) '+' num2str(structure.inv{invIdx}.elements.M.row(j)-1) ']'])]; %#ok
                    info.flops.add = info.flops.add+1*length(indices.bdAndDp.map{i});
                end
                code = [code sprintf(['; /* Element (' num2str(structure.inv{invIdx}.elements.M.row(j)) ',' num2str(structure.inv{invIdx}.elements.M.col(j)) ') */' '\n'])]; %#ok
            end
        else
            code = [code, sprintf([indent '/* Nothing to do */' '\n'])]; %#ok
        end
        % Close if
        if indices.bdAndDp.n > 1
            code = [code sprintf([options.indent.code options.indent.generic '}' '\n'])]; %#ok
        end
    end
    % Compute inverse
    for i=1:indices.inv.n
        % Construct check
        if indices.inv.n > 1
            check = generateCheck(indices.inv.map{i});
            code = [code sprintf([options.indent.code options.indent.generic]) check sprintf([' /* Unique inverse #' num2str(i) ' */ ' '\n'])]; %#ok
            indent = [options.indent.code options.indent.generic options.indent.generic];
        else
            indent = [options.indent.code options.indent.generic];
        end
        if ~isempty(structure.inv{i})
            code = [code, sprintf([indent '/* Compute inverse ' names.Mi ' of ' names.M ' */' '\n' ...
                                   indent structure.inv{i}.names.fun '(' names.Mi ', ' names.M ');' '\n'])]; %#ok
        else
            code = [code, sprintf([indent '/* Nothing to do */' '\n'])]; %#ok
        end
        % Close if
        if indices.inv.n > 1
            code = [code sprintf([options.indent.code options.indent.generic '}' '\n'])]; %#ok
        end
    end
        
    %% Construct m1 and m2
    code = [code sprintf([options.indent.code options.indent.generic '/** Construct ' names.m1 ' and ' names.m2 ', where' '\n' ...
                          options.indent.code options.indent.generic '     ' names.m1 ' = -diag({S_i}_{i=1}^L)*' names.dp '_k(L), where ' names.dp '_k(L) is a matrix with the rows of ' names.dp '_k for which u_k has lower bounds and' '\n' ...
                          options.indent.code options.indent.generic '                      with S_i = 1/(1+' localNames.s '_{k,i}) if i-th lower bound of u has no corresponding upper bound' '\n' ...
                          options.indent.code options.indent.generic '                               = ' localNames.s '_{k,L+j}/(' localNames.s '_{k,i}*' localNames.s '_{k,L+j} + ' localNames.s '_{k,i} + ' localNames.s '_{k,L+j})' ...
                                                                           'if the j-th upper bound of u corresponds to the i-th lower bound ("corresponds" means that they bound the same component of u)' '\n' ...
                          options.indent.code options.indent.generic '     ' names.m2 ' = diag({S_i}_{i=1}^U)*' names.dp '_k(U), where ' names.dp '_k(U) is a matrix with the rows of ' names.dp '_k for which u_k has upper bounds and' '\n' ...
                          options.indent.code options.indent.generic '                     with S_i = 1/(1+' localNames.s '_{k,L+i}) if i-th upper bound of u has no corresponding lower bound' '\n' ...
                          options.indent.code options.indent.generic '                              = ' localNames.s '_{k,j}/(' localNames.s '_{k,L+i}*' localNames.s '_{k,j} + ' localNames.s '_{k,L+i} + ' localNames.s '_{k,j})' ...
                                                                           'if the j-th lower bound of u corresponds to the i-th upper bound ("corresponds" means that they bound the same component of u)' '**/' '\n'])];
	% Structure of m1/m2 only depends on lb/ub and dp, however the values of both depend on lb, ub and dp
    for i=1:indices.bdAndDp.n
        % Construct check
        if indices.bdAndDp.n > 1
            check = generateCheck(indices.bdAndDp.map{i});
            code = [code sprintf([options.indent.code options.indent.generic]) check sprintf([' /* Unique bounds and dp combination #' num2str(i) ' */ ' '\n'])]; %#ok
            indent = [options.indent.code options.indent.generic options.indent.generic];
        else
            indent = [options.indent.code options.indent.generic];
        end
        k = indices.bdAndDp.k(i); % First stage of i-th unique combination
        dpIdx = indices.dp.i(k);
        % Iterate over non-zero elements of m1
        m1Idx = indices.lbAndDp.i(k);
        if ~isempty(structure.m1{m1Idx})
            if structure.m1{m1Idx}.elements.M.num > 0
                code = [code, sprintf([indent '/* Construct ' names.m1 ' */' '\n'])]; %#ok
            end
            for j=1:structure.m1{m1Idx}.elements.M.num
                % Get index of considered row of dp (component index of u)
                idx = find(cumsum(options.bounds.lb(:,k)) == structure.m1{m1Idx}.elements.M.row(j), 1, 'first');
                dpElIdx = structure.dp{dpIdx}.elements.M.data.indices(structure.dp{dpIdx}.elements.M.data.row == idx & structure.dp{dpIdx}.elements.M.data.col == structure.m1{m1Idx}.elements.M.col(j));
                % If there is also an upper upper bound corresponding to the current row of dp
                if options.bounds.ub(idx,k)
                    code = [code sprintf([indent names.m1 '[' num2str(j-1) '] = -' localNames.s '[' num2str(dims.lb(k)) '+' num2str(sum(options.bounds.ub(1:idx,k))-1) ']' ... 
                                                                              '*' localNames.tmp '[' num2str(sum(options.bounds.lb(1:idx,k) | options.bounds.ub(1:idx,k))-1) ']' ... 
                                                                              '*' names.dp '[' num2str(dpElIdx-1) ']; ' ...
                                                                              '/* Element #' num2str(j) ' of ' names.m1 ' (' num2str(structure.m1{m1Idx}.elements.M.row(j)) ',' num2str(structure.m1{m1Idx}.elements.M.col(j)) ') */' '\n'])]; %#ok
                    info.flops.mul = info.flops.mul+2*length(indices.lbAndDp.map{m1Idx});
                else
                    code = [code sprintf([indent names.m1 '[' num2str(j-1) '] = -' localNames.tmp '[' num2str(sum(options.bounds.lb(1:idx,k) | options.bounds.ub(1:idx,k))-1) ']' ... 
                                                                              '*' names.dp '[' num2str(dpElIdx-1) ']; ' ...
                                                                              '/* Element #' num2str(j) ' of ' names.m1 ' (' num2str(structure.m1{m1Idx}.elements.M.row(j)) ',' num2str(structure.m1{m1Idx}.elements.M.col(j)) ') */' '\n'])]; %#ok
                    info.flops.mul = info.flops.mul+1*length(indices.lbAndDp.map{m1Idx});
                end
            end
        else
            code = [code, sprintf([indent '/* Nothing to do */' '\n'])]; %#ok
        end
        % Iterate over non-zero elements of m2
        m2Idx = indices.ubAndDp.i(k);
        if ~isempty(structure.m2{m2Idx})
            if structure.m1{m2Idx}.elements.M.num > 0
                code = [code, sprintf([indent '/* Construct ' names.m2 ' */' '\n'])]; %#ok
            end
            for j=1:structure.m2{m2Idx}.elements.M.num
                % Get index of considered row of dp (component index of u)
                idx = find(cumsum(options.bounds.ub(:,k)) == structure.m2{m2Idx}.elements.M.row(j), 1, 'first');
                dpElIdx = structure.dp{dpIdx}.elements.M.data.indices(structure.dp{dpIdx}.elements.M.data.row == idx & structure.dp{dpIdx}.elements.M.data.col == structure.m2{m2Idx}.elements.M.col(j));
                % If there is also an upper upper bound corresponding to the current row of dp
                if options.bounds.lb(idx,k)
                    code = [code sprintf([indent names.m2 '[' num2str(j-1) '] = ' localNames.s '[' num2str(sum(options.bounds.lb(1:idx,k))-1) ']' ... 
                                                                              '*' localNames.tmp '[' num2str(sum(options.bounds.ub(1:idx,k) | options.bounds.lb(1:idx,k))-1) ']' ... 
                                                                              '*' names.dp '[' num2str(dpElIdx-1) ']; ' ...
                                                                              '/* Element #' num2str(j) ' of ' names.m2 ' (' num2str(structure.m2{m2Idx}.elements.M.row(j)) ',' num2str(structure.m2{m2Idx}.elements.M.col(j)) ') */' '\n'])]; %#ok
                    info.flops.mul = info.flops.mul+2*length(indices.ubAndDp.map{m2Idx});
                else
                    code = [code sprintf([indent names.m2 '[' num2str(j-1) '] = ' localNames.tmp '[' num2str(sum(options.bounds.ub(1:idx,k) | options.bounds.lb(1:idx,k))-1) ']' ... 
                                                                              '*' names.dp '[' num2str(dpElIdx-1) ']; ' ...
                                                                              '/* Element #' num2str(j) ' of ' names.m2 ' (' num2str(structure.m2{m2Idx}.elements.M.row(j)) ',' num2str(structure.m2{m2Idx}.elements.M.col(j)) ') */' '\n'])]; %#ok
                    info.flops.mul = info.flops.mul+1*length(indices.ubAndDp.map{m2Idx});
                end
            end
        else
            code = [code, sprintf([indent '/* Nothing to do */' '\n'])]; %#ok
        end
        % Close if check was added
        if indices.bdAndDp.n > 1
            code = [code sprintf([options.indent.code options.indent.generic '}' '\n'])]; %#ok
        end
    end
    
    %% Construct D1, D2
    % TODO add docu
    code = [code sprintf([options.indent.code options.indent.generic '/** Construct ' names.D1 ' and ' names.D2 ', where' ' **/' '\n'])];
    for i=1:indices.bounds.n
        % Construct check
        if indices.bounds.n > 1
            check = generateCheck(indices.bounds.map{i});
            code = [code sprintf([options.indent.code options.indent.generic]) check sprintf([' /* Unique bounds combination #' num2str(i) ' */ ' '\n'])]; %#ok
            indent = [options.indent.code options.indent.generic options.indent.generic];
        else
            indent = [options.indent.code options.indent.generic];
        end
        k = indices.bounds.k(i); % First stage of i-th unique combination
        % Iterate over non-zero elements of D1
        if structure.D1{i}.elements.M.num > 0
            code = [code, sprintf([indent '/* Construct ' names.D1 ' */' '\n'])]; %#ok
        end
        for j=1:structure.D1{i}.elements.M.num
            % "Diagonal" part of D1
            if structure.D1{i}.elements.M.row(j) <= dims.lb(k)
                % Sanity check
                if ~options.bounds.lb(indices.lb.j(i,structure.D1{i}.elements.M.row(j)),k)
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                idx = indices.lb.j(i,structure.D1{i}.elements.M.row(j)); % Index of component of u corresponding to row(j)-th lower bound
                code = [code sprintf([indent names.D1 '[' num2str(j-1) '] = ' localNames.tmp '[' num2str(sum(options.bounds.lb(1:idx,k) | options.bounds.ub(1:idx,k))-1) ']'])]; %#ok
                % If the row(j)-th lower bound has a corresponding upper bound
                if options.bounds.ub(idx,k)
                    code = [code sprintf(['*(1.0 + ' localNames.s '[' num2str(dims.lb(k)) '+' num2str(sum(options.bounds.ub(1:idx,k))-1) '])'])]; %#ok
                    info.flops.add = info.flops.add + 1*length(indices.bounds.map{i});
                    info.flops.mul = info.flops.mul + 1*length(indices.bounds.map{i});
                end
                code = [code sprintf(['; ' '/* Element #' num2str(j) ' of ' names.D1 ' (' num2str(structure.D1{i}.elements.M.row(j)) ',' num2str(structure.D1{i}.elements.M.col(j)) ')' ...
                                           ', corresponds to lower bound #' num2str(structure.D1{i}.elements.M.row(j))])]; %#ok
                if options.bounds.ub(idx,k)
                    code = [code sprintf([' and upper bound #' num2str(sum(options.bounds.ub(1:idx,k)))])]; %#ok
                end
                code = [code sprintf([', component (' num2str(idx) ') */' '\n'])]; %#ok
            % "Off-Diagonal" part of D1
            elseif indices.lb.j(i,structure.D1{i}.elements.M.col(j)) == indices.ub.j(i,structure.D1{i}.elements.M.row(j)-dims.lb(k)) % Sanity check, should always be true
                idx = indices.lb.j(i,structure.D1{i}.elements.M.col(j));
                if sum(options.bounds.lb(1:idx,k)) ~= structure.D1{i}.elements.M.col(j) || sum(options.bounds.ub(1:idx,k)) ~= structure.D1{i}.elements.M.row(j)-dims.lb(k)
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code sprintf([indent names.D1 '[' num2str(j-1) '] = ' localNames.tmp '[' num2str(sum(options.bounds.lb(1:idx,k) | options.bounds.ub(1:idx,k))-1) ']; ' ...
                                      '/* Element #' num2str(j) ' of ' names.D1 ' (' num2str(structure.D1{i}.elements.M.row(j)) ',' num2str(structure.D1{i}.elements.M.col(j)) ')' ...
                                      ', corresponds to lower bound #' num2str(sum(options.bounds.lb(1:idx,k))) ' and upper bound #' num2str(sum(options.bounds.ub(1:idx,k))) ' */' '\n'])]; %#ok
            else
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
        end
        % Iterate over non-zero elements of D2
        if structure.D2{i}.elements.M.num > 0
            code = [code, sprintf([indent '/* Construct ' names.D2 ' */' '\n'])]; %#ok
        end
        for j=1:structure.D2{i}.elements.M.num
            % "Diagonal" part of D2
            if structure.D2{i}.elements.M.row(j) > dims.lb(k)
                % Sanity check
                if ~options.bounds.ub(indices.ub.j(i,structure.D2{i}.elements.M.row(j)-dims.lb(k)),k)
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                idx = indices.ub.j(i,structure.D2{i}.elements.M.row(j)-dims.lb(k)); % Index of component of u corresponding to (row(j)-dims.lb)-th upper bound
                code = [code sprintf([indent names.D2 '[' num2str(j-1) '] = ' localNames.tmp '[' num2str(sum(options.bounds.ub(1:idx,k) | options.bounds.lb(1:idx,k))-1) ']'])]; %#ok
                % If the row(j)-th lower bound has a corresponding upper bound
                if options.bounds.lb(idx,k)
                    code = [code sprintf(['*(1.0 + ' localNames.s '[' num2str(sum(options.bounds.lb(1:idx,k))-1) '])'])]; %#ok
                    info.flops.add = info.flops.add + 1*length(indices.bounds.map{i});
                    info.flops.mul = info.flops.mul + 1*length(indices.bounds.map{i});
                end
                code = [code sprintf(['; ' '/* Element #' num2str(j) ' of ' names.D2 ' (' num2str(structure.D2{i}.elements.M.row(j)) ',' num2str(structure.D2{i}.elements.M.col(j)) ')' ...
                                           ', corresponds to upper bound #' num2str(structure.D2{i}.elements.M.row(j)-dims.lb(k))])]; %#ok
                if options.bounds.lb(idx,k)
                    code = [code sprintf([' and lower bound #' num2str(sum(options.bounds.lb(1:idx,k)))])]; %#ok
                end
                code = [code sprintf([', component (' num2str(idx) ') */' '\n'])]; %#ok
            % "Off-Diagonal" part of D2
            elseif indices.lb.j(i,structure.D2{i}.elements.M.row(j)) == indices.ub.j(i,structure.D1{i}.elements.M.col(j)) % Sanity check, should always be true
                idx = indices.ub.j(i,structure.D2{i}.elements.M.col(j));
                if sum(options.bounds.lb(1:idx,k)) ~= structure.D2{i}.elements.M.row(j) || sum(options.bounds.ub(1:idx,k)) ~= structure.D2{i}.elements.M.col(j)
                    throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
                end
                code = [code sprintf([indent names.D2 '[' num2str(j-1) '] = ' localNames.tmp '[' num2str(sum(options.bounds.lb(1:idx,k) | options.bounds.ub(1:idx,k))-1) ']; ' ...
                                      '/* Element #' num2str(j) ' of ' names.D2 ' (' num2str(structure.D2{i}.elements.M.row(j)) ',' num2str(structure.D2{i}.elements.M.col(j)) ')' ...
                                      ', corresponds to lower bound #' num2str(sum(options.bounds.lb(1:idx,k))) ' and upper bound #' num2str(sum(options.bounds.ub(1:idx,k))) ' */' '\n'])]; %#ok
            else
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
        end
    end
    
    if indices.bounds.n > 1
        code = [code sprintf([options.indent.code options.indent.generic '}' '\n'])];
    end
    
    %% Close function
    code = [code sprintf([options.indent.code '}' '\n'])];
end

%%
%
function [indices, structure] = extractUnique(dp, dt, dims, options)
    % indices.something.map is a cell array where each cell contains the stage indices for which "something" has the same structure
    % 

    %% Identify unique of lower bounds
    lbIndices = nan(1,options.N);
    idx = 1;
    while any(isnan(lbIndices))
        k = find(isnan(lbIndices), 1, 'first'); % Get index of first non-assigned combination of bounds
        lbIndices(all(repmat(options.bounds.lb(:,k),1,options.N) == options.bounds.lb,1)) = idx;
        idx = idx+1;
    end
    indices.lb.n = length(unique(lbIndices));
    indices.lb.map = cell(1,indices.lb.n);
    for i=1:indices.lb.n
        indices.lb.map{i} = find(lbIndices == i); % Indices of stages belonging to the i-th unique combination of lower/upper bounds
    end
    indices.lb.k = @(i)(indices.lb.map{i}(1));
    indices.lb.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.lb.map)));
    indices.lb.j = @(i,idx)(find(cumsum(options.bounds.lb(:,indices.lb.k(i))) == idx,1,'first')); % Component index j of idx-th lower bound (for i-th unique bound)

    %% Identify unique of upper bounds
    ubIndices = nan(1,options.N);
    idx = 1;
    while any(isnan(ubIndices))
        k = find(isnan(ubIndices), 1, 'first'); % Get index of first non-assigned combination of bounds
        ubIndices(all(repmat(options.bounds.ub(:,k),1,options.N) == options.bounds.ub,1)) = idx;
        idx = idx+1;
    end
    indices.ub.n = length(unique(ubIndices));
    indices.ub.map = cell(1,indices.ub.n);
    for i=1:indices.ub.n
        indices.ub.map{i} = find(ubIndices == i); % Indices of stages belonging to the i-th unique combination of lower/upper bounds
    end
    indices.ub.k = @(i)(indices.ub.map{i}(1));
    indices.ub.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.ub.map)));
    indices.ub.j = @(i,idx)(find(cumsum(options.bounds.ub(:,indices.ub.k(i))) == idx,1,'first')); % Component index j of idx-th upper bound (for i-th unique bound)
    
    %% Identify unique combinations of lower/upper bounds
    boundIndices = nan(1,options.N);    
    idx = 1;
    while any(isnan(boundIndices))
        k = find(isnan(boundIndices), 1, 'first'); % Get index of first non-assigned combination of bounds
        boundIndices(all(repmat(options.bounds.lb(:,k),1,options.N) == options.bounds.lb,1) & all(repmat(options.bounds.ub(:,k),1,options.N) == options.bounds.ub,1)) = idx;
        if any(options.bounds.lb(:,k) | options.bounds.ub(:,k))
            D = zeros(dims.lb(k),dims.ub(k));
            for i=1:dims.nu(k)
                if options.bounds.lb(i,k) && options.bounds.ub(i,k)
                    D(sum(options.bounds.lb(1:i,k)),sum(options.bounds.ub(1:i,k))) = 1;
                end
            end
            structureInfo = falcopt.detectMatrixStructure([eye(dims.lb(k)); D'], 'structure', 'sparse');
            structure.D1{idx} = structureInfo;
            structureInfo = falcopt.detectMatrixStructure([D; eye(dims.ub(k))], 'structure', 'sparse');
            structure.D2{idx} = structureInfo;
        else
            structure.D1{idx} = [];
            structure.D2{idx} = [];
        end
        idx = idx+1;
    end
    indices.bounds.n = length(unique(boundIndices));
    indices.bounds.map = cell(1,indices.bounds.n);
    for i=1:indices.bounds.n
        indices.bounds.map{i} = find(boundIndices == i); % Indices of stages belonging to the i-th unique combination of lower/upper bounds
    end
    indices.bounds.k = @(i)(indices.bounds.map{i}(1));
    indices.bounds.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.bounds.map)));

    %% Identify unique dp's
    dpIndices = nan(1,options.N);
    idx = 1;
    while any(isnan(dpIndices))
        k = find(isnan(dpIndices), 1, 'first'); % Get index of first dp for which the structure has not been generated, yet
        if dims.np(k) > 0
            structureInfo = falcopt.detectMatrixStructure(dp{k}, 'structure', options.structure.dp);
            structure.dp{idx} = structureInfo;
        else
            structure.dp{idx} = [];
        end
        I = cellfun(@(c)(all(size(dp{k}) == size(c)) && all(all((dp{k}~=0) == (c~=0)))), dp) & isnan(dpIndices);
        dpIndices(I) = idx;
        idx = idx+1;
    end
    indices.dp.n = length(unique(dpIndices));
    indices.dp.map = cell(1,indices.dp.n);
    for i=1:indices.dp.n
        indices.dp.map{i} = find(dpIndices == i); % Indices of stages belonging to the i-th unique dp
    end
    indices.dp.k = @(i)(indices.dp.map{i}(1));
    indices.dp.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.dp.map)));

    %% Identify and store information about unique combinations of bounds and dp
    indices.bdAndDp.map = indices.dp.map;
    idx = 1;
    while true
        % Find unique bound indices correpsonding to each stage for one unique dp
        bdIndices = zeros(1,length(indices.bdAndDp.map{idx}));
        for i=1:length(indices.bdAndDp.map{idx})
            bdIndices(i) = find(cellfun(@(c)(any(c==indices.bdAndDp.map{idx}(i))), indices.bounds.map));
        end
        ubdIndices = unique(bdIndices);
        maps = cell(1,length(ubdIndices));
        for i=1:length(ubdIndices)
            maps{i} = indices.bdAndDp.map{idx}(bdIndices == ubdIndices(i));
        end
        indices.bdAndDp.map = [indices.bdAndDp.map(1:idx-1) maps indices.bdAndDp.map(idx+1:end)];
        idx = idx+length(ubdIndices);
        if idx > length(indices.bdAndDp.map)
            break;
        end
    end
    indices.bdAndDp.n = length(indices.bdAndDp.map);
    indices.bdAndDp.k = @(i)(indices.bdAndDp.map{i}(1));
    indices.bdAndDp.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.bdAndDp.map)));

    %% Identify and store information about unique combinations of lower bounds and dp
    indices.lbAndDp.map = indices.dp.map;
    idx = 1;
    while true
        % Find unique bound indices correpsonding to each stage for one unique dp
        bdIndices = zeros(1,length(indices.lbAndDp.map{idx}));
        for i=1:length(indices.lbAndDp.map{idx})
            bdIndices(i) = find(cellfun(@(c)(any(c==indices.lbAndDp.map{idx}(i))), indices.lb.map));
        end
        ubdIndices = unique(bdIndices);
        maps = cell(1,length(ubdIndices));
        for i=1:length(ubdIndices)
            maps{i} = indices.lbAndDp.map{idx}(bdIndices == ubdIndices(i));
        end
        for i=idx:idx+length(ubdIndices)-1
            k = maps{i-idx+1}(1);
            if dims.np(k) > 0 && any(options.bounds.lb(:,k))
                structureInfo = falcopt.detectMatrixStructure(dp{k}(options.bounds.lb(:,k),:), 'structure', options.structure.dp);
                structure.m1{i} = structureInfo;
            else
                structure.m1{i} = [];
            end
        end
        indices.lbAndDp.map = [indices.lbAndDp.map(1:idx-1) maps indices.lbAndDp.map(idx+1:end)];
        idx = idx+length(ubdIndices);
        if idx > length(indices.lbAndDp.map)
            break;
        end
    end
    indices.lbAndDp.n = length(indices.lbAndDp.map);
    indices.lbAndDp.k = @(i)(indices.lbAndDp.map{i}(1));
    indices.lbAndDp.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.lbAndDp.map)));

    %% Identify and store information about unique combinations of upper bounds and dp
    indices.ubAndDp.map = indices.dp.map;
    idx = 1;
    while true
        % Find unique bound indices correpsonding to each stage for one unique dp
        bdIndices = zeros(1,length(indices.ubAndDp.map{idx}));
        for i=1:length(indices.ubAndDp.map{idx})
            bdIndices(i) = find(cellfun(@(c)(any(c==indices.ubAndDp.map{idx}(i))), indices.ub.map));
        end
        ubdIndices = unique(bdIndices);
        maps = cell(1,length(ubdIndices));
        for i=1:length(ubdIndices)
            maps{i} = indices.ubAndDp.map{idx}(bdIndices == ubdIndices(i));
        end
        for i=idx:idx+length(ubdIndices)-1
            k = maps{i-idx+1}(1);
            if dims.np(k) > 0 && any(options.bounds.ub(:,k))
                structureInfo = falcopt.detectMatrixStructure(dp{k}(options.bounds.ub(:,k),:), 'structure', options.structure.dp);
                structure.m2{i} = structureInfo;
            else
                structure.m2{i} = [];
            end
        end
        indices.ubAndDp.map = [indices.ubAndDp.map(1:idx-1) maps indices.ubAndDp.map(idx+1:end)];
        idx = idx+length(ubdIndices);
        if idx > length(indices.ubAndDp.map)
            break;
        end
    end
    indices.ubAndDp.n = length(indices.ubAndDp.map);
    indices.ubAndDp.k = @(i)(indices.ubAndDp.map{i}(1));
    indices.ubAndDp.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.ubAndDp.map)));
    
    %% Identify unique dt's
    dtIndices = nan(1,options.N);
    idx = 1;
    while any(isnan(dtIndices))
        k = find(isnan(dtIndices), 1, 'first'); % Get index of first dt for which the structure has not been generated, yet
        if dims.nt(k) > 0
            structureInfo = falcopt.detectMatrixStructure(dt{k}, 'structure', options.structure.dt);
            structure.dt{idx} = structureInfo;
        else
            structure.dt{idx} = [];
        end
        I = cellfun(@(c)(all(size(dt{k}) == size(c)) && all(all((dt{k}~=0) == (c~=0)))), dt) & isnan(dtIndices);
        dtIndices(I) = idx;
        idx = idx+1;
    end
    indices.dt.n = length(unique(dtIndices));
    indices.dt.map = cell(1,indices.dt.n);
    for i=1:indices.dt.n
        indices.dt.map{i} = find(dtIndices == i); % Indices of stages belonging to the i-th unique dt
    end
    indices.dt.k = @(i)(indices.dt.map{i}(1));
    indices.dt.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.dt.map)));
    
    %% Identify and store information about unique combinations of bounds and dt
    indices.bdAndDt.map = indices.dt.map;
    idx = 1;
    while true
        % Find unique bound and dp indices correpsonding to each stage for one unique dt
        I = zeros(1,length(indices.bdAndDt.map{idx}));
        for i=1:length(indices.bdAndDt.map{idx})
            I(i) = find(cellfun(@(c)(any(c==indices.bdAndDt.map{idx}(i))), indices.bounds.map));
        end
        uI = unique(I);
        maps = cell(1,length(uI));
        for i=1:length(uI)
            maps{i} = indices.bdAndDt.map{idx}(I == uI(i));
        end
        indices.bdAndDt.map = [indices.bdAndDt.map(1:idx-1) maps indices.bdAndDt.map(idx+1:end)];
        % Generate structure
        % TODO generate full structure information, not just basics
        for i=idx:idx+length(uI)-1
            k = indices.bdAndDt.map{i}(1);
            dtIdx = indices.dt.i(k);
            if ~isempty(dt{dtIdx}) && size(structure.dt{dtIdx}.structure.M.mat,2) > 1
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
            % D1
            D1Structure = structure.D1{indices.bounds.i(k)};
            if ~isempty(D1Structure) && ~isempty(dt{dtIdx})
                structure.bdAndDt.D1{i}.structure.M.type = 'indexed';
                mat = D1Structure.structure.M.mat';
                mat(mat ~= 0) = 1:D1Structure.structure.M.num; % Indexed elements
                mat = mat';
                mat = mat(:,structure.dt{dtIdx}.structure.M.mat(options.bounds.lb(:,k),:) ~= 0); % Remove columns depending on structure of dt
                % Add zero-columns for components of u for which dt is non-zero but no lower bound is specified
                structure.bdAndDt.D1{i}.structure.M.mat = zeros(dims.lb(k)+dims.ub(k),structure.dt{dtIdx}.structure.M.num);
                structure.bdAndDt.D1{i}.structure.M.mat(:,options.bounds.lb(structure.dt{dtIdx}.structure.M.mat ~= 0,k)) = mat;
                structure.bdAndDt.D1{i}.structure.M.num = sum(sum(structure.bdAndDt.D1{i}.structure.M.mat ~= 0));
            else
                structure.bdAndDt.D1{i} = [];
            end
            % D2
            D2Structure = structure.D2{indices.bounds.i(k)};
            if ~isempty(D2Structure) && ~isempty(dt{dtIdx})
                structure.bdAndDt.D2{i}.structure.M.type = 'indexed';
                mat = D2Structure.structure.M.mat';
                mat(mat ~= 0) = 1:D2Structure.structure.M.num; % Indexed elements
                mat = mat';
                mat = mat(:,structure.dt{dtIdx}.structure.M.mat(options.bounds.ub(:,k),:) ~= 0); % Remove columns depending on structure of dt
                % Add zero-columns for components of u for which dt is non-zero but no lower bound is specified
                structure.bdAndDt.D2{i}.structure.M.mat = zeros(dims.lb(k)+dims.ub(k),structure.dt{dtIdx}.structure.M.num);
                structure.bdAndDt.D2{i}.structure.M.mat(:,options.bounds.ub(structure.dt{dtIdx}.structure.M.mat ~= 0,k)) = mat;
                structure.bdAndDt.D2{i}.structure.M.num = sum(sum(structure.bdAndDt.D2{i}.structure.M.mat ~= 0));
            else
                structure.bdAndDt.D2{i} = [];
            end
            % dtL
            if ~isempty(structure.dt{dtIdx}) && ~isempty(dt{dtIdx})
                structure.bdAndDt.dtL{i}.structure.M.type = 'indexed';
                mat = structure.dt{dtIdx}.structure.M.mat';
                mat(mat ~= 0) = 1:structure.dt{dtIdx}.structure.M.num; % Indexed elements
                structure.bdAndDt.dtL{i}.structure.M.mat = mat(:,options.bounds.lb(:,k));
                structure.bdAndDt.dtL{i}.structure.M.num = sum(sum(structure.bdAndDt.dtL{i}.structure.M.mat ~= 0));
            else
                structure.bdAndDt.dtL{i} = [];
            end
            % dtU
            if ~isempty(structure.dt{dtIdx}) && ~isempty(dt{dtIdx})
                structure.bdAndDt.dtU{i}.structure.M.type = 'indexed';
                mat = structure.dt{dtIdx}.structure.M.mat';
                mat(mat ~= 0) = 1:structure.dt{dtIdx}.structure.M.num; % Indexed elements
                structure.bdAndDt.dtU{i}.structure.M.mat = mat(:,options.bounds.ub(:,k) ~= 0); % Remove columns depending on structure of dt and upper bounds
                structure.bdAndDt.dtU{i}.structure.M.num = sum(sum(structure.bdAndDt.dtU{i}.structure.M.mat ~= 0));
            else
                structure.bdAndDt.dtU{i} = [];
            end
        end
        idx = idx+length(uI);
        if idx > length(indices.bdAndDt.map)
            break;
        end
    end
    indices.bdAndDt.n = length(indices.bdAndDt.map);
    indices.bdAndDt.k = @(i)(indices.bdAndDt.map{i}(1));
    indices.bdAndDt.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.bdAndDt.map)));
    
    %% Identify and store information about unique combinations of dp and dt
    indices.dpAndDt.map = indices.dt.map;
    idx = 1;
    while true
        % Find unique dp indices corresponding to each stage for one unique dt
        I = zeros(1,length(indices.dpAndDt.map{idx}));
        for i=1:length(indices.dpAndDt.map{idx})
            I(i) = find(cellfun(@(c)(any(c==indices.dpAndDt.map{idx}(i))), indices.dp.map));
        end
        uI = unique(I);
        maps = cell(1,length(uI));
        for i=1:length(uI)
            maps{i} = indices.dpAndDt.map{idx}(I == uI(i));
        end
        indices.dpAndDt.map = [indices.dpAndDt.map(1:idx-1) maps indices.dpAndDt.map(idx+1:end)];
        % Generate structure
        for i=idx:idx+length(uI)-1
            k = maps{i-idx+1}(1);
            dpIdx = indices.dp.i(k);
            if dims.np(k) > 0 && dims.nt(k) > 0
                structure.dpAndDt{i}.structure.M.type = 'indexed';
                mat = structure.dp{dpIdx}.structure.M.mat';
                mat(mat ~= 0) = 1:structure.dp{dpIdx}.structure.M.num; % Indexed elements
                structure.dpAndDt{i}.structure.M.mat = mat(:,structure.dt{indices.dt.i(k)}.structure.M.mat ~= 0); % Only select columns for non-zero elements of dt
                structure.dpAndDt{i}.structure.M.num = sum(sum(structure.dpAndDt{i}.structure.M.mat ~= 0));
            else
                structure.dpAndDt{i} = [];
            end
        end
        idx = idx+length(uI);
        if idx > length(indices.bdAndDt.map)
            break;
        end
    end
    indices.dpAndDt.n = length(indices.dpAndDt.map);
    indices.dpAndDt.k = @(i)(indices.dpAndDt.map{i}(1));
    indices.dpAndDt.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.dpAndDt.map)));
    
    %% Identify and store information about unique combinations of bounds, dp and dt
    indices.bdAndDpAndDt.map = indices.dt.map;
    idx = 1;
    while true
        % Find unique bound and dp indices corresponding to each stage for one unique dt
        I = zeros(1,length(indices.bdAndDpAndDt.map{idx}));
        for i=1:length(indices.bdAndDpAndDt.map{idx})
            I(i) = find(cellfun(@(c)(any(c==indices.bdAndDpAndDt.map{idx}(i))), indices.bdAndDp.map));
        end
        uI = unique(I);
        maps = cell(1,length(uI));
        for i=1:length(uI)
            maps{i} = indices.bdAndDpAndDt.map{idx}(I == uI(i));
        end
        indices.bdAndDpAndDt.map = [indices.bdAndDpAndDt.map(1:idx-1) maps indices.bdAndDpAndDt.map(idx+1:end)];
        % Generate structure
        % TODO generate full structure information, not just basics
        for i=idx:idx+length(uI)-1
            k = indices.bdAndDpAndDt.map{i}(1);
            dtIdx = indices.dt.i(k);
            if ~isempty(structure.dt{dtIdx}) && size(structure.dt{dtIdx}.structure.M.mat,2) > 1
                throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
            end
            m1Structure = structure.m1{indices.lbAndDp.i(k)};
            if ~isempty(m1Structure) && dims.nt(k) > 0
                structure.bdAndDpAndDt.lb{i}.structure.M.type = 'indexed';
                mat = m1Structure.structure.M.mat';
                mat(mat ~= 0) = 1:m1Structure.structure.M.num; % Indexed elements
                mat = mat';
                mat = mat(structure.dt{dtIdx}.structure.M.mat(options.bounds.lb(:,k),:) ~= 0,:); % Remove rows depending on structure of dt
                % Add zero-rows for components of u for which dt is non-zero but no lower bound is specified
                structure.bdAndDpAndDt.lb{i}.structure.M.mat = zeros(structure.dt{dtIdx}.structure.M.num,dims.np(k));
                structure.bdAndDpAndDt.lb{i}.structure.M.mat(options.bounds.lb(structure.dt{dtIdx}.structure.M.mat ~= 0,k),:) = mat;
                structure.bdAndDpAndDt.lb{i}.structure.M.num = sum(sum(structure.bdAndDpAndDt.lb{i}.structure.M.mat ~= 0));
            else
                structure.bdAndDpAndDt.lb{i} = [];
            end
            m2Structure = structure.m2{indices.ubAndDp.i(k)};
            if ~isempty(m2Structure) && dims.nt(k) > 0
                structure.bdAndDpAndDt.ub{i}.structure.M.type = 'indexed';
                mat = m2Structure.structure.M.mat';
                mat(mat ~= 0) = 1:m2Structure.structure.M.num; % Indexed elements
                mat = mat';
                mat = mat(structure.dt{dtIdx}.structure.M.mat(options.bounds.ub(:,k),:) ~= 0,:); % Remove rows depending on structure of dt
                % Add zero-rows for components of u for which dt is non-zero but no upper bound is specified
                structure.bdAndDpAndDt.ub{i}.structure.M.mat = zeros(structure.dt{dtIdx}.structure.M.num,dims.np(k));
                structure.bdAndDpAndDt.ub{i}.structure.M.mat(options.bounds.ub(structure.dt{dtIdx}.structure.M.mat ~= 0,k),:) = mat;
                structure.bdAndDpAndDt.ub{i}.structure.M.num = sum(sum(structure.bdAndDpAndDt.ub{i}.structure.M.mat ~= 0));
            else
                structure.bdAndDpAndDt.ub{i} = [];
            end
        end
        idx = idx+length(uI);
        if idx > length(indices.bdAndDpAndDt.map)
            break;
        end
    end
    indices.bdAndDpAndDt.n = length(indices.bdAndDpAndDt.map);
    indices.bdAndDpAndDt.k = @(i)(indices.bdAndDpAndDt.map{i}(1));
    indices.bdAndDpAndDt.i = @(k)(find(cellfun(@(c)(any(c==k)),indices.bdAndDpAndDt.map)));

end

function check = generateCheck(K)
    ranges = reshape([K([diff(K)-1 ~= 0 false]); K([false diff(K)-1 ~= 0])],[],1)';
    if length(K) > 1 && K(1)+1 == K(2)
        ranges = [K(1) ranges];
    elseif length(K) == 1
        ranges = [K(1) K(1)];
    end
    if mod(length(ranges),2) == 1
        ranges = [ranges K(end)];
    end
    ranges = reshape(ranges,2,length(ranges)/2);
    rangeStrs = cell(1,size(ranges,2));
    for j=1:size(ranges,2)
        if ranges(1,j) == ranges(2,j)
            rangeStrs{j} = sprintf(['k == ' num2str(ranges(2,j)-1) ]);
        elseif ranges(1,j) < ranges(2,j)
            rangeStrs{j} = sprintf([num2str(ranges(1,j)-1) ' <= k <= ' num2str(ranges(2,j)-1)]);
        else
            throw(MException('falcopt:generateConstraintInv:ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
        end
    end
    check = sprintf(['if(' strjoin(rangeStrs, ' || ') ') {']);
end
