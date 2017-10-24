%% generateData Generate static data for a vector or matrix,
% additionally generates code to transform from and to MATLAB
% 
% [data, code, info] = generateData(M)
% or
% [data, code, info] = generateData(M, ...) [with options]
% 
% Returns a string 'data' containing the static data and
%  returns an info struct containing information about the number of static
%  data points in 'size' and the used names in 'names' 
%
% the following options are available:
%  .names     - A struct with function and variable names to be used. Default is
%                struct('M', 'M', 'prefix', '', 'toMATLAB', 'MtoMATLAB', 'fromMATLAB', 'MfromMATLAB')
%  .structure - A string that determines the structure of the matrix or vector. Default is 'unique'.
%                the possibilities are:
%                'dense' for no structure (even zero elements are stored),
%                'sparse' for ommitted zero elements (only non-zero elements of non-static matrices are considered)
%                'unique' for ommitted zero and reccurring elements (recurring elements of non-static matrices are assumed to be the same)
%                'ordered' for ommitted zero and reccurring elements including permutation of elements (elements of non-static matrices are stored in ascending order of their values)
%                'indexed' for ommitted zero and recurring elements including permutation of elements (elements of non-static matrices are stored according to their values as indices starting from 1)
%  .noones    - A boolean, if true, then 1 and -1 are not stored and
%                info.ones contains the structure. Default: false.
%  .type      - A string that indicates the data type, either 'single' or 'double'.
%                Default is 'double'.
%  .precision - A string that determines the precision of computation and stored data in the generated code. 
%                Needs to be either 'single' or 'double'. Default: '' (unspecified, will try to extract from type);
%  .types     - A string that indicates the data type, can be 'float', 'double' or a custom type.
%                Needs to be compatible with the floating-point type of .precision. Default: double
%  .transpose - A boolean, if true then the matrix is transposed.
%  .symmetric - A boolean, if true then the matrix is considered symmetric and only half is stored.
%                Default: false.
%                NOTE: Matrices defined as symmetric are always symmetrized via (M+M')/2
%  .indent    - Indentation to be used in data and code generation
%  .static    - Whether the 'static' keyword should be used for the data. Default: true
%  .const     - Whether the 'const' keyword should be used for the data. Default: true
%  .verbose   - Level of procedural output of this function
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
function [data, code, info] = generateData(varargin)
    functionId = 'generateData';
    
    defaultNames = struct('M', 'M', ...
                          'prefix', '', ...
                          'toDense', 'toDense', ...
                          'fromDense', 'fromDense');
    structures = {'dense', 'sparse', 'unique', 'ordered', 'indexed'};
    indentTypes = {'generic', 'data', 'code'};
    defaultStructure = 'unique';
    p = inputParser;
    p.CaseSensitive = true;
    p.addRequired('M', @isnumeric);
    p.addParameter('names', defaultNames, @isstruct);
    p.addParameter('structure', defaultStructure, @(x)((ischar(x) && any(strcmp(x, structures)))));
    p.addParameter('noones', false, @islogical);
    p.addParameter('stripped', false, @islogical);
    p.addParameter('type', 'double', @ischar);
    p.addParameter('precision', '', @(s)(ischar(s) && any(strcmp(s, {'single', 'double', 'integer'}))));
    p.addParameter('transpose', false, @islogical);
    p.addParameter('symmetric', false, @islogical);
    p.addParameter('indent', '    ', @(x)(ischar(x) || (isstruct(x) && isfield(x, 'generic') && all(cellfun(@(y)(~isfield(x, y) || ischar(x.(y))), indentTypes)))));
    p.addParameter('static', true, @islogical);
    p.addParameter('const', true, @islogical);
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && x >= 0));
    p.parse(varargin{:});
    options = p.Results;
    M = options.M;
    names = options.names;
    
    if options.verbose >= 2
        fprintf('Code generation started...\n');
        fprintf('. processing options\n');
    end
    
    dims = struct();
    
    %% Store and check dimension of M
    % Store dimensions
    dims.n = size(M,1);
    dims.m = size(M,2);
    if min(dims.m, dims.n) == 0
        throw(MException([functionId ':InvalidDimension'], 'The input matrix is empty (it has a zero-dimension).'));
    end
    
    %% Check names
    fields = fieldnames(defaultNames);
    for f=1:length(fields)
        % Set defaults
        if ~isfield(names, fields{f})
            names.(fields{f}) = defaultNames.(fields{f});
        end
    end
    for f=1:length(fields)
        % Add prefix
        if ~strcmp(fields{1}, 'prefix')
            names.(fields{f}) = [names.prefix names.(fields{f})];
        end
    end
    names.prefix = ''; % Remove prefix, since it has been incorporated
    info.names = names;
    
    %% Check optional properties
    % Ensure and exploit symmetry
    if options.symmetric
        warning([functionId ':MissingImplemementation'], 'The special treatment of symmetric matrices has not yet been implemented. Ignoring...');
    end
    if options.symmetric && norm(M - (M+M')/2, inf) > eps
        warning([functionId ':InvalidData'], ['The matrix ' names.M ' is not symmetric (max. error: ' num2str(norm(M - (M+M')/2, inf)) '). Symmetrizing...']);
    end
    if options.symmetric && options.transpose
        options.transpose = false;
    end
    if options.symmetric
        M = (M+M')/2;
    end
    
    % Precision (due to legacy usage)
    if isempty(options.precision)
        warning([functionId ':MissingPrecision'], 'Precision was not specified. Will try to determine based on supplied type. NOTE: this feature is depricated and will be removed in future versions.');
        switch options.type
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
        fprintf(['Generating data for ' names.M '()\n']);
    end
    
    % No-ones
    if options.noones && strcmp(options.structure, 'dense')
        throw(MException([functionId ':InvalidParameters'], 'The parameter ''noones'' cannot be combined with ''dense'' structure.'));
    end
    
    %% Generate data
    if options.verbose >= 2
        fprintf('. processing data\n');
    end
    % info.structure contains data about the structure of M
    % .stored.mat - The order in which the elements of the matrix are stored - in matrix form
    if options.noones && ~strcmp(options.structure, 'dense')
        info.ones.mat = zeros(dims.n, dims.m);
        [c, r, v] = find(abs(M)' == 1); % Find -1 and 1-elements
        info.ones.num = length(v);
        info.ones.sign = sign(M((c-1)*dims.n+r));
        if ~options.transpose
            info.ones.mat((c-1)*dims.n+r) = info.ones.sign;
            info.ones.rows = r(:)';
            info.ones.cols = c(:)';
        else
            info.ones.mat((r-1)*dims.m+c) = info.ones.sign;
            info.ones.rows = c(:)';
            info.ones.cols = r(:)';
        end
        M(abs(M) == 1) = 0; % Set to zero, to be ignored
    else
        info.ones.num = 0;
        info.ones.sign = [];
        info.ones.mat = zeros(dims.n, dims.m);
        info.ones.rows = [];
        info.ones.cols = [];
    end
    switch options.structure
        case 'dense'
            c = repmat((1:dims.m)',dims.n,1); cols = c;
            r = kron((1:dims.n)',ones(dims.m,1)); rows = r;
            v = subsref(M', struct('type', '()', 'subs', {{':'}})); values = v;
            indices = 1:length(v);
        case 'sparse'
            [c, r, v] = find(M'); % Find non-zero elements, transposed to switch from column-major to row-major
            indices = 1:length(v);
            values = v;
            cols = c;
            rows = r;
        case 'unique'
            [c, r, v] = find(M'); % Find non-zero elements, transposed to switch from column-major to row-major
            [values, J, indices] = unique(v, 'stable');
            cols = c(J);
            rows = r(J);
        case 'ordered'
            [c, r, v] = find(M'); % Find non-zero elements, transposed to switch from column-major to row-major
            [values, J, indices] = unique(v, 'sorted');
            cols = c(J);
            rows = r(J);
        case 'indexed'
            [c, r, v] = find(M'); % Find non-zero elements, transposed to switch from column-major to row-major
            [values, J, I] = unique(v, 'sorted');
            indices = values(I);
            cols = c(J);
            rows = r(J);
    end
    % .stored.values
    info.structure.stored.values = values(:)';
    % .stored.rows
    info.structure.stored.rows = rows(:)';
    % .stored.cols
    info.structure.stored.cols = cols(:)';
    % .stored.num
    info.structure.stored.num = length(info.structure.stored.values);
    % .stored.mat
    info.structure.stored.mat = zeros(dims.n, dims.m);
    info.structure.stored.mat((c-1)*dims.n+r) = indices;
    
    [c, r, v] = find(info.structure.stored.mat');
    if ~options.transpose
        % .access.rows
        info.structure.access.rows = r(:)';
        % .access.cols
        info.structure.access.cols = c(:)';
        % .access.indices
        info.structure.access.indices = v(:)';
        % .access.num
        info.structure.access.num = length(v);
        % .access.mat
        info.structure.access.mat = info.structure.stored.mat;
    else
        % .access.rows
        info.structure.access.rows = c(:)'; % Inverted rows and columns to transpose
        % .access.cols
        info.structure.access.cols = r(:)';
        % .access.indices
        info.structure.access.indices = v(:)';
        % .access.num
        info.structure.access.num = length(v);
        % .access.mat
        info.structure.access.mat = info.structure.stored.mat';
    end
    
    % Bookkeeping of number of stored data points
    info.size.data = length(info.structure.stored.values);
    
    if info.structure.stored.num == 0 && (~options.noones || info.ones.num == 0)
        warning([functionId ':ZeroMatrix'], ['The selected structure is not ''dense'' and the matrix ' names.M ' has only zero-elements.']);
    end
    
    options.lineBreakThreshold = 5;
    if options.verbose >= 2
        fprintf('. generating static data\n');
    end
    data = '';
    def = '';
    if info.structure.stored.num > 0
        % Definition
        if options.static && options.const
            def = sprintf('static const ');
        elseif options.static
            def = sprintf('static ');
        elseif options.const
            def = sprintf('const ');
        end
        def = [def, sprintf([options.type ' ' names.M '[' num2str(info.structure.stored.num) ']'])];
        % Data
        data = '{';
        if strcmp(options.structure, 'ordered')
            if info.structure.stored.num > options.lineBreakThreshold % If there is more than one row with stored data
                data = [data, sprintf(['\n' options.indent.data options.indent.generic])];
            end
            data = [data, sprintf(falcopt.internal.vec2strjoin(info.structure.stored.values(1:min(options.lineBreakThreshold,info.structure.stored.num)), ', ', 'precision', options.precision))];
            for i=2:ceil(info.structure.stored.num/options.lineBreakThreshold) % Iterate over lines of stored data
                data = [data, sprintf(falcopt.internal.vec2strjoin(info.structure.stored.values((i-1)*options.lineBreakThreshold+1:min(i*options.lineBreakThreshold,info.structure.stored.num)), ', ', 'precision', options.precision))]; %#ok
            end
        else
            if sum(sum(info.structure.stored.mat,2) ~= 0,1) > 1 % If there is more than one row with stored data
                data = [data, sprintf(['\n' options.indent.data options.indent.generic])];
            end
            bFirstRow = true; % TODO: can be solved more elegantly
            for i=1:dims.n % Iterate over rows
                I = info.structure.stored.rows == i; % Stored elements of i-th row
                if any(I)
                    if bFirstRow
                        bFirstRow = false;
                    else
                        data = [data, sprintf([', \n' options.indent.data options.indent.generic])]; %#ok
                    end
                    data = [data, sprintf(falcopt.internal.vec2strjoin(info.structure.stored.values(I), ', ', 'precision', options.precision))]; %#ok
                end
            end
        end
        data = [data, sprintf(' }')];
    end
    if options.stripped
        code = def;
    else
        data = [def ' = ' data ';'];
    end
    info.size = info.structure.stored.num;
    
    %% Generate helper functions
    if ~options.stripped
        if options.verbose >= 2
            fprintf('. generating helper functions\n');
        end
        % Transformation to dense, column-major format
        if options.verbose >= 3
            fprintf(['.. generating ' names.toDense '()\n']);
        end
        % TODO: add documentation
        code.toDense = sprintf([options.indent.code 'void ' names.toDense '(const ' options.type '* ' names.M ', double* ' names.M 't) {' '\n']);
        code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'unsigned int i,j;' '\n'])];
        if info.structure.stored.num < dims.n*dims.m
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic '/* Initialize to zero */' '\n'])];
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'for(i=0; i<%i; i++) { /* Iterate over columns */' '\n'], dims.m)];
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic options.indent.generic 'for(j=0; j<%i; j++) { ' names.M 't[i*%i+j] = 0.0; }' '\n'], dims.n, dims.n)];
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic '}' '\n'])];
        end
        if options.noones && info.ones.num > 0
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic '/* Initialize -1''s and 1''s */' '\n'])];
            if ~options.transpose
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int ones_rows[%i] = {' falcopt.internal.vec2strjoin(info.ones.rows-1, ', ') '};' '\n'], info.ones.num)];
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int ones_cols[%i] = {' falcopt.internal.vec2strjoin(info.ones.cols-1, ', ') '};' '\n'], info.ones.num)];
            else % Negate transpose
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int ones_rows[%i] = {' falcopt.internal.vec2strjoin(info.ones.cols-1, ', ') '};' '\n'], info.ones.num)];
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int ones_cols[%i] = {' falcopt.internal.vec2strjoin(info.ones.rows-1, ', ') '};' '\n'], info.ones.num)];
            end
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const double ones_sign[%i] = {' falcopt.internal.vec2strjoin(info.ones.sign, ', ', 'precision', options.precision) '};' '\n'], info.ones.num)];
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'for(i=0; i<%i; i++) { /* Iterate over 1-elements */' '\n'], info.ones.num)];
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic options.indent.generic names.M 't[cols[i]*%i+rows[i]] = ones_sign[i];' '\n'], dims.n)];
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic '}' '\n'])];
        end
        if info.structure.stored.num > 0
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic '/* Copy non-zero values */' '\n'])];
            if ~options.transpose
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int rows[%i] = {' falcopt.internal.vec2strjoin(info.structure.access.rows-1, ', ') '};' '\n'], info.structure.access.num)];
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int cols[%i] = {' falcopt.internal.vec2strjoin(info.structure.access.cols-1, ', ') '};' '\n'], info.structure.access.num)];
            else % Negate transpose
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int rows[%i] = {' falcopt.internal.vec2strjoin(info.structure.access.cols-1, ', ') '};' '\n'], info.structure.access.num)];
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int cols[%i] = {' falcopt.internal.vec2strjoin(info.structure.access.rows-1, ', ') '};' '\n'], info.structure.access.num)];
            end
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'const unsigned int indices[%i] = {' falcopt.internal.vec2strjoin(info.structure.access.indices-1, ', ') '};' '\n'], info.structure.access.num)];
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic 'for(i=0; i<%i; i++) { /* Iterate over non-zero elements */' '\n'], info.structure.access.num)];
            if strcmp(options.type, 'double')
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic options.indent.generic names.M 't[cols[i]*%i+rows[i]] = ' names.M '[indices[i]];' '\n'], dims.n)];
            else
                code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic options.indent.generic names.M 't[cols[i]*%i+rows[i]] = (double)' names.M '[indices[i]];' '\n'], dims.n)];
            end
            code.toDense = [code.toDense, sprintf([options.indent.code options.indent.generic '}' '\n'])];
        end
        code.toDense = [code.toDense, sprintf([options.indent.code '}'])];
        % Transformation from dense, column-major format
        if options.verbose >= 3
            fprintf(['.. generating ' names.fromDense '()\n']);
        end
        % TODO: add documentation
        code.fromDense = '';
        if info.structure.stored.num > 0
            code.fromDense = [code.fromDense, sprintf([options.indent.code 'void ' names.fromDense '(const double* ' names.M ', ' options.type '* ' names.M 't) {' '\n'])];
            code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
            if ~options.transpose
                code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic 'const unsigned int rows[%i] = {' falcopt.internal.vec2strjoin(info.structure.stored.rows-1, ', ') '};' '\n'], info.structure.stored.num)];
                code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic 'const unsigned int cols[%i] = {' falcopt.internal.vec2strjoin(info.structure.stored.cols-1, ', ') '};' '\n'], info.structure.stored.num)];
            else % Negate transpose
                code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic 'const unsigned int rows[%i] = {' falcopt.internal.vec2strjoin(info.structure.stored.cols-1, ', ') '};' '\n'], info.structure.stored.num)];
                code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic 'const unsigned int cols[%i] = {' falcopt.internal.vec2strjoin(info.structure.stored.rows-1, ', ') '};' '\n'], info.structure.stored.num)];
            end
            code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic 'for(i=0; i<%i; i++) { /* Iterate over non-zero elements */' '\n'], info.structure.stored.num)];
            if strcmp(options.type, 'double')
                code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic options.indent.generic names.M 't[i] = ' names.M '[cols[i]*%i+rows[i]];' '\n'], dims.n)];
            else
                code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic options.indent.generic names.M 't[i] = (' options.type ')' names.M '[cols[i]*%i+rows[i]];' '\n'], dims.n)];
            end
            code.fromDense = [code.fromDense, sprintf([options.indent.code options.indent.generic '}' '\n'])];
            code.fromDense = [code.fromDense, sprintf([options.indent.code '}'])];
        end

        %% Generate additional utility functions (in MATLAB)
        if options.verbose >= 2
            fprintf('. generating utils\n');
        end
        info.utils.rand = @(s,a,b)(reshape(sum( ...
                                               (   repmat(info.structure.stored.mat(:),1,info.structure.stored.num) ...
                                                == repmat(1:info.structure.stored.num,dims.n*dims.m,1) ...
                                               ).*repmat((b-a)*rand(s,1,info.structure.stored.num)+a*ones(1,info.structure.stored.num),dims.n*dims.m,1) ...
                                            ,2), dims.n, dims.m) ...
                                    + info.ones.mat);
    end
    
    if options.verbose >= 2
        fprintf('...done\n');
    end
    if options.verbose >= 1
        fprintf('Data successfully generated\n');
    end

end

