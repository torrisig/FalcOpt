%% generateInverse Generate code for a C-function for 
% matrix inverse of the following form:
%  Mi = M^-1 (or Mi = M'^-1, if the option 'transpose' is enabled)
% 
% [code, info] = generateInverse(M)
% or
% [code, info] = generateInverse(M, ...) [with options]
% 
% Returns a string 'code' containing the code, and
%  returns an info struct containing information about the number of static
%  data points in 'size', the number of FLOPS in 'flops' and the used names in 'names' 
%
% the following options are available:
%  .names     - A struct with function and variable names to be used. Default is
%                struct('fun', 'invert', 'M', 'M', 'Mi', 'Mi', 'prefix', '')
%  .structure - A string that determines the structure of M and Mi. Default is 'sparse'.
%                only affects how the structure of M is interpreted and utilized.
%                the possibilities are:
%                'dense' for no structure,
%                'sparse' for ommitted zero elements (only non-zero elements of M are considered)
%                'unique' for ommitted zero and reccurring elements (recurring elements of M are assumed to be the same)
%                'ordered' for ommitted zero and reccurring elements including permutation of elements (elements of M are stored in ascending order of their values)
%                'indexed' for ommitted zero and recurring elements including permutation of element (elements of M are stored according to their values as indices starting from 1)
%               M and Mi can be treated differently if it is a struct e.g. struct('M', 'dense', 'Mi' 'sparse').
%  .precision - A string that determines the precision of computation and stored data in the generated code. 
%                Needs to be either 'single' or 'double'. Default: '' (unspecified, will try to extract from type);
%  .types     - A string or struct. If a string then it indicates the data type, either 'float' or 'double',
%                if a struct it also indicates the function type. Default is struct('fun', 'static inline void', 'data', 'double')
%  .symmetric - A boolean, if true then M is considered symmetric and code is generate to exploit this. Default: false.
%  .transpose - A boolean, if true then M is transposed before being inverted. Default: false.
%  .indent    - Indentation to be used in code generation
%  .inline    - Inline keyword to be uised. Default is inline
%  .verbose   - Level of procedural output of this function
%  .test      - Level of tests performed
%  .epsFactor - Factor of machine-precision that is considered. Default is 100
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
function [code, info] = generateInverse(varargin)
    defaultNames = struct('fun', 'invert', ...
                          'M', 'M', ...
                          'Mi', 'Mi', ...
                          'prefix', '');
    structures = {'dense', 'sparse', 'unique', 'ordered', 'indexed'};
    defaultStructure = struct('M', 'sparse', 'Mi', 'sparse');
    p = inputParser;
    p.CaseSensitive = true;
    p.addRequired('M', @(x)(isnumeric(x)));
    p.addParameter('names', defaultNames, @isstruct);
    p.addParameter('structure', defaultStructure, @(x)((ischar(x) && any(strcmp(x, structures))) || (isstruct(x) && all(isfield(x, {'M', 'Mi'})) && ischar(x.M) && any(strcmp(x.M, structures)) && ischar(x.Mi) && any(strcmp(x.Mi, structures)))));
    p.addParameter('precision', '', @(s)(ischar(s) && any(strcmp(s, {'single', 'double'}))));
    p.addParameter('types', 'double', @(x)(ischar(x) || (isstruct(x) && any(isfield(x, {'fun', 'data'})) && (~isfield(x, 'fun') || ischar(x.fun)) && (~isfield(x, 'data') || ischar(x.data)))));
    p.addParameter('symmetric', false, @islogical);
    p.addParameter('transpose', false, @islogical);
    p.addParameter('indent', '    ', @(x)(ischar(x) || (isstruct(x) && all(isfield(x, {'generic', 'data', 'code'})) && ischar(x.generic) && ischar(x.data) && ischar(x.code))));
    p.addParameter('inline', 'inline', @ischar);
    p.addParameter('verbose', 0, @(x)(isnumeric(x) && x >= 0));
    p.addParameter('test', 0, @(x)(isnumeric(x) && x >= 0));
    p.addParameter('epsFactor', 100, @(x)(isnumeric(x) && x > 0));
    p.parse(varargin{:});
    options = p.Results;
    M = options.M;
    names = options.names;
    
    if options.verbose >= 2
        fprintf('Code generation started...\n');
        fprintf('. processing options\n');
    end
    
    if any(options.symmetric)
        warning('falcopt:MissingImplemementation:SymmetricMatrix', 'The special treatment of symmetric matrices has not yet been implemented. Ignoring...');
    end
    
    %% Check dimension of M
    dims.n = size(M,1);
    if dims.n ~= size(M,2)
        throw(MException('falcopt:generateInverse:InvalidDimension', 'The matrix needs to be square.'));
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
    if ~isstruct(options.structure)
        options.structure = struct('M', options.structure, 'Mi', options.structure);
    end
        
    % Ensure and exploit symmetry
    if options.symmetric
        % Set transpose flag to zero if also symmetric
        if options.symmetric && options.transpose
            options.transpose = false;
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
        options.types.fun = ['static ' options.inline ' void'];
    end
    if strcmp(options.types.data, 'single')
        warning('falcopt:generateInverse:InvalidType', 'The type "single" is depricated and will not be allowed in a future version.');
        options.types.data = 'float';
    end
    % Precision (due to legacy usage)
    if isempty(options.precision)
        warning('falcopt:generateInverse:MissingPrecision', 'Precision was not specified. Will try to determine based on supplied type. NOTE: this feature is depricated and will be removed in future versions.');
        switch options.types.data
            case 'float'
                options.precision = 'single';
            case 'single'
                options.precision = 'single';
            case 'double'
                options.precision = 'double';
            otherwise
                warning('falcopt:generateInverse:InvalidPrecision', ['Precision was not specified, could not infer proper precision from given type "' options.types.data '". Setting to "double".']);
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
    
    %% Pre-process data
    if options.verbose >= 2
        fprintf('. processing structure\n');
    end
    structureInfo = falcopt.detectMatrixStructure(M, 'structure', options.structure.M, 'symmetric', options.symmetric, 'verbose', options.verbose+1);
    info.structure = structureInfo.structure;
    elements = structureInfo.elements;
    % Indices of solution matrix Mi
    Mi = zeros(dims.n);
    for i=1:elements.M.blocks.num
        Mi((elements.M.col(elements.M.blocks.indices{i})-1)*dims.n+elements.M.row(elements.M.blocks.indices{i})) = i;
        if elements.M.blocks.size(i) == 2
            N = Mi(unique(elements.M.row(elements.M.blocks.indices{i}),'stable'),unique(elements.M.col(elements.M.blocks.indices{i}),'stable'));
            Mi(unique(elements.M.row(elements.M.blocks.indices{i}),'stable'),unique(elements.M.col(elements.M.blocks.indices{i}),'stable')) = [N(4) N(3); N(2) N(1)];
        end
    end
    Mi = Mi';
    Mi(Mi~=0) = 1:sum(sum(Mi~=0));
    elements.Mi.indices = Mi';
    [c, r, ~] = find(elements.Mi.indices'); % Find non-zero elements, transposed to switch from column-major to row-major
    elements.Mi.row = r;
    elements.Mi.col = c;
    elements.Mi.num = length(elements.Mi.row);
    switch options.structure.Mi
        case 'dense'
        case 'sparse'
            elements.Mi.data.row = elements.Mi.row;
            elements.Mi.data.col = elements.Mi.col;
            elements.Mi.data.num = elements.Mi.num;
            elements.Mi.data.indices = elements.Mi.indices;
        case 'unique'
            throw(MException('falcopt:generateInverse:MissingImplementation', 'Structure ''unique'' for Mi has not yet been implemented. Use ''sparse''.'));
        case 'ordered'
            throw(MException('falcopt:generateInverse:MissingImplementation', 'Structure ''ordered'' for Mi has not yet been implemented. Use ''sparse'''.'));
        case 'indexed'
            throw(MException('falcopt:generateInverse:MissingImplementation', 'Structure ''indexed'' for Mi has not yet been implemented. Use ''sparse'''.'));
    end
    % TODO: Check whether there are singular blocks

    info.elements = elements;
    % Store structure of Mi in info
    info.structure.Mi.type = options.structure.Mi;
    info.structure.Mi.mat = elements.Mi.indices;
    info.structure.Mi.num = sum(sum(elements.Mi.indices~=0));
    if strcmp(options.structure.Mi, 'dense') || strcmp(options.structure.Mi, 'sparse')
        % Replace all non-zero elements with 1
        info.structure.Mi.mat(info.structure.Mi.mat ~= 0) = 1;
    elseif strcmp(options.structure.Mi, 'ordered')
        % TODO?
    elseif strcmp(options.structure.Mi, 'indexed')
        % TODO
    end
    
    %% Generate code
    if options.verbose >= 2
        fprintf('. generating code\n');
    end
    info.flops.add = 0;
    info.flops.mul = 0;
    info.flops.inv = 0;
    % Documentation
    StructStr = mat2str(info.structure.M.mat);
    StructStr = regexprep(StructStr, '([; \[ ])(0)([; \] ])', '$1*$3');
    StructStr = regexprep(StructStr, '([; \[ ])(0)([; \] ])', '$1*$3');
    code = sprintf([options.indent.code '/** ' '\n' ...
                    options.indent.code ' * @brief Inverts the matrix ' names.M ' and stores it in ' names.Mi '.' '\n' ...
                    options.indent.code ' * @param ' names.M ' a %ix%i matrix with %i blocks, %i non-zero elements and the following structure:' '\n' ...
                    options.indent.code ' *         ' names.M ' = ' sprintf(strrep(StructStr, ';', ['\n' options.indent.code ' *          ' repmat(' ',1,length([names.M ' = ']))])) ';' '\n' ...
                    options.indent.code ' */' '\n'], dims.n, dims.n, elements.M.blocks.num, info.structure.M.num);
    % Function header
    code = [code, sprintf([options.indent.code options.types.fun ' ' names.fun '(' options.types.data '* ' names.Mi ', const ' options.types.data '* ' names.M ') {' '\n'])];
    if any(elements.M.blocks.size > 1)
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' det_i; /* Temporary variable to store inverse of determinant */' '\n'])];
    end
    if any(elements.M.blocks.size == 2)
        code = [code, sprintf([options.indent.code options.indent.generic options.types.data ' temp; /* Temporary variable */' '\n'])];
    end
    if any(elements.M.blocks.size == 4)
        code = [code, sprintf([options.indent.code options.indent.generic 'unsigned int i;' '\n'])];
    end
    % Code (iterate over blocks)
    for i=1:elements.M.blocks.num
        code = [code, sprintf([options.indent.code options.indent.generic '/* Block #%i (%ix%i) */' '\n'], i, elements.M.blocks.size(i), elements.M.blocks.size(i))]; %#ok
        warning('falcopt:generateInverse:InvalidData', 'Indexing of Mi is not correct, please verify.');
        switch elements.M.blocks.size(i)
            case 1
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = 1.0/' names.M '[%i];' '\n'], elements.Mi.indices(elements.M.row(elements.M.blocks.indices{i}),elements.M.col(elements.M.blocks.indices{i}))-1, elements.M.blocks.indices{i}(1)-1)]; %#ok
                info.flops.inv = info.flops.inv+1;
            case 2
                [~,~,Jr] = unique(elements.M.row(elements.M.blocks.indices{i}), 'stable');
                [~,~,Jc] = unique(elements.M.col(elements.M.blocks.indices{i}), 'stable');
                I = (Jr-1)*elements.M.blocks.size(i)+Jc;
                % Compute determinant
                if length(I) == 4
                    code = [code, sprintf([options.indent.code options.indent.generic 'det_i  = 1.0/(' names.M '[%i]*' names.M '[%i] - ' names.M '[%i]*' names.M '[%i]); /* Compute determinant */' '\n'], elements.M.blocks.indices{i}(I==1)-1, elements.M.blocks.indices{i}(I==4)-1, elements.M.blocks.indices{i}(I==2)-1, elements.M.blocks.indices{i}(I==3)-1)]; %#ok
                    info.flops.inv = info.flops.inv+1;
                    info.flops.mul = info.flops.mul+2;
                    info.flops.add = info.flops.add+1;
                elseif any(I==1) && any(I==4)
                    code = [code, sprintf([options.indent.code options.indent.generic 'det_i  = 1.0/(' names.M '[%i]*' names.M '[%i]); /* Compute determinant */' '\n'], elements.M.blocks.indices{i}(I==1)-1, elements.M.blocks.indices{i}(I==4)-1)]; %#ok
                    info.flops.inv = info.flops.inv+1;
                    info.flops.mul = info.flops.mul+1;
                elseif any(I==2) && any(I==3)
                    code = [code, sprintf([options.indent.code options.indent.generic 'det_i  = -1.0/(' names.M '[%i]*' names.M '[%i]); /* Compute determinant */' '\n'], elements.M.blocks.indices{i}(I==2)-1, elements.M.blocks.indices{i}(I==3)-1)]; %#ok
                    info.flops.inv = info.flops.inv+1;
                    info.flops.mul = info.flops.mul+1;
                else
                    throw(MException('falcopt:generateInverse:InvalidData', 'Matrix is always singular.'));
                end
                % Compute entries of inverse
                idx = 1;
                cols = unique(elements.M.col(elements.M.blocks.indices{i}), 'stable');
                rows = unique(elements.M.row(elements.M.blocks.indices{i}), 'stable');
                if any(I==4) % First element
                    if any(I==1) % Store first element in temporary variable
                        code = [code, sprintf([options.indent.code options.indent.generic 'temp = ' names.M '[%i];' '\n'], elements.M.blocks.indices{i}(I==1)-1)];
                    end
                    code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = det_i*' names.M '[%i];' '\n'], elements.Mi.indices(rows(1), cols(1))-1, elements.M.blocks.indices{i}(I==4)-1)];
                    info.flops.mul = info.flops.mul+1;
                else
                    %code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = 0.0;' '\n'], elements.M.indices(rows(1), cols(1))-1)];
                end
                if any(I==2) % Second element
                    code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = -det_i*' names.M '[%i];' '\n'], elements.Mi.indices(rows(1), cols(2))-1, elements.M.blocks.indices{i}(I==2)-1)];
                    info.flops.mul = info.flops.mul+1;
                else
                    %code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = 0.0;' '\n'], elements.M.indices(rows(1), cols(2))-1)];
                end
                if any(I==3) % Third element
                    code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = -det_i*' names.M '[%i];' '\n'], elements.Mi.indices(rows(2), cols(1))-1, elements.M.blocks.indices{i}(I==3)-1)];
                    info.flops.mul = info.flops.mul+1;
                else
                    %code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = 0.0;' '\n'], elements.M.indices(rows(2), cols(1))-1)];
                end
                if any(I==1) % Last element
                    if any(I==4) % Use stored element
                        code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = det_i*temp;' '\n'], elements.Mi.indices(rows(2), cols(2))-1)];
                        info.flops.mul = info.flops.mul+1;
                    else
                        code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = det_i*' names.M '[%i];' '\n'], elements.Mi.indices(rows(2), cols(2))-1, elements.M.blocks.indices{i}(I==1)-1)];
                        info.flops.mul = info.flops.mul+1;
                    end
                else
                    %code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[%i] = 0.0;' '\n'], elements.M.indices(rows(2), cols(2))-1)];
                end
            case 4
                [~,~,Jr] = unique(elements.M.row(elements.M.blocks.indices{i}), 'stable');
                [~,~,Jc] = unique(elements.M.col(elements.M.blocks.indices{i}), 'stable');
                I = (Jr-1)*elements.M.blocks.size(i)+Jc;
                % Compute entries of inverse
                cols = unique(elements.M.col(elements.M.blocks.indices{i}), 'stable');
                rows = unique(elements.M.row(elements.M.blocks.indices{i}), 'stable');
                % TODO deal with non-dense inverse and non-dense M!
                mult2str = @(idx)([names.M '[' num2str(elements.M.blocks.indices{i}(I==idx)-1) ']']);
                % (1,1)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = ' strjoin(arrayfun(mult2str, [6,11,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [6,12,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [10,7,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [10,8,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [14,7,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [14,8,11], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (1,2)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(2))-1) '] = -' strjoin(arrayfun(mult2str, [2,11,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [2,12,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [10,3,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [10,4,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [14,3,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [14,4,11], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (1,3)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(3))-1) '] = ' strjoin(arrayfun(mult2str, [2,7,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [2,8,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [6,3,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [6,4,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [14,3,8], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [14,4,7], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (1,4)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(4))-1) '] = -' strjoin(arrayfun(mult2str, [2,7,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [2,8,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [6,3,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [6,4,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [10,3,8], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [10,4,7], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (2,1)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(2), cols(1))-1) '] = -' strjoin(arrayfun(mult2str, [5,11,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,12,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,7,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,8,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,7,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,8,11], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (2,2)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(2), cols(2))-1) '] = ' strjoin(arrayfun(mult2str, [1,11,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [1,12,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,3,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,4,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,3,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,4,11], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (2,3)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(2), cols(3))-1) '] = -' strjoin(arrayfun(mult2str, [1,7,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [1,8,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,3,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [5,4,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,3,8], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,4,7], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (2,4)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(2), cols(4))-1) '] = ' strjoin(arrayfun(mult2str, [1,7,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [1,8,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [5,3,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,4,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,3,8], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,4,7], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (3,1)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(3), cols(1))-1) '] = ' strjoin(arrayfun(mult2str, [5,10,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [5,12,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,6,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,8,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,6,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,8,10], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (3,2)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(3), cols(2))-1) '] = -' strjoin(arrayfun(mult2str, [1,10,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [1,12,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,2,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,4,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,2,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,4,10], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (3,3)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(3), cols(3))-1) '] = ' strjoin(arrayfun(mult2str, [1,6,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [1,8,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [5,2,16], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,4,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,2,8], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,4,6], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (3,4)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(3), cols(4))-1) '] = -' strjoin(arrayfun(mult2str, [1,6,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [1,8,10], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,2,12], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [5,4,10], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,2,8], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,4,6], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (4,1)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(4), cols(1))-1) '] = -' strjoin(arrayfun(mult2str, [5,10,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,11,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,6,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,7,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,6,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,7,10], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (4,2)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(4), cols(2))-1) '] = ' strjoin(arrayfun(mult2str, [1,10,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [1,11,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,2,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,3,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,2,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,3,10], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (4,3)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(4), cols(3))-1) '] = -' strjoin(arrayfun(mult2str, [1,6,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [1,7,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,2,15], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [5,3,14], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [13,2,7], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [13,3,6], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % (4,2)
                code = [code, sprintf([options.indent.code options.indent.generic names.Mi '[' num2str(elements.Mi.indices(rows(4), cols(4))-1) '] = ' strjoin(arrayfun(mult2str, [1,6,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [1,7,10], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [5,2,11], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [5,3,10], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '+' strjoin(arrayfun(mult2str, [9,2,7], 'UniformOutput', false),'*') '\n'])];
                code = [code, sprintf([options.indent.code options.indent.generic repmat(' ',1,length([names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] = '])) '-' strjoin(arrayfun(mult2str, [9,3,6], 'UniformOutput', false),'*') ';' '\n'])];
                info.flops.mul = info.flops.mul+12;
                info.flops.add = info.flops.add+5;
                % Compute inverse determinant
                code = [code, sprintf([options.indent.code options.indent.generic 'det_i = 1.0/(' names.M '[' num2str(elements.M.blocks.indices{i}(I==1)-1) ']*' names.Mi '[' num2str(elements.Mi.indices(rows(1), cols(1))-1) '] + ' names.M '[' num2str(elements.M.blocks.indices{i}(I==2)-1) ']*' names.Mi '[' num2str(elements.Mi.indices(rows(2), cols(1))-1) '] + ' names.M '[' num2str(elements.M.blocks.indices{i}(I==3)-1) ']*' names.Mi '[' num2str(elements.Mi.indices(rows(3), cols(1))-1) '] + ' names.M '[' num2str(elements.M.blocks.indices{i}(I==4)-1) ']*' names.Mi '[' num2str(elements.Mi.indices(rows(4), cols(1))-1) ']);' '\n'])];
                info.flops.mul = info.flops.mul+4;
                info.flops.add = info.flops.add+3;
                info.flops.inv = info.flops.inv+1;
                code = [code, sprintf([options.indent.code options.indent.generic 'for(i=0; i<16; i++) { ' names.Mi '[i] = ' names.Mi '[i]*det_i; }' '\n'])];
            otherwise
                throw(MException('falcopt:generateInverse:MissingImplementation', ['There is a block of size ' num2str(max(elements.M.blocks.size)) 'x' num2str(max(elements.M.blocks.size)) '. generateInverse() has only been implemented for blocks of size 1x1, 2x2 and 4x4.']));
        end
    end
    code = [code, sprintf([options.indent.code '}'])];
    
    if options.test > 0
        if options.verbose >= 2
            fprintf('. generating MEX file\n');
        end
        %% Generate MEX file
        filename = [tempname '.c'];
        f = fopen(filename, 'w+');
        fprintf(f, ['/* Temporary MEX file for testing ' names.fun '() */' '\n' ...
                    '\n' ...
                    '#include "mex.h"' '\n' '\n']);
        % Code generated function
        fprintf(f, ['\n' '/***************************' '\n' ...
                         ' * Code-generated Function *' '\n' ...
                         ' ***************************/' '\n']);
        fprintf(f, [code '\n']);
        % Auxiliary functions
        fprintf(f, ['\n' '/***********************' '\n' ...
                         ' * Auxiliary Functions *' '\n' ...
                         ' ***********************/' '\n']);
        fprintf(f, ['\n'  options.types.fun ' transform_' names.M '(const double* M, ' options.types.data '* Mt) {' '\n']);
        if strcmp(options.structure.M, 'dense')
            % TODO
%             fprintf(f, [options.indent.generic 'const unsigned int rows[%i] = {' falcopt.internal.vec2strjoin(elements.M.data.row-1, ', ') '};' '\n'], elements.M.data.num);
%             fprintf(f, [options.indent.generic 'const unsigned int cols[%i] = {' falcopt.internal.vec2strjoin(elements.M.data.col-1, ', ') '};' '\n'], elements.M.data.num);
%             fprintf(f, [options.indent.generic 'const unsigned int indices[%i] = {' falcopt.internal.vec2strjoin(elements.M.data.indices-1, ', ') '};' '\n'], elements.M.data.num);
%             fprintf(f, [options.indent.generic 'unsigned int i;' '\n' ...
%                         options.indent.generic 'for(i=0; i<%i; i++) {' '\n' ...
%                         options.indent.generic options.indent.generic 'Mt[indices[i]] = (' options.types.data ')(M[cols[i]*%i+rows[i]]);' '\n' ...
%                         options.indent.generic '}' '\n'], elements.M.data.num(i), dims.n);
        else
            fprintf(f, [options.indent.generic 'const unsigned int rows[%i] = {' falcopt.internal.vec2strjoin(elements.M.data.row-1, ', ') '};' '\n'], elements.M.data.num);
            fprintf(f, [options.indent.generic 'const unsigned int cols[%i] = {' falcopt.internal.vec2strjoin(elements.M.data.col-1, ', ') '};' '\n'], elements.M.data.num);
            fprintf(f, [options.indent.generic 'unsigned int i;' '\n' ...
                        options.indent.generic 'for(i=0; i<%i; i++) {' '\n' ...
                        options.indent.generic options.indent.generic 'Mt[i] = (' options.types.data ')(M[cols[i]*%i+rows[i]]);' '\n' ...
                        options.indent.generic '}' '\n'], elements.M.data.num, dims.n);
        end
        fprintf(f, ['}' '\n']);
        fprintf(f, ['\n'  options.types.fun ' transform_' names.Mi '_toOutput(const ' options.types.data '* M, double* Mt) {' '\n']);
        if strcmp(options.structure.Mi, 'dense')
            % TODO
        else
            fprintf(f, [options.indent.generic 'const unsigned int rows[%i] = {' falcopt.internal.vec2strjoin(elements.Mi.data.row-1, ', ') '};' '\n'], elements.Mi.data.num);
            fprintf(f, [options.indent.generic 'const unsigned int cols[%i] = {' falcopt.internal.vec2strjoin(elements.Mi.data.col-1, ', ') '};' '\n'], elements.Mi.data.num);
            fprintf(f, [options.indent.generic 'unsigned int i;' '\n' ...
                        options.indent.generic 'for(i=0; i<%i; i++) {' '\n' ...
                        options.indent.generic options.indent.generic 'Mt[cols[i]*%i+rows[i]] = (' options.types.data ')(M[i]);' '\n' ...
                        options.indent.generic '}' '\n'], elements.Mi.data.num, dims.n);
        end
        fprintf(f, ['}' '\n']);
        fprintf(f, ['\n' '/****************' '\n' ...
                         ' * MEX Function *' '\n' ...
                         ' ****************/' '\n']);
        fprintf(f, ['\n' 'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {' '\n']);
        fprintf(f, [options.indent.generic '/* Processing in- and outputs */' '\n']);
        fprintf(f, [options.indent.generic 'if(nrhs != 1) { mexErrMsgIdAndTxt("' names.fun ':InvalidInput", "Requires 1 input."); }' '\n']);
        % Return matrix Mi
        fprintf(f, [options.indent.generic 'if(nlhs != 1) { mexErrMsgIdAndTxt("' names.fun ':InvalidInput", "Requires 1 output."); }' '\n']);
        fprintf(f, [options.indent.generic 'plhs[0] = mxCreateDoubleMatrix(%i,%i,mxREAL); /* Output matrix Mi */' '\n'], dims.n, dims.n);
        fprintf(f, [options.indent.generic options.types.data ' ' names.Mi '[%i];' '\n'], elements.Mi.data.num);
        % Matrix M
        fprintf(f, [options.indent.generic 'if(!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetM(prhs[0]) != %i || mxGetN(prhs[0]) != %i) { mexErrMsgIdAndTxt("' names.fun ':InvalidDimension", "Matrix ' names.M ' has invalid dimension."); }' '\n'], dims.n, dims.n);
        if strcmp(options.structure.M, 'dense')
            fprintf(f, [options.indent.generic options.types.data ' ' names.M '[%i]; '], dims.n^2);
        else
            fprintf(f, [options.indent.generic options.types.data ' ' names.M '[%i]; '], elements.M.data.num);
        end
        fprintf(f, ['double* ' names.M '_in = mxGetPr(prhs[0]);' '\n']);
        fprintf(f, [options.indent.generic 'transform_' names.M '(' names.M '_in, ' names.M ');' '\n']);
        % Call function
        fprintf(f, [options.indent.generic '/* Call function */' '\n']);
        fprintf(f, [options.indent.generic names.fun '(' names.Mi ', ' names .M ');\n']);
        % Outputs
        fprintf(f, [options.indent.generic '/* Prepare outputs */' '\n']);
        fprintf(f, [options.indent.generic 'transform_' names.Mi '_toOutput(' names.Mi ', mxGetPr(plhs[0]));' '\n']);
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
        for i=1:floor(options.test)
            % Generate data
            values = randn(elements.M.data.num);
            Ms = zeros(dims.n);
            if options.transpose
                if strcmp(options.structure.M, 'dense')
                    % TODO
                else
                    % TODO
                end
            else
                if strcmp(options.structure.M, 'dense')
                    Ms(elements.M.data.indices) = values;
                else
                    Ms((elements.M.data.col-1)*dims.n+elements.M.data.row) = values(elements.M.data.indices);
                end
            end
            % Compute Mi using code-generated function
            eval(['Mis = ' names.fun '_mex(Ms);']);
            % Compute correct value
            Mi = inv(Ms);
            if isnan(norm(Mi-Mis,inf)) || norm(Mi-Mis,inf) > options.epsFactor*eps(options.precision)
                errors = errors + 1;
                maxError = max(maxError, norm(Mi-Mis,inf));
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
        if options.test == 0 || errors == 0
            fprintf('Code successfully generated\n');
        else
            fprintf('Code generated, testing revealed %i%% errors, the maximum error is %1.5g\n', floor(100*errors/floor(options.test)), maxError);
        end
    end
    if options.test > 0 && errors == options.test
        throw(MException('falcopt:generatedInverse:FaultyImplementation', 'All tests failed. Likely bug in implementation.'));
    end

end

