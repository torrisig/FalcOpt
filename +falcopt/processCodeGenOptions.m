%% processCodeGenOptions
%

% Copyright (c) 2018, ETH Zurich, Automatic Control Laboratory 
%                    Damian Frick <falcopt@damianfrick.com>
%                    Giampaolo Torrisi <giampaolo.torrisi@gmail.com>
%                    Tommaso Robbiani <tommasro@student.ethz.ch>
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
function options = processCodeGenOptions(options, defaultNames, indentTypes)
    functionId = 'generateCode';
    
    %% Check dynamics
    if options.nw == 0 && nargin(options.dynamics) ~= 2
        throw(MException([functionId ':InvalidDynamics'], 'The dynamics function must have two inputs (x,u).'));
    end
    if options.nw > 0 && nargin(options.dynamics) ~= 3
        throw(MException([functionId ':InvalidDynamics'], 'The dynamics function must have three inputs (x,u,w), since nw > 0.'));
    end
    
    %% Check objective
    % Objective matrices for quadratic objective
    if isfield(options.objective,'Q')
        if size(options.objective.Q,1) ~= options.nx || size(options.objective.Q,2) ~= options.nx
            throw(MException([functionId ':InvalidObjective'], 'The quadratic state stage objective matrix objective.Q must be of dimension nx times nx.'));
        end
    else
        options.objective.Q = [];
    end
    if isfield(options.objective,'R')
        if size(options.objective.R,1) ~= options.nu || size(options.objective.R,2) ~= options.nu
            throw(MException([functionId ':InvalidObjective'], 'The quadratic input stage objective matrix objective.R must be of dimension nu times nu.'));
        end
    else
        options.objective.R = [];
    end
    if isfield(options.objective,'P')
        if size(options.objective.P,1) ~= options.nx || size(options.objective.P,2) ~= options.nx
            throw(MException([functionId ':InvalidObjective'], 'The quadratic final state objective matrix objective.P must be of dimension nx times nx.'));
        end
    else
        options.objective.P = [];
    end
    % Nonlinear objective function
    if isfield(options.objective, 'nonlinear')
        if ~isa(options.objective.nonlinear, 'function_handle')
            throw(MException([functionId ':InvalidObjective'], 'The nonlinear objective function objective.nonlinear must be a function handle.'));
        elseif nargin(options.objective.nonlinear) ~= 2
            throw(MException([functionId ':InvalidObjective'], 'The nonlinear objective function objective.nonlinear must have two inputs (x,u).'));
        end
    end
    if isfield(options.objective, 'nonlinearN')
        if ~isa(options.objective.nonlinear, 'function_handle')
            throw(MException([functionId ':InvalidObjective'], 'The nonlinear objective function objective.nonlinearN must be a function handle.'));
        elseif nargin(options.objective.nonlinear) ~= 1
            throw(MException([functionId ':InvalidObjective'], 'The nonlinear objective function objective.nonlinearN must have one input (x).'));
        end
    end
    % Reference tracking
    if isfield(options.objective,'trackReference')
        if ~islogical(options.objective.trackReference) || numel(options.objective.trackReference) ~= 1
            throw(MException([functionId ':InvalidObjective'], 'The objective option objective.trackReference must be either true or false).'));
        end
    else
        options.objective.trackReference = false;
    end
    
    %% Check constraints
    % Take care of legacy option names
    if ~isempty(options.box_lowerBound)
        if isempty(options.lb)
            options.lb = options.box_lowerBound;
            warning([functionId ':Depricated'], 'The option ''box_lowerBound'' has ben renamed to ''lb'' and is depricated and will be removed in future versions. Please switch to using ''lb''.');
        else
            warning([functionId ':InconsistentOption'], 'The option ''box_lowerBound'' has ben renamed to ''lb'' and is depricated, however both options are used. Will be using ''lb''.');
        end
    end
    options = rmfield(options, 'box_lowerBound');
    if ~isempty(options.box_upperBound)
        if isempty(options.ub)
            options.ub = options.box_upperBound;
            warning([functionId ':Depricated'], 'The option ''box_upperBound'' has ben renamed to ''ub'' and is depricated and will be removed in future versions. Please switch to using ''ub''.');
        else
            warning([functionId ':InconsistentOption'], 'The option ''box_upperBound'' has ben renamed to ''ub'' and is depricated, however both options are used. Will be using ''ub''.');
        end
    end
    options = rmfield(options, 'box_upperBound');
    if ~isempty(options.constraints_handle)
        if isempty(options.constraints)
            options.constraints = options.constraints_handle;
            warning([functionId ':Depricated'], 'The option ''constraints_handle'' has ben renamed to ''constraints'' and is depricated and will be removed in future versions. Please switch to using ''constraints''.');
        else
            warning([functionId ':InconsistentOption'], 'The option ''constraints_handle'' has ben renamed to ''constraints'' and is depricated, however both options are used. Will be using ''constraints''.');
        end
    end
    options = rmfield(options, 'constraints_handle');
    % lower bounds on input
    if ~isempty(options.lb)
        if all(size(options.lb) == [1,options.nu])
            options.lb = options.lb'; % Make sure it is a column vector
        elseif any(size(options.lb) ~= [options.nu,1]) && any(size(options.lb) ~= [options.nu,options.N])
            throw(MException([functionId ':InvalidBounds'], 'lb must be either empty (no lower bounds), of dimension [nu,1] (time invariant lower bounds) or [nu, N] (stage dependent lower bounds).'));
        end
        % Transform vector to matrix
        if numel(options.lb) == options.nu
            options.lb = repmat(options.lb,1,options.N);
        end
    else
        options.lb = -Inf(options.nu,options.N);
    end
    % upper bounds on input
    if ~isempty(options.ub)
        if all(size(options.ub) == [1,options.nu])
            options.ub = options.ub'; % Make sure it is a column vector
        elseif any(size(options.ub) ~= [options.nu,1]) && any(size(options.ub) ~= [options.nu,options.N])
            throw(MException([functionId ':InvalidBounds'], 'ub must be either empty (no upper bounds), of dimension [nu,1] (time invariant upper bounds) or [nu, N] (stage dependent upper bounds).'));
        end
        % Transform vector to matrix
        if numel(options.ub) == options.nu
            options.ub = repmat(options.ub,1,options.N);
        end
    else
        options.ub = Inf(options.nu,options.N);
    end
    % check feasibility of bounds
    for k=1:options.N
        for i=1:options.nu
            if options.lb(i,k) > options.ub(i,k)
                throw(MException([functionId ':InfeasibleBounds'], ['Bound must be feasible, lb > ub at stage ' num2str(k) '.']));
            elseif options.lb(i,k) == options.ub(i,k)
                throw(MException([functionId ':MissingImplementation'], 'The case of lb = ub is not implemented yet. Please make sure that lb < ub.'));
            end
        end
    end
    % nonlinear input stage constraint
    if ~isempty(options.constraints)
        if ~iscell(options.constraints)
            options.constraints = {options.constraints};
        end
        if all(size(options.constraints) == [1,1])
            options.constraints = repmat(options.constraints,1,options.N);
        end    
        if ~all(size(options.constraints) == [1,options.N])
            if all(size(options.constraints) == [options.N,1])
                options.constraints = options.constraints';
            else
                throw(MException([functionId ':InvalidConstraints'], 'The constraints must be either empty, a function handle, or a cell of function handles of length N.'));
            end
        end
        for k=1:options.N
            if ~isa(options.constraints{k}, 'function_handle')
                throw(MException([functionId ':InvalidConstraints'], 'The constraints must be either empty, a function handle, or a cell of function handles of length N.'));
            end
        end
    end
    % nonlinear input stage constraint dimension
    if isempty(options.nn)
        options.nn = 0;
    end
    if iscell(options.nn)
        options.nn = cell2mat(options.nn);
        warning([functionId ':Depricated'], 'The constraint dimension nn is a cell. This is depricated and will be removed in future versions. Please make it a vector.');
    end
    if all(size(options.nn) == [1,1])
        options.nn = options.nn*ones(1,options.N);
    end
    if ~all(size(options.nn) == [1,options.N])
        if all(size(options.nn) == [options.N,1])
            options.nn = options.nn';
        else
            throw(MException([functionId ':InvalidConstraintDimensions'], 'The constraint dimensions must be either empty, a non-negative integer, or a vector of non-negative integers of length N.'));
        end
    end
    if any(options.nn < 0 | mod(options.nn,1) ~= 0)
        throw(MException([functionId ':InvalidConstraintDimensions'], 'The constraint dimensions must be either empty, a non-negative integer, or a vector of non-negative integers of length N.'));
    end
    if any(options.nn == 0 & ~isempty(options.constraints))
        warning([functionId ':InconsistentConstraint'], 'There are non-empty constraints with zero constraint dimension. This may be an error, please check!');
    end
    if ~isempty(options.constraints)
        for k=1:options.N
            if any(size(options.constraints{k}(ones(options.nu,1))) ~= [options.nn(k),1])
                throw(MException([functionId ':InvalidConstraints'], 'The constraints must return a column vector of length nn(k).'));
            end
        end
    end
    % Terminal and contractive are exclusive
    if (options.terminal && options.contractive)
        throw(MException([functionId ':InvalidConstraint'], 'Cannot have both terminal and contractive constraint.'));
    end
    % Check whether there are at least some constraints
    if all(all(isinf(options.lb))) && all(all(isinf(options.ub))) && isempty(options.constraints) && ~options.contractive && ~options.terminal
        throw(MException([functionId ':MissingImplementation'], 'Completely unconstrained problems are not implemented yet. Please contact us if you need this feature.'));
    end

    %% Check code generation options
    % Name TODO: ensure name leads to valid filenames!
    % Check names
    fields = fieldnames(defaultNames);
    for f=1:length(fields)
        % Set defaults
        if ~isfield(options.names, fields{f})
            options.names.(fields{f}) = defaultNames.(fields{f});
        end
    end
    if ~any(strcmp('prefix', fields))
        options.names.prefix = [options.name '_'];
    end
    for f=1:length(fields)
        % Add prefix
        if ~strcmp(fields{1}, 'prefix')
            options.names.(fields{f}) = [options.names.prefix options.names.(fields{f})];
        end
    end
    options.names.fun = options.name;
    % Target director for code generation
    if options.gendir(end) == '/' || options.gendir(end) == '\'
        options.gendir = options.gendir(1:end-1);
    end
    if ~isdir(options.gendir)
        [success, ~, ~] = mkdir(options.gendir);
        if ~success
            throw(MException([functionId ':InvalidDirectory'], ['The target directory ''gendir'': ' options.gendir ' for code generation is invalid/does not exist and could not be created.']));
        end
    end
    % Types
    switch options.precision
        case 'double'
            options.real = 'double';
            options.auxFunctions = struct('sqrt', 'sqrt', ...
                                          'max', 'my_fmax', ...
                                          'min', 'my_fmin', ...
                                          'abs', 'fabs');
        case 'single'
            options.real = 'float';
            options.auxFunctions = struct('sqrt', 'sqrtf', ...
                                          'max', 'my_fmaxf', ...
                                          'min', 'my_fminf', ...
                                          'abs', 'fabsf');
    end
    % Gradients
    switch options.gradients
        case 'matlab'
            if isempty(ver('symbolic'))
                throw(MException([functionId ':MissingToolbox'], '''gradients'' using ''matlab'' requires the Symbolic Math Toolbox. Install it to continue'));
            end
        case 'casadi'
            usedFunctions = {which('casadi.SX.sym'), ...
                             which('casadi.Function'), ...
                             which('is_constant'), ...
                             which('sparsify'), ...
                             which('densify'), ...
                             which('is_equal'), ...
                             which('to_double')};
            if(any(cellfun(@isempty,usedFunctions)))
                throw(MException([functionId ':MissingToolbox'], '''gradients'' using ''casadi'' requires CasADi. Check the CasADi version or install it. Tested with CasADi v3.1.1'));
            end
        case 'manual'
            if isempty(ver('symbolic'))
                throw(MException([functionId ':MissingToolbox'], '''gradients'' using ''manual'' requires the Symbolic Math Toolbox. Install it to continue'));
            end
            % Take care of legacy option names
            if ~isfield(options.external, 'manual')
                if ~isempty(options.external_jacobian_x)
                    options.external.manual.jacobian_x = options.external_jacobian_x;
                    warning([functionId ':Depricated'], '''gradients'' using ''manual'' now uses the ''external'' option. The use of ''external_jacobian_x'' is depricated and will be removed in future versions.');
                end
                if ~isempty(options.external_jacobian_u)
                    options.external.manual.jacobian_u = options.external_jacobian_u;
                    warning([functionId ':Depricated'], '''gradients'' using ''manual'' now uses the ''external'' option. The use of ''external_jacobian_u'' is depricated and will be removed in future versions.');
                end
                if ~isempty(options.external_jacobian_n)
                    options.external.manual.jacobian_n = options.external_jacobian_n;
                    warning([functionId ':Depricated'], '''gradients'' using ''manual'' now uses the ''external'' option. The use of ''external_jacobian_n'' is depricated and will be removed in future versions.');
                end
            end
            if ~isfield(options.external, 'manual') || isempty(options.external.manual.jacobian_x) || isempty(options.external.manual.jacobian_u)
                throw(MException([functionId ':MissingParameter'], '''gradients'' using ''manual'' requires external.manual.jacobian_x and external.manual.jacobian_x to be provided.'));
            end
            % jacobian_x and jacobian_u
            if ~isa(options.external.manual.jacobian_x, 'function_handle') || ~isa(options.external.manual.jacobian_u, 'function_handle')
                throw(MException([functionId ':InvalidParameter'], 'external.manual.jacobian_x and external.manual.jacobian_x must be function handles.'));
            end
            if options.nw == 0
                if nargin(options.external.manual.jacobian_x) ~= 2 || nargin(options.external.manual.jacobian_u) ~= 2
                    throw(MException([functionId ':InvalidParameter'], 'external.manual.jacobian_x and external.manual.jacobian_x must have two inputs (x,u)'));
                end
            else 
                if nargin(options.external.manual.jacobian_x) ~= 3 || nargin(options.external.manual.jacobian_u) ~= 3
                    throw(MException([functionId ':InvalidParameter'], 'external.manual.jacobian_x and external.manual.jacobian_x must have three inputs (x,u,w)'));
                end
            end
            % jacobian_n
            if all(options.nn == 0) && ~isempty(options.external.manual.jacobian_n)
                warning([functionId ':InconsistentParameter'], 'external.manual.jacobian_n is non-empty but nonlinear constraint dimensions are all zero. This may be an error, please check!');
                options.external.manual.jacobian_n = {};
            end
            if any(options.nn > 0)
                if isempty(options.external.manual.jacobian_n)
                    throw(MException([functionId ':MissingParameter'], '''gradients'' using ''manual'' requires external.manual.jacobian_n to be provided.'));
                end
                if ~iscell(options.external.manual.jacobian_n)
                    options.external.manual.jacobian_n = {options.external.manual.jacobian_n};
                end
                if all(size(options.external.manual.jacobian_n) == [1 1])
                    options.external.manual.jacobian_n = repmat(options.external.manual.jacobian_n, 1, options.N);
                end
                if ~all(size(options.external.manual.jacobian_n) == [1,options.N])
                    if all(size(options.external.manual.jacobian_n) == [options.N,1])
                        options.external.manual.jacobian_n = options.external.manual.jacobian_n';
                    else
                        throw(MException([functionId ':InvalidParameterDimensions'], 'external.manual.jacobian_n must be either empty, a single function handle, or a cell of length N of function handles.'));
                    end
                end
                for k=1:options.N
                    if ~isa(options.external.manual.jacobian_n{k},'function_handle')
                        throw(MException([functionId ':InvalidParameterDimensions'], 'external.manual.jacobian_n must be either empty, a single function handle, or a cell of length N of function handles.'));
                    end
                    if options.nn(k) > 0 && nargin(options.external.manual.jacobian_n{k}) ~= 1
                        throw(MException([functionId ':InvalidParameter'], 'external.manual.jacobian_n must have one input u'));
                    end
                end
            end
        case 'ccode'
            % Take care of legacy option names
            if ~isfield(options.external, 'ccode')
                if ~isempty(options.Jac_x_struct)
                    options.external.ccode.jacobian_x.structure = options.Jac_x_struct;
                    warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now uses the ''external'' option. The use of ''Jac_x_struct'' is depricated and will be removed in future versions.');
                end
                if ~isempty(options.Jac_x_static)
                    options.external.ccode.jacobian_x.static = options.Jac_x_static;
                    warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now uses the ''external'' option. The use of ''Jac_x_static'' is depricated and will be removed in future versions.');
                end
                if ~isempty(options.Jac_u_struct)
                    options.external.ccode.jacobian_u.structure = options.Jac_u_struct;
                    warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now uses the ''external'' option. The use of ''Jac_u_struct'' is depricated and will be removed in future versions.');
                end
                if ~isempty(options.Jac_u_static)
                    options.external.ccode.jacobian_u.static = options.Jac_u_static;
                    warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now uses the ''external'' option. The use of ''Jac_u_static'' is depricated and will be removed in future versions.');
                end
                if ~isempty(options.Jac_n_struct)
                    options.external.ccode.jacobian_n.structure = options.Jac_n_struct;
                    warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now uses the ''external'' option. The use of ''Jac_n_struct'' is depricated and will be removed in future versions.');
                end
                if ~isempty(options.K_n)
                    options.external.ccode.K_n = options.K_n;
                    warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now uses the ''external'' option. The use of ''K_n'' is depricated and will be removed in future versions.');
                end
            end
            % external.ccode.c
            if ~isfield(options.external, 'ccode') || ~isfield(options.external.ccode, 'c')
                warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now needs that the path to the .c file for the external functions is specified explicitly. Using ''./external_functions.c''. This will be strictly required in future versions.');
                options.external.ccode.c = './external_functions.c';
            end
            % external.ccode.h
            if ~isfield(options.external, 'ccode') || ~isfield(options.external.ccode, 'h')
                warning([functionId ':Depricated'], '''gradients'' using ''ccode'' now needs that the path to the .h file for the external functions is specified explicitly. Using ''./external_functions.h''. This will be strictly required in future versions.');
                options.external.ccode.h = './external_functions.c';
            end
            % jacobian_x
            if ~isfield(options.external, 'ccode') || ~isfield(options.external.ccode, 'jacobian_x') || ~isfield(options.external.ccode.jacobian_x, 'structure')
                options.external.ccode.jacobian_x.structure = reshape(1:options.nx*options.nx,options.nx,options.nx)';
            end
            if ~isfield(options.external, 'ccode') || ~isfield(options.external.ccode, 'jacobian_x') || ~isfield(options.external.ccode.jacobian_x, 'static')
                options.external.ccode.jacobian_x.static = true;
            end
            if any(size(options.external.ccode.jacobian_x.structure) ~= [options.nx,options.nx])
                throw(MException([functionId ':InvalidParameter'], ['external.ccode.jacobian_x.structure must be of dimensions nx times nx, where nx = ' num2str(options.nx) '.']));
            end
            options.Jac_x_struct = options.external.ccode.jacobian_x.structure;
            options.Jac_x_static = options.external.ccode.jacobian_x.static;
            % jacobian_u
            if ~isfield(options.external.ccode, 'jacobian_u') || ~isfield(options.external.ccode.jacobian_u, 'structure')
                options.external.ccode.jacobian_u.structure = reshape(1:options.nx*options.nu,options.nu,options.nx)';
            end
            if ~isfield(options.external.ccode, 'jacobian_u') || ~isfield(options.external.ccode.jacobian_u, 'static')
                options.external.ccode.jacobian_u.static = true;
            end
            if any(size(options.external.ccode.jacobian_u.structure) ~= [options.nx,options.nu])
                throw(MException([functionId ':InvalidParameter'], ['external.ccode.jacobian_u.structure must be of dimensions nx times nu, where nx = ' num2str(options.nx) ', nu = ' num2str(options.nu) '.']));
            end
            options.Jac_u_struct = options.external.ccode.jacobian_u.structure;
            options.Jac_u_static = options.external.ccode.jacobian_u.static;
            % jacobian_n
            if any(options.nn > 0)
                if ~isfield(options.external.ccode, 'jacobian_n')
                    throw(MException([functionId ':MissingParameter'], '''gradients'' using ''ccode'' requires external.manual.jacobian_n to be provided.'));
                end
                if ~iscell(options.external.ccode.jacobian_n.structure)
                    options.external.ccode.jacobian_n.structure = {options.external.ccode.jacobian_n.structure};
                end
                if all(size(options.external.ccode.jacobian_n.structure) == [1, 1])
                    options.external.ccode.jacobian_n.structure = repmat(options.external.ccode.jacobian_n.structure, 1, options.N);
                end
                if any(size(options.external.ccode.jacobian_n.structure) ~= [1, options.N])
                    if all(size(options.external.ccode.jacobian_n.structure) == [options.N, 1])
                        options.external.ccode.jacobian_n.structure = options.external.ccode.jacobian_n.structure';
                    else
                        throw(MException([functionId ':InvalidParameterDimensions'], 'external.ccode.jacobian_n.structure must be either empty, a matrix, or a cell of length N of matrices.'));
                    end
                end
                for k=1:options.N
                    if ~isnumeric(options.external.ccode.jacobian_n.structure{k}) || any(size(options.external.ccode.jacobian_n.structure{k}) ~= [options.nu, options.nn(k)])
                        throw(MException([functionId ':InvalidParameter'], 'Each element of external.ccode.jacobian_n.structure must be of dimensions nu times nn.'));
                    end
                end
            else
                if isfield(options.external.ccode.jacobian_n.structure)
                    warning([functionId ':InconsistentParameter'], 'external.ccode.jacobian_n.structure is non-empty but nonlinear constraint dimensions are all zero. This may be an error, please check!');
                    options.external.ccode.jacobian_n.structure = {};
                end
            end
            options.Jac_n_struct = options.external.ccode.jacobian_n.structure;
            % K_n TODO: check more
            options.K_n = options.external.ccode.K_n;
    end
    % Construct K_n unless provided with 'gradients', 'ccode'
    if ~strcmp(options.gradients,'ccode')
        options.K_n = falcopt.internal.findEqual(options.constraints);
    end
    % Some warnings
    if ~strcmp(options.gradients,'manual') && isfield(options.external, 'manual')
        warning([functionId ':InconsistentParameter'], ['''gradients'' using ''' options.gradients ''', therefore the options external.manual will not be considered.']);
        options.external = rmfield(options.external, 'manual');
    end
    if ~strcmp(options.gradients,'ccode') && isfield(options.external, 'ccode')
        warning([functionId ':InconsistentParameter'], ['''gradients'' using ''' options.gradients ''', therefore the options external.ccode will not be considered.']);
        options.external = rmfield(options.external, 'ccode');
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
    
    %% Check algorithm parameters
    % variable step size
    if isfield(options.variable_stepSize,'active') && ~islogical(options.variable_stepSize.active)
        throw(MException([functionId ':InvalidParameter'], 'variable_stepSize.active must be a boolean.'));
    elseif ~isfield(options.variable_stepSize,'active') % Default
        options.variable_stepSize.active = true;
    end
    if options.variable_stepSize.active
        % steady_state_state
        if isfield(options.variable_stepSize,'steady_state_state')
            if ~isnumeric(options.variable_stepSize.steady_state_state)
                throw(MException([functionId ':InvalidParameter'], ['variable_stepSize.steady_state_state must be a vector of length nx, where nx = ' num2str(options.nx) '.']));
            end
            if all(size(options.variable_stepSize.steady_state_state) == [1,options.nx])
                options.variable_stepSize.steady_state_state = options.variable_stepSize.steady_state_state';
            end
            if any(size(options.variable_stepSize.steady_state_state) ~= [options.nx,1])
                throw(MException([functionId ':InvalidParameterDimension'], ['variable_stepSize.steady_state_state must be a vector of length nx, where nx = ' num2str(options.nx) '.']));
            end
        else % default value .steady_state_state
            options.variable_stepSize.steady_state_state = zeros(options.nx,1);
        end
        % steady_state_input
        if isfield(options.variable_stepSize,'steady_state_input')
            if ~isnumeric(options.variable_stepSize.steady_state_input)
                throw(MException([functionId ':InvalidParameter'], ['variable_stepSize.steady_state_input must be a vector of length nu or a matrix of dimensionsn nu times N, where nu = ' num2str(options.nu) ', N = ' num2str(options.N) '.']));
            end
            if all(size(options.variable_stepSize.steady_state_input) == [1,options.nu])
                options.variable_stepSize.steady_state_input = options.variable_stepSize.steady_state_input';
            end
            if all(size(options.variable_stepSize.steady_state_input) == [options.nu,1])
                options.variable_stepSize.steady_state_input = repmat(options.variable_stepSize.steady_state_input,1,options.N);
            end
            if any(size(options.variable_stepSize.steady_state_input) ~= [options.nu,options.N])
                throw(MException([functionId ':InvalidParameterDimension'], ['variable_stepSize.steady_state_input must be a vector of length nx, where nx = ' num2str(options.nx) '.']));
            end
        else % default value .steady_state_input
            options.variable_stepSize.steady_state_input = zeros(options.nu,options.N);
        end
        % increase_threshold
        if isfield(options.variable_stepSize,'increase_threshold')
            if ~isnumeric(options.variable_stepSize.increase_threshold) || numel(options.variable_stepSize.increase_threshold) ~= 1 || options.variable_stepSize.increase_threshold <= 0
                throw(MException([functionId ':InvalidParameter'], 'variable_stepSize.increase_threshold must be a positive scalar.'));
            end
        else % default value .increase_threshold
            options.variable_stepSize.increase_threshold = 0.75;
        end
        % increase_coeff
        if isfield(options.variable_stepSize, 'increase_coeff')
            if ~isnumeric(options.variable_stepSize.increase_coeff) || numel(options.variable_stepSize.increase_coeff) ~= 1 || options.variable_stepSize.increase_coeff <= 1
                throw(MException([functionId ':InvalidParameter'], 'variable_stepSize.increase_coeff must be a scalar larger than 1.'));
            end
        else % default value .increase_coeff
            options.variable_stepSize.increase_coeff = 1/0.75;
        end
        % decrease_threshold
        if isfield(options.variable_stepSize,'decrease_threshold')
            if ~isnumeric(options.variable_stepSize.decrease_threshold) || numel(options.variable_stepSize.decrease_threshold) ~= 1 || options.variable_stepSize.decrease_threshold <= 0
                throw(MException([functionId ':InvalidParameter'], 'variable_stepSize.decrease_threshold must be a positive scalar.'));
            end
        else % default value .increase_threshold
            options.variable_stepSize.decrease_threshold = 0.25;
        end
        % decrease_coeff
        if isfield(options.variable_stepSize, 'decrease_coeff')
            if ~isnumeric(options.variable_stepSize.decrease_coeff) || numel(options.variable_stepSize.decrease_coeff) ~= 1 || options.variable_stepSize.decrease_coeff <= 0 || options.variable_stepSize.decrease_coeff >= 1
                throw(MException([functionId ':InvalidParameter'], 'variable_stepSize.decrease_coeff must be a scalar between 0 and 1.'));
            end
        else % default value .decrease_coeff
            options.variable_stepSize.decrease_coeff = 0.75;
        end
        % alpha_max
        if isfield(options.variable_stepSize,'alpha_max')
            if ~isnumeric(options.variable_stepSize.alpha_max) || numel(options.variable_stepSize.alpha_max) ~= 1 || options.variable_stepSize.alpha_max <= 0
                throw(MException([functionId ':InvalidParameter'], 'variable_stepSize.alpha_max must be a positive scalar.'));
            end
        elseif strcmp(options.gradients,'casadi') % automatically compute .alpha_max
            options.variable_stepSize.alpha_max = falcopt.computeStepSize(options.variable_stepSize.steady_state_state, ...
                                                                  options.variable_stepSize.steady_state_input, options);
            % check value of computed .alpha_max                                            
            if isinf(options.variable_stepSize.alpha_max)
                throw(MException([functionId ':InvalidParameter'], ['Error while computing ''variable_stepSize.alpha_max'', consider to manually specify the value of alpha with' ...
                                                                    '''variable_stepSize.alpha_max''']));
            end
        else
            throw(MException([functionId ':InvalidParameter'], ['If option variable_stepSize.active = ''true'' and option gradients ~= ''casadi'',' ...
                                                                'then the option ''variable_stepSize.alpha_max'' must be specified. ' ...
                                                                'Otherwise set the option gradients = ''casadi''']));
        end
        % alpha_min
        if isfield(options.variable_stepSize,'alpha_min')
            if ~isnumeric(options.variable_stepSize.alpha_min) || numel(options.variable_stepSize.alpha_min) ~= 1 || ...
               options.variable_stepSize.alpha_min <= 0 || options.variable_stepSize.alpha_min > options.variable_stepSize.alpha_max
                throw(MException([functionId ':InvalidParameter'], 'variable_stepSize.alpha_min must be a positive scalar, smaller or equal than alpha_max.'));
            end
        else % default value .alpha_min 
            options.variable_stepSize.alpha_min = 0.1*options.variable_stepSize.alpha_max;
        end
    elseif ~isfield(options.variable_stepSize, 'alpha_max')
        % default value .alpha_max in case of constant step size
        options.variable_stepSize.alpha_max = 0.4;
    end
end
