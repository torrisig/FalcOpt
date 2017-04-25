% GENERATECODE Generate the algorithm .c file
%
% [code, data, optCode, libr, info] = generateCode(N,nx,nu,nc,Q,p,R,...);
% or
% [code, data, optCode, libr, info] = generateCode(N,nx,nu,nc,Q,p,R,, 'par1', val1, 'par2', val2, ...) [with options]
%
% where required inputs are:
% N     - Prediction horizon
% nx    - Number of states
% nu    - Number of inputs
% nc    - Number of constraints (per stage)
% Q     - Weight matrix for states (cost)
% P     - Weight matrix for terminal states (cost)
% R     - Weight matrix for inputs (cost)
%
% Outputs:
% - code:  string containing the generated code
% - data: string containing the static data
% - optCode: string containing the generate code for the algorithm only
% - libr: string containing all the included libraries
% - info: a struct containing the number of FLOPS
% Standard use:
% fprintf(file,[libr '\n' data '\n' code_spec '\n' code_alg '\n' optCode])
%
% The following options are available:
%   .nw                 - Known disturbance dimension
%   .trackReference     - A boolean. Track a desired time-varying reference. Default: false
%   .stepSize           - Step size alpha. Default: 0.3
%   .parLs              - Armijo line search step parameter. Default: 0.3
%   .maxIt              - Max number of iterations. Default: 1000
%   .maxItLs            - Max number of line search iterations. Default: 10
%   .eps                - Tolerance. Defualt: 1e-6
%   .tolLs              - Line search min progress. Default: 1e-4
%   .precision          - 'double'(default) or 'single'
%   .indent             - Indentation to be used in code generation.
%                         Default: '\t'
%   .inline             - Inline keyword to be used. Default: 'inline'
%   .gendir             - Path of the .c file folder
%   .verbose            - Level of procedural output of this function.
%                         Default: 0
%   .test               - Level of tests performed. Default: 0
%   .debug              - Level of debug. Default: 1
%   .Jac_x_static       - Boolean, false(default): Jacobian_x is dynamic
%   .Jac_u_static       - Boolean, false(default): Jacobian_u is dynamic
%   .Jac_x_struct       - Matrix containing the structure of the Jacobian_x
%   .Jac_u_struct       - Matrix containing the structure of the Jacobian_u
%   .Jac_g_struct       - Matrix containing the structure of the Jacobian_g
%   .Jac_m_struct       - Matrix containing the structure of the Jacobian_m
%   .Jac_n_struct       - Matrix containing the structure of the Jacobian_n
%   .merit_function     - 1,2(default) or Inf
%   .contractive        - A boolean. Default: 'false'
%   .terminal           - A boolean. Default: 'false'
%   .constraints_handle - A function handle to provided nonlinear constraint function
%   .gradients          - Different way of automatic function differentiation/generation
%                         Can be: 'casadi','matlab','manual' or 'ccode'. Default: 'casadi'
%   .external_jacobian_x- Function handle to provide jacobian_x (only with options gradients = 'manual')
%   .external_jacobian_u- Function handle to provide jacobian_u (only with options gradients = 'manual')
%   .external_jacobian_n- Function handle to provide jacobian_n (only with options gradients = 'manual')
%
% user supplied functions: model_mpc(real* x, real* u, real* xp);
%                          Jacobian_x(real*x, real* u, real* F);
%                          Jacobian_u(real* x, real* u, real* G);
%                          build_g(real* u, real* gpsl);
%                          build_Dg(real* u, real* Dg);
%                          build_inv(real* u, real* sl, real* inv);


function [info] = generateCode(varargin)
indentTypes = {'generic', 'data', 'code'};
p = inputParser;
p.addRequired('dynamics', @(x)isa(x,'function_handle'));
p.addRequired('N', @isnumeric); % prediction horizon
p.addRequired('nx', @isnumeric); % state dimension
p.addRequired('nu', @isnumeric); % input dimension
p.addRequired('Q', @(x) (isnumeric(x) && min(eig(x))>= 0)); % weight on the state (cost)
p.addRequired('P', @(x) (isnumeric(x) && min(eig(x))>= 0)); % weight on the terminal state (cost)
p.addRequired('R', @(x) (isnumeric(x) && min(eig(x))>= 0)); % weight on the input (cost)
p.addParameter('nn', 0, @(x) (isnumeric(x) || iscell(x))); % number of nonlinear constraints (per stage)
p.addParameter('nw', 0, @(x) (isnumeric(x) && x>=0)); % known disturbance dimension
p.addParameter('trackReference',false,@islogical);
p.addParameter('stepSize', 0.3, @(x)(isnumeric(x) && x > 0)); % step size alpha
p.addParameter('parLs', 0.3, @(x)(isnumeric(x) && ( (x > 0) && (x < 1)) )); % Armijo line search step parameter
p.addParameter('maxIt', 4000, @(x)(isnumeric(x) && x>0 && mod(x,1) == 0)); % max number of iter
p.addParameter('maxItLs', 10, @(x)(isnumeric(x) && x>0 && mod(x,1) == 0)); % max number of line search iterations
p.addParameter('eps', 1e-6, @(x)(isnumeric(x) && x >= eps)); % tolerance
p.addParameter('tolLs', 1e-4, @(x)(isnumeric(x) && x >= eps)); % line search min progress
p.addParameter('precision', 'double', @(x)( strcmp(x,'double')|| strcmp(x,'single') ));
p.addParameter('indent', struct('code', '', 'data', '', 'generic', '\t'), @(x)(ischar(x) || (isstruct(x) && isfield(x, 'generic') && all(cellfun(@(y)(~isfield(x, y) || ischar(x.(y))), indentTypes)))));
p.addParameter('inline', '', @ischar);
p.addParameter('name','my_code', @ischar);
p.addParameter('build_MEX',true, @islogical);
p.addParameter('compile', true, @islogical);
p.addParameter('gendir', '', @ischar);
p.addParameter('verbose', 0, @isnumeric);
p.addParameter('test', 0, @isnumeric);
p.addParameter('debug', 1, @(x)(isnumeric(x) && x >= 0));
p.addParameter('Jac_x_static',false,@islogical);
p.addParameter('Jac_u_static',false,@islogical);
p.addParameter('Jac_x_struct',Inf,@isnumeric);
p.addParameter('Jac_u_struct',Inf,@isnumeric);
p.addParameter('Jac_g_struct',Inf,@isnumeric);
p.addParameter('Jac_m_struct',Inf,@isnumeric);
p.addParameter('Jac_n_struct',[],@isnumeric);
p.addParameter('amu_struct',Inf,@iscell);
p.addParameter('umb_struct',Inf,@iscell);
p.addParameter('merit_function', 2, @(x)( (x == 1)|| (x == 2) ) || ((x == Inf) || (x == 0)));
p.addParameter('contractive', false, @islogical);
p.addParameter('terminal', false, @islogical);
p.addParameter('forceGradient',true,@islogical); % This parameter is to be deleted when converter example works without build_tbi()
p.addParameter('box_lowerBound', [], @isnumeric);
p.addParameter('box_upperBound', [], @isnumeric);
p.addParameter('constraints_handle', Inf, @(x) (iscell(x) || isa(x,'function_handle')));
p.addParameter('gradients', 'casadi', @(x)(ischar(x)));
p.addParameter('external_jacobian_x', [], @(x)isa(x,'function_handle'));
p.addParameter('external_jacobian_u', [], @(x)isa(x,'function_handle'));
p.addParameter('external_jacobian_n', [], @(x)isa(x,'function_handle'));
p.parse(varargin{:});
o = p.Results;

%% Processing options
% Dimensions
if size(o.Q,1) ~= o.nx
    error('Q matrix must be of dimension nx');
end
if size(o.R,1) ~= o.nu
    error('R matrix must be of dimension nu');
end
if size(o.P,1) ~= o.nx
    error('P matrix must be of dimension nx');
end
%     if length(o.nn) == 1
%         o.nc = repmat(o.nc, 1, o.N);
%     elseif length(o.nc) ~= o.N
%         error('Invalid nc');
%     end
%     if length(unique(o.nc)) > 1
%         error('For now only same nc are supported.');
%     end

if (max(max(o.Jac_x_struct)) == Inf)
    o.Jac_x_struct = reshape(1:o.nx*o.nx,o.nx,o.nx)';
end
if (max(max(o.Jac_u_struct)) == Inf)
    o.Jac_u_struct = reshape(1:o.nx*o.nu,o.nu,o.nx)';
end

if (~isempty(o.Jac_n_struct))&&(max(size( o.Jac_n_struct) ~= [o.nu,o.nn]))
    error_string = sprintf(['Jac_n_struct must be of size [%d, %d]' '\n'], o.nu,o.nn);
    error(error_string);
end

%     if (max(max(o.Jac_m_struct)) == Inf) % OBSOLETE - TO be deleted
%         o.Jac_m_struct = ones( o.nn, o.nn);
%     end
%     if (max(max(o.Jac_g_struct)) == Inf)
%         o.Jac_g_struct = reshape(1:o.nn*o.nu,o.nn,o.nu)';
%     end
if min(size(o.Jac_x_struct)== [o.nx,o.nx]) == 0
    error('Jacobian x structure must be [nx,nx]');
end
if min(size(o.Jac_u_struct)== [o.nx,o.nu]) == 0
    error('Jacobian u structure must be [nx,nu]');
end
%     if min(size(o.Jac_m_struct)== [o.nn,o.nn]) == 0 % OBSOLETE - TO be deleted
%         error('Jacobian M structure must be [nc,nc]');
%     end
%     if min(size(o.Jac_g_struct)== [o.nu,o.nn]) == 0
%         error('Jacobian constraint structure must be [nu,nc]');
%     end

if (o.terminal && o.contractive)
    error('Cannot have both contractive and terminal constraint');
end

%% lower and upper bounds
one_lowerBound = false;
if ( ~isempty(o.box_lowerBound))
    % transform o.box_lowerBound in a matrix of dims (o.nu,o.N), with -Inf where no bound is given
    % flag one_lowerBound = true if only one lower bound
    if( any(size(o.box_lowerBound)~=[1,o.nu])&& any(size(o.box_lowerBound)~=[o.nu,1])&& ...
            any(size(o.box_lowerBound)~=[o.nu, o.N])&& any(size(o.box_lowerBound)~=[o.N,o.nu]))
        error('box_lowerBound must be either empty (no lower bounds), of dimension [nu,1] (time invariant lower bounds) or [nu, N] (stage dependent lower bounds)')
    end
    if min(size(o.box_lowerBound)==[o.nu,1])
        one_lowerBound = true;
        o.box_lowerBound = repmat(o.box_lowerBound,1,o.N);
    end
    if min(size(o.box_lowerBound)==[1,o.nu])
        one_lowerBound = true;
        o.box_lowerBound = repmat(o.box_lowerBound',1,o.N);
    end
    if min(size(o.box_lowerBound)==[o.N,o.nu])
        o.box_lowerBound = o.box_lowerBound';
    end
else
    one_lowerBound = true;
    o.box_lowerBound = -Inf(o.nu,o.N);
end

% if max(max(o.box_lowerBound)) ~= -Inf
%     o.lb_flag = true;
% else
%     o.lb_flag = false;
% end

one_upperBound = false;
if ( ~isempty(o.box_upperBound))
    % transform o.box_upperBound in a matrix of dims (o.nu,o.N), with +Inf where no bound is given
    % flag one_upperBound = true if only one upper bound
    if( any(size(o.box_upperBound)~=[1,o.nu])&& any(size(o.box_upperBound)~=[o.nu,1])&& ...
            any(size(o.box_upperBound)~=[o.nu, o.N])&& any(size(o.box_upperBound)~=[o.N,o.nu]))
        error('box_upperBound must be either empty (no upper bounds), of dimension [nu,1] (time invariant upper bounds) or [nu, N] (stage dependent upper bounds)')
    end
    if min(size(o.box_upperBound)==[o.nu,1])
        one_upperBound = true;
        o.box_upperBound = repmat(o.box_upperBound,1,o.N);
    end
    if min(size(o.box_upperBound)==[1,o.nu])
        one_upperBound = true;
        o.box_upperBound = repmat(o.box_upperBound',1,o.N);
    end
    if min(size(o.box_upperBound)==[o.N,o.nu])
        o.box_upperBound = o.box_upperBound';
    end
else
    one_upperBound = true;
    o.box_upperBound = +Inf(o.nu,o.N);
end

% check feasibility of the constraints
for jj=1:o.N
    for ii=1:o.nu
        if o.box_lowerBound(ii,jj) > o.box_upperBound(ii,jj)
            if (one_lowerBound)&&(one_upperBound)
                error('box_lowerBound > box_upperBound');
            else
                error('box_lowerBound > box_upperBound at stage %i', jj);
            end
        elseif o.box_lowerBound(ii,jj) == o.box_upperBound(ii,jj)
            error('The case of box_lowerBound = box_upperBound is not implemented yet. Please make sure that box_lowerBound < box_upperBound');
        end
    end
end

if isequal(o.gradients,'ccode')
    o.K_amu = o.amu_struct;
    o.K_umb = o.umb_struct;
    o.K_lb = [];
    o.K_ub = [];
else
    o.K_amu = detect_structure( o.box_lowerBound,o );
    cons_lb = ~isinf(o.box_lowerBound);
    o.K_lb = detect_structure( cons_lb(:,sum(cons_lb,1) >= 1), o );
    o.K_umb = detect_structure( o.box_upperBound,o );
    cons_ub = ~isinf(o.box_upperBound);
    o.K_ub = detect_structure( cons_ub(:,sum(cons_ub,1) >= 1), o );

end

% TO BE DELETED
%     if ( ~isempty(o.box_constraints))
%         if( any(size(o.box_constraints)~=[2,o.nu])&& any(size(o.box_constraints)~=[o.nu,2])&& ...
%                 any(size(o.box_constraints)~=[o.nu*o.N,2])&& any(size(o.box_constraints)~=[2,o.nu*o.N]))
%             error('box_constraints must be of dimension [2,nu] or [2,N*nu]')
%         elseif ( all(size(o.box_constraints)==[2,o.nu])||all(size(o.box_constraints)==[2,o.nu*o.N]))
%                 o.box_constraints = o.box_constraints';
%         end
%         for i= 1:o.nu
%             if( o.box_constraints(i,1) >= o.box_constraints(i,2))
%                 error('Box constraint %i is infeasible',i);
%             end
%         end
%         if( all(size(o.box_constraints)==[o.nu,2]) )
%             o.box_constraints = repmat(o.box_constraints,o.N,1);
%         end
%     end

%% nonlinear constraint handle

if (any(size(o.constraints_handle)~= [1,o.N])&&any(size(o.constraints_handle)~= [o.N,1]))&&...
        (isempty(o.constraints_handle)&&any(size(o.constraints_handle)~= [1,1]))
    error('constraints_handle must be either empty, or a function handle, or a cell of function handles of dimension [N,1]');
end
if length(o.constraints_handle) >1
    for ii=1:o.N
        if ~isa(o.constraints_handle{ii},'function_handle')
            error('constraints_handle{%i} is not a function handle', ii);
        end
    end
end
if min(size(o.constraints_handle)== [o.N,1])
    o.constraints_handle = o.constraints_handle';
elseif min(size(o.constraints_handle)== [1,1])
    if iscell(o.constraints_handle)
        o.constraints_handle = repmat(o.constraints_handle,1,o.N);
    else
        o.constraints_handle = repmat({o.constraints_handle},1,o.N);
    end
    o.K_n = {1:o.N};
end

if min(size(o.external_jacobian_n)== [o.N,1])&& isequal(o.gradients,'manual')
    o.external_jacobian_n = o.external_jacobian_n';
elseif min(size(o.external_jacobian_n)== [1,1])&& isequal(o.gradients,'manual')
    if iscell(o.external_jacobian_n)
        o.external_jacobian_n = repmat(o.external_jacobian_n,1,o.N);
    else
        o.external_jacobian_n = repmat({o.external_jacobian_n},1,o.N);
    end
    o.K_n = {1:o.N};
end

if  ~isempty(o.constraints_handle)
    if ~isfield(o,'K_n')
        o.K_n = detect_different_NLconstraints(o);
    end
else
    o.K_n = [];
end

% check dims of the function handles
if ((any(size(o.nn)~= [1,o.N]))&&(any(size(o.nn)~= [o.N,1])))&&(any(size(o.nn)~= [1,1]))
    error('nn can either be a scalar (same nonlinear constraint for every stage) or a vector of dims [N,1] (variable nonlinear constraint)');
elseif min(size(o.nn) == [o.N,1])
    o.nn = o.nn';
elseif min(size(o.nn) == [1,1])
    if iscell(o.nn)
        o.nn = repmat(o.nn, 1,o.N);
    else
        o.nn = repmat({o.nn}, 1,o.N);
    end
end

% check dims of the function handles
for ii = 1:o.N
    test_u1 = ones(1,o.nu);
    test_f1 = o.constraints_handle{ii}(test_u1);
    
    test_u2 = ones(o.nu,1);
    test_f2 = o.constraints_handle{ii}(test_u2);
    
    if ((size(test_f1, 1) > 1) && (size(test_f1, 2) > 1))&&((size(test_f2, 1) > 1) && (size(test_f2, 2) > 1))
        error('The output of the nonlinear constraint function handles must be a vector');
    end
    if any(test_f1(:) ~= Inf)&&(any(size(test_f1) ~= [o.nn{ii},1]) && any(size(test_f1) ~= [1,o.nn{ii}]))&&...
            (any(size(test_f2) ~= [o.nn{ii},1]) && any(size(test_f2) ~= [1,o.nn{ii}]))
        error('The output of the nonlinear constraint function handle %i is not of length %d',ii,o.nn{ii});
    end
end

% check dimension of o.jac_n_struct if porvided
if ~isempty(o.Jac_n_struct)
    if ~iscell(o.Jac_n_struct)
        o.Jac_n_struct = repmat({o.Jac_n_struct}, 1,o.N);
    end
end

% check external function and 'manual' gradients compatibility
if ~strcmp(o.gradients,'manual')
    if ~isempty(o.external_jacobian_x)
        warning('If gradients options is not manual the parameter external_jacobian_x will not be considered');
    elseif ~isempty(o.external_jacobian_u)
        warning('If gradients options is not manual the parameter external_jacobian_u will not be considered');
    elseif ~isempty(o.external_jacobian_n)
        warning('If gradients options is not manual the parameter external_jacobian_n will not be considered');
    end
else
    if isempty(o.external_jacobian_x)
        error('external_jacobian_x option must be provided');
    elseif isempty(o.external_jacobian_u)
        error('external_jacobian_u option must be provided');
    elseif isempty(o.external_jacobian_n)
        error('external_jacobian_n option must be provided');
    end
    if( o.nw > 0)
        if( nargin(o.external_jacobian_x)~= 3)
            error('Inputs of external_jacobian_x function handle must be x,u,w');
        end
        if( nargin(o.external_jacobian_u)~= 3)
            error('Inputs of external_jacobian_u function handle must be x,u,w');
        end
    else
        if( nargin(o.external_jacobian_x)~= 2)
            error('Inputs of external_jacobian_x function handle must be x,u');
        end
        if( nargin(o.external_jacobian_u)~= 2)
            error('Inputs of external_jacobian_u function handle must be x,u');
        end
    end
    if( o.nn{1}>0)
        for i= 1:length(o.external_jacobian_n)
            if( nargin(o.external_jacobian_n{i})~= 1)
             error('Inputs of external_jacobian_n function handle must be u');
            end
        end
    end
end

%% check version of toolbox for automatic generation
switch o.gradients
    case 'casadi'
        import casadi.*
        a = {which('SX.sym'),which('is_constant'),which('sparsity'),...
            which('densify'),which('is_equal'),which('to_double')};
        if(any(cellfun(@isempty,a)))
            error('Check the CasaDi version or install it. Program tested with CasaDi v3.1.1')
        end
    case 'matlab'
        v = ver('symbolic');
        if( isempty(v))
            error('Symbolic Math Toolbox not installed. Install it to continue')
        end
end

% Indentation
if ~isstruct(o.indent)
    indent = o.indent;
    o.indent = struct();
    for i=1:length(indentTypes)
        o.indent.(indentTypes{i}) = indent;
    end
end
for i=1:length(indentTypes)
    if ~isfield(o.indent, indentTypes{i})
        o.indent.(indentTypes{i}) = o.indent.generic;
    end
    o.indent.(indentTypes{i}) = o.indent.(indentTypes{i})(:)'; % Make sure is row vector
end



%     if (o.contractive || o.terminal)
%         %sum(o.nc) = o.N*o.nc + 1;
%         o.nc(end+1) = 1;
%     else
%         o.nc(end+1) = 0;
%     end

% Types
switch o.precision
    case 'double'
        o.sqrt = 'sqrt';
        o.max = 'my_fmax';
        o.min = 'my_fmin';
        o.abs = 'fabs';
        o.real = 'double';
    case 'single'
        o.sqrt = 'sqrtf';
        o.max = 'my_fmaxf';
        o.min = 'my_fminf';
        o.abs = 'fabsf';
        o.real = 'float';
end

alpha = o.stepSize;
alpha_eps = alpha*o.eps;
alpha_eps2 = alpha_eps* o.eps;
alpha2_eps = alpha*alpha_eps;
alpha2_eps2 = alpha2_eps*o.eps;

%% Setup info
info.flops.it = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0, 'div', 0, 'casadi', 0);
info.flops.ls = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0, 'div', 0, 'casadi', 0);
info.flops.once = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0, 'div', 0, 'casadi', 0);
info.src = [];
info.header = [];

%% Generating code
libr = sprintf(['#include "math.h"' '\n']);
code = '';
data = '';
optCode = '';

%% Generating max and min

code = [code, sprintf([ o.indent.code o.real ' ' o.max '(' o.real ' x, ' o.real ' y){' '\n' ...
    o.indent.code o.indent.generic 'if (x < y) { return y; }' '\n' ...
    o.indent.code o.indent.generic 'else { return x; }' '\n' ...
    o.indent.code '}' '\n\n'])];
code = [code, sprintf([ o.indent.code o.real ' ' o.min '(' o.real ' x, ' o.real ' y){' '\n' ...
    o.indent.code o.indent.generic 'if (x > y) { return y; }' '\n' ...
    o.indent.code o.indent.generic 'else { return x; }' '\n' ...
    o.indent.code '}' '\n\n'])];

%% CasaDi model_mpc, Jacobian_u and Jacobian_x automatic generation
switch o.gradients
    case 'casadi'
        [d,c,i] = casadi_jacobians( o, o.dynamics,'model_mpc','','x', 'u'); % generate casadi functions
        code = [code, c];
        data = [data, d];
        
        o.Jac_x_static = i.in_F.static;                         % 1/0 if static data
        o.Jac_x_struct = i.in_F.struct.structure.stored.mat;    % structure of jacobian_x
        o.Jac_u_static = i.in_G.static;                         % 1/0 if static data
        o.Jac_u_struct = i.in_G.struct.structure.stored.mat;    % structure of jacobian_u
        
        info.flops.once.casadi = info.flops.once.casadi + i.y.flops*(o.N+1);    % flops model
        info.flops.ls.casadi = info.flops.ls.casadi + i.y.flops*(o.N+1);        % flops model
        info.flops.it.casadi = info.flops.it.casadi + i.in_F.flops*(o.N);       % flops jacobian_x
        info.flops.it.casadi = info.flops.it.casadi + i.in_G.flops*(o.N);       % flops jacobian_u
        info.src = [info.src; i.src];                                           % external src file
        info.header = [info.header; i.header];                                  % external header file
        libr = [libr, sprintf(['#include "' i.header{:} '"' '\n'])];            % add .h file to libr
    case 'matlab'
        [d,c,i] = matlab_jacobians(o,o.dynamics,'model_mpc','x','u');          % generate function and jacobians
        data = [data, d];
        code = [code, c];
        
        o.Jac_x_static = i.in_F.static;
        o.Jac_x_struct = i.in_F.struct.structure.mat;                        % structure of jacobian_x
        o.Jac_u_static = i.in_G.static;
        o.Jac_u_struct = i.in_G.struct.structure.mat;                        % structure of jacobian_u
        
        info.flops.once = falcopt.addFlops( info.flops.once, falcopt.multFlops(i.y.flops,o.N+1));  % flops model_mpc
        info.flops.ls = falcopt.addFlops( info.flops.ls, falcopt.multFlops(i.y.flops,o.N+1));      % flops model_mpc
        info.flops.it = falcopt.addFlops( info.flops.it, falcopt.multFlops(i.in_F.flops,o.N));     % flops jacobian_x
        info.flops.it = falcopt.addFlops( info.flops.it, falcopt.multFlops(i.in_G.flops,o.N));     % flops jacobian_u
    case 'manual'
        
        % generate external function (matlab)
        [ d, c, i] = external_jacobians(o);
        data = [data, d];
        code = [code, c];
        
        o.Jac_x_static = i.in_F.static;
        o.Jac_x_struct = i.in_F.struct.structure.mat;                        % structure of jacobian_x
        o.Jac_u_static = i.in_G.static;
        o.Jac_u_struct = i.in_G.struct.structure.mat;                        % structure of jacobian_u
        
        info.flops.once = falcopt.addFlops( info.flops.once, falcopt.multFlops(i.y.flops,o.N+1));  % flops model_mpc
        info.flops.ls = falcopt.addFlops( info.flops.ls, falcopt.multFlops(i.y.flops,o.N+1));      % flops model_mpc
        info.flops.it = falcopt.addFlops( info.flops.it, falcopt.multFlops(i.in_F.flops,o.N));     % flops jacobian_x
        info.flops.it = falcopt.addFlops( info.flops.it, falcopt.multFlops(i.in_G.flops,o.N));     % flops jacobian_u
        
    case 'ccode'
        % all functions will be provided in external_functions.c file
        path = cd;
        %path = [path sprintf(['\\' o.gendir])];
        if( exist(fullfile(path, 'external_functions.c'), 'file'))
            info.src = [info.src; {sprintf('external_functions.c')}]; 
            % warning string construction
            string = ['Make sure you provide the following c functions:\n' ...
                'void model_mpc(const double* x, const double* u, const double* v, double* xp)\n'];
            if( o.nw>1)
                string = [string, 'void Jacobian_u(const double* x, const double* u, const double* q, double* G)\n' ...
                    'void Jacobian_x(const double* x, const double* u, const double* q, double* F)\n' ];
            else
                string = [ string, 'void Jacobian_u(const double* x, const double* u, double* G)\n' ...
                    'void Jacobian_x(const double* x, const double* u, double* F)\n'];
            end
            string = [string,  '\n potentially also functions: \n void build_n(const double* u, const unsigned int k, double* n))\n'];
            string = [string, 'void build_Dn(const double* u, const unsigned int k, double* Dn)\n'];
            string = [string, 'void build_amu(const double* u, const unsigned int k, double* r)\n'];
            string = [string, 'void build_umb(const double* u, const unsigned int k, double* r)\n'];
            string = [string,'no check of dimensions, the structure of matrices must be provided in an "ordered" manner: \n'...
                'Jac_x_struct, Jac_u_struct and Jac_n_struct parameters must be provided \n' ];
            warning(sprintf(string));
            
            % add .h file to library
            if( exist(fullfile(path, 'external_functions.h'), 'file'))
                if ~isempty(o.gendir)
                    libr = [libr, sprintf('#include "../external_functions.h" \n')];
                else
                    libr = [libr, sprintf('#include "external_functions.h" \n')];
                end
                info.header = [info.header; {'external_function.h'}];
            else
                if ~isempty(o.gendir)
                    libr = [libr, sprintf('#include "../external_functions.c" \n')];
                else
                    libr = [libr, sprintf('#include "external_functions.c" \n')];
                end
            end
            data = [ data, sprintf('/*static data for jacobian_u*/\nstatic double G[%d];\n',max(o.Jac_u_struct(:)))];
            data = [ data, sprintf('/*static data for jacobian_x*/\nstatic double F[%d];\n',max(o.Jac_x_struct(:)))];
        else
            error(sprintf(['file external_functions.c not found in directory ' path]));
        end
    otherwise
        error(' "gradients" option can be: casadi, matlab, manual or ccode');
end

%% constraints automatic generation
switch o.gradients
    case 'casadi'
        % generate build_g, build_Dg, build_Dn functions
        [ d, c, i] = generate_n_and_Dn( o,'casadi');
        code = [code, c];
        data = [data, d];
        o.nc = i.nc;  
        for jj=1:length(o.K_n)
            %             o.Jac_g_struct = i.in_Dg_g.struct.structure.stored.mat;           % structure of Dg
            o.Jac_n_struct{jj} = i.in_Dn_n{jj}.struct.structure.stored.mat;     % structure of Dn
            %             o.Jac_g_static = i.in_Dg_g.static;
                                                          % structure of constraints
            info.flops.once.casadi = info.flops.once.casadi + i.in_n{jj}.flops*length(o.K_n{jj});     % flops build_n
            info.flops.ls.casadi = info.flops.ls.casadi + i.in_n{jj}.flops*length(o.K_n{jj});         % flops build_n
            info.flops.it.casadi = info.flops.it.casadi + i.in_Dn_n{jj}.flops*length(o.K_n{jj});      % flops build_Dn
        end
        try
            info.src = [info.src; i.src];                                           % external src file
            info.header = [info.header; i.header];                                  % external header file
            libr = [libr, sprintf(['#include "' i.header{:} '"' '\n'])];
        end
        
        if ~o.forceGradient
            o.Jac_m_struct = i.Jac_m_struct;    % structure of inverse matrix
        end
    case {'matlab','manual'}
        % generate build_g, build_Dg, build_Dn functions
        [ d, c, i] = generate_n_and_Dn( o, 'matlab');
        code = [code, c];
        data = [data, d];
        
        %o.Jac_g_struct = i.in_Dg_g.struct.structure.mat;        % structure of Dg
        %o.Jac_g_static = i.in_Dg_g.static;
        o.nc = i.nc;
        
        for jj=1:length(o.K_n)
            o.Jac_n_struct{jj} = i.in_Dn_n{jj}.struct.structure.mat;        % structure of Dn
            
            info.flops.once = falcopt.addFlops( info.flops.once, falcopt.multFlops(i.in_n{jj}.flops,length(o.K_n{jj})));  % flops build_n
            info.flops.ls = falcopt.addFlops( info.flops.ls, falcopt.multFlops(i.in_n{jj}.flops,length(o.K_n{jj})));      % flops build_n
            info.flops.it = falcopt.addFlops( info.flops.it, falcopt.multFlops(i.in_Dn_n{jj}.flops,length(o.K_n{jj})));   % flops build_Dn
        end
        if ~o.forceGradient
            o.Jac_m_struct = i.Jac_m_struct;
        end
    case 'ccode'
        % build o.nc (constraints structure) % ToDo TOmmaso
        o.nc = cell2mat(o.nn);
        
        
        if( o.terminal || o.contractive)
            o.nc = [o.nc 1];
        else
            o.nc = [o.nc 0];
        end
        
        for i=1:o.N 
            lb = sum(~isinf(o.box_lowerBound(:,i)));
            ub = sum(~isinf(o.box_upperBound(:,i)));
            o.nc(i) = o.nc(i)+lb+ub;
        end
    end   
        
        
    % find the structure of the current Jac_n
    o.Jac_n_struct_hor = cell(1,o.N);
    
    for ii= 1:o.N
        flag_exit = false;
        for jj=1:length(o.K_n)
            for kk=1:length(o.K_n{jj})
                if ii == o.K_n{jj}(kk)
                    o.Jac_n_struct_hor{ii} = o.Jac_n_struct{jj};
                    flag_exit = true;
                    break;
                end
            end
            if flag_exit
                break;
            end
        end
    end


%     else % case external .c functions
%         % build o.nc (constraints structure)
%         box = o.box_constraints;
%         o.nc = repmat(o.nn,1,o.N);
%         if( o.terminal || o.contractive)
%             o.nc = [o.nc 1];
%         else
%             o.nc = [o.nc 0];
%         end
%         ind = 1;
%
%         for( i=1:o.nu:length(box)) % not tested ToDo Tommaso
%             nu = sum(~isinf(box(i:i+o.nu-1,:)));
%             o.nc(ind) = o.nc(ind)+nu(1)+nu(2);
%             ind = ind+1;
%         end
%     end

% check patterns of stage constraints
o.K_nc = detect_structure_constraints( o.nc,o );


%% Initialization main algorithm

c_w_dec = argument_w(o,true);
c_tr_dec = argument_def(o,true);
c_contr_dec = argument_contr_value(o,true);
if o.debug > 1
    optCode = [optCode, sprintf(['int proposed_algorithm(const ' o.real '* x0, ' o.real '* u' c_w_dec c_tr_dec c_contr_dec ...
        ', ' o.real '* x, ' o.real '* fval, unsigned int* iter, unsigned int* iter_ls){' '\n' '\n'])];
else
    optCode = [optCode, sprintf(['int proposed_algorithm(const ' o.real '* x0, ' o.real '* u' c_w_dec c_tr_dec c_contr_dec ...
        ', unsigned int* iter, unsigned int* iter_ls){' '\n' '\n'])];
end


%% Declaration of variables

optCode = [optCode, sprintf(['\t' 'unsigned int conditions_f = 1, conditions_x = 1,' '\n' ...
    '\t\t' 'conditions_n = 1, cond_err = 1,' '\n',...
    '\t\t' 'it = 0, it_ls = 0, ii = 0, jj = 0;' '\n'...
    '\t' 'int cond = -2;' '\n'])];

if o.merit_function == 0
    optCode = [optCode, sprintf(['\t' 'unsigned int flag_ini = 0;' '\n'])];
else
    optCode = [optCode, sprintf(['\t' 'unsigned int reset_rho = 0;' '\n'])];
end

optCode = [optCode, sprintf(['\t' o.real ' J= 0.0, alpha = ' falcopt.num2str(alpha, o.precision) ', Jt = 0.0,' '\n'])];
if o.debug > 1
    optCode = [optCode, sprintf(['\t\t' 'xp[%d],' '\n'],o.N*o.nx)];
else
    optCode = [optCode, sprintf(['\t\t' 'x[%i], xp[%d],' '\n'],o.N*o.nx,o.N*o.nx)];
end
optCode = [optCode, sprintf(['\t\t' 'dot_J[%d], du[%d], up[%d],' '\n'],o.N*o.nu,o.N*o.nu,o.N*o.nu)];
optCode = [optCode, sprintf(['\t\t' 'gps[%d], gpsp[%d],' '\n'],sum(o.nc),sum(o.nc))];
optCode = [optCode, sprintf(['\t\t' 'sl[%d], sl_sqr[%d], dsl[%d], slp[%d], slp_sqr[%d], muG[%d],' '\n'],...
    sum(o.nc),sum(o.nc),sum(o.nc),sum(o.nc),sum(o.nc),sum(o.nc))];

if o.merit_function == 0
    optCode = [optCode, sprintf(['\t\t' 'mu[%d], dm[%d], mup[%d],' '\n' ...
        '\t\t' 'dsl_sqr = 0.0, gps_sqr = 0.0, dm_sqr = 0.0, rho_hat_tmp = 0.0,' '\n'],sum(o.nc),sum(o.nc),sum(o.nc))];
end

if (o.contractive || o.terminal)
    optCode = [optCode, sprintf(['\t\t' 'psi_N[1], psi_Np[1], dot_psi_N[' num2str(o.N*o.nu) '],' '\n'])];
end

optCode = [optCode, sprintf(['\t\t' 'rho = 0.0, rho_hat = 0.0,' '\n'])];
optCode = [optCode, sprintf(['\t\t' 'phi0 = 0.0, phit = 0.0, phi0_dot = 0.0,' '\n'])];
optCode = [optCode, sprintf(['\t\t' 't = 0.0, t_u = 0.0,' '\n'])];
optCode = [optCode, sprintf(['\t\t' 'du_sqr = 0.0, u_sqr = 0.0;' '\n' '\n'])];

% possibility to define a variable in the future that reuses
% computationally expensive expressions

%% initialization of the algorithm

[c, d] = generate_forward_simulation(o);
% this function generates det_x
code = [code, c];
data = [data, d];
if o.nw >0
    optCode = [optCode, sprintf(['\t' 'det_x(x0,u,w,x);' '\n'])];
else
    optCode = [optCode, sprintf(['\t' 'det_x(x0,u,x);' '\n'])];
end

[c, d, in] = generate_objective_gradient_oracle(o);
% this function generates det_J_and_dot_J and det_J
code = [code, c];
data = [data, d];
info.flops.it = falcopt.addFlops(info.flops.it, in.flops.it);
info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops.ls);

c_w_use = argument_w(o,false);
c_tr_use = argument_def(o,false);
c_psi_use = argument_def_internal_psi(o,false);
c_psi_dot_use = argument_def_internal_psi_dot(o,false);
c_contr_use = argument_contr_value(o,false);

if (o.contractive || o.terminal)
    optCode = [optCode, sprintf(['\t' 'det_J_and_dot_J(x0, u, x' c_w_use c_tr_use ', &J, dot_J' c_psi_use c_psi_dot_use ');' '\n\n'])];
end

%     for jj=1:length(o.K_nc)
%         [d, c, in] = falcopt.generateMVMult({eye(o.nc(1)), 0.5*eye(o.nc(1))}, ...
%         'names', struct('fun', ['half_sum_nc_' num2str(jj) ], 'M', {{'I', 'hI'}},...
%         'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', real, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
%
%         if ~isempty(d)
%             data = [data, d, sprintf('\n')];
%         end
%         code = [code, c, sprintf('\n\n')];
%         in.flops = falcopt.multFlops( in.flops, length(o.K_nc{jj}));
%         info.flops.once = falcopt.addFlops( info.flops.once, in.flops); % in initialize_slack
%         info.flops.ls = falcopt.addFlops( info.flops.ls, in.flops); % in build_gpsl
%     end

[c, d, in] = generate_slack_initialization(o);

% this function initializes the slack variables
code = [code, c];
data = [data, d];
info.flops.once = falcopt.addFlops(info.flops.once, in.flops);

optCode = [optCode, sprintf(['\t' 'initialize_slack(u' c_psi_use c_contr_use ', sl, sl_sqr, gps' ');' '\n\n'])];





%% here we loop over the iterations it

optCode = [optCode, sprintf(['\t' 'for (it=0; it < %i; it++) {' '\n'], o.maxIt)];

if (o.contractive || o.terminal)
    optCode = [optCode, sprintf(['\t\t' 'if (it > 0) {' '\n'...
        '\t\t\t' 'det_J_and_dot_J(x0, u, x' c_w_use c_tr_use ', &J, dot_J' c_psi_use c_psi_dot_use '); }' '\n\n'])];
else
    optCode = [optCode, sprintf(['\t\t' 'det_J_and_dot_J(x0, u, x' c_w_use c_tr_use ', &J, dot_J' c_psi_use c_psi_dot_use ');' '\n\n'])];
end



% generate functions: "build_sqr_Nnc" and "build_gpsl"
[c, d, in] = generate_build_gps(o);

code = [code, c];
data = [data, d];
info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);

[c, d, in] = generate_dot_product_Nnu(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.multFlops( in.flops, 2); % called twice
info.flops.it = falcopt.addFlops( info.flops.it, in.flops);

[c, d, in] = generate_gradient_step(o);
code = [code, c];
data = [data, d];
info.flops.it = falcopt.addFlops(info.flops.it, in.flops);

optCode = [optCode, sprintf(['\t\t' 'gradient_step(dot_J, u, sl, sl_sqr, gps' c_psi_dot_use ', du, dsl, muG);' '\n\n'])];
    
% optCode = [optCode, sprintf(['if (it==2) {\n copy_Nnu(u,muG+30); return;}; \n'])]; % DEBUG

% optCode = [optCode, sprintf(['copy_Nnu(u,muG); \n return; \n'])]; % DEBUG

optCode = [optCode, sprintf(['\t\t' 'dot_product_Nnu(&du_sqr, du, du);' '\n',...
    '\t\t' 'dot_product_Nnu(&u_sqr, u, u);' '\n'])];

[c, d, in] = generate_dot_product_Nnc(o);
code = [code, c];
data = [data, d];
if o.merit_function == 0
    optCode = [optCode, sprintf(['\t\t' 'dot_product_Nnc(&gps_sqr, gps, gps);' '\n',...
        '\t\t' 'dot_product_Nnc(&dsl_sqr, dsl, dsl);' '\n'])];
    in.flops = falcopt.multFlops( in.flops, 3);
    info.flops.it = falcopt.addFlops( info.flops.it, in.flops);
end

if o.merit_function == 0
    optCode = [optCode, sprintf(['\t\t' 'if (it==0)' '\n',...
        '\t\t\t' 'copy_Nnc(mu, muG);' '\n'])];
    info.flops.it.comp = info.flops.it.comp +1;
    
    [c, d, in] = generate_difference(o);
    % generate diff_Nnc
    code = [code, c];
    data = [data, d];
    optCode = [optCode, sprintf(['\t\t' 'diff_Nnc(dm,muG,mu);' '\n\n'])];
    info.flops.it = falcopt.addFlops(info.flops.it, in.flops);
    
    [c, d, in] = generate_det_phi(o);
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);
    info.flops.it = falcopt.addFlops(info.flops.it, in.flops); % inside condition
    
    [c, d, in] = generate_det_dot_phi(o);
    code = [code, c];
    data = [data, d];
    info.flops.it = falcopt.addFlops(info.flops.it, in.flops);
    info.flops.it = falcopt.addFlops(info.flops.it, in.flops); % inside condition
    
    [c, d, in] = generate_conditions_rho(o);
    code = [code, c];
    data = [data, d];
    info.flops.it = falcopt.addFlops(info.flops.it, in.flops);
    
else
    
    [c,d,in] = generate_Lagrangian_oracles_lp(o);
    % generate det_dot_phi, det_phi, condition_rho
    code = [code, c];
    data = [data, d];
    info.flops.it = falcopt.addFlops(info.flops.it, in.flops.it);
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops.ls);
end

%% check conditions on the penalty parameter

if o.merit_function == 0
    
    
    optCode = [optCode, sprintf(['\t\t' 'det_dot_phi (du,dot_J, rho, gps, mu, dm, &phi0_dot);' '\n'])];
    
    
    optCode = [optCode, sprintf(['\t\t' 'if (conditions_rho_PM_simpler(phi0_dot,du_sqr,dsl_sqr,alpha)==0){' '\n',...
        '\t\t\t' 'dot_product_Nnc(&dm_sqr,dm,dm);' '\n',...
        '\t\t\t' 'rho_hat_tmp = dm_sqr / gps_sqr;' '\n'])];
    optCode = [optCode, sprintf(['\t\t\t' 'rho_hat = 2.0 * ' o.sqrt '(rho_hat_tmp);' '\n'])];
    optCode = [optCode, sprintf(['\t\t\t' 'rho = ' o.max '(2.0*rho,rho_hat);' '\n'])];
    optCode = [optCode, sprintf(['\t\t\t' 'flag_ini = 1;' '\n'])];
    optCode = [optCode, sprintf(['\t\t\t' 'det_dot_phi (du,dot_J, rho, gps, mu, dm, &phi0_dot);' '\n'])];
    optCode = [optCode, sprintf(['\t\t' '}' '\n'])];
    
    
    optCode = [optCode, sprintf(['\t\t' 'if ((flag_ini == 1)||(it == 0)){' '\n'])];
    optCode = [optCode, sprintf(['\t\t\t' 'det_phi (J, gps, mu, rho, &phi0);' '\n'])];
    optCode = [optCode, sprintf(['\t\t\t' 'flag_ini = 0;' '\n'])];
    optCode = [optCode, sprintf(['\t\t' '}' '\n\n'])];
    
    info.flops.it.comp = info.flops.it.comp + 4;  % 3+ 1 in fmax
    info.flops.it.mul = info.flops.it.mul +3;
    info.flops.it.sqrt = info.flops.it.sqrt +1;
    info.flops.it.div = info.flops.it.div +1;
    
else
    
    optCode = [optCode, sprintf(['\t\t' 'update_rho(muG, &rho, &rho_hat);' '\n',...
        '\t\t' 'det_phi (J, gps, rho, &phi0);' '\n',...
        '\t\t' 'det_dot_phi (du,dot_J, rho, gps, &phi0_dot);' '\n\n'])];
    
end

%% start line search

if o.merit_function == 0
    
    optCode = [optCode, sprintf(['\t\t' 'if (phi0_dot <= ' falcopt.num2str(-alpha_eps2, o.precision) ') {' '\n',...
        '\t\t\t' 't = 1.0;' '\n',...
        '\t\t\t' 't_u = 1.0;' '\n',...
        '\t\t\t' 'for (it_ls = 0; it_ls<%d; it_ls++) {' '\n'], o.maxItLs)];
    info.flops.it.comp = info.flops.it.comp +1;
    
    [c,d,in] = generate_weighted_sum_nc(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);
    
    [c,d,in] = generate_weighted_sum_nu(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);
    if o.merit_function == 0
        info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops); % Used twice
    end
    
    [c,d,in] = generate_quadratic_interpolation(o);
    
    code = [code, c, sprintf('\n\n')];
    data = [data, d];
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);
    
    
    optCode = [optCode, sprintf(['\t\t\t\t' 'weighted_sum_Nnu(up,du,u,&t);' '\n',...
        '\t\t\t\t' 'weighted_sum_Nnc(slp,dsl,sl,&t);' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\t\t\t\t' 'weighted_sum_Nnc(mup,dm,mu,&t);' '\n\n'])];
    end
    
    if o.nw > 0
        optCode = [optCode, sprintf([...
            '\t\t\t\t' 'det_x(x0,up,w,xp);' '\n'])];
    else
        optCode = [optCode, sprintf([...
            '\t\t\t\t' 'det_x(x0,up,xp);' '\n'])];
    end
    
    c_psip_use = argument_def_internal_psi_plus(o,false);
    
    
    optCode = [optCode, sprintf(['\t\t\t\t' 'det_J(x0, up, xp' c_tr_use ', &Jt' c_psip_use ');' '\n\n'])];
    
    optCode = [optCode, sprintf([...
        '\t\t\t\t' 'build_sqr_Nnc(slp, slp_sqr);' '\n',...
        '\t\t\t\t' 'build_gpsl(up' c_psip_use c_contr_use ',slp_sqr, gpsp);' '\n'])];
    
%     optCode = [optCode, sprintf(['if (it==1) {copy_Nnu(u,gpsp + 30); \n return;}; \n'])]; % DEBUG

    
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\t\t\t\t' 'det_phi (Jt,gpsp,mup,rho,&phit);' '\n'])];
    else
        optCode = [optCode, sprintf(['\t\t\t\t' 'det_phi (Jt,gpsp,rho,&phit);' '\n'])];
    end
    
    
    optCode = [optCode, sprintf(['\t\t\t\t' 'if (phit - phi0 <= ' falcopt.num2str(o.parLs, o.precision) '*t*phi0_dot)' '\n',...
        '\t\t\t\t\t' 'break;' '\n',...
        '\t\t\t\t' 'else {' '\n',...
        '\t\t\t\t\t' 't_u = t;' '\n',...
        '\t\t\t\t\t' 't = quadratic_interp (phi0, phi0_dot, t_u, phit);' '\n',...
        '\t\t\t\t' '}' '\n',...
        '\t\t\t' '}' '\n',...
        '\t\t\t' 'if (t_u <= ' falcopt.num2str(o.tolLs, o.precision) ')' '\n',...
        '\t\t\t\t' 'cond_err = 0;' '\n',...
        '\t\t' '}' '\n',...
        '\t\t' 'else' '\n',...
        '\t\t\t' 'conditions_f = 0;' '\n\n'])];
    info.flops.ls.comp = info.flops.ls.comp +1;
    info.flops.ls.mul = info.flops.ls.mul +2;
    info.flops.it.comp = info.flops.it.comp +1;
    
else
    %% recompute rho
    
    optCode = [optCode, sprintf(['\t\t' 'if (phi0_dot <= ' falcopt.num2str(-alpha_eps2, o.precision) ') {' '\n',...
        '\t\t\t' 't = 1.0;' '\n',...
        '\t\t\t' 't_u = 1.0;' '\n',...
        '\t\t\t' 'for (it_ls = 0; it_ls<%d; it_ls++) {' '\n'], o.maxItLs)];
    info.flops.it.comp = info.flops.it.comp +1;
    
    [c,d,in] = generate_weighted_sum_nc(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);
    
    [c,d,in] = generate_weighted_sum_nu(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);
    if o.merit_function == 0
        info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops); % Used twice
    end
    
    [c,d,in] = generate_quadratic_interpolation(o);
    
    code = [code, c, sprintf('\n\n')];
    data = [data, d];
    info.flops.ls = falcopt.addFlops(info.flops.ls, in.flops);
    
    
    optCode = [optCode, sprintf(['\t\t\t\t' 'weighted_sum_Nnu(up,du,u,&t);' '\n',...
        '\t\t\t\t' 'weighted_sum_Nnc(slp,dsl,sl,&t);' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\t\t\t\t' 'weighted_sum_Nnc(mup,dm,mu,&t);' '\n\n'])];
    end
    
    if o.nw > 0
        optCode = [optCode, sprintf([...
            '\t\t\t\t' 'det_x(x0,up,w,xp);' '\n'])];
    else
        optCode = [optCode, sprintf([...
            '\t\t\t\t' 'det_x(x0,up,xp);' '\n'])];
    end
    
    c_psip_use = argument_def_internal_psi_plus(o,false);
    
    
    optCode = [optCode, sprintf(['\t\t\t\t' 'det_J(x0, up, xp' c_tr_use ', &Jt' c_psip_use ');' '\n\n'])];
    
    optCode = [optCode, sprintf([...
        '\t\t\t\t' 'build_sqr_Nnc(slp, slp_sqr);' '\n',...
        '\t\t\t\t' 'build_gpsl(up' c_psip_use c_contr_use ',slp_sqr, gpsp);' '\n'])];
    
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\t\t\t\t' 'det_phi (Jt,gpsp,mup,rho,&phit);' '\n'])];
    else
        optCode = [optCode, sprintf(['\t\t\t\t' 'det_phi (Jt,gpsp,rho,&phit);' '\n'])];
    end
    
    optCode = [optCode, sprintf(['\t\t\t\t' 'if (phit - phi0 <= ' falcopt.num2str(o.parLs, o.precision) '*t*phi0_dot)' '\n',...
        '\t\t\t\t\t' 'break;' '\n',...
        '\t\t\t\t' 'else {' '\n',...
        '\t\t\t\t\t' 't_u = t;' '\n',...
        '\t\t\t\t\t' 't = quadratic_interp (phi0, phi0_dot, t_u, phit);' '\n',...
        '\t\t\t\t' '}' '\n',...
        '\t\t\t\t' 'if (t_u <= ' falcopt.num2str(o.tolLs, o.precision) ')' '\n'...
        '\t\t\t\t\t' 'reset_rho = 1;' '\n'...
        '\t\t\t' '}' '\n'])];
    
    optCode = [optCode, sprintf([...
        '\t\t\t' '/* Reset rho to rho_hat if rho is overly big */' '\n'...
        '\t\t\t' 'if (reset_rho == 1){' '\n'...
        '\t\t\t\t' 'reset_rho = 0;' '\n'...
        '\t\t\t\t' 'rho = rho_hat;' '\n'...
        '\t\t\t\t' 'det_phi (J, gps, rho, &phi0);' '\n',...
        '\t\t\t\t' 'det_dot_phi (du,dot_J, rho, gps, &phi0_dot);' '\n\n',...
        ...
        '\t\t\t\t' 't = 1.0;' '\n',...
        '\t\t\t\t' 't_u = 1.0;' '\n',...
        '\t\t\t\t' 'for (it_ls = 0; it_ls<%d; it_ls++) {' '\n'], o.maxItLs)];
    info.flops.it.comp = info.flops.it.comp +1;
    
    optCode = [optCode, sprintf(['\t\t\t\t\t' 'weighted_sum_Nnu(up,du,u,&t);' '\n',...
        '\t\t\t\t\t' 'weighted_sum_Nnc(slp,dsl,sl,&t);' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\t\t\t\t' 'weighted_sum_Nnc(mup,dm,mu,&t);' '\n\n'])];
    end
    
    if o.nw > 0
        optCode = [optCode, sprintf([...
            '\t\t\t\t\t' 'det_x(x0,up,w,xp);' '\n'])];
    else
        optCode = [optCode, sprintf([...
            '\t\t\t\t\t' 'det_x(x0,up,xp);' '\n'])];
    end
    
    c_psip_use = argument_def_internal_psi_plus(o,false);
    
    
    optCode = [optCode, sprintf(['\t\t\t\t\t' 'det_J(x0, up, xp' c_tr_use ', &Jt' c_psip_use ');' '\n\n'])];
    
    optCode = [optCode, sprintf([...
        '\t\t\t\t\t' 'build_sqr_Nnc(slp, slp_sqr);' '\n',...
        '\t\t\t\t\t' 'build_gpsl(up' c_psip_use c_contr_use ',slp_sqr, gpsp);' '\n'])];
    
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\t\t\t\t\t' 'det_phi (Jt,gpsp,mup,rho,&phit);' '\n'])];
    else
        optCode = [optCode, sprintf(['\t\t\t\t\t' 'det_phi (Jt,gpsp,rho,&phit);' '\n'])];
    end
    
    optCode = [optCode, sprintf(['\t\t\t\t\t' 'if (phit - phi0 <= ' falcopt.num2str(o.parLs, o.precision) '*t*phi0_dot)' '\n',...
        '\t\t\t\t\t\t' 'break;' '\n',...
        '\t\t\t\t\t' 'else {' '\n',...
        '\t\t\t\t\t\t' 't_u = t;' '\n',...
        '\t\t\t\t\t\t' 't = quadratic_interp (phi0, phi0_dot, t_u, phit);' '\n',...
        '\t\t\t\t\t' '}' '\n',...
        '\t\t\t\t\t' 'if (t_u <= ' falcopt.num2str(o.tolLs, o.precision) ')' '\n',...
        '\t\t\t\t\t\t' 'cond_err = 0;' '\n',...
        '\t\t\t\t' '}' '\n',...
        '\t\t\t' '}' '\n',...
        '\t\t' '}' '\n',...
        '\t\t' 'else' '\n',...
        '\t\t\t' 'conditions_f = 0;' '\n\n'])];
    info.flops.ls.comp = info.flops.ls.comp +1;
    info.flops.ls.mul = info.flops.ls.mul +2;
    info.flops.it.comp = info.flops.it.comp +1;
    
end

%% the code continues

optCode = [optCode, sprintf(['\t\t' 'iter_ls[it] = it_ls+1;','\n',...
    '\t\t' 'if (it_ls== %d)' '\n',...
    '\t\t\t' 'cond_err = 0;' '\n\n'],o.maxItLs)];
info.flops.once.add = info.flops.once.add +1;
%         optCode = [optCode, sprintf(['if (it==1)

%      optCode = [optCode, sprintf(['if (cond_err == 0){' '\n' ...
%          '\t' 'u[0] = rho;\n' ...
%          '\t' 'u[1] = phi0_dot; \n' ...
%          '\t' 'return cond;}' '\n'])];

%% update u, slack and dual variables

if o.merit_function == 0
    optCode = [optCode, sprintf([ '\t\t' 'copy_Nnc(mu,mup);' '\n'])];
end

optCode = [optCode, sprintf(['\t\t' 'copy_Nnu(u,up);' '\n',...
    '\t\t' 'copy_Nnc(sl,slp);' '\n',...
    '\t\t' 'copy_Nnc(sl_sqr,slp_sqr);' '\n',...
    '\t\t' 'copy_Nnc(gps,gpsp);' '\n',...
    '\t\t' 'copy_Nnx(x,xp);' '\n',...
    '\t\t' 'J = Jt;' '\n',...
    '\t\t' 'phi0 = phit;' '\n\n'])];

%% check convergence

[c, d, in] = generate_compute_max(o);
info.flops.it = falcopt.addFlops(info.flops.it, in.flops);
% compute_max_Nnc

code = [code, c, sprintf('\n\n')];
data = [data, d];

optCode = [optCode,sprintf(['\t\t' 'if ((du_sqr >= ' falcopt.num2str(alpha2_eps2, o.precision) ')||(compute_max_Nnc(gpsp) >= ' falcopt.num2str(o.eps, o.precision) '))' '\n'])];
optCode = [optCode,sprintf(['\t\t\t' 'conditions_x = 1;' '\n',...
    '\t\t' 'else' '\n',...
    '\t\t\t' 'conditions_x = 0;' '\n'])];


optCode = [optCode,sprintf(['\t\t' 'if (it == %d)' '\n'],o.maxIt-1)];
optCode = [optCode,sprintf(['\t\t\t' 'conditions_n = 0;' '\n'])];


optCode = [optCode,sprintf(['\t\t' 'if(!((conditions_f && cond_err) && (conditions_n && conditions_x)))' '\n'])];
optCode = [optCode,sprintf(['\t\t\t' 'break;' '\n\n',...
    '\t' '}' '\n\n'])];
info.flops.it.comp = info.flops.it.comp +3;

optCode = [optCode,sprintf(['\t' '(*iter) = it+1;' '\n',...
    '\t' 'if (conditions_f == 0)' '\n',...
    '\t\t' 'cond = 2;' '\n',...
    '\t' 'else {' '\n',...
    '\t\t' 'if (conditions_x == 0)' '\n',...
    '\t\t\t' 'cond = 1;' '\n',...
    '\t\t' 'else {' '\n',...
    '\t\t\t' 'if (conditions_n == 0)' '\n',...
    '\t\t\t\t' 'cond = 0;' '\n',...
    '\t\t\t' 'else' '\n',...
    '\t\t\t\t' 'cond = -1;' '\n'...
    '\t\t' '}' '\n'...
    '\t' '}' '\n\n'])];
if o.debug > 1
    optCode = [optCode, sprintf(['\t' '*fval = J;' '\n'])];
end
optCode = [optCode, sprintf(['\t' 'return cond;' '\n'])];
optCode = [optCode, sprintf('} \n\n')];
info.flops.once.add = info.flops.once.add +1;
info.flops.once.comp = info.flops.once.comp +3;

%% generate MEX

if o.build_MEX
    mexName = o.name;
    mexCode = falcopt.generateMEX(o.N, struct('x',o.nx,'u',o.nu,'w',o.nw), 'names',...
        struct('fun','proposed_algorithm','mex',mexName), 'ref', o.trackReference,...
        'terminalContraction', o.terminal || o.contractive, 'indexContraction',...
        o.contractive, 'type', o.real,'maxIt',o.maxIt, 'timing', true, 'debug', o.debug,...
        'verbose', o.verbose, 'indent', o.indent, 'inline', o.inline, 'precision', o.precision);
    final_code = [libr '\n' data '\n'  code '\n' optCode '\n' mexCode];
else
    final_code = [libr '\n' data '\n'  code '\n' optCode];
end

if ~isempty(o.gendir)
    if (o.gendir(end) == '/')||(o.gendir(end) == '\')
        o.gendir = o.gendir(1:end-1);
    end
    file_folder = o.gendir;
    if exist(file_folder, 'dir')~=7
        mkdir(file_folder);
    end
else
    file_folder = '';
end


filename = [file_folder '/' o.name '.c'];
ext_file = [];
for k = 1:length(info.src)
    ext_file = [ext_file, ' ', info.src{k}];
end
f = fopen(filename, 'w+');
fprintf(f, final_code);
fclose(f);


% if ~isempty(o.gendir)
%     cd ..
% end

if o.compile
    if o.build_MEX
        compile = ['mex ' filename ext_file ' -v -output ' mexName];
    else
        throw(MException('activate:build_MEX', 'The option build_MEX must be set to true to compile the code'));
    end
    disp(compile);
    eval(compile);
end



end

function [ data, code, info] = casadi_jacobians(o,f,fName,staticName, varargin)
% given function f, jacobians are automatically computed using casadi and c-code
% is generated
%
% Inputs:
% - o: options
% - f: function handle
% - staticName: name given at C variable
% - varargin: can be 'x', 'u' or 'xu' and indicate what jacobian should be generated
%
% Outputs:
% - data
% - code
% - sxfcn: casadi functions
% - info: containg structure of function and flops
import casadi.*
x = SX.sym('x',o.nx);
u = SX.sym('u',o.nu);
w = SX.sym('w',o.nw);

data = [];
code = [];
sxfcn = {};


% function
if( nargin(f) == 3)
    y = f(x,u,w);
elseif( nargin(f) == 2)
    y = f(x,u);
else
    if( strcmp(fName,'model_mpc'))
        error('Check number of dynamics inputs, can be (x,u) or (x,u,w)');
    else
        error(['Check number of ' fName ' inputs, can be (x,u) or (x,u,w)']);
    end
end

if( isnumeric(y))   % generate static data if is not symbolic
    name.M = staticName;
    [d, ~, in_y] = falcopt.generateData(y, 'names', name, ...
        'type', o.real, 'precision', o.precision, 'structure', 'unique', 'noones', false, 'indent', o.indent, ...
        'static', true, 'const', true, 'verbose', o.verbose);
    data = [data, d, sprintf('\n')];
    in_y.static = 1;
    info.y = in_y;
    info.y.flops = 0;
else
    [~,in_y] = falcopt.casadi2struct( y);
    y_f = Function([fName '_casadi'],{x,u,w},{in_y.stored.values});
    sxfcn{1} = y_f;
    try
        info.y.flops =  y_f.getAlgorithmSize(); %flops
    catch
        warning('Cannot use casadi.Function.getAlgoirthmSize()');
    end
    info.y.static = 0;
    
    % wrapper function
    code = [code, sprintf('/* Dynamics of the system */\n')];
    if( o.nw > 0)
        code = [code, sprintf(['void ' fName '( const ' o.real '* x, const ' o.real '* u, const ' o.real '* v, ' o.real '* xp){' '\n\n'])];
        code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],3)];
        code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)]; %TODO Tommaso check parameter
        code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\t' 'in[2] = v;\n\n'])];
    else
        code = [code, sprintf(['void ' fName '( const ' o.real '* x, const ' o.real '* u, ' o.real '* xp){' '\n\n'])];
        code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],2)];
        code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)]; %TODO Tommaso check parameter
        code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\n'])];
    end
    
    code = [code, sprintf(['\t' fName '_casadi( in, &xp, &iw, w, mem);\t/* external generated casadi function*/\n}\n\n'])];
    if( ~isempty(staticName))
        data = [data, sprintf(['\t' 'static ' o.real ' ' staticName '[%i];' '\n'], in_y.stored.num)];
    end
    info.y.structure = in_y;
end

% jacobians
for k = 1:length(varargin)
    switch varargin{k}
        case 'x'
            jac = jacobian( y, x);
            static_name = 'F';
            struct_name = 'in_F';
            J_name = 'jacobian_x_casadi';
        case 'u'
            jac = jacobian( y, u);
            static_name = 'G';
            struct_name = 'in_G';
            J_name = 'jacobian_u_casadi';
        otherwise
            error('invalid variable name');
    end
    [info.(struct_name).static,in] = falcopt.casadi2struct(jac);
    
    if( info.(struct_name).static) %is jac matrix constant ?
        name.M = static_name;
        data = [data, sprintf('\t/* Static data for jacobian w.r.t %c variable */\n',varargin{k})];
        [d, ~, in_d] = falcopt.generateData(full(DM(jac)), 'names', name, ...
            'type', o.real, 'precision', o.precision, 'structure', 'unique', 'noones', false, 'indent', o.indent, ...
            'static', true, 'const', true, 'verbose', o.verbose);
        info.(struct_name).struct = in_d;
        info.(struct_name).flops = 0;
        if ~isempty(d)
            data = [data, d, sprintf('\n')];
        end
        
    else
        jac = Function(J_name,{x,u,w},{in.stored.values});
        sxfcn{length(sxfcn)+1} =  jac;
        
        % wrapper jacobian
        code = [code, sprintf('/* System dynmics jacobian w.r.t %c variable */\n',varargin{k})];
        if( o.nw > 0)
            code = [code, sprintf(['void Jacobian_' varargin{k} '( const ' o.real '* x, const ' o.real '* u, const ' o.real '* v, ' o.real '* res){' '\n\n'])];
            code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],3)];
            code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)];
            code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\t' 'in[2] = v;\n\n' ...
                        '\t' J_name '( in, &res, &iw, w, mem);\t/* external generated casadi function*/\n}\n\n']) ];
        else
            code = [code, sprintf(['void Jacobian_' varargin{k} '( const ' o.real '* x, const ' o.real '* u, ' o.real '* res){' '\n\n'])];
            code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],2)];
            code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)];
            code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\n' ...
                        '\t' J_name '( in, &res, &iw, w, mem);\t/* external generated casadi function*/\n}\n\n']) ];
        end
        
        data = [data, sprintf('\t/* Static data for jacobian w.r.t %c variable */\n',varargin{k})];
        data = [data, sprintf(['\t' 'static ' o.real ' ' static_name '[%i];' '\n'], in.stored.num)];
        info.(struct_name).struct.structure = in;
        try
            info.(struct_name).flops =   jac.getAlgorithmSize(); %flops
        catch
            warning('Cannot use casadi.Function.getAlgoirthmSize()');
        end
    end
end

% generate .c file with casadi functions
[info.src,info.header] = generate_casadi_c(o, 'casadi_fcn', sxfcn);

end

function K = detect_structure ( ind, o )
% ind is a matrix of bools of dims (nu, N)

n = size(ind,2);
indices = nan (1,n);
i = 1;
K = {};

while any(isnan(indices))
    k = find(isnan(indices),1,'first');
    I = all(repmat(ind(:,k),1,n) == ind ,1);
    A = find(I);
    indices(A) = 0;
    if any(sum(~isinf(ind(:,A)),1) >= 1)
        K{i} = A(sum(~isinf(ind(:,A)),1) >= 1); %#ok
    end
    i = i+1;
end

% if any(size(K))==0
%     K = {};
% end

end

function K = detect_different_NLconstraints(o)

indices = nan (1,o.N);
i = 1;
while any(isnan(indices))
    K{i} = [];
    k_all = find(isnan(indices));
    k = k_all(1);
    %k = find(isnan(indices),1,'first');
    for j = k_all
        if isequal(o.constraints_handle{k},o.constraints_handle{j})
            K{i} = [K{i},j];
        end
    end
    indices(K{i}) = 0;
    i = i+1;
end

end

function K = detect_structure_constraints( ind ,o )
% ind is a vector of bools of dims (1, N+1)

ind = ind(1:o.N);

indices = nan (1,o.N);
i = 1;
while any(isnan(indices))
    k = find(isnan(indices),1,'first');
    I = repmat(ind(:,k),1,o.N) == ind;
    K{i} = find(I);
    indices(K{i}) = 0;
    i = i+1;
end

end

function [src,header] = generate_casadi_c( o, fileName, fcns)
% given functions 'fcns' generate c code using casadi in the give 'fileName'
% returns name of .c src file and .h header file
import casadi.*

opts = struct('with_header',true,'real_t',o.real); % options for code generation




C = CodeGenerator ( fileName,opts );

for k = 1:length(fcns)
    if( ~isempty(fcns{k}))
        C.add(fcns{k})
    end
end


if ~isempty(o.gendir)
    % create folder
    if (o.gendir(end) == '/')||(o.gendir(end) == '\')
        o.gendir = o.gendir(1:end-1);
    end
    file_folder = o.gendir;
    if exist(file_folder, 'dir')~=7
        mkdir(file_folder);
    end
    cd (sprintf(o.gendir));
    C.generate( ) ;
    cd ..;
    src = { sprintf([o.gendir '/%s.c'],fileName)};
    header = { sprintf('%s.h',fileName)};
else
    C.generate();
    src = { sprintf('%s.c',fileName)};
    header = { sprintf('%s.h',fileName)};
end
end

function [ data, code, info] = matlab_jacobians(o,f,fName,varargin)
% generate c code from provided function f and automatically generate jacobians
% using matlab symbolic toolbox
code = [];
data = [];
info.y = [];
info.in_F = [];
info.in_G = [];
x = sym('x',[o.nx,1]);
u = sym('u',[o.nu,1]);
w = sym('w',[o.nw,1]);


if( nargin(f) == 3)
    y = f(x,u,w);
elseif( nargin(f) == 2)
    y = f(x,u);
else
    if( strcmp(fName,'model_mpc'))
        error('Check number of dynamics inputs, can be (x,u) or (x,u,w)');
    else
        error(['Check number of ' fName ' inputs, can be (x,u) or (x,u,w)']);
    end
end

% function
[d, in_y] = falcopt.fcn2struct( y,o,'name','xp');

code = [ code, sprintf('/*Dynamics of the system*/\n')];
if( o.nw > 0)
    code = [code, sprintf(['void ' fName '(const double* x, const double* u, const double* w, double* xp){\n\n'])];
else
    code = [code, sprintf(['void ' fName '(const double* x, const double* u, double* xp){\n\n'])];
end
code = [ code, d, sprintf('}\n\n')];
info.y.flops = in_y.flops;


% jacobians
for k = 1:length(varargin)
    switch varargin{k}
        case 'x'
            jac = jacobian( y, x);
            static_name = 'F';
            struct_name = 'in_F';
            J_name = 'Jacobian_x';
        case 'u'
            jac = jacobian( y, u);
            static_name = 'G';
            struct_name = 'in_G';
            J_name = 'Jacobian_u';
        otherwise
            error('invalid variable name');
    end
    [d, i] = falcopt.fcn2struct( jac,o,'name', static_name);
    info.(struct_name).static = i.static;
    data = [ data, sprintf([o.indent.data '/*Static data for Jacobian w.r.t %c*/\n'],varargin{k})];
    if( info.(struct_name).static)
        data = [data, sprintf([o.indent.data 'static const ' o.real ' ' static_name '[%i] = {\n'], i.structure.num)];
        data = [data, d, sprintf([o.indent.data '};\n\n'])];
    else
        data = [data, sprintf([o.indent.data 'static ' o.real ' ' static_name '[%i];\n'], i.structure.num)];
        code = [ code, sprintf([o.indent.data '/*Jacobian w.r.t %c*/\n'],varargin{k})];
        if( o.nw > 0)
            code = [code, sprintf([o.indent.code 'static ' o.inline ' void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, const ' o.real ' * w, ' o.real ' * ' static_name ') {\n\n'])];
        else
            code = [code, sprintf([o.indent.code 'static ' o.inline ' void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, ' o.real ' * ' static_name ') {\n\n'])];
        end
        code = [code, d, sprintf([o.indent.code '}\n\n'])];
    end
    info.(struct_name).struct.structure = i.structure;
    info.(struct_name).flops = i.flops;
end
end


function [ data, code, info] = external_jacobians(o, varargin)
% generate c-code from provided matlab function for model_mpc, jacobian_x, jacobian_u
code = [];
data = [];
info.y = [];
info.in_F = [];
info.in_G = [];
x = sym('x',[o.nx,1]);
u = sym('u',[o.nu,1]);
w = sym('w',[o.nw,1]);

if( o.nw>0)
    y = o.dynamics(x,u,w);
    jac_u = o.external_jacobian_u(x,u,w);
    jac_x = o.external_jacobian_x(x,u,w);
else
    y = o.dynamics(x,u);
    jac_u = o.external_jacobian_u(x,u);
    jac_x = o.external_jacobian_x(x,u);
end

% dynamics of system
[d, in_y] = falcopt.fcn2struct( y,o,'name','xp');
code = [ code, sprintf('/*Dynamics of the system*/\n')];
if( o.nw > 0)
    code = [code, sprintf(['void model_mpc(const double* x, const double* u, const double* w, double* xp){\n\n'])];
else
    code = [code, sprintf(['void model_mpc(const double* x, const double* u, double* xp){\n\n'])];
end
code = [ code, d, sprintf('}\n\n')];
info.y.flops = in_y.flops;

% jacobians
for k = 1:2
    if(k==1)
        jac = jac_x;
        static_name = 'F';
        struct_name = 'in_F';
        J_name = 'Jacobian_x';
    else
        jac = jac_u;
        static_name = 'G';
        struct_name = 'in_G';
        J_name = 'Jacobian_u';
    end
    [d, i] = falcopt.fcn2struct( jac,o,'name', static_name);
    info.(struct_name).static = i.static;
    data = [ data, sprintf(['\t/*Static data for ' J_name '*/\n'])];
    if( info.(struct_name).static)
        data = [data, sprintf(['\tstatic const ' o.real ' ' static_name '[%i] = {\n'], i.structure.num)];
        data = [data, d, sprintf(['};\n\n'])];
    else
        data = [data, sprintf(['\tstatic ' o.real ' ' static_name '[%i];\n'], i.structure.num)];
        code = [ code, sprintf(['/* ' J_name '*/\n'])];
        if( o.nw > 0)
            code = [code, sprintf(['static inline void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, const ' o.real ' * w, ' o.real ' * ' static_name ') {\n\n'])];
        else
            code = [code, sprintf(['static inline void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, ' o.real ' * ' static_name ') {\n\n'])];
        end
        code = [code, d, sprintf(['}\n\n'])];
    end
    info.(struct_name).struct.structure = i.structure;
    info.(struct_name).flops = i.flops;
end
end

function [ data, code, info] = generate_n_and_Dn(o,grad)
code = [];
data = [];
info.nc = [];
info.in_amu = {};
info.in_umb = {};
info.in_n = {};
info.in_Dn_n = {};
sxfcn = {};

if ( strcmp(o.gradients,'casadi'))
    import casadi.*
end

% check if there are box constraints
% if( isempty(o.box_constraints))
%     ibox = 0;
% else
%     ibox = 1;
% end

if ~isempty(o.K_n)
    nl_con = 1;
else
    nl_con = 0;
end
%nl_con = ~isempty(o.K_n); % check if there are nonlinear constraints n(u)

% if( isa(o.constraints_handle,'function_handle'))
%     nl_con = 1;
% else
%     nl_con = 0;
% end

% build nonlinear constraints n(u)
switch grad
    case 'casadi'
        z = SX.sym('u',[o.nu,1]);
    case {'matlab','manual'}
        z = sym('u',[o.nu,1],'real');
end
if( nl_con)
    for jj=1:length(o.K_n)
        if( nargin(o.constraints_handle{jj}) == 1)
            try
                n{jj} = o.constraints_handle{jj}(z);
            catch
                try
                    n{jj} = o.constraints_handle{jj}(z');
                catch
                    error('inequality constraints can depend only on inputs' );
                end
            end
        else
            error('inequality constraints can depend only on input u' );
        end
    end
end



%% generate a - u functions
for ii=1:length(o.K_amu)
    cons_Bound = o.box_lowerBound(:,o.K_amu{ii}(1));
    cons_BoundTemp = diag(~isinf(cons_Bound));
    new_Identity = double(cons_BoundTemp(any(cons_BoundTemp,2),:));
    cons_Bound(isinf(cons_Bound)) = 0;
    new_Bound = new_Identity*cons_Bound;
    [d, c, in] = falcopt.generateMVMult(-new_Identity, new_Bound, ...
        'names', struct('fun', ['build_amu_' num2str(ii)], 'M', ['minusI_' num2str(ii)],...
        'm', ['umin_' num2str(ii)], 'v', 'u'), 'types', o.real, 'verbose', o.verbose,...
        'structure', 'ordered', 'test', o.test, 'inline', o.inline, 'indent', o.indent,...
        'precision', o.precision);
    
    info.amu{ii}.flops = in.flops;
    
    if ~isempty(d)
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    
end

% generate a selector that chooses which a - u to use
code = [code, sprintf(['\n' o.indent.code o.inline ' void build_amu(const ' o.real '* u, const unsigned int k, ' o.real '* amu){' '\n\n'])];
for ii=1:length(o.K_amu)
    check = generateCheck_custom(o.K_amu{ii}, [o.indent.code o.indent.generic], 'k');
    code = [code, check, sprintf('\n')]; %#ok
    code = [code, sprintf([o.indent.code o.indent.generic o.indent.generic 'build_amu_' num2str(ii) '(&amu[0], &u[0]);' '\n'])]; %#ok
    code = [code, sprintf([o.indent.code o.indent.generic '}' '\n'])]; %#ok
end
code = [code, sprintf([o.indent.code '}' '\n'])];

%% generate u - b functions
for ii=1:length(o.K_umb)
    cons_Bound = o.box_upperBound(:,o.K_umb{ii}(1));
    cons_BoundTemp = diag(~isinf(cons_Bound));
    new_Identity = double(cons_BoundTemp(any(cons_BoundTemp,2),:));
    cons_Bound(isinf(cons_Bound)) = 0;
    new_Bound = new_Identity*cons_Bound;
    [d, c, in] = falcopt.generateMVMult(new_Identity, -new_Bound, ...
        'names', struct('fun', ['build_umb_' num2str(ii)], 'M', ['I_' num2str(ii)],...
        'm', ['minus_umax_' num2str(ii)], 'v', 'u'), 'types', o.real, 'verbose', o.verbose,...
        'structure', 'ordered', 'test', o.test, 'inline', o.inline, 'indent', o.indent,...
        'precision', o.precision);
    
    info.umb{ii}.flops = in.flops;
    
    if ~isempty(d)
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
end

% generate a selector that chooses which u - b to use
code = [code, sprintf(['\n' o.indent.code o.inline ' void build_umb(const ' o.real '* u, const unsigned int k, ' o.real '* umb){' '\n\n'])];
for ii=1:length(o.K_umb)
    check = generateCheck_custom(o.K_umb{ii}, [o.indent.code '\t'], 'k');
    code = [code, check, sprintf('\n')]; %#ok
    code = [code, sprintf([o.indent.code '\t\t' 'build_umb_' num2str(ii) '(&umb[0], &u[0]);' '\n'])]; %#ok
    code = [code, sprintf([o.indent.code '\t' '}' '\n'])]; %#ok
end
code = [code, sprintf([o.indent.code '}' '\n'])];

% TO BE DELETED do not generate build_g anymore
% add variable n
% g = [];
% for jj = 1:o.N
%     g_add = [];
%     for ii = 1:o.nu
%         if ~isinf(o.box_lowerBound(ii,jj))
%             g_add = [g_add; o.box_lowerBound(ii,jj) - z];
%         end
%         if ~isinf(o.box_upperBound(ii,jj))
%             g_add = [g_add; z-o.box_upperBound(ii,jj)];
%         end
%         g = [g;g_add;n{jj}];
%     end
% end

% TO BE DELETED
% box = o.box_constraints;
% g = n;
% if( ibox)
%     for i = size(box,1):-1:1
%         if( i>o.nu)
%             g = [ box(i,1)-z(i-o.nu); z(i-o.nu)-box(i,2); g];
%         else
%             g = [ box(i,1)-z(i); z(i)-box(i,2); g];
%         end
%     end
% end

% build o.nc (constraints structure)
info.nc = cell2mat(o.nn);


if( o.terminal || o.contractive)
    info.nc = [info.nc 1];
else
    info.nc = [info.nc 0];
end

for i=1:o.N % not tested ToDo Tommaso
    lb = sum(~isinf(o.box_lowerBound(:,i)));
    ub = sum(~isinf(o.box_upperBound(:,i)));
    info.nc(i) = info.nc(i)+lb+ub;
end



% generate build_n and build_Dn
if( nl_con)
    switch grad
        case 'casadi' % use casadi generation
        
        %     % build_Dg
        %     Dg = jacobian( g, z);
        %     Dg = transpose(Dg);
        %     [info.in_Dg_g.static,i] = falcopt.casadi2struct(Dg);
        %     Dg_f = Function('build_Dg_casadi',{z},{i.stored.values});
        %     sxfcn{length(sxfcn)+1} =  Dg_f;
        %
        %     % wrapper build_Dg
        %     code = [code, sprintf('/* Jacobian_u of constraints */\n')];
        %     code = [code, sprintf(['void build_Dg( const ' o.real '* u, ' o.real '* Dg_fun){' '\n\n'])];
        %     code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],1)];
        %     code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)];
        %     code = [code, sprintf(['\tin[0] = u;\n\n' ...
        %         '\tbuild_Dg_casadi( in, &Dg_fun, &iw, w, mem);\t/* external casadi generated function*/\n}\n\n']) ];
        %     info.in_Dg_g.struct.structure = i;
        %     try
        %         info.in_Dg_g.flops =   Dg_f.getAlgorithmSize(); %flops
        %     catch
        %         warning('Cannot use casadi.Function.getAlgoirthmSize()');
        %     end
        sxfcn = {};
        
        if( nl_con)  % if exist n(u)
            for jj=1:length(o.K_n)
                
                % build_n
                [~,in_n{jj}] = falcopt.casadi2struct( n{jj});
                
                y_f{jj} = Function(['build_n_' num2str(jj) '_casadi'],{z},{in_n{jj}.stored.values});
                
                sxfcn{end + 1} = y_f{jj};
                try
                    info.in_n{jj}.flops =  y_f{jj}.getAlgorithmSize(); %flops
                catch
                    warning('Cannot use casadi.Function.getAlgoirthmSize()');
                end
                info.in_n{jj}.static = 0;
                info.in_n{jj}.struct.structure = in_n{jj};
                
                % wrapper function
                code = [code, sprintf(['/* Constraints evaluation*/' '\n'])];
                code = [code, sprintf(['void build_n_%d( const ' o.real '* z,'  o.real '* n){' '\n\n'],jj)];
                code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],1)];
                code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)];
                code = [code, sprintf(['\tin[0] = z;\n\n' ...
                    '\t' 'build_n_%d_casadi( in, &n, &iw, w, mem);\t/* external casadi generated function */\n}\n\n'],jj) ];
                
%                 % wrapper build_Dg
%                 code = [code, sprintf('/* Jacobian_u of constraints */\n')];
%                 code = [code, sprintf(['void build_Dg( const ' o.real '* u, ' o.real '* Dg_fun){' '\n\n'])];
%                 code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],1)];
%                 code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)];
%                 code = [code, sprintf(['\tin[0] = u;\n\n' ...
%                     '\tbuild_Dg_casadi( in, &Dg_fun, &iw, w, mem);\t/* external casadi generated function*/\n}\n\n']) ];
                
                % build_Dn
                Dn{jj} = jacobian( n{jj}, z);
                Dn{jj} = transpose(Dn{jj});
                [info.in_Dn_n{jj}.static,i] = falcopt.casadi2struct(Dn{jj});
                
                if( info.in_Dn_n{jj}.static) %is jac matrix constant ?
                    name.M = sprintf(['build_Dn_%d_casadi'],jj);
                    [d, ~, in] = falcopt.generateData(full(DM(Dn{jj})), 'names', name, ...
                        'type', o.real, 'precision', o.precision, 'structure', 'unique', 'noones', false, 'indent', o.indent, ...
                        'static', true, 'const', true, 'verbose', o.verbose);
                    info.in_Dn_n{jj}.struct = in;
                    info.in_Dn_n{jj}.flops = 0;
                    if ~isempty(d)
                        data = [data, d, sprintf('\n')];
                    end
                else
                    Dn_f{jj} = Function(sprintf(['build_Dn_%d_casadi'],jj),{z},{i.stored.values});
                    sxfcn{end + 1} =  Dn_f{jj};
                    
                    % wrapper build_Dn
                    code = [code, sprintf('/* Jacobian_u of nonlinear constraints */\n')];
                    code = [code, sprintf(['void build_Dn_%d( const ' o.real '* u, ' o.real '* Dn_fun){' '\n\n'],jj)];
                    code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],1)];
                    code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)];
                    code = [code, sprintf(['\tin[0] = u;\n\n' ...
                        '\t' 'build_Dn_%d_casadi( in, &Dn_fun, &iw, w, mem);\t/* external casadi generated function*/\n}\n\n'],jj) ];
                    info.in_Dn_n{jj}.struct.structure = i;
                    try
                        info.in_Dn_n{jj}.flops =   Dn_f{jj}.getAlgorithmSize(); %flops
                    catch
                        warning('Cannot use casadi.Function.getAlgoirthmSize()');
                    end
                end
            end
            
        end
        
        
        case {'manual','matlab'} % use matlab jacobian generation
        
        %
        %     %build_Dg
        %     Dg = jacobian( g,z);
        %     Dg = Dg';
        %     [c,i] = falcopt.fcn2struct( Dg, o, 'name', 'Dg_fun');
        %     code = [code, sprintf('/* Jacobian_u of constraints matrix*/\n')];
        %     code = [code, sprintf(['void build_Dg( const ' real '* u, ' real '* Dg_fun){' '\n'])];
        %     code = [code, c, sprintf('}\n\n')];
        %     info.in_Dg_g.struct.structure = i.structure;
        %     info.in_Dg_g.flops =   i.flops; %flops
        %     info.in_Dg_g.static = 0;
        %     sxfcn = {};
        
        
        
        for jj=1:length(o.K_n)
            
            [c, in_n{jj}] = falcopt.fcn2struct( n{jj},o,'name', 'n' );
            code = [code, sprintf('/* Constraints evaluation*/ \n')];
            code = [code, sprintf(['void build_n_' num2str(jj) '( const ' o.real '* u,'  o.real '* n){' '\n'])];
            code = [code, c, sprintf('}\n\n')];
            
            info.in_n{jj}.flops =  in_n{jj}.flops; %flops
            info.in_n{jj}.static = 0;
            info.in_n{jj}.struct.structure = in_n{jj}.structure;
            
            %build_Dn
            if(  isequal(o.gradients,'matlab'))
                Dn{jj} = jacobian( n{jj},z);
                Dn{jj} = Dn{jj}';
            else
                try
                    Dn{jj} = o.external_jacobian_n{jj}(z);
                catch
                    error('o.external_jacobian_n function not found. Or invalid output: it should return a struct of size o.N');
                end
            end
            [c,i] = falcopt.fcn2struct( Dn{jj}, o, 'name', 'Dn_fun');
            code = [code, sprintf('/* Jacobian_u of nonlinear constraints*/\n')];
            code = [code, sprintf(['void build_Dn_' num2str(jj) '( const ' o.real '* u, ' o.real '* Dn_fun){' '\n'])];
            code = [code, c, sprintf('}\n\n')];
            info.in_Dn_n{jj}.struct.structure = i.structure;
            info.in_Dn_n{jj}.flops =   i.flops; %flops
            info.in_Dn_n{jj}.static = 0;
        end
        
        
        
        sxfcn = {};
    end
    % generate a selector that chooses which n to use
    code = [code, sprintf(['\n' o.indent.code 'void build_n(const ' o.real '* u, const unsigned int k, ' o.real '* n){' '\n\n'])];
    for ii=1:length(o.K_n)
        check = generateCheck_custom(o.K_n{ii}, [o.indent.code '\t'], 'k');
        code = [code, check, sprintf('\n')]; %#ok
        code = [code, sprintf([o.indent.code '\t\t' 'build_n_' num2str(ii) '(&u[0], &n[0]);' '\n'])]; %#ok
        code = [code, sprintf([o.indent.code '\t' '}' '\n'])]; %#ok
    end
    code = [code, sprintf([o.indent.code '}' '\n'])];
    
    % generate a selector that chooses which Dn to use
    code = [code, sprintf(['\n' o.indent.code 'void build_Dn(const ' o.real '* u, const unsigned int k, ' o.real '* Dn){' '\n\n'])];
    for ii=1:length(o.K_n)
        check = generateCheck_custom(o.K_n{ii}, [o.indent.code '\t'], 'k');
        code = [code, check, sprintf('\n')]; %#ok
        code = [code, sprintf([o.indent.code '\t\t' 'build_Dn_' num2str(ii) '(&u[0], &Dn[0]);' '\n'])]; %#ok
        code = [code, sprintf([o.indent.code '\t' '}' '\n'])]; %#ok
    end
    code = [code, sprintf([o.indent.code '}' '\n'])];
end


% build inverse: ToDo eventually delete this part
if( ~o.forceGradient)
    if( isequal(o.gradients,'casadi'))
        % build_tbi
        
        %sl = SX.sym('sl',o.nu);
        sl = SX.sym('sl',size(Dg,2));
        tbi = transpose(Dg)*Dg + diag(sl)*diag(sl);
        [info.in_tbi.static,i] = falcopt.casadi2struct(tbi);
        
        % the inversion does not work with this variant
        %tbi_f = Function('build_tbi_casadi',{z,sl},{i.stored.values});
        tbi_f = Function('build_tbi_casadi',{z,sl},{tbi});
        sxfcn{length(sxfcn)+1} =  tbi_f;
        
        % wrapper build_tbi
        code = [code, sprintf(['void build_tbi( const ' o.real '* u, const ' o.real '* sl, ' o.real '* A){' '\n\n'])];
        code = [code, sprintf(['\t' 'const ' o.real ' *in[%i];\n'],2)];
        code = [code, sprintf(['\t' o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)];
        code = [code, sprintf(['\tin[0] = u;\n\tin[1] = sl;\n\n' ...
            '\tbuild_tbi_casadi( in, &A, &iw, w, mem);\n}\n']) ];
        info.in_tbi.struct.structure = i;
        try
            info.in_tbi.flops =   tbi_f.getAlgorithmSize(); %flops
        catch
            warning('Cannot use casadi.Function.getAlgoirthmSize()');
        end
        
        
        % build inverse
        in_M.structure.stored.mat = i.stored.mat;
        info.Jac_m_struct = i.stored.mat;
        c = falcopt.generateInverse(in_M.structure.stored.mat, 'names', struct('fun', 'build_inv'), 'symmetric', true, 'indent', o.indent, 'inline', o.inline, 'types', o.real, 'precision', o.precision, 'test', o.test, 'verbose', o.verbose);
        code = [code, c];
        
        
    elseif ( isequal(o.gradients,'matlab'))
        sl = sym('sl',[size(Dg,2),1],'real');
        tbi = Dg'*Dg + diag(sl)*diag(sl);
        [c,i] = falcopt.fcn2struct( tbi, o, 'name', 'A', 'structure' ,'sparse');
        code = [code, sprintf(['void build_tbi( const ' o.real '* u, const ' o.real '* sl, ' o.real '* A){' '\n\n'])];
        code = [code, c, sprintf('\n}\n')];
        
        % build inverse
        [c,i] = falcopt.fcn2struct( tbi, o, 'name', 'A');
        in_M.structure.stored.mat = i.structure.mat;
        info.Jac_m_struct = i.structure.mat;
        c = falcopt.generateInverse(in_M.structure.stored.mat, 'names', struct('fun', 'build_inv'), 'symmetric', true, 'indent', o.indent, 'inline', o.inline, 'types', o.real, 'precision', o.precision, 'test', o.test, 'verbose', o.verbose);
        code = [code, c];
    end
end

% generate .c file
if( isequal(o.gradients,'casadi')&& ~isempty(sxfcn))
    [info.src,info.header] = generate_casadi_c(o, 'constraints', sxfcn);
end
end

function [c] = argument_w(o,declaration)

if declaration
    str_real = sprintf(['const ' o.real '* ']);
else
    str_real = [];
end

c = [];

if o.nw > 0
    c = [c, sprintf([', ' str_real 'w'])];
end

end

function c = argument_contr_value(o,declaration)

if declaration
    str_real = sprintf(['const ' o.real ' ']);
else
    str_real = [];
end

c = [];

if (o.contractive || o.terminal)
    c = [c, sprintf([', ' str_real 'c_contr'])];
end
end

function [c] = argument_def(o,declaration)

if declaration
    str_real = sprintf(['const ' o.real '* ']);
    str_int  = sprintf('const unsigned int ');
else
    str_real = [];
    str_int = [];
end

c = [];

if o.trackReference
    c = [c, sprintf([', ' str_real 'xref, ' str_real 'uref'])];
end
if o.contractive
    c = [c, sprintf([', ' str_int 'ind'])];
end
end

function [c] = argument_def_internal_psi(o,declaration)

if declaration
    str_real = sprintf(['const ' o.real '* ']);
else
    str_real = [];
end

c = [];

if (o.contractive || o.terminal)
    c = [c, sprintf([', ' str_real 'psi_N'])];
end
end

function [c] = argument_def_internal_psi_noconst(o,declaration)

if declaration
    str_real = sprintf([o.real '* ']);
else
    str_real = [];
end

c = [];

if (o.contractive || o.terminal)
    c = [c, sprintf([', ' str_real 'psi_N'])];
end
end

function [c] = argument_def_internal_psi_plus(o,declaration)

if declaration
    str_real = sprintf(['const ' o.real '* ']);
else
    str_real = [];
end

c = [];

if (o.contractive || o.terminal)
    c = [c, sprintf([', ' str_real 'psi_Np'])];
end
end

function [c] = argument_def_internal_psi_dot_noconst(o,declaration)

if declaration
    str_real = sprintf([o.real '* ']);
else
    str_real = [];
end

c = [];

if (o.contractive || o.terminal)
    c = [c, sprintf([', ' str_real 'dot_psi_N'])];
end
end

function [c] = argument_def_internal_psi_dot(o,declaration)

if declaration
    str_real = sprintf(['const ' o.real '* ']);
else
    str_real = [];
end

c = [];

if (o.contractive || o.terminal)
    c = [c, sprintf([', ' str_real 'dot_psi_N'])];
end
end


function [code, data] = generate_forward_simulation(o)

nx = o.nx;
nu = o.nu;
nw = o.nw;
N = o.N;

code = [];
data = [];

if o.nw >0
    code = [code, sprintf([o.inline ' void det_x (const ' o.real '* x0, const ' o.real '* u, const ' o.real '* w, ' o.real '* x){' '\n\n'])];
else
    code = [code, sprintf([o.inline ' void det_x (const ' o.real '* x0, const ' o.real '* u, ' o.real '* x){' '\n\n'])];
end

code = [code, sprintf(['\t' 'unsigned int ii = 0;' '\n\n'])];

if o.nw > 0
    code = [code, sprintf(['\t' 'model_mpc(x0,u,w,x);' '\n'])];
else
    code = [code, sprintf(['\t' 'model_mpc(x0,u,x);' '\n'])];
end

code = [code, sprintf(['\t' 'for (ii = 1;ii < %d; ++ii)' '\n'],N)];

if o.nw > 0
    code = [code, sprintf(['\t\t' 'model_mpc(x + (ii-1)* %d, u + ii* %d, w+ ii* %d, x + ii* %d);' '\n\n'],nx,nu,nw,nx)];
else
    code = [code, sprintf(['\t\t' 'model_mpc(x + (ii-1)* %d, u + ii* %d, x + ii* %d);' '\n\n'],nx,nu,nx)];
end
code = [code, sprintf(['}' '\n'])];

end

function [code, data, info] = generate_objective_gradient_oracle(o)

nx = o.nx;
nu = o.nu;
nw = o.nw;
N = o.N;
trackRef = o.trackReference;

code = [];
data = [];
info.flops.it = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
info.flops.ls = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

%it generates: dot_product_nx_nx, Qmul, Pmul, Rmul, dot_product_nu_nu,
%               product_and_sum_nu
[c, d, in] = generate_auxiliary_functions(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.multFlops(in.flops, N);
info.flops.it = falcopt.addFlops(info.flops.it,in.flops);
info.flops.ls = falcopt.addFlops(info.flops.ls,in.flops);


[c, d, in] = generate_product_and_sum_nx(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.multFlops(in.flops, N-1);
info.flops.it = falcopt.addFlops(info.flops.it,in.flops);
info.flops.ls = falcopt.addFlops(info.flops.ls,in.flops);


if trackRef
    [c, d, in] = generate_diffXU(o);
    code = [code, c];
    data = [data, d];
    in.flops = falcopt.multFlops(in.flops, N);
    info.flops.it = falcopt.addFlops(info.flops.it,in.flops);
    info.flops.ls = falcopt.addFlops(info.flops.ls,in.flops);
end


%% det_J_and_dot_J
[c,d] = generate_copy(o);
% generate copy_Nnc copy_Nnu copy_Nnx
code = [code, c];
data = [data, d];

c_w_dec = argument_w(o,true);
c_tr_dec = argument_def(o,true);
c_psi_dec = argument_def_internal_psi_noconst(o,true);
c_psi_dot_dec = argument_def_internal_psi_dot_noconst(o,true);

code = [code, sprintf([o.inline ' void det_J_and_dot_J(const ' o.real '* x0, const ' o.real ...
    '* u, const ' o.real '* x' c_w_dec c_tr_dec ', ' ...
    o.real '* J, ' o.real '* dot_J' c_psi_dec c_psi_dot_dec '){' '\n\n',...
    ])];

code = [code, sprintf(['\t' o.real ' Px[%d], Qx[%d], mem_tmp2[%d], ' '\n'...
    '\t\t' 'Ru[%d], tmp_x = 0.0, tmp_u = 0.0;' '\n'],...
    nx,nx,nx,nu)];

if o.contractive
    code = [code, sprintf(['\t' 'int index = 0;' '\n'])];
end

if trackRef
    code = [code, sprintf(['\t' o.real ' dx[%d], du[%d];' '\n'],...
        nx,nu)];
end
if o.contractive
    code = [code, sprintf(['\t' o.real ' Px_contr[%d], mem_tmp_contr[%d], tmp_contr = 0.0;' '\n'], nx, nx)];
    if trackRef
        code = [code, sprintf(['\t' o.real ' dx_contr[%d];' '\n'], nx)];
    end
elseif o.terminal
    code = [code, sprintf(['\t' o.real ' Px_contr[%d], mem_tmp_contr[%d];' '\n'], nx, nx)];
end

code = [code, sprintf(['\t' 'unsigned int ii = 0;' '\n\n'])];

if trackRef
    code = [code, sprintf(['\t' 'diffX(dx, x + %d, xref + %d);' '\n'],(N-1)*nx, (N-1)*nx)];
    code = [code, sprintf(['\t' 'Pmul(Px, dx);' '\n'])];
    code = [code, sprintf(['\t' 'dot_product_nx_nx(&tmp_x,Px, dx);' '\n\n'])];
    if o.contractive
        code = [code, sprintf(['\t' 'diffX(dx_contr, x + (ind-1)*%d, xref + (ind-1)*%d);' '\n'],nx, nx)];
        code = [code, sprintf(['\t' 'Pmul(Px_contr, dx_contr);' '\n'...
            '\t' 'dot_product_nx_nx(&tmp_contr,Px_contr,dx_contr);' '\n'...
            '\t' '(*psi_N) = 0.5*tmp_contr;' '\n'])];
    end
else
    code = [code, sprintf(['\t' 'Pmul(Px, x + %d);' '\n'], (N-1)*nx)];
    code = [code, sprintf(['\t' 'dot_product_nx_nx(&tmp_x,Px, x + %d);' '\n\n'],(N-1)*nx)];
    if o.contractive
        code = [code, sprintf(['\t' 'Pmul(Px_contr, x+ (ind - 1)*%d);' '\n'...
            '\t' 'dot_product_nx_nx(&tmp_contr,Px_contr,x+ (ind - 1)*%d);' '\n'...
            '\t' '(*psi_N) = 0.5*tmp_contr;' '\n'],nx,nx)];
    end
end

if o.terminal
    code = [code, sprintf(['\t' '(*psi_N) = 0.5*tmp_x;' '\n',...
        '\t' 'copy_Nnu(Px_contr,Px);' '\n'])];
    
end

if trackRef
    code = [code, sprintf(['\t' 'diffU(du, u + %d, uref + %d);' '\n'],(N-1)*nu, (N-1)*nu)];
    code = [code, sprintf(['\t' 'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf(['\t' 'dot_product_nu_nu(&tmp_u,Ru,du);' '\n'])];
else
    code = [code, sprintf(['\t' 'Rmul(Ru, u + %d);' '\n'], (N-1)*nu)];
    code = [code, sprintf(['\t' 'dot_product_nu_nu(&tmp_u,Ru, u + %d);' '\n'],(N-1)*nu)];
end

code = [code, sprintf(['\t' '(*J) = 0.5*(tmp_x + tmp_u);' '\n'])];


info.flops.it.mul = info.flops.it.mul+1;
info.flops.it.add = info.flops.it.add+1;


if o.nw > 0
    if ~o.Jac_u_static
        code = [code, sprintf(['\t' 'Jacobian_u(x + %d,u + %d, w + %d, G);' '\n'],(N-2)*nx,(N-1)*nu,(N-1)*nw)];
    end
else
    if ~o.Jac_u_static
        code = [code, sprintf(['\t' 'Jacobian_u(x + %d,u + %d, G);' '\n'],(N-2)*nx,(N-1)*nu)];
    end
end

code = [code, sprintf(['\t' 'product_and_sum_nu(dot_J + %d,Px, Ru, G);' '\n\n'], (N-1)*nu)];

if o.terminal
    code = [code, sprintf(['\t' 'product_contr_nu(&dot_psi_N[%d], Px, G);' '\n\n'], (N-1)*nu)];
elseif o.contractive
    code = [code, sprintf(['\t' 'if (ind == %d) ' '\n'...
        '\t\t' 'product_contr_nu(&dot_psi_N[%d], Px_contr, G);' '\n'], N, (N-1)*nu)];
    code = [code, sprintf(['\t' 'else' '\n'...
        '\t\t' 'set_zero_nu(&dot_psi_N[%d]);' '\n\n'], (N-1)*nu)];
end



if o.contractive
    code = [code, sprintf(['\t' 'if (ind == %d) ' '\n'...
        '\t\t' 'product_contr_nu(&dot_psi_N[%d], Px_contr, G);' '\n'], N, (N-1)*nu)];
    code = [code, sprintf(['\t' 'else' '\n'...
        '\t\t' 'set_zero_nu(&dot_psi_N[%d]);' '\n\n'], (N-1)*nu)];
end

code = [code, sprintf(['\t' 'for (ii=%d; ii-->0; ) {' '\n'],N-1)];

if o.nw > 0
    if ~o.Jac_x_static
        code = [code, sprintf(['\t\t' 'Jacobian_x(x + ii*%d,u + (ii+1)*%d, w + (ii+1)*%d, F);' '\n'],nx,nu,nw)];
        info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
        
    end
else
    if ~o.Jac_x_static
        code = [code, sprintf(['\t\t' 'Jacobian_x(x + ii*%d,u + (ii+1)*%d, F);' '\n'],nx,nu)];
        info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
    end
end

if o.nw > 0
    if ~o.Jac_u_static
        code = [code, sprintf(['\t\t' 'if (ii==0)' '\n',...
            '\t\t\t' 'Jacobian_u(x0,u,w,G);' '\n',...
            '\t\t' 'else' '\n',...
            '\t\t\t' 'Jacobian_u(x + (ii-1)*%d, u + ii*%d, w + ii*%d, G);' '\n\n'],nx,nu,nw)];
        info.flops.it.mul = info.flops.it.mul+ 3*(N-1);
        info.flops.it.comp = info.flops.it.comp+ (N-1);
    end
else
    if ~o.Jac_u_static
        code = [code, sprintf(['\t\t' 'if (ii==0)' '\n',...
            '\t\t\t' 'Jacobian_u(x0,u,G);' '\n',...
            '\t\t' 'else' '\n',...
            '\t\t\t' 'Jacobian_u(x + (ii-1)*%d, u + ii*%d, G);' '\n\n'],nx,nu)];
        info.flops.it.mul = info.flops.it.mul+ 3*(N-1);
        info.flops.it.comp = info.flops.it.comp+ (N-1);
    end
end


if trackRef
    code = [code, sprintf(['\t\t' 'diffX(dx, x + ii*%d, xref + ii*%d);' '\n'],nx, nx)];
    info.flops.it.mul = info.flops.it.mul+ 3*(N-1);
    code = [code, sprintf(['\t\t' 'Qmul(Qx, dx);' '\n'])];
    code = [code, sprintf(['\t\t' 'dot_product_nx_nx(&tmp_x,Qx,dx);' '\n'])];
else
    code = [code, sprintf(['\t\t' 'Qmul(Qx, x + ii*%d);' '\n'], nx)];
    code = [code, sprintf(['\t\t' 'dot_product_nx_nx(&tmp_x,Qx, x + ii*%d);' '\n'],nx)];
    info.flops.it.mul = info.flops.it.mul+ (nx+nx)*(N-1);
end

if trackRef
    code = [code, sprintf(['\t\t' 'diffU(du, u + ii*%d, uref + ii*%d);' '\n'],nu, nu)];
    info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
    code = [code, sprintf(['\t\t' 'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf(['\t\t' 'dot_product_nu_nu(&tmp_u,Ru,du);' '\n'])];
else
    code = [code, sprintf(['\t\t' 'Rmul(Ru, u + ii*%d);' '\n'], nu)];
    code = [code, sprintf(['\t\t' 'dot_product_nu_nu(&tmp_u,Ru, u + ii*%d);' '\n'],nu)];
    info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
end

code = [code, sprintf(['\t\t' '(*J) += .50*(tmp_x + tmp_u);' '\n\n'])];
info.flop.it.mul = info.flops.it.mul+1*(N-1);
info.flops.it.add = info.flops.it.add+2*(N-1);



if o.terminal
    code = [code, sprintf(['\t\t' 'product_contr_nx(mem_tmp_contr, Px_contr, F);' '\n',...
        '\t\t' 'copy_nx(Px_contr,mem_tmp_contr);' '\n',...
        '\t\t' 'product_contr_nu(&dot_psi_N[ii*%d], Px_contr, G);' '\n\n'], nu)];
end

if o.contractive
    code = [code, sprintf(['\t\t' 'index = ind - ii - 1;' '\n'])];
    
    code = [code, sprintf(['\t\t' 'if ( index >= 0){' '\n',...
        '\t\t\t' 'if (index > 0){' '\n',...
        '\t\t\t\t' 'product_contr_nx(mem_tmp_contr, Px_contr, F);' '\n',...
        '\t\t\t\t' 'copy_nx(Px_contr,mem_tmp_contr);' '\n',...
        '\t\t\t' '}' '\n',...
        '\t\t\t' 'product_contr_nu(&dot_psi_N[ii*%d], Px_contr, G);' '\n'...
        '\t\t' '}' '\n'], nu)];
    code = [code, sprintf(['\t\t' 'else' '\n',...
        '\t\t\t' 'set_zero_nu(&dot_psi_N[ii*%d]);' '\n\n'],nu)];
    
end

code = [code, sprintf(['\t\t' 'product_and_sum_nx(mem_tmp2, Px, Qx, F);' '\n',...
    '\t\t' 'copy_nx(Px,mem_tmp2);' '\n',...
    '\t\t' 'product_and_sum_nu(dot_J + ii*%d, Px, Ru, G);' '\n',...
    '\t' '}' '\n',...
    '}' '\n\n'],nu)];
info.flops.it.mul = info.flops.it.mul+ 1*(N-1);

%% det_J


code = [code, sprintf([o.inline ' void det_J(const ' o.real '* x0, const ' o.real ...
    '* u, const ' o.real '* x' c_tr_dec ', ' ...
    o.real '* J' c_psi_dec '){' '\n\n',...
    ])];

code = [code, sprintf(['\t' o.real ' Qx[%d], Ru[%d], tmp_x = 0.0, tmp_u = 0.0;' '\n'],...
    nx,nu)];
if trackRef
    code = [code, sprintf(['\t' o.real ' dx[%d], du[%d];' '\n'],...
        nx,nu)];
end
if o.contractive
    code = [code, sprintf(['\t' o.real ' Px_contr[%d], tmp_contr = 0.0;' '\n'], nx)];
    if trackRef
        code = [code, sprintf(['\t' o.real ' dx_contr[%d];' '\n'], nx)];
    end
end

code = [code, sprintf(['\t' 'unsigned int ii = 0;' '\n\n'])];

if trackRef
    code = [code, sprintf(['\t' 'diffX(dx, x + %d, xref + %d);' '\n'],(N-1)*nx, (N-1)*nx)];
    code = [code, sprintf(['\t' 'Pmul(Qx, dx);' '\n'])];
    code = [code, sprintf(['\t' 'dot_product_nx_nx(&tmp_x,Qx, dx);' '\n'])];
    if o.contractive
        code = [code, sprintf(['\t' 'diffX(dx_contr, x + (ind - 1)*%d, xref + (ind-1)*%d);' '\n'],nx, nx)];
        code = [code, sprintf(['\t' 'Pmul(Px_contr, dx_contr);' '\n'...
            '\t' 'dot_product_nx_nx(&tmp_contr,Px_contr,dx_contr);' '\n',...
            '\t' '(*psi_N) = 0.5*tmp_contr;' '\n'])];
    end
else
    code = [code, sprintf(['\t' 'Pmul(Qx, x + %d);' '\n'], (N-1)*nx)];
    code = [code, sprintf(['\t' 'dot_product_nx_nx(&tmp_x,Qx, x + %d);' '\n'],(N-1)*nx)];
    if o.contractive
        code = [code, sprintf(['\t' 'Pmul(Px_contr, x+ (ind - 1)*%d);' '\n'...
            '\t' 'dot_product_nx_nx(&tmp_contr,Px_contr,x+ (ind - 1)*%d);' '\n',...
            '\t' '(*psi_N) = 0.5*tmp_contr;' '\n'],nx,nx)];
    end
end

if o.terminal
    code = [code, sprintf(['\t' '(*psi_N) = 0.5*tmp_x;' '\n'])];
end

if trackRef
    code = [code, sprintf(['\t' 'diffU(du, u + %d, uref + %d);' '\n'],(N-1)*nu, (N-1)*nu)];
    code = [code, sprintf(['\t' 'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf(['\t' 'dot_product_nu_nu(&tmp_u,Ru,du);' '\n\n'])];
else
    code = [code, sprintf(['\t' 'Rmul(Ru, u + %d);' '\n'], (N-1)*nu)];
    code = [code, sprintf(['\t' 'dot_product_nu_nu(&tmp_u,Ru, u + %d);' '\n\n'],(N-1)*nu)];
end

code = [code, sprintf(['\t' '(*J) = 0.5*(tmp_x + tmp_u);' '\n'])];
info.flops.ls.mul = info.flops.ls.mul+1;
info.flops.ls.add = info.flops.ls.add+1;

code = [code, sprintf(['\t' 'for (ii=%d; ii-->0; ) {' '\n\n'],N-1)];

if trackRef
    code = [code, sprintf(['\t\t' 'diffX(dx, x + ii*%d, xref + ii*%d);' '\n'],nx, nx)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
    code = [code, sprintf(['\t\t' 'Qmul(Qx, dx);' '\n'])];
    code = [code, sprintf(['\t\t' 'dot_product_nx_nx(&tmp_x,Qx,dx);' '\n\n'])];
else
    code = [code, sprintf(['\t\t' 'Qmul(Qx, x + ii*%d);' '\n'], nx)];
    code = [code, sprintf(['\t\t' 'dot_product_nx_nx(&tmp_x,Qx, x + ii*%d);' '\n\n'],nx)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
end

if trackRef
    code = [code, sprintf(['\t\t' 'diffU(du, u + ii*%d, uref + ii*%d);' '\n'],nu, nu)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
    code = [code, sprintf(['\t\t' 'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf(['\t\t' 'dot_product_nu_nu(&tmp_u,Ru,du);' '\n\n'])];
else
    code = [code, sprintf(['\t\t' 'Rmul(Ru, u + ii*%d);' '\n'], nu)];
    code = [code, sprintf(['\t\t' 'dot_product_nu_nu(&tmp_u,Ru, u + ii*%d);' '\n\n'],nu)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
end


code = [code, sprintf(['\t\t' '(*J) += .5*(tmp_x + tmp_u);' '\n'...
    '\t' '}' '\n'...
    '}' '\n\n'])];

info.flops.ls.mul = info.flops.ls.mul+ 1*(N-1);
info.flops.ls.add = info.flops.ls.add+ 2*(N-1);



end

function [code, data, info] = generate_auxiliary_functions(o)

code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nx = o.nx;
nu = o.nu;
nc = o.nc(1);
Q = o.Q;
R = o.R;
P = o.P;




[d, c, in] = falcopt.generateMVMult(Q, ...
    'names', struct('fun', 'Qmul', 'M', 'Q', 'v', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);

[d, c] = falcopt.generateMVMult(P, ...
    'names', struct('fun', 'Pmul', 'M', 'P', 'v', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
% flops of Pmul are counted together with Qmul


if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

[d, c, in] = falcopt.generateMVMult(R, ...
    'names', struct('fun', 'Rmul', 'M', 'R', 'v', 'du'), 'types', o.real,...
    'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);

[d, c, in] = falcopt.generateMVMult(ones(1,nx), ...
    'names', struct('fun', 'dot_product_nx_nx', 'M', 'R',...
    'v', 'du'), 'static', false, 'types', o.real, 'verbose', o.verbose, 'precision', o.precision, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);

[d, c, in] = falcopt.generateMVMult(ones(1,nu), ...
    'names', struct('fun', 'dot_product_nu_nu', 'M', 'R',...
    'v', 'du'), 'static', false, 'types', o.real, 'verbose', o.verbose, 'precision', o.precision, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);


[d, c, in] = falcopt.generateMVMult({o.Jac_u_struct,eye(nu)}, ...
    'names', struct('fun', 'product_and_sum_nu', 'M', {{'A', 'I'}},...
    'v', {{'u1', 'u2'}}), 'static', [false,true], 'structure', 'ordered', 'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);


% [d, c] = falcopt.generateMVMult({eye(nc),eye(nc)}, ... % To Be deleted
%     'names', struct('fun', 'sum_nc', 'M', {{'I1', 'I2'}},...
%     'v', {{'u1', 'u2'}}), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
% 
% if ~isempty(d)
%     data = [data, d, sprintf('\n')];
% end
% code = [code, c, sprintf('\n\n')];

end

function [code, data, info] = generate_product_and_sum_nx(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nx = o.nx;
nu = o.nu;

if o.contractive
    code =  [code, sprintf(['void set_zero_nu (' o.real '* x){\n\n'])];
    
    for jj=0:nu-1
        code =  [code, sprintf(['\t' 'x[%d] = 0.0;' '\n'],jj)];
    end
    
    code =  [code, sprintf(['\n\n' '}' '\n\n'])];
end

if o.contractive|| o.terminal
    [~, c, in] = falcopt.generateMVMult({o.Jac_x_struct}, ...
        'names', struct('fun', 'product_contr_nx', 'M', {{'A'}},...
        'v', {{'u'}}), 'static', false, 'structure', 'ordered',...
        'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test,...
        'inline', o.inline, 'indent', o.indent);
    
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
    
    [~, c, in] = falcopt.generateMVMult({o.Jac_u_struct}, ...
        'names', struct('fun', 'product_contr_nu', 'M', {{'A'}},...
        'v', {{'u'}}), 'static', false, 'structure', 'ordered',...
        'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test,...
        'inline', o.inline, 'indent', o.indent);
    
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
    
end

[d, c, in] = falcopt.generateMVMult({o.Jac_x_struct,eye(nx)}, ...
    'names', struct('fun', 'product_and_sum_nx', 'M', {{'A', 'I'}},...
    'v', {{'u1', 'u2'}}), 'static', [false,true], 'structure', 'ordered', 'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);

end

function [code, data, info] = generate_diffXU(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nx = o.nx;
nu = o.nu;

[d, c, in] = falcopt.generateMVMult({eye(nx), -eye(nx)}, ...
    'names', struct('fun', 'diffX', 'M', {{'I', 'mI'}}, 'v', {{'x', 'xref'}}, 'r', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);

[d, c, in] = falcopt.generateMVMult({eye(nu), -eye(nu)}, ...
    'names', struct('fun', 'diffU', 'M', {{'I', 'mI'}}, 'v', {{'u', 'uref'}}, 'r', 'du'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);
end

function [code, data, info, optCode] = generate_slack_initialization(o)

code = [];
data = [];
optCode = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nc = o.nc(1);
nu = o.nu;
N = o.N;

c_psi_dec = argument_def_internal_psi(o,true);
c_contr_dec = argument_contr_value(o,true);

code = [code, sprintf([o.inline ' void build_sl_slsqr('])];
if ~isempty(o.K_lb)
    code = [code, sprintf(['const ' o.real '* amu, '])];
end
if ~isempty(o.K_ub)
    code = [code, sprintf(['const ' o.real '* umb, '])];
end
if ~isempty(o.K_n)
    code = [code, sprintf(['const ' o.real '* n, '])];
end
code = [code,  sprintf(['const unsigned int nc, const unsigned int na, '...
    'const unsigned int nb, ' o.real '* sl, ' o.real '* sl_sqr, ' o.real '* gps){' '\n\n'...
    '\t' 'unsigned int jj;' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code,  sprintf(['\t' 'for (jj=0;jj< na;jj++){' '\n',...
        '\t\t' 'sl_sqr[jj] = ' o.max '(1.0, -2.0* amu[jj]);' '\n',...
        '\t\t' 'gps[jj]= amu[jj] + 0.5*sl_sqr[jj];' '\n',...
        '\t\t' 'sl[jj] = ' o.sqrt '(sl_sqr[jj]);' '\n',...
        '\t' '}' '\n'])];
end
if ~isempty(o.K_ub)
    code = [code,  sprintf(['\t' 'for (jj=na;jj< nb;jj++){' '\n'...
        '\t\t' 'sl_sqr[jj] = ' o.max '(1.0, -2.0* umb[jj-na]);' '\n',...
        '\t\t' 'gps[jj]= umb[jj-na] + 0.5*sl_sqr[jj];' '\n',...
        '\t\t' 'sl[jj] = ' o.sqrt '(sl_sqr[jj]);' '\n',...
        '\t' '}' '\n'])];
end
if ~isempty(o.K_n)
    code = [code,  sprintf(['\t' 'for (jj= nb;jj< nc;jj++) {' '\n'...
        '\t\t' 'sl_sqr[jj] = ' o.max '(1.0, -2.0* n[jj - nb]);' '\n',...
        '\t\t' 'gps[jj]= n[jj - nb] + 0.5*sl_sqr[jj];' '\n',...
        '\t\t' 'sl[jj] = ' o.sqrt '(sl_sqr[jj]);' '\n',...
        '\t' '}' '\n'])];
end

code = [code,  sprintf(['}' '\n\n'])];

code = [code, sprintf([o.inline ' void initialize_slack( const ' o.real '* u' c_psi_dec c_contr_dec ', ' o.real '* sl, ' o.real '* sl_sqr, ' o.real '* gps){' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code, sprintf(['\t'  o.real ' amu[%d];' '\n'], max(sum(~isinf( o.box_lowerBound))) )];
end
if ~isempty(o.K_ub)
    code = [code, sprintf(['\t'  o.real ' umb[%d];' '\n'], max(sum(~isinf( o.box_upperBound))) )];
end
if ~isempty(o.K_n)
    code = [code, sprintf(['\t'  o.real ' n[%d];' '\n'], max(cell2mat(o.nn)) )];
end

if (o.contractive || o.terminal)
    code = [code, sprintf(['\t' o.real ' g_contr = 0.0;' '\n'])];
end



info.flops.mul = info.flops.mul+ N*3; % mul inside only first cycle % NOT CLEAR, ToDo


lb = sum(~isinf(o.box_lowerBound),1);
ub = sum(~isinf(o.box_upperBound),1);

for k = 1:o.N
    code = [code, sprintf(['\n'...
        '\t' '/* Unrolling the for loop: iteration %i of %i */' '\n'], k-1, o.N - 1)];   %#ok
    if ~isempty(o.K_lb)
        code = [code, sprintf(['\t' 'build_amu(&u[%i], %i, &amu[0]);' '\n'],o.nu*(k-1),k-1)];%#ok
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf(['\t' 'build_umb(&u[%i],%i,&umb[0]);' '\n'],o.nu*(k-1),k-1)];%#ok
    end
    if ~isempty(o.K_n)
        code = [code, sprintf(['\t' 'build_n(&u[%i],%i,&n[0]);' '\n'],o.nu*(k-1),k-1)];%#ok
    end
    code = [code, sprintf(['\t' 'build_sl_slsqr('])]; %#ok
    if ~isempty(o.K_lb)
        code = [code, sprintf([' &amu[0], '])]; %#ok
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf([' &umb[0], '])]; %#ok
    end
    if ~isempty(o.K_n)
        code = [code, sprintf(' &n[0], ')]; %#ok
    end
    code = [code,  sprintf(['%i, %i, %i, &sl[%i], &sl_sqr[%i], &gps[%i]);' '\n'],...
        o.nc(k), lb(:,k), lb(:,k) + ub(:,k), sum(o.nc(1:k-1)),...
        sum(o.nc(1:k-1)),sum(o.nc(1:k-1)))]; %#ok
end

% Obsolete, To be Deleted
% code = [code, sprintf(['\t\t' 'for (jj=0;jj<%d;jj++){' '\n'],nc)];
%
% code = [code, sprintf(['\t\t\t' 'sl_sqr[ii*%d + jj] = ' o.max '(1.0, -2.0*g[jj]);' '\n',...
%     '\t\t\t' 'sl[ii*%d + jj] = ' o.sqrt '(sl_sqr[ii*%d + jj]);' '\n',...
%     '\t\t' '}' '\n',...
%     '\t\t' 'half_sum_nc(&gps[ii*%d], &g[0], &sl_sqr[ii*%d]);' '\n',...
%     '\t' '}' '\n\n'],nc,nc,nc,nc,nc)];

if (o.contractive || o.terminal)
    code = [code, sprintf(['\t' 'g_contr = *psi_N - c_contr;' '\n'...
        '\t' 'sl_sqr[%d] = ' o.max '(1.0, -2.0*g_contr);' '\n',...
        '\t' 'sl[%d] = ' o.sqrt '(sl_sqr[%d]);' '\n'...
        '\t' 'gps[%d] = g_contr + 0.5*sl_sqr[%d];' '\n'],...
        sum(o.nc) - 1, sum(o.nc) - 1, sum(o.nc) - 1, sum(o.nc)-1, sum(o.nc)-1)];
    info.flops.add = info.flops.add + 1;
end

code = [code, sprintf(['}' '\n\n'])];


info.flops.mul = info.flops.mul+ 4*sum(o.nc); % mult inside the two cycle
info.flops.sqrt = info.flops.sqrt +sum(o.nc);
info.flops.comp = info.flops.comp +sum(o.nc);   % max function


% the flops from function "build_g" are conuted in generatedConverter

end

function [code, data, info] = generate_gradient_step(o)

code = [];
data = [];

nc = o.nc(1);
nu = o.nu;
N = o.N;
alpha = o.stepSize;

c_psi_dot_dec = argument_def_internal_psi_dot(o,true);


info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);



%% gradient_step

% lower bounds
for jj = 1:length(o.K_lb) % if o.K_lb is empty, this loop is not executed
    struct_mat = diag(~isinf(o.box_lowerBound(:,o.K_lb{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    [d, c, in] = falcopt.generateMVMult({1/alpha*eye(sum(~isinf(o.box_lowerBound(:,o.K_lb{jj}(1))))), struct_mat_cut}, ...
        'names', struct('fun', ['build_vNnc_lb_' num2str(jj)], 'M', {{ ['alpha_inv_lb_' num2str(jj)], ['I_lb_' num2str(jj)]}},...
        'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
end
if ~isempty(o.K_lb)
    code = [code, sprintf([ o.inline ' void build_vNnc_lb(const ' o.real '* gps, '...
        'const ' o.real '* dot_J, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_lb)
        code = [code, generateCheck_custom(o.K_lb{jj}, '\t', 'k'),...
            sprintf(['\n',...
            '\t\t' 'build_vNnc_lb_%i(&res[0], &gps[0], &dot_J[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])];
    end
    code = [code, sprintf(['}' '\n'])];
end

% upper bounds
for jj = 1:length(o.K_ub) % if o.K_ub is empty, this loop is not executed
    struct_mat = diag(~isinf(o.box_upperBound(:,o.K_ub{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    [d, c, in] = falcopt.generateMVMult({1/alpha*eye(sum(~isinf(o.box_upperBound(:,o.K_ub{jj}(1))))), -struct_mat_cut}, ...
        'names', struct('fun', ['build_vNnc_ub_' num2str(jj)], 'M', {{ ['alpha_inv_ub_' num2str(jj)], ['I_ub_' num2str(jj)]}},...
        'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
end

if ~isempty(o.K_ub)
    code = [code, sprintf([ o.inline ' void build_vNnc_ub(const ' o.real '* gps, '...
        'const ' o.real '* dot_J, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_ub)
        code = [code, generateCheck_custom(o.K_ub{jj}, '\t', 'k'),...
            sprintf(['\n',...
            '\t\t' 'build_vNnc_ub_%i(&res[0], &gps[0], &dot_J[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])];
    end
    code = [code, sprintf(['}' '\n'])];
end

% nonlinear constraints
for jj = 1:length(o.K_n) % if o.K_n is empty, this loop is not executed
    [d, c, in] = falcopt.generateMVMult(o.Jac_n_struct{jj}, ...
        'names', struct('fun', ['Dntop_times_dotJ_' num2str(jj)], 'M',...
        {{ ['Dn_' num2str(jj)]}},...
        'v', {{'x'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'transpose', true, 'static', false,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
    
    [d, c, in] = falcopt.generateMVMult({1/alpha*eye(o.nn{o.K_n{jj}(1)}), -eye(o.nn{o.K_n{jj}(1)})}, ...
        'names', struct('fun', ['build_vNnc_n_' num2str(jj)], 'M', {{ ['alpha_inv_n_' num2str(jj)], ['temp_n' num2str(jj)]}},...
        'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'static', [true,true],...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
end

if ~isempty(o.K_n)
    code = [code, sprintf([ o.inline ' void Dntop_times_dotJ_n(const ' o.real '* Dn, '...
        'const ' o.real '* x, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_n)
        code = [code, generateCheck_custom(o.K_n{jj}, '\t', 'k'),...
            sprintf(['\n',...
            '\t\t' 'Dntop_times_dotJ_%i(&res[0], &x[0], &Dn[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])];
    end
    code = [code, sprintf(['}' '\n'])];
    
    code = [code, sprintf([ o.inline ' void build_vNnc_n(const ' o.real '* gps, '...
        'const ' o.real '* temp_n, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_n)
        code = [code, generateCheck_custom(o.K_n{jj}, '\t', 'k'),...
            sprintf(['\n',...
            '\t\t' 'build_vNnc_n_%i(&res[0], &gps[0], &temp_n[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])];
    end
    code = [code, sprintf(['}' '\n'])];
end


% [d, c, in] = falcopt.generateMVMult({o.Jac_g_struct}, ... % TO BE DELETED
%     'names', struct('fun', 'product_Dgtop_DJ', 'M', {{'A'}},...
%     'v', {{'x'}}, 'r', 'z'), 'indent', o.indent, 'types', o.real,'precision', o.precision, 'structure', 'ordered', 'verbose', o.verbose,...
%     'test', o.test, 'inline', o.inline, 'indent', o.indent, 'static', false, 'transpose', true);
% if ~isempty(d)
%     data = [data, d, sprintf('\n')];
% end
% code = [code, c, sprintf('\n\n')];
% info.flops = falcopt.addFlops(info.flops, in.flops);
%
% [d, c, in] = falcopt.generateMVMult({1/alpha*eye(nc), -eye(nc)}, ...
%     'names', struct('fun', 'scale_sub_nc', 'M', {{'alpha_inv', 'mI'}},...
%     'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
% if ~isempty(d)
%     data = [data, d, sprintf('\n')];
% end
% code = [code, c, sprintf('\n\n')];
% info.flops = falcopt.addFlops(info.flops, in.flops);


if (o.contractive || o.terminal)
    
    %     ... % TO BE DELETED
    %     [d, c, in] = falcopt.generateMVMult({eye(nu), o.Jac_g_struct, eye(nu)},
    %         'names', struct('fun', 'product_matrix_nu', 'M', {{'I', 'B', 'N'}},...
    %         'v', {{'x1', 'x2', 'x3'}}, 'r', 'z'),'structure', 'ordered', 'types', o.real, 'precision', o.precision, 'static', [true,false,false],...
    %         'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    %     [d, c, in] = falcopt.generateMVMult({eye(nu), eye(nu)}, ...
    %         'names', struct('fun', 'sum_terminal', 'M', {{'I', 'N'}},...
    %         'v', {{'x1', 'x2'}}, 'r', 'z'),'structure', 'ordered', 'types', o.real,...
    %         'precision', o.precision, 'static', [true,false],...
    %         'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    
    
    [d, c, in] = falcopt.generateMVMult({eye(nu)}, ...
        'names', struct('fun', 'sum_terminal', 'M', {{'I'}},...
        'v', {{'x1'}}, 'r', 'z'),'structure', 'ordered', 'types', o.real,...
        'precision', o.precision, 'static', false,...
        'add',true,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
    
end
    
    %     [d, c, in] = falcopt.generateMVMult({eye(nu), o.Jac_g_struct}, ... %TO BE DELETED
    %         'names', struct('fun', 'product_matrix_nu', 'M', {{'I', 'B'}},...
    %         'v', {{'x1', 'x2'}}, 'r', 'z'),'structure', 'ordered', 'types', o.real, 'precision', o.precision, 'static', [true,false],...
    %         'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    
    % lower bounds
for jj = 1:length(o.K_lb) % if o.K_lb is empty, this loop is not executed
    struct_mat = diag(~isinf(o.box_lowerBound(:,o.K_lb{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    [d, c, in] = falcopt.generateMVMult( - struct_mat_cut', ...
        'names', struct('fun', ['minus_Ina_muG_' num2str(jj)], 'M', {{['minus_Ina_' num2str(jj)]}},...
        'v', {{'x1'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
end
if ~isempty(o.K_lb)
    code = [code, sprintf([ o.inline ' void minus_Ina_muG(const ' o.real '* muG, '...
        'const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_lb)
        code = [code, generateCheck_custom(o.K_lb{jj}, '\t', 'k'),...
            sprintf(['\n',...
            '\t\t' 'minus_Ina_muG_%i(&res[0], &muG[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])];
    end
    code = [code, sprintf(['}' '\n'])];
end

%upper bounds
for jj = 1:length(o.K_ub) % if o.K_ub is empty, this loop is not executed
    struct_mat = diag(~isinf(o.box_upperBound(:,o.K_ub{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    [d, c, in] = falcopt.generateMVMult( struct_mat_cut', ...
        'names', struct('fun', ['Inb_muG_' num2str(jj)], 'M', {{['Inb_' num2str(jj)]}},...
        'v', {{'x1'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
end
if ~isempty(o.K_ub)
    code = [code, sprintf([ o.inline ' void Inb_muG(const ' o.real '* muG, '...
        'const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_ub)
        code = [code, generateCheck_custom(o.K_ub{jj}, '\t', 'k'),...
            sprintf(['\n',...
            '\t\t' 'Inb_muG_%i(&res[0], &muG[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])];
    end
    code = [code, sprintf(['}' '\n'])];
end

% nonlinear constraints
for jj = 1:length(o.K_n) % if o.K_n is empty, this loop is not executed
    [d, c, in] = falcopt.generateMVMult(o.Jac_n_struct{jj}, ...
        'names', struct('fun', ['Dn_times_muG_' num2str(jj)], 'M',...
        {{ ['Dn_' num2str(jj)]}},...
        'v', {{'x'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'transpose', false, 'static', false,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);

end

if ~isempty(o.K_n)
    code = [code, sprintf([ o.inline ' void Dn_times_muG_n(const ' o.real '* Dn, '...
        'const ' o.real '* x, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_n)
        code = [code, generateCheck_custom(o.K_n{jj}, '\t', 'k'),...
            sprintf(['\n',...
            '\t\t' 'Dn_times_muG_%i(&res[0], &x[0], &Dn[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])];
    end
    code = [code, sprintf(['}' '\n'])];
end

nr_constr = double(~isempty(o.K_lb)) + double(~isempty(o.K_ub)) + double(~isempty(o.K_n)) + 1;
[d, c, in] = falcopt.generateMVMult(repmat({eye(o.nu)},1,nr_constr), ...
    'names', struct('fun', 'sum_nr_constr'),'structure', 'ordered', 'types', o.real,...
    'precision', o.precision, 'static', true(1,nr_constr),...
    'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);
    
[d, c, in] = falcopt.generateMVMult(-alpha*eye(nu), ...
    'names', struct('fun', 'minus_scale_nu', 'M', {{'m_alpha'}},...
    'v', {{'x'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);



% [d,c, in] = falcopt.generateMVMult(o.Jac_m_struct,... % TO BE DELETED once ForceGradient is active
%     'names',struct('fun', 'product_nc_nc', 'M', {{'B'}},...
%     'v', {{'x'}}, 'r', 'r'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent, ...
%     'static',false,...
%     'symmetric',true);
% if ~isempty(d)
%     data = [data, d, sprintf('\n')];
% end
% code = [code, c, sprintf('\n\n')];
% info.flops = falcopt.addFlops(info.flops, in.flops);

for jj=1:length(o.K_nc)
    [~, c, in] = falcopt.generateMVMult(-alpha*eye(o.nc(o.K_nc{jj}(1))), ...
        'names', struct('fun', ['minus_scale_nc_' num2str(jj)], 'M', {{'m_alpha'}},...
        'v', {{'x'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision, 'static', true,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    %     if ~isempty(d)
    %         data = [data, d, sprintf('\n')];
    %     end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.addFlops(info.flops, in.flops);
    
    code = [code, sprintf([o.inline ' void product_matlab_nc_' num2str(jj) '(const ' o.real '* x, const '...
        o.real '* y, ' o.real '* z){' '\n\n',...
        '\t' 'unsigned int ii=0;' '\n\n',...
        '\t' 'for (ii=0;ii<%d;ii++)' '\n',...
        '\t\t' 'z[ii] = x[ii]*y[ii];' '\n\n',...
        '}' '\n'],o.nc(o.K_nc{jj}(1)) )];
    info.flops.mul = info.flops.mul+ o.nc(o.K_nc{jj}(1));
end
if ~isempty(o.K_nc)
    code = [code, sprintf([o.inline ' void minus_scale_nc(const ' o.real ...
        '* x, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj = 1:length(o.K_nc)
        check = generateCheck_custom(o.K_nc{jj}, '\t', 'k' );
        code = [code, sprintf([check, '\n' '\t\t' 'minus_scale_nc_%i( &res[0], &x[0]);' '\n'], jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])]; %#ok
    end
    code = [code, sprintf(['}' '\n'])];
    
    code = [code, sprintf([o.inline ' void product_matlab_nc(const ' o.real ...
        '* x, const ' o.real '* y, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj = 1:length(o.K_nc)
        check = generateCheck_custom(o.K_nc{jj}, '\t', 'k' );
        code = [code, sprintf([check, '\n' '\t\t' 'product_matlab_nc_%i(&x[0], &y[0], &res[0]);' '\n'], jj)]; %#ok
        code = [code, sprintf(['\t' '}' '\n'])]; %#ok
    end
    code = [code, sprintf(['}' '\n'])];
end

% TODO incroporate Dn, dt
if (o.contractive || o.terminal)
        ter_struct = repmat({ones(nu,1)},1,o.N); % Activate when ready
elseif o.forceGradient
        ter_struct = repmat([],1,o.N); % Activate when ready
end

[c, in] = falcopt.generateConstraintInv(o.Jac_n_struct_hor, ter_struct, 'N', N, 'types', o.real,'precision', o.precision, ...
    'indent', o.indent, 'inline', o.inline, 'verbose', max(0,o.verbose-1), 'test', o.test);
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);

% Gradient step function
code = [code, sprintf([o.inline ' void gradient_step(const ' o.real '* dot_J, const ' o.real '* u, const ' o.real '* sl,' '\n',...
    '\t' 'const ' o.real '* sl_sqr, const ' o.real '* gps' c_psi_dot_dec ', ' o.real '* du, ' o.real '* dsl, ' o.real '* muG){' '\n\n'])];

if (o.forceGradient) % This if is to be deleted
    
    if ~isempty(o.Jac_n_struct)
        
        % determine required size of Dn
        size_Dn = 0;
        for jj=1:length(o.K_n)
            current_size_Dn = max(max(o.Jac_n_struct{o.K_n{jj}(1)}))*length(o.K_n{jj});
            size_Dn = size_Dn + current_size_Dn;
        end
        
        code = [code, sprintf(['\t' o.real ' Dn[%d];' '\n'],size_Dn)];
    end
    
    code = [code, sprintf(['\t' o.real ' v_Nnc[%d], tmp_nu[%d], tmp_nc_m[%d], tmp_contr = 0.0;' '\n'],...
        sum(o.nc), nu, max(o.nc))];
    
    if ~isempty(o.K_lb)
        code = [code, sprintf(['\t' o.real ' temp_lb[%i];' '\n'],o.nu)];
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf(['\t' o.real ' temp_ub[%i];' '\n'],o.nu)];
    end
    if ~isempty(o.K_n)
        code = [code, sprintf(['\t' o.real ' temp_n[%i], temp_n2[%i];' '\n'], max(cell2mat(o.nn)), o.nu)];
    end
    
    
    Dn_need = 0;
    Dn_need_vec = zeros(1,o.N);
    for ii = 1:o.N
        
        code = [code, sprintf(['\n' '\t' '/* loop unrolling: step %i of %i */' '\n'],ii-1,o.N-1)];
        lb_size = 0;
        ub_size = 0;
        if ~isempty(o.K_lb)
            code = [code, sprintf(['\t' 'build_vNnc_lb(&gps[%i], &dot_J[%i], %i, &v_Nnc[%i]);' '\n'],...
                sum(o.nc(1:ii-1)), (ii-1)*o.nu, ii-1, sum(o.nc(1:ii-1)))];
            lb_size = sum(~isinf(o.box_lowerBound(:,ii)));
        end
        if ~isempty(o.K_ub)
            code = [code, sprintf(['\t' 'build_vNnc_ub(&gps[%i], &dot_J[%i], %i, &v_Nnc[%i]);' '\n'],...
                sum(o.nc(1:ii-1)) + lb_size, (ii-1)*o.nu, ii-1, sum(o.nc(1:ii-1)) + lb_size )];
            ub_size = sum(~isinf(o.box_upperBound(:,ii)));
        end
        if ~isempty(o.K_n)
            code = [code, sprintf(['\t' 'build_Dn(&u[%d], %i, &Dn[%d]);' '\n'], (ii-1)*nu, ii-1, Dn_need)];
            
            code = [code, sprintf(['\t' 'Dntop_times_dotJ_n(&Dn[%i], &dot_J[%i], %i, &temp_n[0]);' '\n'],...
                Dn_need, (ii-1)*o.nu, ii-1)];
            
            code = [code, sprintf(['\t' 'build_vNnc_n(&gps[%i], &temp_n[0], %i, &v_Nnc[%i]);' '\n'],...
                sum(o.nc(1:ii-1)) + lb_size + ub_size, ii-1, sum(o.nc(1:ii-1)) + lb_size + ub_size )];
            
            
            Dn_need = Dn_need + max(max(o.Jac_n_struct_hor{ii}));
            Dn_need_vec(ii) = Dn_need;

        end
    end
    
    %         To be deleted
    %         code = [code, sprintf(['\t\t' 'build_Dg(&u[ii*%d],Dg);' '\n'],nu)];
    %         code = [code, sprintf(['\t\t' 'product_Dgtop_DJ(Dg_dot_J, &dot_J[ii*%d], Dg);' '\n'],nu)];
    %         code = [code, sprintf(['\t\t' 'scale_sub_nc(&v_Nnc[ii*%d], &gps[ii*%d], Dg_dot_J); ' '\n'],nc,nc)];
    if (o.contractive || o.terminal)
        code = [code, sprintf(['\t' 'dot_product_Nnu(&tmp_contr,dot_psi_N,dot_J);' '\n',...
            '\t' 'v_Nnc[%d] = ' falcopt.num2str(1/alpha, o.precision) ' * gps[%d] - tmp_contr;' '\n'],sum(o.nc)-1, sum(o.nc)-1)];
    end
    
    % TODO: incorporate Dn when available
    code = [code, sprintf(['\n' '\t' 'solveConstraintSystem(&muG[0], '])];
    if ~isempty(o.K_n)
        code = [code, sprintf(['&Dn[0], '])];
    end
    if (o.contractive || o.terminal)
        code = [code, sprintf(['&dot_psi_N[0], '])];
    end
    code = [code, sprintf(['&v_Nnc[0], &sl_sqr[0]);' '\n'])];
    
    for ii = 1:o.N
        code = [code, sprintf(['\n' '\t' '/* loop unrolling: step %i of %i */' '\n'],ii-1,o.N-1)];
        lb_size = 0;
        ub_size = 0;
        if ~isempty(o.K_lb)
            code = [code, sprintf(['\t' 'minus_Ina_muG(&muG[%i], %i, &temp_lb[0]);' '\n'],...
                sum(o.nc(1:ii-1)), ii-1)];
            lb_size = sum(~isinf(o.box_lowerBound(:,ii)));
        end
        if ~isempty(o.K_ub)
            code = [code, sprintf(['\t' 'Inb_muG(&muG[%i], %i, &temp_ub[0]);' '\n'],...
                sum(o.nc(1:ii-1)) + lb_size, ii-1)];
            ub_size = sum(~isinf(o.box_upperBound(:,ii)));
        end
        if ~isempty(o.K_n)
            code = [code, sprintf(['\t' 'Dn_times_muG_n(&Dn[%i], &muG[%i], %i, &temp_n2[0]);' '\n'],...
                Dn_need_vec(ii) - Dn_need_vec(1), sum(o.nc(1:ii-1)) + lb_size + ub_size, ii-1)];
        end

        code = [code, sprintf(['\t' 'sum_nr_constr(&tmp_nu[0]'])];

        if ~isempty(o.K_lb)
            code = [code, ', &temp_lb[0]'];
        end
        if ~isempty(o.K_ub)
            code = [code, ', &temp_ub[0]'];
        end
        if ~isempty(o.K_n)
            code = [code, ', &temp_n2[0]'];
        end
        code = [code, sprintf([', &dot_J[%i]);' '\n'], nu*(ii-1))];

        if (o.contractive || o.terminal)
            code = [code, sprintf(['\t' 'sum_terminal(&tmp_nu[0], &dot_psi_N[%i], &muG[%i]);'...
                '\n'], nu*(ii-1), sum(o.nc) - 1)];
        end

        code = [code, sprintf(['\t' 'minus_scale_nu(&du[%i], &tmp_nu[0]);'...
            '\n'],(ii-1)*nu)];

        code = [code, sprintf(['\t' 'product_matlab_nc(&sl[%d], &muG[%d], %i, &tmp_nc_m[0]);'...
            '\n'],sum(o.nc(1:ii-1)), sum(o.nc(1:ii-1)), ii-1 )];
        code = [code, sprintf(['\t' 'minus_scale_nc(&tmp_nc_m[0], %i, &dsl[%d]);'...
            '\n'], ii-1, sum(o.nc(1:ii-1)) )];
    end
    if (o.contractive || o.terminal)
        code = [code, sprintf(['\t' 'dsl[%d] = ' falcopt.num2str(-alpha, o.precision) '* sl[%d] * muG[%d];' '\n'], sum(o.nc) - 1, sum(o.nc) - 1, sum(o.nc) - 1)];
    end
    
    code = [code, sprintf(['}' '\n\n'])];
    
    
else % the following has to be deleted
    
    code = [code, sprintf(['\t' o.real ' M_tbi[%d], M[%d], Dg[%d], tmp_nc[%d], tmp_nu[%d], Dg_dot_J[%d], ' ...
        'tmp_nc_m[%d];' '\n'], nc*nc, nc*nc, nc*nu, nc, nu, nc, nc)];
    
    code = [code, sprintf(['\t' 'unsigned int ii= 0;' '\n\n',...
        '\t' 'for (ii= 0;ii<%d;ii++){' '\n'],N)];
    code = [code, sprintf(['\t\t' 'build_Dg(&u[ii*%d],Dg);' '\n'],nu)];
    
    code = [code, sprintf(['\t\t' 'build_tbi(&u[ii*%d], &sl[ii*%d], M_tbi);' '\n'],nu,nc)];
    code = [code, sprintf(['\t\t' 'build_inv(M, M_tbi); ' '\n'])];
    
    code = [code, sprintf(['\t\t' 'product_Dgtop_DJ(Dg_dot_J, &dot_J[ii*%d], Dg);'...
        '\n'],nc)];
    
    code = [code, sprintf(['\t\t' 'scale_sub_nc(tmp_nc, &gps[ii*%d], Dg_dot_J); '...
        '\n'],nc)];
    
    code = [code, sprintf(['\t\t' 'product_nc_nc(&muG[ii*%d], tmp_nc, M);'...
        '\n'],nc)];
    code = [code, sprintf(['\t\t' 'product_matrix_nu(tmp_nu, &dot_J[ii*%d], &muG[ii*%d], Dg);'...
        '\n'],nu,nc)];
    code = [code, sprintf(['\t\t' 'minus_scale_nu(&du[ii*%d], tmp_nu);'...
        '\n'],nu)];
    code = [code, sprintf(['\t\t' 'product_matlab_nc(&sl[ii*%d], &muG[ii*%d], tmp_nc_m);'...
        '\n'],nc,nc)];
    code = [code, sprintf(['\t\t' 'minus_scale_nu(&dsl[ii*%d], tmp_nc_m);'...
        '\n'],nc)];
    code = [code, sprintf(['\t' '}' '\n'])];
    code = [code, sprintf(['}' '\n\n'])];
    
end

% multiply the for N the flops ( inside for-cycle)
info.flops = falcopt.multFlops( info.flops, N);

% product_matrix_nc(a,x,B,y, z)
% a is a static scalar, gpy is a nc dynamic vector, Dg is a nu*nc dynamic
% mat and dot_J is nu dynamic
% z = a*x - B'*y

% [d,c] = falcopt.generateMVMult(ones(nu,nc),...
%     'names',struct('fun', 'Btmul', 'M', {{'B'}},...
%     'v', {{'x'}}, 'r', 'r'),...
%     'static',false,...
%     'transpose',true, 'verbose', o.verbose, 'test', o.test);


end

function [code, data, info] = generate_build_gps(o)
code = [];
data = [];

nc = o.nc(1);
nu = o.nu;
N = o.N;

info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

c_psi_dec = argument_def_internal_psi(o,true);
c_contr_dec = argument_contr_value(o,true);

%% build_gpsl

code = [code, sprintf([o.inline ' void build_gpsl_lowLevel('])];
if ~isempty(o.K_lb)
    code = [code, sprintf(['const ' o.real '* amu, '])];
end
if ~isempty(o.K_ub)
    code = [code, sprintf(['const ' o.real '* umb, '])];
end
if ~isempty(o.K_n)
    code = [code, sprintf(['const ' o.real '* n, '])];
end
code = [code,  sprintf(['const unsigned int nc, const unsigned int na, '...
    'const unsigned int nb, const ' o.real '* sl_sqr, ' o.real '* gps){' '\n\n'...
    '\t' 'unsigned int jj;' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code,  sprintf(['\t' 'for (jj=0;jj< na;jj++){' '\n',...
        '\t\t' 'gps[jj]= amu[jj] + 0.5*sl_sqr[jj];' '\n',...
        '\t' '}' '\n'])];
end
if ~isempty(o.K_ub)
    code = [code,  sprintf(['\t' 'for (jj=na;jj< nb;jj++){' '\n'...
        '\t\t' 'gps[jj]= umb[jj-na] + 0.5*sl_sqr[jj];' '\n',...
        '\t' '}' '\n'])];
end
if ~isempty(o.K_n)
    code = [code,  sprintf(['\t' 'for (jj= nb;jj< nc;jj++) {' '\n'...
        '\t\t' 'gps[jj]= n[jj - nb] + 0.5*sl_sqr[jj];' '\n',...
        '\t' '}' '\n'])];
end
code = [code,  sprintf(['}' '\n\n'])];


code = [code, sprintf([o.inline ' void build_gpsl(const ' o.real '* u' c_psi_dec ...
    c_contr_dec ', const ' o.real '* sl_sqr, ' o.real '* gps){' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code, sprintf(['\t'  o.real ' amu[%d];' '\n'], max(sum(~isinf( o.box_lowerBound))) )];
end
if ~isempty(o.K_ub)
    code = [code, sprintf(['\t'  o.real ' umb[%d];' '\n'], max(sum(~isinf( o.box_upperBound))) )];
end
if ~isempty(o.K_n)
    code = [code, sprintf(['\t'  o.real ' n[%d];' '\n'], max(cell2mat(o.nn)) )];
end

if (o.contractive || o.terminal)
    code = [code, sprintf(['\t' o.real ' g_contr = 0.0;' '\n'])];
end


info.flops.mul = info.flops.mul+ N*3; % mul inside only first cycle % NOT CLEAR, ToDo


lb = sum(~isinf(o.box_lowerBound),1);
ub = sum(~isinf(o.box_upperBound),1);




for k = 1:o.N
    code = [code, sprintf(['\n',...
        '\t' '/* Unrolling the for loop: iteration %i of %i */' '\n'], k-1, o.N -1)];
    if ~isempty(o.K_lb)
        code = [code, sprintf(['\t' 'build_amu(&u[%i],%i,&amu[0]);' '\n'], (k-1)*o.nu, k-1)];
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf(['\t' 'build_umb(&u[%i],%i,&umb[0]);' '\n'], (k-1)*o.nu, k-1)];
    end
    if ~isempty(o.K_n)
        code = [code, sprintf(['\t' 'build_n(&u[%i],%i,&n[0]);' '\n'], (k-1)*o.nu, k-1)];
    end
    
    code = [code, sprintf(['\t' 'build_gpsl_lowLevel('])]; %#ok
    if ~isempty(o.K_lb)
        code = [code, sprintf([' &amu[0], '])]; %#ok
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf([' &umb[0], '])]; %#ok
    end
    if ~isempty(o.K_n)
        code = [code, sprintf(' &n[0], ')]; %#ok
    end
    code = [code,  sprintf(['%i, %i, %i, &sl_sqr[%i], &gps[%i]);' '\n'],...
        o.nc(k), lb(:,k), lb(:,k) + ub(:,k), ...
        sum(o.nc(1:k-1)),sum(o.nc(1:k-1)))]; %#ok
end

code = [code, sprintf('\n')];

% code = [code, sprintf(['\t' 'for (ii=0;ii<%d;ii++){' '\n',... % TBD
%     '\t\t' 'build_g(u + ii*%d, g);' '\n',...
%     '\t\t' 'half_sum_nc(gps+ ii*%d, g, sl_sqr + ii*%d);' '\n'...
%     '\t' '}' '\n'],N,nu,nc,nc)];

if (o.contractive || o.terminal)
    code = [code, sprintf(['\t' 'g_contr = *psi_N - c_contr;' '\n'...
        '\t' 'gps[%d] = g_contr + 0.5*sl_sqr[%d];' '\n'],...
        sum(o.nc)-1, sum(o.nc)-1)];
    info.flops.add = info.flops.add + 1;
end


code = [code, sprintf(['}' '\n'])];


info.flops.mul = info.flops.mul+ 3*N;
% The flops of "build_g" are conuted outside this function (in generateConverter)

%% build_sqr

code = [code, sprintf([o.inline ' void build_sqr_Nnc( const ' o.real ' *x, ' o.real ' *x_sqr){' '\n\n'])];

code = [code, sprintf(['\t' 'unsigned int ii = 0;' '\n\n'])];

code = [code, sprintf(['\t' 'for (ii=0;ii<%d;ii++)' '\n',...
    '\t\t' 'x_sqr[ii] = x[ii]*x[ii];' '\n',...
    '}' '\n\n'],sum(o.nc))];

info.flops.mul = info.flops.mul+ sum(o.nc);
end

function [code, data, info] = generate_dot_product_Nnu(o)
% Damian please use your function to generate products

code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

[d, c, in] = falcopt.generateMVMult(ones(1,o.N*o.nu), ...
    'names', struct('fun', 'dot_product_Nnu', 'M', 'u',...
    'v', 'du'), 'static', false, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);
end

function [code, data, info] = generate_dot_product_Nnc(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

[d, c, in] = falcopt.generateMVMult(ones(1,sum(o.nc)), ...
    'names', struct('fun', 'dot_product_Nnc', 'M', 'u',...
    'v', 'du'), 'static', false, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);
end

function [code, data] = generate_copy(o)
% copy_nx

nx = o.nx;
nu = o.nu;
nc = o.nc(1);
N = o.N;

code = [];
data = [];

[d, c] = falcopt.generateMVMult({eye(nx)}, ...
    'names', struct('fun', 'copy_nx', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

[d, c] = falcopt.generateMVMult({eye(N*nx)}, ...
    'names', struct('fun', 'copy_Nnx', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

[d, c] = falcopt.generateMVMult({eye(sum(o.nc))}, ...
    'names', struct('fun', 'copy_Nnc', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

[d, c] = falcopt.generateMVMult({eye(N*nu)}, ...
    'names', struct('fun', 'copy_Nnu', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

end

function check = generateCheck_custom(K, indent, string)

n = length(K);


first = 1;
ind = 1;
K = [K,Inf];

for i = 1:n
    if (K(i) + 1 ~= K(i+1))
        
        last = i;
        if first ~= last
            Out{ind} = [K(first),K(last)]; %#ok
        else
            Out{ind} = K(first);  %#ok
        end
        ind = ind + 1;
        
        first = i+1;
    end
end
rangeStrs = cell(1,length(Out));

for ind = 1:length(Out)
    
    ranges = Out{ind};
    
    if length(ranges) == 1
        rangeStrs{ind} = sprintf([string ' == ' num2str(ranges-1) ]);
    else
        rangeStrs{ind} = sprintf([num2str(ranges(1)-1) ' <= ' string ' <= ' num2str(ranges(2)-1)]);
    end
    
end

check = sprintf([indent 'if(' strjoin(rangeStrs, ' || ') ') {']);
end

% function check = generateCheck(K, indent, string) % TO BE DELETED
% ranges = reshape([K([diff(K)-1 ~= 0 false]); K([false diff(K)-1 ~= 0])],[],1)';
% if length(K) > 1 && K(1)+1 == K(2)
%     ranges = [K(1) ranges]; %#ok
% elseif length(K) == 1
%     ranges = [K(1) K(1)];
% end
% if mod(length(ranges),2) == 1
%     ranges = [ranges K(end)]; %#ok
% end
% ranges = reshape(ranges,2,length(ranges)/2);
% rangeStrs = cell(1,size(ranges,2));
% for j=1:size(ranges,2)
%     if ranges(1,j) == ranges(2,j)
%         rangeStrs{j} = sprintf([string ' == ' num2str(ranges(2,j)-1) ]);
%     elseif ranges(1,j) < ranges(2,j)
%         rangeStrs{j} = sprintf([num2str(ranges(1,j)-1) ' <= ' string ' <= ' num2str(ranges(2,j)-1)]);
%     else
%         throw(MException('ImplementationError', 'Something went wrong. This is a bug, please report to the project owner.'));
%     end
% end
% check = sprintf([indent 'if(' strjoin(rangeStrs, ' || ') ') {']);
% end

function [code, data, info, optCode] = generate_difference(o)
% generate  diff_Nnc

data = [];
code = [];

nc = o.nc(1);
N = o.N;

info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

[d, c, in] = falcopt.generateMVMult({eye(sum(o.nc)), -eye(sum(o.nc))}, ...
    'names', struct('fun', 'diff_Nnc', 'M', {{'I', 'mI'}},...
    'v', {{'x', 'xref'}}, 'r', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.addFlops(info.flops, in.flops);

optCode = [];

end


function [code, data, info] = generate_det_phi(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nc = o.nc(1);
N = o.N;

code = [code, sprintf([o.inline ' void det_phi (const ' o.real ' J, const ' o.real '* gps, const ',...
    o.real '* mu, const ' o.real ' rho, ' o.real '* phi){' '\n\n'])];

code = [code, sprintf(['\t' o.real ' tmp[%d], pr = 0.0;' '\n'],sum(o.nc))];
code = [code, sprintf(['\t' 'unsigned int ii = 0;' '\n\n'])];

code = [code, sprintf(['\t' 'for (ii=%d; ii--; )' '\n'], sum(o.nc))];
code = [code, sprintf(['\t\t' 'tmp[ii] = mu[ii] + 0.50*rho*gps[ii];' '\n'], sum(o.nc))];
info.flops.mul = info.flops.mul+ 2*sum(o.nc);
info.flops.add = info.flops.add+ sum(o.nc);

code = [code, sprintf(['\t' 'dot_product_Nnc(&pr, gps,tmp);' '\n\n'])];
info.flops.add = info.flops.add+ (sum(o.nc)-1); % addition of dot_product_Nnc
info.flops.mul = info.flops.mul+ sum(o.nc); % multiplication of dot_product_Nnc

code = [code, sprintf(['\t' '(*phi) = J + pr;' '\n'])];
info.flops.add = info.flops.add+ sum(o.nc);
code = [code, sprintf(['}' '\n\n'])];

end

function [code, data, info] = generate_det_dot_phi(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nc = o.nc(1);
N = o.N;

code = [code, sprintf([o.inline ' void det_dot_phi (const ' o.real '* du, const ' o.real '* DJ, const ' ...
    o.real ' rho, const '  o.real '* gps, ' '\n',...
    '\t const ' o.real '* mu, const ' o.real '* dm, ' o.real '* dot_phi){' '\n\n'])];

code = [code, sprintf(['\t' o.real ' tmp_prod[%d], prod_1 = 0.0, prod_2 = 0.0;' '\n'],sum(o.nc))];
code = [code, sprintf(['\t' 'unsigned int ii = 0;' '\n\n'])];

code = [code, sprintf(['\t' 'for (ii=%d; ii--; )' '\n'], sum(o.nc))];
code = [code, sprintf(['\t\t' 'tmp_prod[ii] = mu[ii] - dm[ii] + rho*gps[ii];' '\n'], sum(o.nc))];
info.flops.mul = info.flops.mul+ sum(o.nc);
info.flops.add = info.flops.add+ 2*sum(o.nc);

code = [code, sprintf(['\t' 'dot_product_Nnu(&prod_1, du, DJ);' '\n\n'])];
code = [code, sprintf(['\t' 'dot_product_Nnc(&prod_2, gps,tmp_prod);' '\n\n'])];

info.flops.add = info.flops.add+ o.N*(o.nu + sum(o.nc) - 1); % additions of dot_product_Nnc and _Nnu
info.flops.mul = info.flops.mul+ o.N*(o.nu + sum(o.nc)); % multiplications of dot_product_Nnc and _Nnu

code = [code, sprintf(['\t' '(*dot_phi) = prod_1 - prod_2;' '\n'])];

info.flops.add = info.flops.add+ 1;
code = [code, sprintf(['}' '\n\n'])];
end

function [ code, data, info] = generate_conditions_rho(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

alpha = o.stepSize;

code = [code, sprintf([o.inline ' int conditions_rho_PM_simpler (const ' o.real ' dot_phi, const ' o.real ' du_sqr, const ' ...
    o.real ' dsl_sqr, const ' o.real ' alpha){' '\n\n'])];
code = [code, sprintf(['\t' 'unsigned int res = 2;' '\n\n'])];

code = [code, sprintf(['\t' 'if (dot_phi <= ' falcopt.num2str(-0.50/alpha, o.precision) '*(du_sqr + dsl_sqr))' '\n'])];
code = [code, sprintf(['\t\t' 'res = 1;' '\n',...
    '\t' 'else' '\n',...
    '\t\t' 'res = 0;' '\n\n',...
    '\t' 'return res;' '\n'...
    '}' '\n\n' ])];
info.flops.add = info.flops.add+ 1;
info.flops.comp = info.flops.comp +1;


end

function [code, data, info] = generate_Lagrangian_oracles_lp(o)

code = [];
data = [];
info.flops.it = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
info.flops.ls = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nc = o.nc(1);
N = o.N;

%% define one norm and inf norm of dimension N*nc

if o.merit_function ~= 2
    
    code = [code, sprintf([o.inline ' ' o.real ' one_norm (const ' o.real '* g){' '\n\n',...
        '\t' 'int ii = 0;' '\n'...
        '\t' o.real ' norm = 0.0;' '\n\n',...
        '\t' 'for (ii = %d; ii-- >0; ){' '\n'],sum(o.nc))];
    
    code = [code, sprintf(['\t\t' 'norm += ' o.abs '(g[ii]);' '\n',...
        '\t' '}' '\n',...
        '\t' 'return norm;' '\n',...
        '}' '\n\n'])];
    
    
    code = [code, sprintf([o.inline ' ' o.real ' inf_norm (const ' o.real '* g){' '\n\n',...
        '\t' 'int ii = 0;' '\n'...
        '\t' o.real ' norm = 0.0;' '\n\n',...
        '\t' 'for (ii = %d; ii-- >0; ){' '\n'],sum(o.nc))];
    
    code = [code, sprintf(['\t\t' 'norm = ' o.max '(norm,' o.abs '(g[ii]));' '\n',...
        '\t' '}' '\n',...
        '\t' 'return norm;' '\n',...
        '}' '\n\n'])];
    
else
    
    code = [code, sprintf([o.inline ' ' o.real ' two_norm (const ' o.real '* g){' '\n\n',...
        '\t' 'int ii = 0;' '\n'...
        '\t' o.real ' norm = 0.0, norm2 = 0.0;' '\n\n',...
        '\t' 'for (ii = %d; ii-- >0; ){' '\n'],sum(o.nc))];
    
    code = [code, sprintf(['\t\t' 'norm2 += g[ii]*g[ii];' '\n',...
        '\t' '}' '\n',...
        '\t' 'norm = sqrt(norm2);' '\n',...
        '\t' 'return norm;' '\n',...
        '}' '\n\n'])];
    
end

%% det_phi
code = [code, sprintf([o.inline ' void det_phi (const ' o.real ' J, const ' o.real '* gps, const ',...
    o.real ' rho, ' o.real '* phi){' '\n\n'])];

if o.merit_function == 1
    code = [code, sprintf(['\t' '(*phi) = J + rho * one_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
elseif o.merit_function == Inf
    code = [code, sprintf(['\t' '(*phi) = J + rho * inf_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
elseif o.merit_function == 2
    code = [code, sprintf(['\t' '(*phi) = J + rho * two_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
end

code = [code, sprintf(['}' '\n\n'])];
info.flops.ls = info.flops.it; %det_phi called in ls too.

%% det_dot_phi


code = [code, sprintf([o.inline ' void det_dot_phi (const ' o.real '* du, const ' o.real '* DJ, const ' ...
    o.real ' rho, const '  o.real '* gps, ' o.real '* dot_phi){' '\n\n'])];

code = [code, sprintf(['\t' o.real ' pr = 0.0;' '\n'])];

code = [code, sprintf(['\t' 'dot_product_Nnu(&pr, du, DJ);' '\n\n'])];
info.flops.it.add = info.flops.it.add+ o.N*(o.nu-1); % addition of dot_product_Nnu
info.flops.it.mul = info.flops.it.mul+ o.N * o.nu; % multiplication of dot_product_Nnu

if o.merit_function == 1
    code = [code, sprintf(['\t' '(*dot_phi) = pr - rho*one_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
elseif o.merit_function == Inf
    code = [code, sprintf(['\t' '(*dot_phi) = pr - rho*inf_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
elseif o.merit_function == 2
    code = [code, sprintf(['\t' '(*dot_phi) = pr - rho*two_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
end

code = [code, sprintf(['}' '\n\n'])];

%% conditions_rho

code = [code, sprintf([o.inline ' void update_rho (const ' o.real '* muG, ' o.real '* rho, ' o.real '* rho_hat){' '\n\n'])];

if o.merit_function == 1
    code = [code, sprintf(['\t' '(*rho_hat) = inf_norm(muG);' '\n'])];
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
elseif o.merit_function == Inf
    code = [code, sprintf(['\t' '(*rho_hat) = one_norm(muG);' '\n'])];
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
elseif o.merit_function == 2
    code = [code, sprintf(['\t' '(*rho_hat) = two_norm(muG);' '\n'])];
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
end

code = [code, sprintf(['\t' '(*rho) = ' o.max '( (*rho), (*rho_hat)); ' '\n'])];
info.flops.it.comp = info.flops.it.comp+ 1;

code = [code, sprintf(['}' '\n\n' ])];

end

function [code, data, info] = generate_weighted_sum_nc(o)
code = '';
data = '';
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

[d, c, in] = falcopt.generateMVMult({eye(sum(o.nc)), eye(sum(o.nc))}, ...
    'names', struct('fun', 'weighted_sum_Nnc', 'M',...
    {'I'}, 'v', {{'x', 'xref'}}, 'r', 'dx'),...
    'static', struct('M',[false,true]),...
    'types', o.real, 'precision', o.precision, 'structure', 'unique', ...
    'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
info.flops = falcopt.addFlops(info.flops, in.flops);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
end

function [code, data, info] = generate_weighted_sum_nu(o)
code = '';
data = '';
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
[d, c, in] = falcopt.generateMVMult({eye(o.N*o.nu), eye(o.N*o.nu)}, ...
    'names', struct('fun', 'weighted_sum_Nnu', 'M',...
    {'I'}, 'v', {{'x', 'xref'}}, 'r', 'dx'),...
    'static', struct('M',[false,true]),...
    'types', o.real, 'precision', o.precision, 'structure', 'unique', ...
    'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
info.flops = falcopt.addFlops(info.flops, in.flops);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

end

function [code, data, info] = generate_quadratic_interpolation(o)
code = '';
data = '';
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0, 'div', 0);

% TODO: make parameter struct
p.tau = [0.01 1-0.01];
code = [code, sprintf([o.indent.code o.inline ' ' o.real ' quadratic_interp (const ' o.real ' f_l, const ', o.real ' g_l, const ' o.real ' t_u, const ' o.real ' f_u) {' '\n' '\n'])];
code = [code, sprintf([o.indent.code o.indent.generic o.real ' t_theo;' '\n' '\n'])];
code = [code, sprintf([o.indent.code o.indent.generic 't_theo = -0.5*(g_l*t_u*t_u)/(f_u - f_l - g_l*t_u);' '\n'])];
info.flops.mul = info.flops.mul+5;
info.flops.add = info.flops.add+2;
info.flops.div = info.flops.div +1;
% code = [code, sprintf([o.indent.code o.indent.generic 'a = 0.1*t_u;'
% '\n'])]; % TBD
% info.flops.mul = info.flops.mul+1;
code = [code, sprintf([o.indent.code o.indent.generic 'return ' o.max '(' o.min '(t_theo,' falcopt.num2str(p.tau(2), o.precision) '*t_u), ' falcopt.num2str(p.tau(1), o.precision) '*t_u);' '\n'])];
info.flops.mul = info.flops.mul+2;
info.flops.add = info.flops.add+1;
info.flops.comp = info.flops.comp+ 3;
code = [code, sprintf([o.indent.code '}'])];
end

function [code, data, info] = generate_compute_max(o)
code = '';
data = '';
code = [code, sprintf([o.indent.code o.inline ' ' o.real ' compute_max_Nnc( const ' o.real '* x) {' '\n' ...
    o.indent.code o.indent.generic 'unsigned int ii = 0;' '\n',...
    o.indent.code o.indent.generic o.real ' m = -100.0;' '\n' '\n'])];
code = [code, sprintf([o.indent.code o.indent.generic 'for (ii=%d; ii--; ) {' '\n'],sum(o.nc))];
code = [code, sprintf([o.indent.code o.indent.generic o.indent.generic 'm = ' o.max '(m,x[ii]);' '\n' ...
    o.indent.code o.indent.generic '}' '\n' ...
    o.indent.code o.indent.generic 'return m;' '\n' ...
    o.indent.code '}'])];
info.flops.comp = sum(o.nc);
end

