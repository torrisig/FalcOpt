%% generateCode Generate an embeddable FalcOpt solver in C
% tailored to your model predictive control (MPC) problem.
%
% Given a nonlinear MPC problem of the following form:
% min_{x,u} sum_{k=0}^{N-1} l(xk,uk) + lN(xN)
%     s.t.  xk+1 = f(xk,uk)        for k=0,...,N-1
%           x0 = initial state
%           uk in [ak,bk]          for k=0,...,N-1
%           nk(uk) <= 0            for k=0,...,N-1
%           0.5*xN'*D*xN <= c
%
% with prediction horizon N, state dimension nx and input dimension nu.
%
% This function can be called as
%
% [info] = generateCode(dynamics, N, nx, nu, objective) [with options]
% or
% [info] = generateCode(dynamics, N, nx, nu, objective, ...) [with options]
%
% it generates C code and a MEX file in the folder specificed with the option
%  'gendir' (default: '../'), where the file name of the main .c file is given by
%  'name' (default: 'my_code').
% it returns a structure containing TODO
%
% The required inputs are:
%  dynamics  - A function handle for the system dynamics of the form 
%               d/dx = f = @(x,u)( ... ), e.g. for linear dynamics
%               d/dx = f = @(x,u)(A*x + B*u).
%               Dynamics can be made parameter-dependent via the option 'nw' (see below).
%  N         - The prediction horizon, a positive integer.
%  nx        - The dimension of the state vector at time k.
%  nu        - The dimension of the input vector at time k.
%  objective - A struct defining the MPC objective function, i.e. defining
%                l(xk,uk) and l(xN).
%               Either struct('Q', Q, 'R' R, 'P', P) for quadratic objective, where
%                .Q a matrix of dimension nx times nx is the state stage objective,
%                .R a matrix of dimension nu times nu is the input stage objective,
%                .P a matrix of dimension nx times nx is the final state objective,
%                 i.e., l(xk,uk) = xk'*Q*xk + uk'*R*uk and lN(xN) = xN'*P*xN.
%               Or struct('nonlinear', Q, 'R' R, 'P', P) for arbitrary nonlinear objective
%                 .nonlinear a function handle for the nonlinear stage cost function
%                   of the form l = @(x,u)( ... ),
%                 .nonlinearN a function handle for the nonlinear final state cost function
%                   of the form lN = @(x)( ... ).
%               Additionally, 
%                .trackReference a boolean to enable reference tracking, to track a
%                  time-varying reference. Default: false.
%
% The following options are available and can be set as key-value pairs, e.g. 'optionName', optionValue
% 0) General settings
%   verbose - Level of procedural output of the code generation. Default: 0.
%   test    - Level of internal numerical tests performed. Default: 0.
%   debug   - Level of debug (relates to number of inputs/outputs of generated functions). Default: 1.
% 1) Problem constraints
%   lb          - Lower bound constraints on the inputs. Either
%                  a vector of length nu, or
%                  a matrix of dimension nu times N for individual bounds for each stage.
%                  Default: [] (no bounds).
%   ub          - Upper bound constraints on the inputs. Either
%                         a vector of length nu, or
%                         a matrix of dimension nu times N for individual bounds for each stage.
%                         Default: [] (no bounds).
%   constraints - Nonlinear input stage constraints nk(uk). Either
%                  a function handle if each stage has the same constraint, i.e., nk(uk)=n(uk), or
%                  a cell of N function handles, with individual constrains for each stage.
%                  The function handles are of form n = @(u)(...). Default: {} (mo constraints).
%   nn          - The dimension of the constraint vector n, or
%                  a vector of length N with dimensions for the onstraint vectors nk, for each k.
%   terminal    - TODO
%   contractive - TODO
% 2) Parameters/known disturbances
%     Parameters entering the dynamics can be added by specifying the folowing option:
%   nw          - Parameter/disturbance dimension. Default: 0.
%     If nw > 0 then the dynamics need to be a function of w, i.e. d/dx = f = @(x,u,w)( ... ),
%      e.g. for linear dynamics d/dx = f = @(x,u,w)(A*x + B*u + C*w).
% 3) Code generation settings
%   name       - Name of the generated solver. Will also be the name of the .c, .mex and simulink file (if generated).
%                 Default: 'my_code'.
%   gendir     - Target directory for generated files. Needs to be a valid (existing) directory. Default: './'.
%   precision  - The precision of the computation and data used. Can be either 'double' or 'single'. Defaukt: 'double'.
%   buildTypes - A cell array of build types to be generated, can include
%                 'mex', 'simulink', 'production'. Default: 'mex'.
%   compile    - Whether the generated code is compiled. Default: true.
%   gradients  - How to do code generation for function evaluation and gradient/Jacobian computation.
%                 Can be 'matlab' using MATLAB code generation,
%                        'casadi' using CasADi (https://github.com/casadi/casadi/wiki),
%                        'manual' where user supplied functions for gradient/Jacobian computation
%                          is translated into C code using MATLAB code generation. See 'external' below.
%                        'ccode', where user supplied C code for gradient/Jacobian computation
%                          is incorporated into the generated code.
%                 Default: 'matlab'.
%   external   - A struct containing needed data for the 'manual' or 'ccode' gradients setting.
%                 for 'manual' .manual.jacobian_x is a function @(x,u)(...) that returns
%                                the Jacobian of the dynamics with respect to x.
%                              .manual.jacobian_u is a function @(x,u)(...) that returns
%                                the Jacobian of the dynamics with respect to u.
%                              .manual.jacobian_n is a function (or cell of functions) @(u)(...) that returns
%                                the Jacobian of n(uk) (or Jacobians of nk(uk)).
%                 for 'ccode' .ccode.c is a file with C code implementing the evaluation of the dynamics and Jacobian computations.
%                             .ccode.h is the corresponding header file.
%                             .ccode.jacobian_x.structure is a Matrix containing the structure of
%                               the Jacobian of the dynamics with respect to x. Default: no structure.
%                             .ccode.jacobian_x.static is a Boolean, which is true if the Jacobian is static,
%                               i.e. if it does not depend on x and u. Default: true.
%                             .ccode.jacobian_u.structure is a Matrix containing the structure of
%                               the Jacobian of the dynamics with respect to u. Default: no structure.
%                             .ccode.jacobian_u.static is a Boolean, which is true if the Jacobian is static,
%                               i.e. if it does not depend on x and u. Default: true.
%                             .ccode.jacobian_n.structure is a Matrix containing the structure of
%                               the Jacobian of n(uk).
%                             .ccode.K_n a cell array containing the vectors of the different stages
%                               to which the nonlinear constraints apply.
%   * Internal code generation settings
%   indent     - Indentation used in code generation. Default: struct('code', '', 'data', '', 'generic', '\t' ).
%   inline     - Inline keyword to be used. Default: ''.
%   names      - A struct used for internal naming. Do not modify unless you know what you are doing.
% 4) Algorithm parameters, tolerance and max. iteration settings
%   eps               - Tolerance on KKT optimality. Default: 1e-3.
%   maxIt             - Max number of iterations. Default: 4000.
%   maxItLs           - Max number of line search iterations. Default: 10
%   variable_stepSize - struct with fields: 
%                        .active a boolean that activates variable step size, via a trust-region procedure
%                          if false is set to constant (alpha_max specified next). Default: true
%                        .alpha_max maximum step size. If not specified, it is computed via second order 
%                          information around the assumed equilibrium x = 0, u = 0 (requires CasADi)
%                        .alpha_min minimum step size. If not specified, is set to 0.1*alpha_max
%                        .steady_state_state if specified, alpha_max is computed around this equilibrium
%                          otherwise around x = 0.
%                        .steady_state_input if specified, alpha_max is computed around this equilibrium
%                          otherwise around u = 0.
%                        .increase_threshold define upper bound for trust-region. Default: 0.75.
%                        .increase_coeff coefficient >1 by which we increase the trust-region. Default: 1.33.
%                        .decrease_threshold lower bound for the trust-region: Default: 0.25.
%                        .decrease_coeff coefficient <1 by which we decrease the trust-region. Default: 0.75.
%   merit_function    - Merit function selection
%                        0: Augmented Lagrangian
%                        1, 2 or Inf: l_1, l_2(default) or l_Inf non-smooth penalty function
%                        Default: 2.
%   parLs             - Armijo line search step parameter. Default: 0.3.
%   tolLs             - Line search min progress. Default: 1e-4.

% Copyright (c) 2017, ETH Zurich, Automatic Control Laboratory 
%                    Giampaolo Torrisi <giampaolo.torrisi@gmail.com>
%                    Tommaso Robbiani <tommasro@student.ethz.ch>
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

function [info] = generateCode(varargin)
    gradientTypes = {'matlab', 'casadi', 'ccode', 'manual'};
    indentTypes = {'generic', 'data', 'code'};
    defaultNames = struct('maximumIterations', 'maximumIterations', ...
                          'maximumLineSearchIterations', 'maximumLineSearchIterations', ...
                          'lineSearchAlphaMax', 'lineSearchAlphaMax', ...
                          'lineSearchAlphaMin', 'lineSearchAlphaMin', ...
                          'lineSearchAlphaTol', 'lineSearchAlphaTol', ...
                          'lineSearchAlphaTolSq', 'lineSearchAlphaTolSq', ...
                          'lineSearchAlphaSqTol', 'lineSearchAlphaSqTol', ...
                          'lineSearchAlphaSqTolSq', 'lineSearchAlphaSqTolSq', ...
                          'KKTOptimalityTolerance', 'KKTOptimalityTolerance', ...
                          'KKTOptimalityToleranceSq', 'KKTOptimalityToleranceSq');
    
    p = inputParser;
    % Required inputs
    p.addRequired('dynamics', @(x)(isa(x,'function_handle')));
    p.addRequired('N', @(x)(isnumeric(x) && numel(x) == 1 && x > 0 && mod(x,1) == 0)); % prediction horizon
    p.addRequired('nx', @(x)(isnumeric(x) && numel(x) == 1 && x > 0 && mod(x,1) == 0)); % state dimension
    p.addRequired('nu', @(x)(isnumeric(x) && numel(x) == 1 && x > 0 && mod(x,1) == 0)); % input dimension
    p.addRequired('objective', @isstruct);
    % Optional inputs: 0) General settings
    p.addParameter('verbose', 0, @isnumeric);
    p.addParameter('test', 0, @isnumeric);
    p.addParameter('debug', 1, @isnumeric);
    % Optional inputs: 1) Problem constraints
    p.addParameter('lb', [], @isnumeric);
    p.addParameter('ub', [], @isnumeric);
    p.addParameter('constraints', {}, @(x) (iscell(x) || isa(x,'function_handle')));
    p.addParameter('nn', 0, @(x) (isnumeric(x) || iscell(x)));
    p.addParameter('terminal', false, @(x)(islogical(x) && numel(x) == 1));
    p.addParameter('contractive', false, @(x)(islogical(x) && numel(x) == 1));
        % Depricated options
        p.addParameter('box_lowerBound', [], @isnumeric); % Kept for legacy purposes, proper option is 'lb'
        p.addParameter('box_upperBound', [], @isnumeric); % Kept for legacy purposes, proper option is 'ub'
        p.addParameter('constraints_handle', [], @(x) (iscell(x) || isa(x,'function_handle'))); % Kept for legacy purposes, proper option is 'constraints'
    % Optional inputs: 2) Parameters/known disturbances
    p.addParameter('nw', 0, @(x) (isnumeric(x) && numel(x) == 1 && x>=0 && mod(x,1) == 0)); % Dimension of parameter/known disturbance
    % Optional inputs: 3) Code generation settings
    p.addParameter('name', 'my_code', @ischar);
    p.addParameter('gendir', './', @ischar);
    p.addParameter('precision', 'double', @(x)( strcmp(x,'double')|| strcmp(x,'single') ));
    p.addParameter('buildTypes', {'mex'}, @(x)(iscell(x) && all(cellfun(@(c)(any(strcmp(c, {'mex', 'simulink', 'production'}))), x))));
    p.addParameter('compile', true, @islogical);
    p.addParameter('gradients', 'matlab', @(x)(ischar(x) && any(strcmp(x, gradientTypes))));
    p.addParameter('external', struct(), @isstruct);
        % Depricated options
        p.addParameter('external_jacobian_x', [], @(x)isa(x,'function_handle')); % Kept for legacy purposes, proper option is 'external' with external.manual.jacobian_x
        p.addParameter('external_jacobian_u', [], @(x)isa(x,'function_handle')); % Kept for legacy purposes, proper option is 'external' with external.manual.jacobian_u
        p.addParameter('external_jacobian_n', [], @(x)isa(x,'function_handle')); % Kept for legacy purposes, proper option is 'external' with external.manual.jacobian_n
        p.addParameter('Jac_x_static',[],@islogical); % Kept for legacy purposes, proper option is 'external' with external.ccode.jacobian_x.static
        p.addParameter('Jac_u_static',[],@islogical); % Kept for legacy purposes, proper option is 'external' with external.ccode.jacobian_u.static
        p.addParameter('Jac_x_struct',[],@isnumeric); % Kept for legacy purposes, proper option is 'external' with external.ccode.jacobian_x.structure
        p.addParameter('Jac_u_struct',[],@isnumeric); % Kept for legacy purposes, proper option is 'external' with external.ccode.jacobian_u.structure
        p.addParameter('Jac_n_struct',[],@isnumeric); % Kept for legacy purposes, proper option is 'external' with external.ccode.jacobian_n.structure
        p.addParameter('K_n',{},@iscell); % Kept for legacy purposes, proper option is 'external' with external.ccode.K_n
    p.addParameter('names', defaultNames, @isstruct); % Internal naming
    p.addParameter('indent', struct('code', '', 'data', '', 'generic', '\t' ), @(x)(ischar(x) || (isstruct(x) && isfield(x, 'generic') && all(cellfun(@(y)(~isfield(x, y) || ischar(x.(y))), indentTypes)))));
    p.addParameter('inline', '', @ischar);
    % Optional inputs: 4) Algorithm parameters
    p.addParameter('eps', 1e-3, @(x)(isnumeric(x) && numel(x) == 1 && x >= eps)); % tolerance
    p.addParameter('maxIt', 4000, @(x)(isnumeric(x) && numel(x) == 1 && x>0 && mod(x,1) == 0)); % max number of iter
    p.addParameter('maxItLs', 10, @(x)(isnumeric(x) && numel(x) == 1 && x>0 && mod(x,1) == 0)); % max number of line search iterations
    p.addParameter('variable_stepSize', struct(), @isstruct);
    p.addParameter('merit_function', 0, @(x)( isnumeric(x) && numel(x) == 1 && (x == 1 || x == 2 || x == Inf || x == 0)));
    p.addParameter('parLs', 0.3, @(x)(isnumeric(x) && numel(x) == 1 && x > 0 && x < 1)); % Armijo line search step parameter
    p.addParameter('tolLs', 1e-4, @(x)(isnumeric(x) && numel(x) == 1 && x >= eps)); % line search min progress
    p.parse(varargin{:});
    o = p.Results;
    
    % TODO this also requires YALMIP, need to check and throw an error early!

    %% Processing options
    o = falcopt.processCodeGenOptions(o, defaultNames, indentTypes);
    
    o.K_amu = falcopt.internal.findNInf(o.lb);
    cons_lb = ~isinf(o.lb);
    o.K_lb = falcopt.internal.findNInf(cons_lb(:,sum(cons_lb,1) >= 1));
    o.K_umb = falcopt.internal.findNInf(o.ub);
    cons_ub = ~isinf(o.ub);
    o.K_ub = falcopt.internal.findNInf(cons_ub(:,sum(cons_ub,1) >= 1));

%% Setup info
info.names = o.names;
info.flops.it = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0, 'div', 0, 'casadi', 0);
info.flops.ls = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0, 'div', 0, 'casadi', 0);
info.flops.once = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0, 'div', 0, 'casadi', 0);
info.src = {};
info.header = {};

%% Generating code
copyright = sprintf(['/* Copyright (c) 2017, ETH Zurich, Automatic Control Laboratory */' '\n' ...
    '/* Main developers: */' '\n' ...
    '/* Giampaolo Torrisi <giampaolo.torrisi@gmail.com> */' '\n' ...
    '/* Damian Frick <falcopt@damianfrick.com> */' '\n' ...
    '/* Tommaso Robbiani <tommasro@student.ethz.ch> */' '\n' ...
    '/* Scientic Contributors: */' '\n' ...
    '/* Sergio Grammatico */' '\n' ...
    '/* Roy S. Smith */' '\n' ...
    '/* Manfred Morari */' '\n']);

libr = sprintf(['#include "math.h"' '\n']);
defs = ''; % Defines
code = '';
data = '';
optCode = '';

%% Generating max and min

code = [code, sprintf([ o.indent.code o.real ' ' o.auxFunctions.max '(' o.real ' x, ' o.real ' y){' '\n' ...
    o.indent.code o.indent.generic 'if (x < y) { return y; }' '\n' ...
    o.indent.code o.indent.generic 'else { return x; }' '\n' ...
    o.indent.code '}' '\n\n'])];
code = [code, sprintf([ o.indent.code o.real ' ' o.auxFunctions.min '(' o.real ' x, ' o.real ' y){' '\n' ...
    o.indent.code o.indent.generic 'if (x > y) { return y; }' '\n' ...
    o.indent.code o.indent.generic 'else { return x; }' '\n' ...
    o.indent.code '}' '\n\n'])];



%% CasaDi model_mpc, Jacobian_u and Jacobian_x automatic generation
switch o.gradients
    case 'casadi'
        [d,c,i] = casadi_jacobians( o, o.dynamics,'model_mpc','jac',{'x','u'}); % generate casadi functions
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
        [d,c,i] = matlab_jacobians(o,o.dynamics,'model_mpc','jac',{'x','u'});   % generate function and jacobians
        data = [data, d];
        code = [code, c];
        
        o.Jac_x_static = i.in_F.static;
        o.Jac_x_struct = i.in_F.struct.structure.mat;                        % structure of jacobian_x
        o.Jac_u_static = i.in_G.static;
        o.Jac_u_struct = i.in_G.struct.structure.mat;                        % structure of jacobian_u
        
        info.flops.once = falcopt.internal.addFlops( info.flops.once, falcopt.internal.multFlops(i.y.flops,o.N+1));  % flops model_mpc
        info.flops.ls = falcopt.internal.addFlops( info.flops.ls, falcopt.internal.multFlops(i.y.flops,o.N+1));      % flops model_mpc
        info.flops.it = falcopt.internal.addFlops( info.flops.it, falcopt.internal.multFlops(i.in_F.flops,o.N));     % flops jacobian_x
        info.flops.it = falcopt.internal.addFlops( info.flops.it, falcopt.internal.multFlops(i.in_G.flops,o.N));     % flops jacobian_u
    case 'manual'
        
        % generate external function (matlab)
        [ d, c, i] = external_jacobians(o);
        data = [data, d];
        code = [code, c];
        
        o.Jac_x_static = i.in_F.static;
        o.Jac_x_struct = i.in_F.struct.structure.mat;                        % structure of jacobian_x
        o.Jac_u_static = i.in_G.static;
        o.Jac_u_struct = i.in_G.struct.structure.mat;                        % structure of jacobian_u
        
        info.flops.once = falcopt.internal.addFlops( info.flops.once, falcopt.internal.multFlops(i.y.flops,o.N+1));  % flops model_mpc
        info.flops.ls = falcopt.internal.addFlops( info.flops.ls, falcopt.internal.multFlops(i.y.flops,o.N+1));      % flops model_mpc
        info.flops.it = falcopt.internal.addFlops( info.flops.it, falcopt.internal.multFlops(i.in_F.flops,o.N));     % flops jacobian_x
        info.flops.it = falcopt.internal.addFlops( info.flops.it, falcopt.internal.multFlops(i.in_G.flops,o.N));     % flops jacobian_u
        
    case 'ccode'
        % all functions will be provided in external_functions.c file
        if( exist(o.external.ccode.c, 'file'))
            filename = strsplit(o.external.ccode.c, {'/','\'});
            filename = filename{end};
            copyfile(o.external.ccode.c,[o.gendir '/' filename]);
            info.src = [info.src; {sprintf(o.external.ccode.c)}];
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
                'external.ccode.jacobian_x, external.ccode.jacobian_u and external.ccode.jacobian_n parameters must be provided \n' ];
            warning(sprintf(string)); %#ok
            
            % add .h file to library
            if(exist(o.external.ccode.h, 'file'))
                filename = strsplit(o.external.ccode.h, {'/','\'});
                filename = filename{end};
                copyfile(o.external.ccode.h,[o.gendir '/' filename]);
                libr = [libr, sprintf(['#include "' filename '" \n'])];
                info.header = [info.header; {filename}];
            else % TODO including a .c file is STRANGE
                error(['File ' o.external.ccode.h ' for gradient/jacobian computation not found.']);
            end
            data = [ data, sprintf([o.indent.data '/*static data for jacobian_u*/' '\n'...
                o.indent.data 'static double G[%d];' '\n'],max(o.Jac_u_struct(:)))];
            data = [ data, sprintf([o.indent.data '/*static data for jacobian_x*/' '\n' ...
                o.indent.data 'static double F[%d];' '\n'],max(o.Jac_x_struct(:)))];
        else
            error(['File ' o.external.ccode.c ' for gradient/jacobvian computation not found.']);
        end
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
            o.Jac_n_struct{jj} = i.in_Dn_n{jj}.struct.structure.stored.mat;     % structure of Dn
            info.flops.once.casadi = info.flops.once.casadi + i.in_n{jj}.flops*length(o.K_n{jj});     % flops build_n
            info.flops.ls.casadi = info.flops.ls.casadi + i.in_n{jj}.flops*length(o.K_n{jj});         % flops build_n
            info.flops.it.casadi = info.flops.it.casadi + i.in_Dn_n{jj}.flops*length(o.K_n{jj});      % flops build_Dn
        end
        try %#ok
            info.src = [info.src; i.src];                                           % external src file
            info.header = [info.header; i.header];                                  % external header file
            libr = [libr, sprintf(['#include "' i.header{:} '"' '\n'])];
        end
        
    case {'matlab','manual'}
        % generate build_g, build_Dg, build_Dn functions
        [ d, c, i] = generate_n_and_Dn( o, 'matlab');
        code = [code, c];
        data = [data, d];
        
        o.nc = i.nc;
        
        for jj=1:length(o.K_n)
            o.Jac_n_struct{jj} = i.in_Dn_n{jj}.struct.structure.mat;        % structure of Dn
            
            info.flops.once = falcopt.internal.addFlops( info.flops.once, falcopt.internal.multFlops(i.in_n{jj}.flops,length(o.K_n{jj})));  % flops build_n
            info.flops.ls = falcopt.internal.addFlops( info.flops.ls, falcopt.internal.multFlops(i.in_n{jj}.flops,length(o.K_n{jj})));      % flops build_n
            info.flops.it = falcopt.internal.addFlops( info.flops.it, falcopt.internal.multFlops(i.in_Dn_n{jj}.flops,length(o.K_n{jj})));   % flops build_Dn
        end
        
    case 'ccode'
        [ d, c, i] = generate_n_and_Dn( o,'ccode');
        code = [code, c];
        data = [data, d];
        o.nc = i.nc;
        
end


% find the structure of the current Jac_n
o.Jac_n_struct_hor = cell(1,o.N);

if isequal(o.gradients,'ccode')
    o.Jac_n_struct_hor = o.Jac_n_struct;
    
else
    
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
end

% check patterns of stage constraints
o.K_nc = falcopt.internal.findNInf(o.nc(1:o.N));


%% Initialization main algorithm

c_w_dec = argument_w(o,true);
c_tr_dec = argument_def(o,true);
c_contr_dec = argument_contr_value(o,true);
optCode = [optCode, sprintf(['/* main function of the algorithm */' '\n'])];
if o.debug == 3
    optCode = [optCode, sprintf(['int proposed_algorithm(const ' o.real '* x0, ' o.real '* u' c_w_dec c_tr_dec c_contr_dec ...
        ', ' o.real '* x, ' o.real '* fval, unsigned int* iter, unsigned int* iter_ls, ' o.real '* optimval, ' o.real '* feasval, ' o.real '* meritval){' '\n' '\n'])];    
elseif o.debug == 2
    optCode = [optCode, sprintf(['int proposed_algorithm(const ' o.real '* x0, ' o.real '* u' c_w_dec c_tr_dec c_contr_dec ...
        ', ' o.real '* x, ' o.real '* fval, unsigned int* iter, unsigned int* iter_ls){' '\n' '\n'])];
else
    optCode = [optCode, sprintf(['int proposed_algorithm(const ' o.real '* x0, ' o.real '* u' c_w_dec c_tr_dec c_contr_dec ...
        ', unsigned int* iter, unsigned int* iter_ls){' '\n' '\n'])];
end


%% Definitions fo parameters for main algorithm
defs = [defs, sprintf(['/** FALCOPT parameters **/' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.maximumIterations) ' (' falcopt.internal.num2str(o.maxIt) ')' ' /* Maximum number of iterations */' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) ' (' falcopt.internal.num2str(o.eps, o.precision) ')' ' /* KKT optimality tolerance */' '\n'])];
defs = [defs, sprintf(['/* Derived parameters DO NOT TOUCH */' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.KKTOptimalityToleranceSq) ' (' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) '*' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) ')' '\n'])];

%% Declaration of variables
optCode = [optCode, sprintf([o.indent.generic  'unsigned int conditions_f = 1, conditions_x = 1,' '\n' ...
    o.indent.generic o.indent.generic 'conditions_n = 1, cond_err = 1,' '\n',...
    o.indent.generic o.indent.generic 'it = 0, it_ls = 0, ii = 0, jj = 0;' '\n'...
    o.indent.generic  'int cond = -2;' '\n'])];

if o.merit_function == 0
    optCode = [optCode, sprintf([o.indent.generic  'unsigned int flag_ini = 0;' '\n'])];
else
    optCode = [optCode, sprintf([o.indent.generic  'unsigned int reset_rho = 0;' '\n'])];
end

optCode = [optCode, sprintf([o.indent.generic o.real ' J= 0.0, Jt = 0.0, constr_viol = 1.0;' '\n'])];

if o.variable_stepSize.active
    data = [data, sprintf('\n/* static data for variable step alpha */ \n')];
    data = [data, sprintf([o.indent.data 'static ' o.real ' alpha = ' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) ';\n'])];
    data = [data, sprintf([o.indent.data 'static ' o.real ' alpha_old = ' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) ';\n'])];
    data = [data, sprintf([o.indent.data 'static ' o.real ' alpha_inverse = 1.0/' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) ';\n'])];
    data = [data, sprintf([o.indent.data 'static const ' o.real ' alpha_max = ' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) ';\n'])];
    data = [data, sprintf([o.indent.data 'static const ' o.real ' alpha_min = ' falcopt.internal.toDefineName(o.names.lineSearchAlphaMin) ';\n'])];
    data = [data, sprintf([o.indent.data 'static const ' o.real ' rat_1 = ' falcopt.internal.num2str(o.variable_stepSize.decrease_threshold,o.precision) ';\n'])];
    data = [data, sprintf([o.indent.data 'static const ' o.real ' rat_2 = ' falcopt.internal.num2str(o.variable_stepSize.increase_threshold,o.precision) ';\n'])];
    data = [data, sprintf([o.indent.data 'static const ' o.real ' gamma_1 = ' falcopt.internal.num2str(o.variable_stepSize.decrease_coeff,o.precision) ';\n'])];
    data = [data, sprintf([o.indent.data 'static const ' o.real ' gamma_2 = ' falcopt.internal.num2str(o.variable_stepSize.increase_coeff,o.precision) ';\n\n'])];

    optCode = [optCode, sprintf([o.indent.generic o.real ' ared = 0.0, pred = 0.0, rat = 0.0' ';\n'])];
end

if o.debug > 1
    optCode = [optCode, sprintf([o.indent.generic o.real ' xp[%d],' '\n'],o.N*o.nx)];
else
    optCode = [optCode, sprintf([o.indent.generic o.real ' x[%i], xp[%d],' '\n'],o.N*o.nx,o.N*o.nx)];
end
optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'dot_J[%d], du[%d], up[%d],' '\n'],o.N*o.nu,o.N*o.nu,o.N*o.nu)];
optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'gps[%d], gpsp[%d],' '\n'],sum(o.nc),sum(o.nc))];
optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'sl[%d], sl_sqr[%d], dsl[%d], slp[%d], slp_sqr[%d], muG[%d],' '\n'],...
    sum(o.nc),sum(o.nc),sum(o.nc),sum(o.nc),sum(o.nc),sum(o.nc))];

if o.merit_function == 0
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'mu[%d], dm[%d], mup[%d],' '\n' ...
        o.indent.generic o.indent.generic 'dsl_sqr = 0.0, gps_sqr = 0.0, dm_sqr = 0.0, rho_hat_tmp = 0.0,' '\n'],sum(o.nc),sum(o.nc),sum(o.nc))];
else
    if o.debug > 2
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'dsl_sqr = 0.0, \n'])];
    end
end


if (o.contractive || o.terminal)
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'psi_N[1], psi_Np[1], dot_psi_N[' num2str(o.N*o.nu) '],' '\n'])];
end

optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'rho = 0.0, rho_hat = 0.0,' '\n'])];
optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'phi0 = 0.0, phit = 0.0, phi0_dot = 0.0,' '\n'])];
optCode = [optCode, sprintf([o.indent.generic o.indent.generic 't = 0.0, t_u = 0.0,' '\n'])];
optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'du_sqr = 0.0, u_sqr = 0.0;' '\n' '\n'])];

% possibility to define a variable in the future that reuses
% computationally expensive expressions

%% initialization of the algorithm

[c, d] = generate_forward_simulation(o);
% this function generates det_x
code = [code, c];
data = [data, d];
optCode = [optCode, sprintf(['/* Initialization of the predicted state x (dimension: N*nx) */' '\n'])];
if o.nw >0
    optCode = [optCode, sprintf([o.indent.generic  'det_x(x0,u,w,x);' '\n'])];
else
    optCode = [optCode, sprintf([o.indent.generic  'det_x(x0,u,x);' '\n'])];
end

if isempty(o.objective.Q)&&isempty(o.objective.R)&&isempty(o.objective.P)
    [c, d, in] = general_objective_gradient_oracle(o);
    info.src = [info.src; in.src];
    info.header = [info.header; in.header];
    if ~isempty(in.header)
        libr = [libr, sprintf(['#include "' in.header{:} '"' '\n'])];
    end
else
    [c, d, in] = generate_objective_gradient_oracle(o);
end

% this function generates det_J_and_dot_J and det_J
code = [code, c];
data = [data, d];
info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops.it);
info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops.ls);

c_w_use = argument_w(o,false);
c_tr_use = argument_def(o,false);
c_psi_use = argument_def_internal_psi(o,false);
c_psi_dot_use = argument_def_internal_psi_dot(o,false);
c_contr_use = argument_contr_value(o,false);

if (o.contractive || o.terminal)
    if o.terminal
        optCode = [optCode, sprintf(['\n' o.indent.generic  '/* Initialization of the cost J, its derivative dot_J,' '\n'...
            'the terminal constraint psi_N = x_N^top P x_N and its derivative dot_psi_N */' '\n'])];
    elseif o.contractive
        optCode = [optCode, sprintf(['\n' o.indent.generic  '/* Initialization of the cost J, its derivative dot_J,' '\n'...
            'the terminal constraint psi_N = x_{\tilde k}^top P x_{\tilde k} and its derivative dot_psi_N */' '\n'])];
    end
    optCode = [optCode, sprintf([o.indent.generic  'det_J_and_dot_J(x0, u, x' c_w_use c_tr_use ', &J, dot_J' c_psi_use c_psi_dot_use ');' '\n\n'])];
end

[c, d, in] = generate_slack_initialization(o);

% this function initializes the slack variables
code = [code, c];
data = [data, d];
info.flops.once = falcopt.internal.addFlops(info.flops.once, in.flops);

optCode = [optCode, sprintf(['\n' o.indent.generic  '/* This function initializes the slack variables sl, its squares sl_sqr ' '\n'...
    'and the function gps = g + 0.5 * diag(sl_sqr) */' '\n'])];
optCode = [optCode, sprintf([o.indent.generic  'initialize_slack(u' c_psi_use c_contr_use ', sl, sl_sqr, gps' ');' '\n\n'])];

%% here we loop over the iterations it

optCode = [optCode, sprintf(['\n' o.indent.generic  '/* Outer loop over the iterates it. The break command interrups the cycle \n'...
    o.indent.generic  'if convergence or error is encountered before termination */' '\n'])];

optCode = [optCode, sprintf([o.indent.generic  'for (it=0; it < ' falcopt.internal.toDefineName(o.names.maximumIterations) '; it++) {' '\n'])];

if (o.contractive || o.terminal)
    if o.terminal
        optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* Computation of the cost J, its derivative dot_J,' '\n'...
            'the terminal constraint psi_N = x_N^top P x_N and its derivative dot_psi_N */' '\n'])];
    elseif o.contractive
        optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* Computation of the cost J, its derivative dot_J,' '\n'...
            'the terminal constraint psi_N = x_{\tilde k}^top P x_{\tilde k} and its derivative dot_psi_N */' '\n'])];
    end
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if (it > 0) {' '\n'...
        o.indent.generic o.indent.generic o.indent.generic 'det_J_and_dot_J(x0, u, x' c_w_use c_tr_use ', &J, dot_J' c_psi_use c_psi_dot_use '); }' '\n\n'])];
else
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* Computation of the cost J and its derivative dot_J */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'det_J_and_dot_J(x0, u, x' c_w_use c_tr_use ', &J, dot_J' c_psi_use c_psi_dot_use ');' '\n\n'])];
end



% generate functions: "build_sqr_Nnc" and "build_gpsl"
[c, d, in] = generate_build_gps(o);

code = [code, c];
data = [data, d];
info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);

[c, d, in] = generate_dot_product_Nnu(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.internal.multFlops( in.flops, 2); % called twice
info.flops.it = falcopt.internal.addFlops( info.flops.it, in.flops);

[c, d, in] = generate_gradient_step(o);
code = [code, c];
data = [data, d];
info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops);

optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* Compute the gradient step, projected onto the linearization of the constraints \n' ...
    'du: primal variable gradient steps, dsl: slack variables, muG: dual variables */' '\n'])];
optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'gradient_step(dot_J, u, sl, sl_sqr, gps' c_psi_dot_use ', du, dsl, muG);' '\n\n'])];

optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'dot_product_Nnu(&du_sqr, du, du);' '\n',...
    o.indent.generic o.indent.generic 'dot_product_Nnu(&u_sqr, u, u);' '\n'])];

[c, d, in] = generate_dot_product_Nnc(o);
code = [code, c];
data = [data, d];
if o.merit_function == 0
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'dot_product_Nnc(&gps_sqr, gps, gps);' '\n',...
        o.indent.generic o.indent.generic 'dot_product_Nnc(&dsl_sqr, dsl, dsl);' '\n'])];
    in.flops = falcopt.internal.multFlops( in.flops, 3);
    info.flops.it = falcopt.internal.addFlops( info.flops.it, in.flops);
end

if o.merit_function == 0
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if (it==0)' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'copy_Nnc(mu, muG);' '\n'])];
    info.flops.it.comp = info.flops.it.comp +1;
    
    [c, d, in] = generate_difference(o);
    % generate diff_Nnc
    code = [code, c];
    data = [data, d];
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* dm = muG - mu; */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'diff_Nnc(dm,muG,mu);' '\n\n'])];
    info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops);
    
    % generate function that computes the merit function
    [c, d, in] = generate_det_phi(o);
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);
    info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops); % inside condition
    
    % generate function that computes the derivative of the merit function
    [c, d, in] = generate_det_dot_phi(o);
    code = [code, c];
    data = [data, d];
    info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops);
    info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops); % inside condition
    
    % generate the condition to update the penalty parameter rho
    [c, d, in] = generate_conditions_rho(o);
    code = [code, c];
    data = [data, d];
    info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops);
    
else
    
    % generate function that computes the merit function and its
    % directional derivative (non-smooth penalty function)
    [c,d,in] = generate_Lagrangian_oracles_lp(o);
    % generate det_dot_phi, det_phi, condition_rho
    code = [code, c];
    data = [data, d];
    info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops.it);
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops.ls);
end

%% check conditions on the penalty parameter

if o.merit_function == 0
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* compute the directional derivative of the merit function, phi0_dot */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'det_dot_phi (du,dot_J, rho, gps, mu, dm, &phi0_dot);' '\n'])];
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* Check the penalty parameter rho via phi0_dot */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if (conditions_rho_PM_simpler(phi0_dot,du_sqr,dsl_sqr'])];
    if o.variable_stepSize.active
        optCode = [optCode, sprintf(',alpha')];
    else
        optCode = [optCode, sprintf([',' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax)])];
    end
    optCode = [optCode, sprintf([')==0) {' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic 'dot_product_Nnc(&dm_sqr,dm,dm);' '\n',...
                                 o.indent.generic o.indent.generic o.indent.generic 'rho_hat_tmp = dm_sqr / gps_sqr;' '\n'])];
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic '/* Update the penalty parameter rho */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic 'rho_hat = 2.0 * ' o.auxFunctions.sqrt '(rho_hat_tmp);' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic 'rho = ' o.auxFunctions.max '(2.0*rho,rho_hat);' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic 'flag_ini = 1;' '\n'])];
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic '/* Recompute phi0_dot */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic 'det_dot_phi (du,dot_J, rho, gps, mu, dm, &phi0_dot);' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic '}' '\n'])];
    
    
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if ((flag_ini == 1)||(it == 0)){' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic 'det_phi (J, gps, mu, rho, &phi0);' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic 'flag_ini = 0;' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic '}' '\n\n'])];
    
    info.flops.it.comp = info.flops.it.comp + 4;  % 3+ 1 in fmax
    info.flops.it.mul = info.flops.it.mul +3;
    info.flops.it.sqrt = info.flops.it.sqrt +1;
    info.flops.it.div = info.flops.it.div +1;
    
else
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* Update the penalty parameter rho \n' ...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'and compute merit function phi0 and its derivative phi0_dot */' '\n\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'update_rho(muG, &rho, &rho_hat);' '\n',...
        o.indent.generic o.indent.generic 'det_phi (J, gps, rho, &phi0);' '\n',...
        o.indent.generic o.indent.generic 'det_dot_phi (du,dot_J, rho, gps, &phi0_dot);' '\n\n'])];
    
end

%% start line search
% Add line search parameters to defines
defs = [defs, sprintf(['/** Line Search Parameters **/' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.maximumLineSearchIterations) ' (' falcopt.internal.num2str(o.maxItLs) ')' ' /* Maximum number of line search iterations */' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) ' (' falcopt.internal.num2str(o.variable_stepSize.alpha_max, o.precision) ')' ' /* Alpha max'])];
if ~o.variable_stepSize.active
    defs = [defs, sprintf(' - DO NOT CHANGE, switch to variable step size to make it adjustable')];
end
defs = [defs, sprintf([' */' '\n'])];
if o.variable_stepSize.active
    defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.lineSearchAlphaMin) ' (' falcopt.internal.num2str(o.variable_stepSize.alpha_min, o.precision) ')' ' /* Alpha min */' '\n'])];
end
defs = [defs, sprintf(['/* Derived quantities DO NOT TOUCH */' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.lineSearchAlphaTol) ' (' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) '*' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) ')' ' /* alphaMax*KKTTol */' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.lineSearchAlphaTolSq) ' (' falcopt.internal.toDefineName(o.names.lineSearchAlphaTol) '*' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) ')' ' /* alphaMax*KKTTol^2 */' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.lineSearchAlphaSqTol) ' (' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) '*' falcopt.internal.toDefineName(o.names.lineSearchAlphaTol) ')' ' /* alphaMax^2*KKTTol */' '\n'])];
defs = [defs, sprintf(['#define ' falcopt.internal.toDefineName(o.names.lineSearchAlphaSqTolSq) ' (' falcopt.internal.toDefineName(o.names.lineSearchAlphaSqTol) '*' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) ')' ' /* alphaMax^2*KKTTol^2 */' '\n'])];

% Generate line search
if o.merit_function == 0
    

    if o.variable_stepSize.active
         optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if (phi0_dot <= -alpha*' falcopt.internal.toDefineName(o.names.KKTOptimalityToleranceSq) ') {' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't_u = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'for (it_ls = 0; it_ls < ' falcopt.internal.toDefineName(o.names.maximumLineSearchIterations) '; it_ls++) {' '\n'])];
    else
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if (phi0_dot <= -' falcopt.internal.toDefineName(o.names.lineSearchAlphaTolSq) ') {' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't_u = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'for (it_ls = 0; it_ls < ' falcopt.internal.toDefineName(o.names.maximumLineSearchIterations) '; it_ls++) {' '\n'])];
    end
    
    info.flops.it.comp = info.flops.it.comp +1;
    
    [c,d,in] = generate_weighted_sum_nc(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);
    
    [c,d,in] = generate_weighted_sum_nu(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);
    if o.merit_function == 0
        info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops); % Used twice
    end
    
    % generate the function that computes the next merit function step
    % size t based on a safeguarded quadratic interpolation
    [c,d,in] = generate_quadratic_interpolation(o);
    
    code = [code, c, sprintf('\n\n')];
    data = [data, d];
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);
    
    % update u and slp
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* up = u + t*du; */'])];
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* slp = sl + t*dsl; */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnu(up,du,u,&t);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnc(slp,dsl,sl,&t);' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* mup = mu + t*dm; */' '\n'])];
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnc(mup,dm,mu,&t);' '\n\n'])];
    end
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the predicted state xp (dim. N*nx) using the new predicted input up (dim. N*nu) */' '\n'])];
    if o.nw > 0
        optCode = [optCode, sprintf([...
            o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_x(x0,up,w,xp);' '\n'])];
    else
        optCode = [optCode, sprintf([...
            o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_x(x0,up,xp);' '\n'])];
    end
    
    c_psip_use = argument_def_internal_psi_plus(o,false);
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the new cost Jt, square the slack variables slp ' '\n'...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'and compute gpsp = gps + 0.5 * slp_sqr */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_J(x0, up, xp' c_tr_use ', &Jt' c_psip_use ');' '\n\n'])];
    
    optCode = [optCode, sprintf([...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'build_sqr_Nnc(slp, slp_sqr);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'build_gpsl(up' c_psip_use c_contr_use ',slp_sqr, gpsp);' '\n'])];
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the merit function phit (function of the step size t) */' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_phi (Jt,gpsp,mup,rho,&phit);' '\n'])];
    else
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_phi (Jt,gpsp,rho,&phit);' '\n'])];
    end
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Check the Armijo condition */' '\n\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'if (phit - phi0 <= ' falcopt.internal.num2str(o.parLs, o.precision) '*t*phi0_dot){' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* step size t accepted. Break the line-search */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'break;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'else {' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't_u = t;' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Reduce the step size t via quadratic interpolation */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't = quadratic_interp (phi0, phi0_dot, t_u, phit);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic '/* if t gets too small, output an error */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'if (t_u <= ' falcopt.internal.num2str(o.tolLs, o.precision) ')' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'cond_err = 0;' '\n',...
        o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic 'else {' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic '/* phi0_dot is sufficiently small and we have converged */' '\n' ...
        o.indent.generic o.indent.generic o.indent.generic 'conditions_f = 0;' '\n' ...
        o.indent.generic o.indent.generic '}' '\n\n'])];
    info.flops.ls.comp = info.flops.ls.comp +1;
    info.flops.ls.mul = info.flops.ls.mul +2;
    info.flops.it.comp = info.flops.it.comp +1;
    
else
    %% recompute rho
    

    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic '/* Check conditions on phi0_dot for convergence and start line-search */' '\n'])];
    if o.variable_stepSize.active
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if (phi0_dot <= -alpha*' falcopt.internal.toDefineName(o.names.KKTOptimalityToleranceSq) ') {' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't_u = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'for (it_ls = 0; it_ls < ' falcopt.internal.toDefineName(o.names.maximumLineSearchIterations) '; it_ls++) {' '\n'])];
    else
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'if (phi0_dot <= -' falcopt.internal.toDefineName(o.names.lineSearchAlphaTolSq) ') {' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 't_u = 1.0;' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'for (it_ls = 0; it_ls < ' falcopt.internal.toDefineName(o.names.maximumLineSearchIterations) '; it_ls++) {' '\n'])];
    end
    
    info.flops.it.comp = info.flops.it.comp +1;
    
    [c,d,in] = generate_weighted_sum_nc(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);
    
    [c,d,in] = generate_weighted_sum_nu(o);
    % generate weighted_sum_Nnu() and weighted_sum_Nnc()
    code = [code, c];
    data = [data, d];
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);
    if o.merit_function == 0
        info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops); % Used twice
    end
    
    [c,d,in] = generate_quadratic_interpolation(o);
    
    code = [code, c, sprintf('\n\n')];
    data = [data, d];
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls, in.flops);
    
    % update u and slp
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* up = u + t*du; */'])];
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* slp = sl + t*dsl; */' '\n\n'])];
    
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnu(up,du,u,&t);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnc(slp,dsl,sl,&t);' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* mup = mu + t*dm; */' '\n\n'])];
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnc(mup,dm,mu,&t);' '\n\n'])];
    end
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the predicted state xp (dim. N*nx)' '\n'...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'using the new predicted input up (dim. N*nu) */' '\n'])];
    if o.nw > 0
        optCode = [optCode, sprintf([...
            o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_x(x0,up,w,xp);' '\n'])];
    else
        optCode = [optCode, sprintf([...
            o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_x(x0,up,xp);' '\n'])];
    end
    
    c_psip_use = argument_def_internal_psi_plus(o,false);
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the new cost Jt, square the slack variables slp' '\n'...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'and compute gpsp = gps + 0.5 * slp_sqr */' '\n\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_J(x0, up, xp' c_tr_use ', &Jt' c_psip_use ');' '\n\n'])];
    
    optCode = [optCode, sprintf([...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'build_sqr_Nnc(slp, slp_sqr);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'build_gpsl(up' c_psip_use c_contr_use ',slp_sqr, gpsp);' '\n'])];
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the merit function phit (function of the step size t) */' '\n\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_phi (Jt,gpsp,mup,rho,&phit);' '\n'])];
    else
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_phi (Jt,gpsp,rho,&phit);' '\n'])];
    end
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Check the Armijo condition */' '\n\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'if (phit - phi0 <= ' falcopt.internal.num2str(o.parLs, o.precision) '*t*phi0_dot){' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* step size t accepted. Break the line-search */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'break;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'else {' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Reduce the step size t via quadratic interpolation */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't_u = t;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't = quadratic_interp (phi0, phi0_dot, t_u, phit);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* if t gets too small, try resetting rho */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'if (t_u <= ' falcopt.internal.num2str(o.tolLs, o.precision) ')' '\n'...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'reset_rho = 1;' '\n'...
        o.indent.generic o.indent.generic o.indent.generic '}' '\n'])];
    
    optCode = [optCode, sprintf([...
        o.indent.generic o.indent.generic o.indent.generic '/* Reset rho to rho_hat if rho is overly big */' '\n'...
        o.indent.generic o.indent.generic o.indent.generic 'if (reset_rho == 1){' '\n'...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'reset_rho = 0;' '\n'...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'rho = rho_hat;' '\n'...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Reset phi0 and phi0_dot */' '\n'...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_phi (J, gps, rho, &phi0);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_dot_phi (du,dot_J, rho, gps, &phi0_dot);' '\n\n',...
        ...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't = 1.0;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't_u = 1.0;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'for (it_ls = 0; it_ls < ' falcopt.internal.toDefineName(o.names.maximumLineSearchIterations) '; it_ls++) {' '\n'])];
    info.flops.it.comp = info.flops.it.comp +1;
    
    % update u and slp
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* up = u + t*du; */'])];
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* slp = sl + t*dsl; */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnu(up,du,u,&t);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnc(slp,dsl,sl,&t);' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'weighted_sum_Nnc(mup,dm,mu,&t);' '\n\n'])];
    end
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the predicted state xp (dim. N*nx)' '\n'...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'using the new predicted input up (dim. N*nu) */' '\n\n'])];
    if o.nw > 0
        optCode = [optCode, sprintf([...
            o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_x(x0,up,w,xp);' '\n'])];
    else
        optCode = [optCode, sprintf([...
            o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_x(x0,up,xp);' '\n'])];
    end
    
    c_psip_use = argument_def_internal_psi_plus(o,false);
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the new cost Jt, square the slack variables slp and compute gpsp = gps + 0.5 * slp_sqr */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_J(x0, up, xp' c_tr_use ', &Jt' c_psip_use ');' '\n\n'])];
    
    optCode = [optCode, sprintf([...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'build_sqr_Nnc(slp, slp_sqr);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'build_gpsl(up' c_psip_use c_contr_use ',slp_sqr, gpsp);' '\n'])];
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Compute the merit function phit (function of the step size t) */' '\n'])];
    if o.merit_function == 0
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_phi (Jt,gpsp,mup,rho,&phit);' '\n'])];
    else
        optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'det_phi (Jt,gpsp,rho,&phit);' '\n'])];
    end
    
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Check the Armijo condition */' '\n\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'if (phit - phi0 <= ' falcopt.internal.num2str(o.parLs, o.precision) '*t*phi0_dot){' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* step size t accepted. Break the line-search */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'break;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'else {' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* Reduce the step size t via quadratic interpolation */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't_u = t;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 't = quadratic_interp (phi0, phi0_dot, t_u, phit);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        '\n' o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic '/* if t gets too small, output an error */' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'if (t_u <= ' falcopt.internal.num2str(o.tolLs, o.precision) ')' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'cond_err = 0;' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic 'else' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'conditions_f = 0;' '\n\n'])];
    info.flops.ls.comp = info.flops.ls.comp +1;
    info.flops.ls.mul = info.flops.ls.mul +2;
    info.flops.it.comp = info.flops.it.comp +1;
    
end

% update variable stepSize
if o.variable_stepSize.active
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic '/* update step size via trust-region procedure */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'ared = (phi0 - phit); \n'...
                o.indent.generic o.indent.generic 'pred = (-t*phi0_dot);\n'...
                o.indent.generic o.indent.generic 'rat = ared/pred ;\n'...
                o.indent.generic o.indent.generic 'if( rat < rat_1)\n'...
                o.indent.generic o.indent.generic o.indent.generic 'alpha = ' o.auxFunctions.max '( alpha_min, alpha*gamma_1);\n'...
                o.indent.generic o.indent.generic 'if( rat > rat_2)\n'...
                o.indent.generic o.indent.generic o.indent.generic 'alpha = ' o.auxFunctions.min '( alpha_max, alpha*gamma_2);\n'...
                o.indent.generic o.indent.generic 'alpha_old = alpha; \n'...
                o.indent.generic o.indent.generic 'alpha_inverse = 1.0/alpha;\n\n'])];
end

%% the code continues

optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'iter_ls[it] = it_ls+1;','\n',...
    o.indent.generic o.indent.generic 'if (it_ls == ' falcopt.internal.toDefineName(o.names.maximumLineSearchIterations) ')' '\n',...
    o.indent.generic o.indent.generic o.indent.generic 'cond_err = 0;' '\n\n'])];
info.flops.once.add = info.flops.once.add +1;
%         optCode = [optCode, sprintf(['if (it==1)

%      optCode = [optCode, sprintf(['if (cond_err == 0){' '\n' ...
%          o.indent.generic  'u[0] = rho;\n' ...
%          o.indent.generic  'u[1] = phi0_dot; \n' ...
%          o.indent.generic  'return cond;}' '\n'])];

%% update u, slack and dual variables

optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* update step */' '\n'])];
if o.merit_function == 0
    optCode = [optCode, sprintf([ o.indent.generic o.indent.generic 'copy_Nnc(mu,mup);' '\n'])];
end

optCode = [optCode, sprintf([o.indent.generic o.indent.generic 'copy_Nnu(u,up);' '\n',...
    o.indent.generic o.indent.generic 'copy_Nnc(sl,slp);' '\n',...
    o.indent.generic o.indent.generic 'copy_Nnc(sl_sqr,slp_sqr);' '\n',...
    o.indent.generic o.indent.generic 'copy_Nnc(gps,gpsp);' '\n',...
    o.indent.generic o.indent.generic 'copy_Nnx(x,xp);' '\n',...
    o.indent.generic o.indent.generic 'J = Jt;' '\n',...
    o.indent.generic o.indent.generic 'phi0 = phit;' '\n\n'])];

%% check convergence

[c, d, in] = generate_compute_max(o);
info.flops.it = falcopt.internal.addFlops(info.flops.it, in.flops);
% compute_max_Nnc

code = [code, c, sprintf('\n\n')];
data = [data, d];


optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* check convergence (update step) */' '\n'])];
optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic 'constr_viol = compute_max_Nnc(gpsp); ' '\n'])];
if o.variable_stepSize.active 
    optCode = [optCode,sprintf([o.indent.generic o.indent.generic 'if ((du_sqr >= alpha_old*alpha_old*' falcopt.internal.toDefineName(o.names.KKTOptimalityToleranceSq) ')||(constr_viol>= ' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) '))' '\n'])];
else
    optCode = [optCode,sprintf([o.indent.generic o.indent.generic 'if ((du_sqr >= ' falcopt.internal.toDefineName(o.names.lineSearchAlphaSqTolSq) ')||(constr_viol >= ' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) '))' '\n'])];
end

optCode = [optCode,sprintf([o.indent.generic o.indent.generic o.indent.generic 'conditions_x = 1;' '\n',...
    o.indent.generic o.indent.generic 'else' '\n',...
    o.indent.generic o.indent.generic 'conditions_x = 0;' '\n'])];



optCode = [optCode,sprintf([o.indent.generic o.indent.generic 'if (it == ' falcopt.internal.toDefineName(o.names.maximumIterations) '-1)' '\n'])];
optCode = [optCode,sprintf([o.indent.generic o.indent.generic o.indent.generic 'conditions_n = 0;' '\n'])];


optCode = [optCode,sprintf([o.indent.generic o.indent.generic 'if(!((conditions_f && cond_err) && (conditions_n && conditions_x)))' '\n'])];
optCode = [optCode,sprintf([o.indent.generic o.indent.generic o.indent.generic 'break;' '\n\n',...
    o.indent.generic  '}' '\n\n'])];
info.flops.it.comp = info.flops.it.comp +3;

if o.debug >2
    optCode = [optCode, sprintf(['\n' o.indent.generic  '/* Assign optimality */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic 'dot_product_Nnc(&dsl_sqr, dsl, dsl);' '\n'])];
    if o.variable_stepSize.active 
        optCode = [optCode, sprintf([o.indent.generic '(*optimval) = (' o.auxFunctions.sqrt '(du_sqr + dsl_sqr))*alpha_inverse; ' '\n'])];
    else 
        optCode = [optCode, sprintf([o.indent.generic '(*optimval) = (' o.auxFunctions.sqrt '(du_sqr + dsl_sqr))/' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) '; \n'])];
    end
    optCode = [optCode, sprintf([o.indent.generic '(*feasval) = constr_viol; ' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic '(*meritval) = phi0_dot; ' '\n'])];
end
    
optCode = [optCode, sprintf(['\n' o.indent.generic  '/* Assign exitflag */' '\n'])];
optCode = [optCode,sprintf([o.indent.generic  '(*iter) = it+1;' '\n',...
    o.indent.generic  'if (conditions_f == 0)' '\n',...
    o.indent.generic o.indent.generic 'cond = 2;' '\n',...
    o.indent.generic  'else {' '\n',...
    o.indent.generic o.indent.generic 'if (conditions_x == 0)' '\n',...
    o.indent.generic o.indent.generic o.indent.generic 'cond = 1;' '\n',...
    o.indent.generic o.indent.generic 'else {' '\n',...
    o.indent.generic o.indent.generic o.indent.generic 'if (conditions_n == 0)' '\n',...
    o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'cond = 0;' '\n',...
    o.indent.generic o.indent.generic o.indent.generic 'else {' '\n'])];
if o.variable_stepSize.active
    optCode = [optCode,sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'if (phi0_dot < -alpha_old* ' falcopt.internal.toDefineName(o.names.KKTOptimalityTolerance) ') \n'])];
else
    optCode = [optCode,sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'if (phi0_dot < -' falcopt.internal.toDefineName(o.names.lineSearchAlphaTol) ') \n'])];
end    
optCode = [optCode,sprintf([o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'cond = -1;' '\n'...
    o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'else' '\n'...
    o.indent.generic o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'cond = -10;' '\n'...
    o.indent.generic o.indent.generic o.indent.generic '}' '\n'...
    o.indent.generic o.indent.generic '}' '\n'...
    o.indent.generic  '}' '\n\n'])];

if o.debug > 1
    optCode = [optCode, sprintf(['\n' o.indent.generic o.indent.generic '/* Output cost function */' '\n'])];
    optCode = [optCode, sprintf([o.indent.generic  '*fval = J;' '\n'])];
end
optCode = [optCode, sprintf([o.indent.generic  'return cond;' '\n'])];
optCode = [optCode, sprintf('} \n\n')];
info.flops.once.add = info.flops.once.add +1;
info.flops.once.comp = info.flops.once.comp +3;

%% Generate build types
% generate MEX
if any(strcmp('mex', o.buildTypes))
    mexName = o.name;
    mex_code = falcopt.generateMEX(o.N, struct('x',o.nx,'u',o.nu,'w',o.nw), ...
        'names', struct('fun','proposed_algorithm','mex',mexName, 'maximumIterations', o.names.maximumIterations), ...
        'ref', o.objective.trackReference,...
        'terminalContraction', o.terminal || o.contractive, 'indexContraction',...
        o.contractive, 'type', o.real,'maxIt',o.maxIt, 'timing', true, 'debug', o.debug,...
        'verbose', o.verbose, 'indent', o.indent, 'inline', o.inline, 'precision', o.precision);
	mex_code = [copyright  '\n' libr '\n' defs '\n' data '\n'  code '\n' optCode '\n' mex_code];
    if length(o.buildTypes) > 1
        name = [o.name '_mex'];
    else
        name = o.name;
    end
    filename = [o.gendir '/' name '.c'];
    ext_file = [];
    for k = 1:length(info.src)
        ext_file = [ext_file, ' ', info.src{k}]; %#ok
    end
    if strcmp(o.gradients,'ccode')
        if any(cellfun('length',regexp(info.src,'external_functions.c')))&&~any(cellfun('length',regexp(info.header,'external_functions.h')))
            % avoid to compile twice external_functions.c
            ext_file = [];
        end
    end
    f = fopen(filename, 'w+');
    fprintf(f, mex_code);
    fclose(f);

    if o.compile
        compile = ['mex ' filename ext_file];
        if o.verbose >= 2
            compile = [compile ' -v'];
        end
        compile = [compile ' -output ' mexName];
        if o.verbose >= 1
            disp(compile);
        end
        eval(compile);
    end
end
% Generate production code
if any(strcmp('production', o.buildTypes))
    production_code = [copyright  '\n' libr '\n' defs '\n' data '\n'  code '\n' optCode];
    filename = [o.gendir '/' o.name];
    ext_file = [];
    for k = 1:length(info.src)
        ext_file = [ext_file, ' ', info.src{k}]; %#ok
    end
    if strcmp(o.gradients,'ccode')
        if any(cellfun('length',regexp(info.src,'external_functions.c')))&&~any(cellfun('length',regexp(info.header,'external_functions.h')))
            % avoid to compile twice external_functions.c
            ext_file = [];
        end
    end
    f = fopen([filename '.c'], 'w+');
    fprintf(f, production_code);
    fclose(f);
end
% generate simulink
if any(strcmp('simulink', o.buildTypes))
    simulink_fcn = falcopt.generateSFunction(o.nx,o.nu,o.nw,o.N, 'maxIt',o.maxIt,...
                            'names', struct('fun','proposed_algorithm', 'maximumIterations', o.names.maximumIterations), ...
                            'name',sprintf(['simulink_' o.name]), 'trackReference',o.objective.trackReference,...
                            'terminal',o.terminal, 'contractive',o.contractive, ...
                            'indent',o.indent.generic, 'precision',o.precision, ...
                            'type', o.real, 'debug', o.debug, ...
                            'maskImage_name','+falcopt\simulinkMask_logo.png');
    simulink_code = [copyright '\n' libr '\n' defs '\n' data '\n'  code '\n' optCode '\n' simulink_fcn];
    name = sprintf([o.gendir '/simulink_' o.name '.c']);
    f = fopen(name, 'w+');
    fprintf(f, simulink_code);
    fclose(f);
    ext_file = [];
    for k = 1:length(info.src)
        ext_file = [ext_file, ' ', info.src{k}]; %#ok
    end
    % compile
    eval(sprintf(['mex ' name ext_file]));
end
end

function [ data, code, info] = casadi_jacobians(varargin)
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
p = inputParser;
p.addRequired('o');
p.addRequired('f',@(x)isa(x,'function_handle'));
p.addRequired('fname',@ischar);
p.addParameter('staticName','',@ischar);
p.addParameter('jac',{},@iscell);
p.addParameter('jac_x',{'F','in_F','Jacobian_x'},@iscell);
p.addParameter('jac_u',{'G','in_G','Jacobian_u'},@iscell);
p.addParameter('fileName','casadi_fcn',@ischar);
p.addParameter('generate_code',true,@islogical); 
p.addParameter('structure','unique',@ischar); 
p.parse(varargin{:});                                   
                                                                        
o = p.Results.o;                                                                
f = p.Results.f;                                                                 
fName = p.Results.fname; 
x = casadi.SX.sym('x',o.nx,1);
u = casadi.SX.sym('u',o.nu,1);
w = casadi.SX.sym('w',o.nw,1);

data = [];
code = [];
info.sxfcn = {};


% function
if( nargin(f) == 3)
    y = f(x,u,w);
elseif( nargin(f) == 2)
    y = f(x,u);
else
    if( strcmp(fName,'model_mpc'))
        error('Check number of dynamics inputs, can be (x,u) or (x,u,w)');
    elseif strcmp(fName,'cost_N')
        y = f(x);
    else
        error(['Check number of ' fName ' inputs, can be (x,u) or (x,u,w)']);
    end
end

if( isnumeric(y))&&(~isempty(p.Results.staticName))   % generate static data if is not symbolic
    name.M = staticName;
    [d, ~, in_y] = falcopt.generateData(y, 'names', name, ...
        'type', o.real, 'precision', o.precision, 'structure', 'unique', 'noones', false, 'indent', o.indent, ...
        'static', true, 'const', true, 'verbose', o.verbose);
    data = [data, d, sprintf('\n')];
    in_y.static = 1;
    info.y = in_y;
    info.y.flops = 0;
else
    [~,in_y] = falcopt.casadi2struct( y, 'structure', p.Results.structure,'errorName',fName);
    y_f = casadi.Function([fName '_casadi'],{x,u,w},{in_y.stored.values});
    info.sxfcn{1} = y_f;
   
    % Get number of FLOPS
    try
        info.y.flops = y_f.getAlgorithmSize();
    catch
        try
            info.y.flops = y_f.n_instructions();
        catch
            warning('Cannot use casadi.Function.getAlgorithmSize() or casadi.Function.n_instructions(), number of FLOPS may too low as a result.');
            info.y.flops = 0;
        end
    end
    

    info.y.static = 0;
    
    % wrapper function
    if strcmp(fName,'model_mpc') 
        code = [code, sprintf('/* Dynamics of the system */\n')];
    end
    if( nargin(f) == 3)
        code = [code, sprintf(['void ' fName '( const ' o.real '* x, const ' o.real '* u, const ' o.real '* v, ' o.real '* xp){' '\n\n'])];
        code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],3)];
        code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],y_f.getWorkSize(),0,0)]; %TODO Tommaso check parameter
        code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\t' 'in[2] = v;\n\n'])];
    elseif( nargin(f) == 2)
        code = [code, sprintf(['void ' fName '( const ' o.real '* x, const ' o.real '* u, ' o.real '* xp){' '\n\n'])];
        code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],2)];
        code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],y_f.getWorkSize(),0,0)]; %TODO Tommaso check parameter
        code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\n'])];
    elseif( nargin(f) == 1)
        code = [code, sprintf(['void ' fName '( const ' o.real '* x, const ' o.real '* xp){' '\n\n'])];
        code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],1)];
        code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],y_f.getWorkSize(),0,0)]; %TODO Tommaso check parameter
        code = [code, sprintf('\tin[0] = x;\n\n')];
    end
    
    code = [code, sprintf([o.indent.generic  fName '_casadi( in, &xp, &iw, w, mem);' o.indent.generic '/* external generated casadi function*/\n}\n\n'])];
    if( isnumeric(y))&&(~isempty(p.Results.staticName))
        data = [data, sprintf([o.indent.generic  'static ' o.real ' ' p.Results.staticName '[%i];' '\n'], in_y.stored.num)];
    end
    info.y.structure = in_y;
end

% jacobians
for k = 1:length(p.Results.jac)
    switch p.Results.jac{k}
        case 'x'
            jac = jacobian( y, x);
            static_name = p.Results.jac_x{1};
            struct_name = p.Results.jac_x{2};
            J_name = p.Results.jac_x{3};
        case 'u'
            jac = jacobian( y, u);
            static_name = p.Results.jac_u{1};
            struct_name = p.Results.jac_u{2};
            J_name = p.Results.jac_u{3};
        otherwise
            error('invalid variable name');
    end
    [info.(struct_name).static,in] = falcopt.casadi2struct(jac,'structure', p.Results.structure,'errorName',fName);
    
    if( info.(struct_name).static) %is jac matrix constant ?
        name.M = static_name;
        data = [data, sprintf([o.indent.data '/* Static data for %s */\n'],J_name)]; %#ok
        if ~isempty(find((full(DM(jac)) ~= 0), 1))
            [d, ~, in_d] = falcopt.generateData(full(DM(jac)), 'names', name, ...
                'type', o.real, 'precision', o.precision, 'structure', 'unique', 'noones', false, 'indent', o.indent, ...
                'static', true, 'const', true, 'verbose', o.verbose);
            info.(struct_name).struct = in_d;
            info.(struct_name).flops = 0;
        else
            d = sprintf([o.indent.data 'static ' o.real ' ' static_name '[%i]={ '],length(full(DM(jac))));
            for j = 1:length(full(DM(jac)))-1
                d = [d, '0.0, ']; %#ok
            end
            d = [d, sprintf(['0.0 };\n'])]; %#ok
            info.(struct_name).flops = 0;
        end
        if ~isempty(d)
            data = [data, d, sprintf('\n')]; %#ok
        end
        info.(struct_name).struct.structure = in;
    else
        jac = casadi.Function(sprintf([J_name '_casadi']),{x,u,w},{in.stored.values});
        info.sxfcn{length(info.sxfcn)+1} =  jac; 
        
        % wrapper jacobian
        if strcmp(fName,'model_mpc')
            code = [code, sprintf('/* System dynmics jacobian w.r.t %c variable */\n',p.Results.jac{k})]; %#ok
        end
        if( nargin(f) == 3)
            code = [code, sprintf(['void ' J_name '( const ' o.real '* x, const ' o.real '* u, const ' o.real '* v, ' o.real '* res){' '\n\n'])]; %#ok
            code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],3)]; %#ok
            code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],jac.getWorkSize(),0,0)]; %#ok
            code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\t' 'in[2] = v;\n\n' ...
                o.indent.generic  J_name '_casadi( in, &res, &iw, w, mem);' o.indent.generic '/* external generated casadi function*/\n}\n\n']) ]; %#ok
        elseif( nargin(f) == 2)
            code = [code, sprintf(['void ' J_name '( const ' o.real '* x, const ' o.real '* u, ' o.real '* res){' '\n\n'])]; %#ok
            code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],2)]; %#ok
            code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],jac.getWorkSize(),0,0)]; %#ok
            code = [code, sprintf(['\tin[0] = x;\n\t' 'in[1] = u;\n\n' ...
                o.indent.generic  J_name '_casadi( in, &res, &iw, w, mem);' o.indent.generic '/* external generated casadi function*/\n}\n\n']) ]; %#ok
        elseif( nargin(f) == 1)
            code = [code, sprintf(['void ' J_name '( const ' o.real '* x, const ' o.real '* res){' '\n\n'])]; %#ok
            code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],1)]; %#ok
            code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],jac.getWorkSize(),0,0)]; %#ok
            code = [code, sprintf(['\tin[0] = x;\n\n' ...
                o.indent.generic  J_name '_casadi( in, &res, &iw, w, mem);' o.indent.generic '/* external generated casadi function*/\n}\n\n']) ]; %#ok
        end
        
        data = [data, sprintf([o.indent.data '/* Static data for %s */\n'],J_name)]; %#ok
        data = [data, sprintf([o.indent.data  'static ' o.real ' ' static_name '[%i];' '\n'], in.stored.num)]; %#ok
        info.(struct_name).struct.structure = in;
        
        % Get number of FLOPS
        try
            info.(struct_name).flops = jac.getAlgorithmSize();
        catch
            try
                info.(struct_name).flops = jac.n_instructions();
            catch
                warning('Cannot use casadi.Function.getAlgorithmSize() or casadi.Function.n_instructions(), number of FLOPS may too low as a result.');
                info.(struct_name).flops = 0;
            end
        end
    end
end

% generate .c file with casadi functions
if p.Results.generate_code
    [info.src,info.header] = generate_casadi_c(o, p.Results.fileName, info.sxfcn);
end

end

function [src,header] = generate_casadi_c(o, fileName, fcns)
    % given functions 'fcns' generate c code using casadi in the give 'fileName'
    % returns name of .c src file and .h header file
    C = casadi.CodeGenerator(fileName, struct('with_header',true,'real_t',o.real));

    for k = 1:length(fcns)
        if( ~isempty(fcns{k}))
            C.add(fcns{k})
        end
    end
    
    % TODO: do this without cd!!!
    currentPath = pwd;
    cd(sprintf(o.gendir));
    C.generate();
    cd(currentPath);
    src = {[o.gendir '/' fileName '.c']};
    header = {[fileName '.h']};
end

function [ data, code, info] = matlab_jacobians(varargin)
% generate c code from provided function f and automatically generate jacobians
% using matlab symbolic toolbox
p = inputParser;
p.addRequired('o');
p.addRequired('f',@(x)isa(x,'function_handle'));
p.addRequired('fname',@ischar);
p.addParameter('staticName','',@ischar);
p.addParameter('jac',{},@iscell);
p.addParameter('jac_x',{'F','in_F','Jacobian_x'},@iscell);
p.addParameter('jac_u',{'G','in_G','Jacobian_u'},@iscell);
p.addParameter('fileName','casadi_fcn',@ischar);
p.addParameter('generate_code',true,@islogical);
p.addParameter('structure','unique',@ischar);
p.parse(varargin{:});

o = p.Results.o;
f = p.Results.f;
fName = p.Results.fname;

code = [];
data = [];
info.y = [];
info.in_F = [];
info.in_G = [];
x = sym('x',[o.nx,1],'real');
u = sym('u',[o.nu,1],'real');
w = sym('w',[o.nw,1],'real');


if( nargin(f) == 3)
    y = f(x,u,w);
elseif( nargin(f) == 2)
    y = f(x,u);
else
    if( strcmp(fName,'model_mpc'))
        error('Check number of dynamics inputs, can be (x,u) or (x,u,w)');
    elseif strcmp( fName,'cost_N')  
        y = f(x);
    else
        error(['Check number of ' fName ' inputs, can be (x,u) or (x,u,w)']);
    end
end

% function
[d, in_y] = falcopt.fcn2struct( y,o,'name','xp','structure',p.Results.structure);
info.y.static = in_y.static;
if info.y.static
   error(['function ' fName ' does not depend on any variables']);
end

if strcmp(fName,'model_mpc')
    code = [ code, sprintf('/*Dynamics of the system*/\n')];
end
if( nargin(f) == 3)
    code = [code, sprintf(['void ' fName '(const ' o.real '* x, const ' o.real '* u, const ' o.real '* w, ' o.real '* xp){\n\n'])];
elseif( nargin(f) == 2)
    code = [code, sprintf(['void ' fName '(const ' o.real '* x, const ' o.real '* u, ' o.real '* xp){\n\n'])];
elseif (nargin(f) == 1)
    code = [code, sprintf(['void ' fName '(const ' o.real '* x, ' o.real '* xp){\n\n'])];
end
code = [ code, d, sprintf('}\n\n')];
info.y.flops = in_y.flops;


% jacobians
for k = 1:length(p.Results.jac)
    switch p.Results.jac{k}
        case 'x'
            jac = jacobian( y, x);
            static_name = p.Results.jac_x{1};
            struct_name = p.Results.jac_x{2};
            J_name = p.Results.jac_x{3};
        case 'u'
            jac = jacobian( y, u);
            static_name = p.Results.jac_u{1};
            struct_name = p.Results.jac_u{2};
            J_name = p.Results.jac_u{3};
        otherwise
            error('invalid variable name');
    end
    [d, i] = falcopt.fcn2struct( jac,o,'name', static_name,'structure',p.Results.structure); 
    info.(struct_name).static = i.static;
    data = [ data, sprintf([o.indent.data '/*Static data for %s*/\n'],J_name)]; %#ok
    if( info.(struct_name).static)
        data = [data, sprintf([o.indent.data 'static const ' o.real ' ' static_name '']),d,sprintf(';\n')]; %#ok
        %data = [data, d, sprintf([o.indent.data '};\n\n'])]; %#ok
    else
        data = [data, sprintf([o.indent.data 'static ' o.real ' ' static_name '[%i];\n'], i.structure.num)]; %#ok
        code = [ code, sprintf([o.indent.data '/*%s*/\n'],J_name)]; %#ok
        if( nargin(f) == 3)
            code = [code, sprintf([o.indent.code 'static ' o.inline ' void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, const ' o.real ' * w, ' o.real ' * ' static_name ') {\n\n'])]; %#ok
        elseif( nargin(f) == 2)
            code = [code, sprintf([o.indent.code 'static ' o.inline ' void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, ' o.real ' * ' static_name ') {\n\n'])]; %#ok
        elseif( nargin(f) == 1)
            code = [code, sprintf([o.indent.code 'static ' o.inline ' void ' J_name '(const ' o.real ' * x, ' o.real ' * ' static_name ') {\n\n'])]; %#ok
        end
        code = [code, d, sprintf([o.indent.code '}\n\n'])]; %#ok
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
        try
            y = o.dynamics(x,u,w);
            jac_u = o.external.manual.jacobian_u(x,u,w);
            jac_x = o.external.manual.jacobian_x(x,u,w);
        catch
            error(['Error in converting output of ''external_jacobian'' functions into sym variables\n'...
                'Try initialize output of these functions as sym. Example: sym(zeros(2,2)))']);
        end
    else
        try
            y = o.dynamics(x,u);
            jac_u = o.external.manual.jacobian_u(x,u);
            jac_x = o.external.manual.jacobian_x(x,u);
        catch
            error(['Error in converting output of ''external_jacobian'' functions into sym variables\n'...
                'Try initialize output of these functions as sym. Example: sym(zeros(2,2)))']);
        end
    end

    % dynamics of system
    [d, in_y] = falcopt.fcn2struct( y,o,'name','xp');
    code = [ code, sprintf('/*Dynamics of the system*/\n')];
    if( o.nw > 0)
        code = [code, sprintf(['void model_mpc(const ' o.real '* x, const ' o.real '* u, const ' o.real '* w, ' o.real '* xp){\n\n'])];
    else
        code = [code, sprintf(['void model_mpc(const ' o.real '* x, const ' o.real '* u, ' o.real '* xp){\n\n'])];
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
        data = [ data, sprintf([o.indent.generic '/*Static data for ' J_name '*/\n'])]; %#ok
        if( info.(struct_name).static)
            data = [data, sprintf(['\tstatic const ' o.real ' ' static_name '[%i] = {\n'], i.structure.num)]; %#ok
            data = [data, d, sprintf(['};\n\n'])]; %#ok
        else
            data = [data, sprintf(['\tstatic ' o.real ' ' static_name '[%i];\n'], i.structure.num)]; %#ok
            code = [ code, sprintf(['/* ' J_name '*/\n'])]; %#ok
            if( o.nw > 0)
                code = [code, sprintf(['static ' o.inline ' void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, const ' o.real ' * w, ' o.real ' * ' static_name ') {\n\n'])]; %#ok
            else
                code = [code, sprintf(['static ' o.inline ' void ' J_name '(const ' o.real ' * x, const ' o.real ' * u, ' o.real ' * ' static_name ') {\n\n'])]; %#ok
            end
            code = [code, d, sprintf(['}\n\n'])]; %#ok
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

% check if there are box constraints
% if( isempty(o.box_constraints))
%     ibox = 0;
% else
%     ibox = 1;
% end

if ~isempty(o.K_n)&& ~isequal(grad,'ccode')
    nl_con = 1;
else
    nl_con = 0;
end
%nl_con = ~isempty(o.K_n); % check if there are nonlinear constraints n(u)

% if( isa(o.constraints,'function_handle'))
%     nl_con = 1;
% else
%     nl_con = 0;
% end

% build nonlinear constraints n(u)
switch grad
    case 'casadi'
        z = casadi.SX.sym('u',o.nu,1);
    case {'matlab','manual'}
        z = sym('u',[o.nu,1],'real');
end
if( nl_con)
    for jj=1:length(o.K_n)
        if( nargin(o.constraints{jj}) == 1)
            try
                n{jj} = o.constraints{jj}(z); %#ok
            catch
                try
                    n{jj} = o.constraints{jj}(z'); %#ok
                catch
                    error('inequality constraints can depend only on inputs' );
                end
                error('inequality constraints must depend on a column vector u (not a row vector)' );
            end
        else
            error('inequality constraints can depend only on input u' );
        end
    end
end



%% generate a - u functions
for ii=1:length(o.K_amu)
    cons_Bound = o.lb(:,o.K_amu{ii}(1));
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
code = [code, sprintf(['\n' o.indent.code o.inline 'void build_amu(const ' o.real '* u, const unsigned int k, ' o.real '* amu){' '\n\n'])];
for ii=1:length(o.K_amu)
    code = [code, sprintf([o.indent.code o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_amu{ii}-1, 'k') ') {' '\n'])]; %#ok
    code = [code, sprintf([o.indent.code o.indent.generic o.indent.generic 'build_amu_' num2str(ii) '(&amu[0], &u[0]);' '\n'])]; %#ok
    code = [code, sprintf([o.indent.code o.indent.generic '}' '\n'])]; %#ok
end
code = [code, sprintf([o.indent.code '}' '\n'])];

%% generate u - b functions
for ii=1:length(o.K_umb)
    cons_Bound = o.ub(:,o.K_umb{ii}(1));
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
code = [code, sprintf(['\n' o.indent.code o.inline 'void build_umb(const ' o.real '* u, const unsigned int k, ' o.real '* umb){' '\n\n'])];
for ii=1:length(o.K_umb)
    code = [code, sprintf([o.indent.code o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_umb{ii}-1, 'k') ') {' '\n'])]; %#ok
    code = [code, sprintf([o.indent.code o.indent.generic o.indent.generic 'build_umb_' num2str(ii) '(&umb[0], &u[0]);' '\n'])]; %#ok
    code = [code, sprintf([o.indent.code o.indent.generic  '}' '\n'])]; %#ok
end
code = [code, sprintf([o.indent.code '}' '\n'])];

% build o.nc (constraints structure)
info.nc = o.nn;


if( o.terminal || o.contractive)
    info.nc = [info.nc 1];
else
    info.nc = [info.nc 0];
end

for i=1:o.N
    lb = sum(~isinf(o.lb(:,i)));
    ub = sum(~isinf(o.ub(:,i)));
    info.nc(i) = info.nc(i)+lb+ub;
end

% generate build_n and build_Dn
if( nl_con)
    switch grad
        case 'casadi' % use casadi generation
            sxfcn = {};
            
            in_n = cell(length(o.K_n));
            y_f = cell(length(o.K_n));
            Dn = cell(length(o.K_n));
            Dn_f = cell(length(o.K_n));
            
            if( nl_con)  % if exist n(u)
                for jj=1:length(o.K_n)
                    
                    % build_n
                    [~,in_n{jj}] = falcopt.casadi2struct( n{jj});
                    
                    y_f{jj} = casadi.Function(['build_n_' num2str(jj) '_casadi'],{z},{in_n{jj}.stored.values});
                    
                    sxfcn{end + 1} = y_f{jj}; %#ok
                    
                    % Get number of FLOPS
                    try
                        info.in_n{jj}.flops =  y_f{jj}.getAlgorithmSize();
                    catch
                        try
                            info.in_n{jj}.flops =  y_f{jj}.n_instructions();
                        catch
                            warning('Cannot use casadi.Function.getAlgorithmSize() or casadi.Function.n_instructions(), number of FLOPS may too low as a result.');
                            info.in_n{jj}.flops = 0;
                        end
                    end
                    
                    info.in_n{jj}.static = 0;
                    info.in_n{jj}.struct.structure = in_n{jj};
                    
                    % wrapper function
                    code = [code, sprintf(['/* Constraints evaluation*/' '\n'])]; %#ok
                    code = [code, sprintf(['void build_n_%d( const ' o.real '* z,'  o.real '* n){' '\n\n'],jj)]; %#ok
                    code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],1)]; %#ok
                    code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)]; %#ok
                    code = [code, sprintf(['\tin[0] = z;\n\n' ...
                        o.indent.generic  'build_n_%d_casadi( in, &n, &iw, w, mem);' o.indent.generic '/* external casadi generated function */\n}\n\n'],jj) ]; %#ok
                    
                    % build_Dn
                    Dn{jj} = jacobian( n{jj}, z);
                    Dn{jj} = transpose(Dn{jj});
                    [info.in_Dn_n{jj}.static,i] = falcopt.casadi2struct(Dn{jj});
                    
                    if( info.in_Dn_n{jj}.static) %is jac matrix constant ?
                        name.M = sprintf('build_Dn_%d_casadi',jj);
                        [d, ~, in] = falcopt.generateData(full(DM(Dn{jj})), 'names', name, ...
                            'type', o.real, 'precision', o.precision, 'structure', 'unique', 'noones', false, 'indent', o.indent, ...
                            'static', true, 'const', true, 'verbose', o.verbose);
                        info.in_Dn_n{jj}.struct = in;
                        info.in_Dn_n{jj}.flops = 0;
                        if ~isempty(d)
                            data = [data, d, sprintf('\n')]; %#ok
                        end
                    else
                        Dn_f{jj} = casadi.Function(sprintf('build_Dn_%d_casadi',jj),{z},{i.stored.values});
                        sxfcn{end + 1} =  Dn_f{jj}; %#ok
                        
                        % wrapper build_Dn
                        code = [code, sprintf('/* Jacobian_u of nonlinear constraints */\n')]; %#ok
                        code = [code, sprintf(['void build_Dn_%d( const ' o.real '* u, ' o.real '* Dn_fun){' '\n\n'],jj)]; %#ok
                        code = [code, sprintf([o.indent.generic  'const ' o.real ' *in[%i];\n'],1)]; %#ok
                        code = [code, sprintf([o.indent.generic  o.real ' w[%i];\n\tint iw = %i;\n\t' 'int mem = %i; \n\n'],3,0,0)]; %#ok
                        code = [code, sprintf(['\tin[0] = u;\n\n' ...
                            o.indent.generic  'build_Dn_%d_casadi( in, &Dn_fun, &iw, w, mem);' o.indent.generic '/* external casadi generated function*/\n}\n\n'],jj) ]; %#ok
                        info.in_Dn_n{jj}.struct.structure = i;
                        
                        % Get number of FLOPS
                        try
                            info.in_Dn_n{jj}.flops = Dn_f{jj}.getAlgorithmSize();
                        catch
                            try
                                info.in_Dn_n{jj}.flops = Dn_f{jj}.n_instructions();
                            catch
                                warning('Cannot use casadi.Function.getAlgorithmSize() or casadi.Function.n_instructions(), number of FLOPS may too low as a result.');
                                info.in_Dn_n{jj}.flops = 0;
                            end
                        end
                    end
                end
                
            end
            
            
        case {'manual','matlab'} % use matlab jacobian generation
            
            for jj=1:length(o.K_n)
                
                in_n = cell(length(o.K_n));
                Dn = cell(length(o.K_n));
                
                [c, in_n{jj}] = falcopt.fcn2struct( n{jj},o,'name', 'n' );
                code = [code, sprintf('/* Constraints evaluation*/ \n')]; %#ok
                code = [code, sprintf(['void build_n_' num2str(jj) '( const ' o.real '* u,'  o.real '* n){' '\n'])]; %#ok
                code = [code, c, sprintf('}\n\n')]; %#ok
                
                info.in_n{jj}.flops =  in_n{jj}.flops; %flops
                info.in_n{jj}.static = 0;
                info.in_n{jj}.struct.structure = in_n{jj}.structure;
                
                %build_Dn
                if(  isequal(o.gradients,'matlab'))
                    Dn{jj} = jacobian( n{jj},z);
                    Dn{jj} = Dn{jj}';
                else
                    try
                        Dn{jj} = o.external.manual.jacobian_n{jj}(z);
                    catch
                        error('o.external.manual.jacobian_n function not found. Or invalid output: it should return a struct of size o.N');
                    end
                end
                [c,i] = falcopt.fcn2struct( Dn{jj}, o, 'name', 'Dn_fun');
                code = [code, sprintf('/* Jacobian_u of nonlinear constraints*/\n')]; %#ok
                code = [code, sprintf(['void build_Dn_' num2str(jj) '( const ' o.real '* u, ' o.real '* Dn_fun){' '\n'])]; %#ok
                code = [code, c, sprintf('}\n\n')]; %#ok
                info.in_Dn_n{jj}.struct.structure = i.structure;
                info.in_Dn_n{jj}.flops =   i.flops; %flops
                info.in_Dn_n{jj}.static = 0;
            end
            
            
            
            sxfcn = {};
    end
    % generate a selector that chooses which n to use
    code = [code, sprintf(['\n' o.indent.code 'void build_n(const ' o.real '* u, const unsigned int k, ' o.real '* n){' '\n\n'])];
    for ii=1:length(o.K_n)
        code = [code, sprintf([o.indent.code o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_n{ii}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        code = [code, sprintf([o.indent.code o.indent.generic o.indent.generic 'build_n_' num2str(ii) '(&u[0], &n[0]);' '\n'])]; %#ok
        code = [code, sprintf([o.indent.code o.indent.generic  '}' '\n'])]; %#ok
    end
    code = [code, sprintf([o.indent.code '}' '\n'])];
    
    % generate a selector that chooses which Dn to use
    code = [code, sprintf(['\n' o.indent.code 'void build_Dn(const ' o.real '* u, const unsigned int k, ' o.real '* Dn){' '\n\n'])];
    for ii=1:length(o.K_n)
        code = [code, sprintf([o.indent.code o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_n{ii}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        code = [code, sprintf([o.indent.code o.indent.generic o.indent.generic 'build_Dn_' num2str(ii) '(&u[0], &Dn[0]);' '\n'])]; %#ok
        code = [code, sprintf([o.indent.code o.indent.generic  '}' '\n'])]; %#ok
    end
    code = [code, sprintf([o.indent.code '}' '\n'])];
end

% generate .c file
if( isequal(o.gradients,'casadi')&& ~isempty(sxfcn))
    [info.src,info.header] = generate_casadi_c(o, 'constraints', sxfcn);
end
end

function [c] = argument_w(o,declaration)
% add or not w as an argument of the functions
% declaration: logic variable, if true it adds the type of variable

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
% consider the additional argument for the constant term introduced by the
% contractive/terminal constraint
% declaration: logic variable, if true it adds the type of variable

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
% consider the additional arguments introduced by the contractive approach
% (index of the contractive constraint) and possible reference tracking
% declaration: logic variable, if true it adds the type of variable

if declaration
    str_real = sprintf(['const ' o.real '* ']);
    str_int  = sprintf('const unsigned int ');
else
    str_real = [];
    str_int = [];
end

c = [];

if o.objective.trackReference
    c = [c, sprintf([', ' str_real 'xref, ' str_real 'uref'])];
end
if o.contractive
    c = [c, sprintf([', ' str_int 'ind'])];
end
end

function [c] = argument_def_internal_psi(o,declaration)
% consider the additional argument introduced by the contractive/terminal
% approach (variable psi_N)
% declaration: logic variable, if true it adds const and the type of variable

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
% consider the additional argument introduced by the contractive/terminal
% approach (variable psi_N)
% declaration: logic variable, if true it adds the type of variable

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
% consider the additional argument introduced by the contractive/terminal
% approach (variable psi_Np)
% declaration: logic variable, if true it adds the type of variable

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
% consider the additional argument introduced by the contractive/terminal
% approach (variable dot_psi_N)
% declaration: logic variable, if true it adds the type of variable

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
% consider the additional argument introduced by the contractive/terminal
% approach (variable dot_psi_N)
% declaration: logic variable, if true it adds const and the type of variable

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
% function that generates the forward simulation 'det_x'

nx = o.nx;
nu = o.nu;
nw = o.nw;
N = o.N;

code = [];
data = [];

code = [code, sprintf(['\n' '/* det_x is a forward simulation of model_mpc with initial state x0 \n'...
    'and sequence of predicted inputs u (of dim. N*nu) */' '\n'])];
if o.nw >0
    code = [code, sprintf([o.inline 'void det_x (const ' o.real '* x0, const ' o.real '* u, const ' o.real '* w, ' o.real '* x){' '\n\n'])];
else
    code = [code, sprintf([o.inline 'void det_x (const ' o.real '* x0, const ' o.real '* u, ' o.real '* x){' '\n\n'])];
end

code = [code, sprintf([o.indent.generic  'unsigned int ii = 0;' '\n\n'])];

if o.nw > 0
    code = [code, sprintf([o.indent.generic  'model_mpc(x0,u,w,x);' '\n'])];
else
    code = [code, sprintf([o.indent.generic  'model_mpc(x0,u,x);' '\n'])];
end

code = [code, sprintf([o.indent.generic  'for (ii = 1;ii < %d; ++ii)' '\n'],N)];

if o.nw > 0
    code = [code, sprintf([o.indent.generic o.indent.generic 'model_mpc(x + (ii-1)* %d, u + ii* %d, w+ ii* %d, x + ii* %d);' '\n\n'],nx,nu,nw,nx)];
else
    code = [code, sprintf([o.indent.generic o.indent.generic 'model_mpc(x + (ii-1)* %d, u + ii* %d, x + ii* %d);' '\n\n'],nx,nu,nx)];
end
code = [code, sprintf(['}' '\n'])];

end

function [code, data, info] = generate_objective_gradient_oracle(o)
% generate the function that computes the gradient steps, along with
% auxiliary functions

nx = o.nx;
nu = o.nu;
nw = o.nw;
N = o.N;
trackRef = o.objective.trackReference;

code = [];
data = [];
info.flops.it = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
info.flops.ls = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

%it generates: dot_product_nx_nx, Qmul, Pmul, Rmul, dot_product_nu_nu,
%               product_and_sum_nu
[c, d, in] = generate_auxiliary_functions(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.internal.multFlops(in.flops, N);
info.flops.it = falcopt.internal.addFlops(info.flops.it,in.flops);
info.flops.ls = falcopt.internal.addFlops(info.flops.ls,in.flops);


[c, d, in] = generate_product_and_sum_nx(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.internal.multFlops(in.flops, N-1);
info.flops.it = falcopt.internal.addFlops(info.flops.it,in.flops);
info.flops.ls = falcopt.internal.addFlops(info.flops.ls,in.flops);


if trackRef
    [c, d, in] = generate_diffXU(o);
    code = [code, c];
    data = [data, d];
    in.flops = falcopt.internal.multFlops(in.flops, N);
    info.flops.it = falcopt.internal.addFlops(info.flops.it,in.flops);
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls,in.flops);
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

code = [code, sprintf([o.indent.generic  o.real ' Px[%d], Qx[%d], mem_tmp2[%d], ' '\n'...
    o.indent.generic o.indent.generic 'Ru[%d], tmp_x = 0.0, tmp_u = 0.0;' '\n'],...
    nx,nx,nx,nu)];

if o.contractive
    code = [code, sprintf([o.indent.generic  'int index = 0;' '\n'])];
end

if trackRef
    code = [code, sprintf([o.indent.generic  o.real ' dx[%d], du[%d];' '\n'],...
        nx,nu)];
end
if o.contractive
    code = [code, sprintf([o.indent.generic  o.real ' Px_contr[%d], mem_tmp_contr[%d], tmp_contr = 0.0;' '\n'], nx, nx)];
    if trackRef
        code = [code, sprintf([o.indent.generic  o.real ' dx_contr[%d];' '\n'], nx)];
    end
elseif o.terminal
    code = [code, sprintf([o.indent.generic  o.real ' Px_contr[%d], mem_tmp_contr[%d];' '\n'], nx, nx)];
end

code = [code, sprintf([o.indent.generic  'unsigned int ii = 0;' '\n\n'])];

code = [code, sprintf(['\n' o.indent.generic  '/* Compute tmp_x (cost associated to last stage) */' '\n'])];
if trackRef
    code = [code, sprintf([o.indent.generic  'diffX(dx, x + %d, xref + %d);' '\n'],(N-1)*nx, (N-1)*nx)];
    code = [code, sprintf([o.indent.generic  'Pmul(Px, dx);' '\n'])];
    code = [code, sprintf([o.indent.generic  'dot_product_nx_nx(&tmp_x,Px, dx);' '\n\n'])];
    if o.contractive
        code = [code, sprintf([o.indent.generic  'diffX(dx_contr, x + (ind-1)*%d, xref + (ind-1)*%d);' '\n'],nx, nx)];
        code = [code, sprintf([o.indent.generic  'Pmul(Px_contr, dx_contr);' '\n'...
            o.indent.generic  'dot_product_nx_nx(&tmp_contr,Px_contr,dx_contr);' '\n'...
            o.indent.generic  '(*psi_N) = 0.5*tmp_contr;' '\n'])];
    end
else
    code = [code, sprintf([o.indent.generic  'Pmul(Px, x + %d);' '\n'], (N-1)*nx)];
    code = [code, sprintf([o.indent.generic  'dot_product_nx_nx(&tmp_x,Px, x + %d);' '\n\n'],(N-1)*nx)];
    if o.contractive
        code = [code, sprintf([o.indent.generic  'Pmul(Px_contr, x+ (ind - 1)*%d);' '\n'...
            o.indent.generic  'dot_product_nx_nx(&tmp_contr,Px_contr,x+ (ind - 1)*%d);' '\n'...
            o.indent.generic  '(*psi_N) = 0.5*tmp_contr;' '\n'],nx,nx)];
    end
end

if o.terminal
    code = [code, sprintf([o.indent.generic  '(*psi_N) = 0.5*tmp_x;' '\n',...
        o.indent.generic  'copy_nx(Px_contr,Px);' '\n'])];
    
end

if trackRef
    code = [code, sprintf([o.indent.generic  'diffU(du, u + %d, uref + %d);' '\n'],(N-1)*nu, (N-1)*nu)];
    code = [code, sprintf([o.indent.generic  'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf([o.indent.generic  'dot_product_nu_nu(&tmp_u,Ru,du);' '\n'])];
else
    code = [code, sprintf([o.indent.generic  'Rmul(Ru, u + %d);' '\n'], (N-1)*nu)];
    code = [code, sprintf([o.indent.generic  'dot_product_nu_nu(&tmp_u,Ru, u + %d);' '\n'],(N-1)*nu)];
end

code = [code, sprintf([o.indent.generic  '(*J) = 0.5*(tmp_x + tmp_u);' '\n'])];


info.flops.it.mul = info.flops.it.mul+1;
info.flops.it.add = info.flops.it.add+1;


if o.nw > 0
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic  'Jacobian_u(x + %d,u + %d, w + %d, G);' '\n'],(N-2)*nx,(N-1)*nu,(N-1)*nw)];
    end
else
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic  'Jacobian_u(x + %d,u + %d, G);' '\n'],(N-2)*nx,(N-1)*nu)];
    end
end

code = [code, sprintf(['\n' o.indent.generic  '/* Start computing dot_J from the bottom */' '\n'])];
code = [code, sprintf([o.indent.generic  'product_and_sum_nu(dot_J + %d,Px, Ru, G);' '\n\n'], (N-1)*nu)];

if o.terminal
    code = [code, sprintf([o.indent.generic  'product_contr_nu(&dot_psi_N[%d], Px, G);' '\n\n'], (N-1)*nu)];
elseif o.contractive
    code = [code, sprintf([o.indent.generic  'if (ind == %d) ' '\n'...
        o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[%d], Px_contr, G);' '\n'], N, (N-1)*nu)];
    code = [code, sprintf([o.indent.generic  'else' '\n'...
        o.indent.generic o.indent.generic 'set_zero_nu(&dot_psi_N[%d]);' '\n\n'], (N-1)*nu)];
end



if o.contractive
    code = [code, sprintf([o.indent.generic  'if (ind == %d) ' '\n'...
        o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[%d], Px_contr, G);' '\n'], N, (N-1)*nu)];
    code = [code, sprintf([o.indent.generic  'else' '\n'...
        o.indent.generic o.indent.generic 'set_zero_nu(&dot_psi_N[%d]);' '\n\n'], (N-1)*nu)];
end

code = [code, sprintf([o.indent.generic  'for (ii=%d; ii-->0; ) {' '\n'],N-1)];

if o.nw > 0
    if ~o.Jac_x_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'Jacobian_x(x + ii*%d,u + (ii+1)*%d, w + (ii+1)*%d, F);' '\n'],nx,nu,nw)];
        info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
        
    end
else
    if ~o.Jac_x_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'Jacobian_x(x + ii*%d,u + (ii+1)*%d, F);' '\n'],nx,nu)];
        info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
    end
end

if o.nw > 0
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'if (ii==0)' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x0,u,w,G);' '\n',...
            o.indent.generic o.indent.generic 'else' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x + (ii-1)*%d, u + ii*%d, w + ii*%d, G);' '\n\n'],nx,nu,nw)];
        info.flops.it.mul = info.flops.it.mul+ 3*(N-1);
        info.flops.it.comp = info.flops.it.comp+ (N-1);
    end
else
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'if (ii==0)' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x0,u,G);' '\n',...
            o.indent.generic o.indent.generic 'else' '\n',...
            o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x + (ii-1)*%d, u + ii*%d, G);' '\n\n'],nx,nu)];
        info.flops.it.mul = info.flops.it.mul+ 3*(N-1);
        info.flops.it.comp = info.flops.it.comp+ (N-1);
    end
end


if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffX(dx, x + ii*%d, xref + ii*%d);' '\n'],nx, nx)];
    info.flops.it.mul = info.flops.it.mul+ 3*(N-1);
    code = [code, sprintf([o.indent.generic o.indent.generic 'Qmul(Qx, dx);' '\n'])];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nx_nx(&tmp_x,Qx,dx);' '\n'])];
else
    code = [code, sprintf([o.indent.generic o.indent.generic 'Qmul(Qx, x + ii*%d);' '\n'], nx)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nx_nx(&tmp_x,Qx, x + ii*%d);' '\n'],nx)];
    info.flops.it.mul = info.flops.it.mul+ (nx+nx)*(N-1);
end

if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffU(du, u + ii*%d, uref + ii*%d);' '\n'],nu, nu)];
    info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
    code = [code, sprintf([o.indent.generic o.indent.generic 'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nu_nu(&tmp_u,Ru,du);' '\n'])];
else
    code = [code, sprintf([o.indent.generic o.indent.generic 'Rmul(Ru, u + ii*%d);' '\n'], nu)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nu_nu(&tmp_u,Ru, u + ii*%d);' '\n'],nu)];
    info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
end

code = [code, sprintf(['\n' o.indent.generic o.indent.generic '/* Increment J */' '\n'])];
code = [code, sprintf([o.indent.generic o.indent.generic '(*J) += .50*(tmp_x + tmp_u);' '\n\n'])];
info.flop.it.mul = info.flops.it.mul+1*(N-1);
info.flops.it.add = info.flops.it.add+2*(N-1);



if o.terminal
    code = [code, sprintf([o.indent.generic o.indent.generic 'product_contr_nx(mem_tmp_contr, Px_contr, F);' '\n',...
        o.indent.generic o.indent.generic 'copy_nx(Px_contr,mem_tmp_contr);' '\n',...
        o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[ii*%d], Px_contr, G);' '\n\n'], nu)];
end

if o.contractive
    code = [code, sprintf([o.indent.generic o.indent.generic 'index = ind - ii - 1;' '\n'])];
    
    code = [code, sprintf([o.indent.generic o.indent.generic 'if ( index >= 0){' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'if (index > 0){' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'product_contr_nx(mem_tmp_contr, Px_contr, F);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'copy_nx(Px_contr,mem_tmp_contr);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[ii*%d], Px_contr, G);' '\n'...
        o.indent.generic o.indent.generic '}' '\n'], nu)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'else' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'set_zero_nu(&dot_psi_N[ii*%d]);' '\n\n'],nu)];
    
end

code = [code, sprintf(['\n' o.indent.generic o.indent.generic '/* Compute dot_J */' '\n'])];
code = [code, sprintf([o.indent.generic o.indent.generic 'product_and_sum_nx(mem_tmp2, Px, Qx, F);' '\n',...
    o.indent.generic o.indent.generic 'copy_nx(Px,mem_tmp2);' '\n',...
    o.indent.generic o.indent.generic 'product_and_sum_nu(dot_J + ii*%d, Px, Ru, G);' '\n',...
    o.indent.generic  '}' '\n',...
    '}' '\n\n'],nu)];
info.flops.it.mul = info.flops.it.mul+ 1*(N-1);

%% det_J


code = [code, sprintf([o.inline ' void det_J(const ' o.real '* x0, const ' o.real ...
    '* u, const ' o.real '* x' c_tr_dec ', ' ...
    o.real '* J' c_psi_dec '){' '\n\n',...
    ])];

code = [code, sprintf([o.indent.generic  o.real ' Qx[%d], Ru[%d], tmp_x = 0.0, tmp_u = 0.0;' '\n'],...
    nx,nu)];
if trackRef
    code = [code, sprintf([o.indent.generic  o.real ' dx[%d], du[%d];' '\n'],...
        nx,nu)];
end
if o.contractive
    code = [code, sprintf([o.indent.generic  o.real ' Px_contr[%d], tmp_contr = 0.0;' '\n'], nx)];
    if trackRef
        code = [code, sprintf([o.indent.generic  o.real ' dx_contr[%d];' '\n'], nx)];
    end
end

code = [code, sprintf([o.indent.generic  'unsigned int ii = 0;' '\n\n'])];

if trackRef
    code = [code, sprintf([o.indent.generic  'diffX(dx, x + %d, xref + %d);' '\n'],(N-1)*nx, (N-1)*nx)];
    code = [code, sprintf([o.indent.generic  'Pmul(Qx, dx);' '\n'])];
    code = [code, sprintf([o.indent.generic  'dot_product_nx_nx(&tmp_x,Qx, dx);' '\n'])];
    if o.contractive
        code = [code, sprintf([o.indent.generic  'diffX(dx_contr, x + (ind - 1)*%d, xref + (ind-1)*%d);' '\n'],nx, nx)];
        code = [code, sprintf([o.indent.generic  'Pmul(Px_contr, dx_contr);' '\n'...
            o.indent.generic  'dot_product_nx_nx(&tmp_contr,Px_contr,dx_contr);' '\n',...
            o.indent.generic  '(*psi_N) = 0.5*tmp_contr;' '\n'])];
    end
else
    code = [code, sprintf([o.indent.generic  'Pmul(Qx, x + %d);' '\n'], (N-1)*nx)];
    code = [code, sprintf([o.indent.generic  'dot_product_nx_nx(&tmp_x,Qx, x + %d);' '\n'],(N-1)*nx)];
    if o.contractive
        code = [code, sprintf([o.indent.generic  'Pmul(Px_contr, x+ (ind - 1)*%d);' '\n'...
            o.indent.generic  'dot_product_nx_nx(&tmp_contr,Px_contr,x+ (ind - 1)*%d);' '\n',...
            o.indent.generic  '(*psi_N) = 0.5*tmp_contr;' '\n'],nx,nx)];
    end
end

if o.terminal
    code = [code, sprintf([o.indent.generic  '(*psi_N) = 0.5*tmp_x;' '\n'])];
end

if trackRef
    code = [code, sprintf([o.indent.generic  'diffU(du, u + %d, uref + %d);' '\n'],(N-1)*nu, (N-1)*nu)];
    code = [code, sprintf([o.indent.generic  'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf([o.indent.generic  'dot_product_nu_nu(&tmp_u,Ru,du);' '\n\n'])];
else
    code = [code, sprintf([o.indent.generic  'Rmul(Ru, u + %d);' '\n'], (N-1)*nu)];
    code = [code, sprintf([o.indent.generic  'dot_product_nu_nu(&tmp_u,Ru, u + %d);' '\n\n'],(N-1)*nu)];
end

code = [code, sprintf([o.indent.generic  '(*J) = 0.5*(tmp_x + tmp_u);' '\n'])];
info.flops.ls.mul = info.flops.ls.mul+1;
info.flops.ls.add = info.flops.ls.add+1;

code = [code, sprintf([o.indent.generic  'for (ii=%d; ii-->0; ) {' '\n\n'],N-1)];

if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffX(dx, x + ii*%d, xref + ii*%d);' '\n'],nx, nx)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
    code = [code, sprintf([o.indent.generic o.indent.generic 'Qmul(Qx, dx);' '\n'])];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nx_nx(&tmp_x,Qx,dx);' '\n\n'])];
else
    code = [code, sprintf([o.indent.generic o.indent.generic 'Qmul(Qx, x + ii*%d);' '\n'], nx)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nx_nx(&tmp_x,Qx, x + ii*%d);' '\n\n'],nx)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
end

if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffU(du, u + ii*%d, uref + ii*%d);' '\n'],nu, nu)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
    code = [code, sprintf([o.indent.generic o.indent.generic 'Rmul(Ru, du);' '\n'])];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nu_nu(&tmp_u,Ru,du);' '\n\n'])];
else
    code = [code, sprintf([o.indent.generic o.indent.generic 'Rmul(Ru, u + ii*%d);' '\n'], nu)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'dot_product_nu_nu(&tmp_u,Ru, u + ii*%d);' '\n\n'],nu)];
    info.flops.ls.mul = info.flops.ls.mul+ 2*(N-1);
end


code = [code, sprintf([o.indent.generic o.indent.generic '(*J) += .5*(tmp_x + tmp_u);' '\n'...
    o.indent.generic  '}' '\n'...
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
Q = o.objective.Q;
P = o.objective.P;
R = o.objective.R;

if ~isempty(Q)
    code = [code, sprintf(['\n' '/* It computes Q*x */' '\n'])];
    [d, c, in] = falcopt.generateMVMult(Q, ...
        'names', struct('fun', 'Qmul', 'M', 'Q', 'v', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end

if ~isempty(P)
    code = [code, sprintf(['\n' '/* It computes P*x */' '\n'])];
    [d, c] = falcopt.generateMVMult(P, ...
        'names', struct('fun', 'Pmul', 'M', 'P', 'v', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    % flops of Pmul are counted together with Qmul
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
end


if ~isempty(R)
    code = [code, sprintf(['\n' '/* It computes R*u */' '\n'])];
    [d, c, in] = falcopt.generateMVMult(R, ...
        'names', struct('fun', 'Rmul', 'M', 'R', 'v', 'du'), 'types', o.real,...
        'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    
    if ~isempty(d)
        data = [data, d, sprintf('\n')];
    end
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end

code = [code, sprintf(['\n' '/* dot product x^top *x */' '\n'])];
[d, c, in] = falcopt.generateMVMult(ones(1,nx), ...
    'names', struct('fun', 'dot_product_nx_nx', 'M', 'R',...
    'v', 'du'), 'static', false, 'types', o.real, 'verbose', o.verbose, 'precision', o.precision, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

code = [code, sprintf(['\n' '/* dot product u^top *u */' '\n'])];
[d, c, in] = falcopt.generateMVMult(ones(1,nu), ...
    'names', struct('fun', 'dot_product_nu_nu', 'M', 'R',...
    'v', 'du'), 'static', false, 'types', o.real, 'verbose', o.verbose, 'precision', o.precision, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

code = [code, sprintf(['\n' '/* It computes G^top * x + u (exploiting structure of G) */' '\n'])];
[d, c, in] = falcopt.generateMVMult({o.Jac_u_struct,eye(nu)}, ...
    'names', struct('fun', 'product_and_sum_nu', 'M', {{'A', 'I'}},...
    'v', {{'u1', 'u2'}}), 'static', [false,true], 'structure', 'ordered', 'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

end

function [code, data, info] = generate_product_and_sum_nx(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nx = o.nx;
nu = o.nu;

if o.contractive
    code = [code, sprintf(['\n' '/* It generates a null vector */' '\n'])];
    code =  [code, sprintf(['void set_zero_nu (' o.real '* x){\n\n'])];
    
    for jj=0:nu-1
        code =  [code, sprintf([o.indent.generic  'x[%d] = 0.0;' '\n'],jj)]; %#ok
    end
    
    code =  [code, sprintf(['\n\n' '}' '\n\n'])];
end

if o.contractive|| o.terminal
    code = [code, sprintf(['\n' '/* It computes F^top * u (exploiting structure of F) */' '\n'])];
    [~, c, in] = falcopt.generateMVMult({o.Jac_x_struct}, ...
        'names', struct('fun', 'product_contr_nx', 'M', {{'A'}},...
        'v', {{'u'}}), 'static', false, 'structure', 'ordered',...
        'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test,...
        'inline', o.inline, 'indent', o.indent);
    
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
    
    code = [code, sprintf(['\n' '/* It computes G^top * u (exploiting structure of G) */' '\n'])];
    [~, c, in] = falcopt.generateMVMult({o.Jac_u_struct}, ...
        'names', struct('fun', 'product_contr_nu', 'M', {{'A'}},...
        'v', {{'u'}}), 'static', false, 'structure', 'ordered',...
        'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test,...
        'inline', o.inline, 'indent', o.indent);
    
    code = [code, c, sprintf('\n\n')];
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
    
end

code = [code, sprintf(['\n' '/* It computes F^top * u1 (exploiting structure of F) + u2 */' '\n'])];
[d, c, in] = falcopt.generateMVMult({o.Jac_x_struct,eye(nx)}, ...
    'names', struct('fun', 'product_and_sum_nx', 'M', {{'A', 'I'}},...
    'v', {{'u1', 'u2'}}), 'static', [false,true], 'structure', 'ordered', 'transpose', true, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

end

function [code, data, info] = generate_diffXU(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

nx = o.nx;
nu = o.nu;

code = [code, sprintf(['\n' '/* It computes dx = x - xref */' '\n'])];
[d, c, in] = falcopt.generateMVMult({eye(nx), -eye(nx)}, ...
    'names', struct('fun', 'diffX', 'M', {{'I', 'mI'}}, 'v', {{'x', 'xref'}}, 'r', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

code = [code, sprintf(['\n' '/* It computes du = u - uref */' '\n'])];
[d, c, in] = falcopt.generateMVMult({eye(nu), -eye(nu)}, ...
    'names', struct('fun', 'diffU', 'M', {{'I', 'mI'}}, 'v', {{'u', 'uref'}}, 'r', 'du'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end

function [code, data, info, optCode] = generate_slack_initialization(o)

code = [];
data = [];
optCode = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

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
    o.indent.generic  'unsigned int jj;' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code,  sprintf([o.indent.generic  'for (jj=0;jj< na;jj++){' '\n',...
        o.indent.generic o.indent.generic 'sl_sqr[jj] = ' o.auxFunctions.max '(1.0, -2.0* amu[jj]);' '\n',...
        o.indent.generic o.indent.generic 'gps[jj]= amu[jj] + 0.5*sl_sqr[jj];' '\n',...
        o.indent.generic o.indent.generic 'sl[jj] = ' o.auxFunctions.sqrt '(sl_sqr[jj]);' '\n',...
        o.indent.generic  '}' '\n'])];
end
if ~isempty(o.K_ub)
    code = [code,  sprintf([o.indent.generic  'for (jj=na;jj< nb;jj++){' '\n'...
        o.indent.generic o.indent.generic 'sl_sqr[jj] = ' o.auxFunctions.max '(1.0, -2.0* umb[jj-na]);' '\n',...
        o.indent.generic o.indent.generic 'gps[jj]= umb[jj-na] + 0.5*sl_sqr[jj];' '\n',...
        o.indent.generic o.indent.generic 'sl[jj] = ' o.auxFunctions.sqrt '(sl_sqr[jj]);' '\n',...
        o.indent.generic  '}' '\n'])];
end
if ~isempty(o.K_n)
    code = [code,  sprintf([o.indent.generic  'for (jj= nb;jj< nc;jj++) {' '\n'...
        o.indent.generic o.indent.generic 'sl_sqr[jj] = ' o.auxFunctions.max '(1.0, -2.0* n[jj - nb]);' '\n',...
        o.indent.generic o.indent.generic 'gps[jj]= n[jj - nb] + 0.5*sl_sqr[jj];' '\n',...
        o.indent.generic o.indent.generic 'sl[jj] = ' o.auxFunctions.sqrt '(sl_sqr[jj]);' '\n',...
        o.indent.generic  '}' '\n'])];
end

code = [code,  sprintf(['}' '\n\n'])];

code = [code, sprintf(['\n' '/* It initialize the slack variables sl and its squares sl_sqr \n'...
    'such that, if possible, gps = g + 0.5*sl_sqr = 0 */' '\n'])];
code = [code, sprintf([o.inline ' void initialize_slack( const ' o.real '* u' c_psi_dec c_contr_dec ', ' o.real '* sl, ' o.real '* sl_sqr, ' o.real '* gps){' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code, sprintf([o.indent.generic   o.real ' amu[%d];' '\n'], max(sum(~isinf( o.lb))) )];
end
if ~isempty(o.K_ub)
    code = [code, sprintf([o.indent.generic   o.real ' umb[%d];' '\n'], max(sum(~isinf( o.ub))) )];
end
if ~isempty(o.K_n)
    code = [code, sprintf([o.indent.generic   o.real ' n[%d];' '\n'], max(o.nn) )];
end

if (o.contractive || o.terminal)
    code = [code, sprintf([o.indent.generic  o.real ' g_contr = 0.0;' '\n'])];
end

info.flops.mul = info.flops.mul+ N*3; % mul inside only first cycle % NOT CLEAR, ToDo


lb = sum(~isinf(o.lb),1);
ub = sum(~isinf(o.ub),1);

for k = 1:o.N
    code = [code, sprintf(['\n'...
        o.indent.generic  '/* Unrolling the for loop: iteration %i of %i */' '\n'], k-1, o.N - 1)];   %#ok
    if ~isempty(o.K_lb)
        code = [code, sprintf([o.indent.generic  'build_amu(&u[%i], %i, &amu[0]);' '\n'],o.nu*(k-1),k-1)];%#ok
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf([o.indent.generic  'build_umb(&u[%i],%i,&umb[0]);' '\n'],o.nu*(k-1),k-1)];%#ok
    end
    if ~isempty(o.K_n)
        code = [code, sprintf([o.indent.generic  'build_n(&u[%i],%i,&n[0]);' '\n'],o.nu*(k-1),k-1)];%#ok
    end
    code = [code, sprintf([o.indent.generic  'build_sl_slsqr('])]; %#ok
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

if (o.contractive || o.terminal)
    code = [code, sprintf([o.indent.generic  'g_contr = *psi_N - c_contr;' '\n'...
        o.indent.generic  'sl_sqr[%d] = ' o.auxFunctions.max '(1.0, -2.0*g_contr);' '\n',...
        o.indent.generic  'sl[%d] = ' o.auxFunctions.sqrt '(sl_sqr[%d]);' '\n'...
        o.indent.generic  'gps[%d] = g_contr + 0.5*sl_sqr[%d];' '\n'],...
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

nu = o.nu;
N = o.N;

c_psi_dot_dec = argument_def_internal_psi_dot(o,true);


info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);



%% gradient_step

% lower bounds
for jj = 1:length(o.K_lb) % if o.K_lb is empty, this loop is not executed
    struct_mat = diag(~isinf(o.lb(:,o.K_lb{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    
    if o.variable_stepSize.active
        M_struct = {eye(sum(~isinf(o.lb(:,o.K_lb{jj}(1))))),struct_mat_cut};
        M_names = {{'alpha_inv', ['I_lb_' num2str(jj)]}};
        M_static = [false,true];
    else
        M_struct = {1/o.variable_stepSize.alpha_max*eye(sum(~isinf(o.lb(:,o.K_lb{jj}(1))))), struct_mat_cut};
        M_names = {{ ['alpha_inv_lb_' num2str(jj)], ['I_lb_' num2str(jj)]}};
        M_static = [true,true];
    end
    
    [d, c, in] = falcopt.generateMVMult(M_struct, ...
        'names', struct('fun', ['build_vNnc_lb_' num2str(jj)], 'M', M_names,...
        'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent,'static',M_static,'structure','ordered');
    if ~isempty(d)
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end
if ~isempty(o.K_lb)
    code = [code, sprintf([ o.inline ' void build_vNnc_lb(const ' o.real '* gps, '...
        'const ' o.real '* dot_J, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_lb)
        code = [code, sprintf([o.indent.code o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_lb{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        if o.variable_stepSize.active
            code = [code, sprintf([o.indent.generic o.indent.generic 'build_vNnc_lb_%i(&res[0], &gps[0], &dot_J[0], &alpha_inverse);' '\n'],jj)];%#ok
        else
            code = [code, sprintf([o.indent.generic o.indent.generic 'build_vNnc_lb_%i(&res[0], &gps[0], &dot_J[0]);' '\n'],jj)]; %#ok
        end
        code = [code, sprintf([o.indent.generic '}' '\n'])];%#ok
    end
    code = [code, sprintf(['}' '\n'])];
end

% upper bounds
for jj = 1:length(o.K_ub) % if o.K_ub is empty, this loop is not executed
    struct_mat = diag(~isinf(o.ub(:,o.K_ub{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    
    if o.variable_stepSize.active
        M_struct = {eye(sum(~isinf(o.ub(:,o.K_ub{jj}(1))))), -struct_mat_cut};
        M_names = {{ 'alpha_inv', ['I_ub_' num2str(jj)]}};
        M_static = [false,true];
    else
        M_struct = {1/o.variable_stepSize.alpha_max*eye(sum(~isinf(o.ub(:,o.K_ub{jj}(1))))), -struct_mat_cut};
        M_names = {{ ['alpha_inv_ub_' num2str(jj)], ['I_ub_' num2str(jj)]}};
        M_static = [true,true];
    end
    
    [d, c, in] = falcopt.generateMVMult(M_struct, ...
        'names', struct('fun', ['build_vNnc_ub_' num2str(jj)], 'M', M_names,...
        'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent,'static',M_static,'structure','ordered');
    if ~isempty(d)
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end

if ~isempty(o.K_ub)
    code = [code, sprintf([ o.inline ' void build_vNnc_ub(const ' o.real '* gps, '...
        'const ' o.real '* dot_J, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_ub)
        code = [code, sprintf([o.indent.code o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_ub{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        if o.variable_stepSize.active
            code = [code, sprintf([o.indent.generic o.indent.generic 'build_vNnc_ub_%i(&res[0], &gps[0], &dot_J[0], &alpha_inverse);' '\n'],jj)]; %#ok
        else
            code = [code, sprintf([o.indent.generic o.indent.generic 'build_vNnc_ub_%i(&res[0], &gps[0], &dot_J[0]);' '\n'],jj)]; %#ok
        end
        code = [code, sprintf([o.indent.generic '}' '\n'])]; %#ok
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
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
    
    if o.variable_stepSize.active
        M_struct = {eye(o.nn(o.K_n{jj}(1))), -eye(o.nn(o.K_n{jj}(1)))};
        M_names = {{ 'alpha_inv', ['temp_n' num2str(jj)]}};
        M_static = [false,true];
    else
        M_struct = {1/o.variable_stepSize.alpha_max*eye(o.nn(o.K_n{jj}(1))), -eye(o.nn(o.K_n{jj}(1)))};
        M_names = {{ ['alpha_inv_n_' num2str(jj)], ['temp_n' num2str(jj)]}};
        M_static = [true,true];
    end
    
    [d, c, in] = falcopt.generateMVMult(M_struct, ...
        'names', struct('fun', ['build_vNnc_n_' num2str(jj)], 'M', M_names,...
        'v', {{'x1', 'x2'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'static', M_static,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent,'structure','ordered');
    if ~isempty(d)
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end

if ~isempty(o.K_n)
    code = [code, sprintf([ o.inline ' void Dntop_times_dotJ_n(const ' o.real '* Dn, '...
        'const ' o.real '* x, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_n)
        code = [code, sprintf([o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_n{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        code = [code, sprintf([o.indent.generic o.indent.generic 'Dntop_times_dotJ_%i(&res[0], &x[0], &Dn[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf([o.indent.generic  '}' '\n'])]; %#ok
    end
    code = [code, sprintf(['}' '\n'])];
    
    code = [code, sprintf([ o.inline ' void build_vNnc_n(const ' o.real '* gps, '...
        'const ' o.real '* temp_n, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_n)
        code = [code, sprintf([o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_n{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        if o.variable_stepSize.active
            code = [code, sprintf([o.indent.generic o.indent.generic 'build_vNnc_n_%i(&res[0], &gps[0], &temp_n[0], &alpha_inverse);' '\n'],jj)];%#ok
        else
            code = [code, sprintf([o.indent.generic o.indent.generic 'build_vNnc_n_%i(&res[0], &gps[0], &temp_n[0]);' '\n'],jj)]; %#ok
        end
        code = [code, sprintf([o.indent.generic '}' '\n'])];%#ok
    end
    code = [code, sprintf(['}' '\n'])];
end

if (o.contractive || o.terminal)
    
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
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
    
end

% lower bounds
for jj = 1:length(o.K_lb) % if o.K_lb is empty, this loop is not executed
    struct_mat = diag(~isinf(o.lb(:,o.K_lb{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    [d, c, in] = falcopt.generateMVMult( - struct_mat_cut', ...
        'names', struct('fun', ['minus_Ina_muG_' num2str(jj)], 'M', {{['minus_Ina_' num2str(jj)]}},...
        'v', {{'x1'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end
if ~isempty(o.K_lb)
    code = [code, sprintf([ o.inline ' void minus_Ina_muG(const ' o.real '* muG, '...
        'const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_lb)
        code = [code, sprintf([o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_lb{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        code = [code, sprintf([o.indent.generic o.indent.generic 'minus_Ina_muG_%i(&res[0], &muG[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf([o.indent.generic  '}' '\n'])]; %#ok
    end
    code = [code, sprintf(['}' '\n'])];
end

%upper bounds
for jj = 1:length(o.K_ub) % if o.K_ub is empty, this loop is not executed
    struct_mat = diag(~isinf(o.ub(:,o.K_ub{jj}(1))));
    struct_mat_cut = double(struct_mat(any(struct_mat,1),:));
    [d, c, in] = falcopt.generateMVMult( struct_mat_cut', ...
        'names', struct('fun', ['Inb_muG_' num2str(jj)], 'M', {{['Inb_' num2str(jj)]}},...
        'v', {{'x1'}}, 'r', 'z'), 'types', o.real, 'precision', o.precision,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
    if ~isempty(d)
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end
if ~isempty(o.K_ub)
    code = [code, sprintf([ o.inline ' void Inb_muG(const ' o.real '* muG, '...
        'const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_ub)
        code = [code, sprintf([o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_ub{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        code = [code, sprintf([o.indent.generic o.indent.generic 'Inb_muG_%i(&res[0], &muG[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf([o.indent.generic  '}' '\n'])]; %#ok
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
        data = [data, d, sprintf('\n')]; %#ok
    end
    code = [code, c, sprintf('\n\n')]; %#ok
    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
    
end

if ~isempty(o.K_n)
    code = [code, sprintf([ o.inline ' void Dn_times_muG_n(const ' o.real '* Dn, '...
        'const ' o.real '* x, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj=1:length(o.K_n)
        code = [code, sprintf([o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_n{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        code = [code, sprintf([o.indent.generic o.indent.generic 'Dn_times_muG_%i(&res[0], &x[0], &Dn[0]);' '\n'],jj)]; %#ok
        code = [code, sprintf([o.indent.generic  '}' '\n'])]; %#ok
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
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

if o.variable_stepSize.active
    M_struct = eye(nu);
    M_static = false;
else
    M_struct = o.variable_stepSize.alpha_max*eye(nu);
    M_static = true;
end

[d, c, in] = falcopt.generateMVMult(M_struct, ...
    'names', struct('fun', 'scale_nu', 'M', {{'alpha_I'}}, 'v', {{'x'}}, 'r', 'z'), 'scale', -1, ...
    'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent,'static',M_static,'structure','ordered');
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

for jj=1:length(o.K_nc)
    
    if o.variable_stepSize.active
        M_struct = eye(o.nc(o.K_nc{jj}(1)));
        M_static = false;
    else
        M_struct = o.variable_stepSize.alpha_max*eye(o.nc(o.K_nc{jj}(1)));
        M_static = true;
    end
    [~, c, in] = falcopt.generateMVMult(M_struct, ...
        'names', struct('fun', ['scale_nc_' num2str(jj)], 'M', {{'alpha_I'}}, 'v', {{'x'}}, 'r', 'z'), 'scale', -1, ...
        'types', o.real, 'precision', o.precision, 'static', M_static,...
        'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent,'structure','ordered');
  
    code = [code, c, sprintf('\n\n')];%#ok

    info.flops = falcopt.internal.addFlops(info.flops, in.flops);
    
    code = [code, sprintf([o.inline ' void product_matlab_nc_' num2str(jj) '(const ' o.real '* x, const '...
        o.real '* y, ' o.real '* z){' '\n\n',...
        o.indent.generic  'unsigned int ii=0;' '\n\n',...
        o.indent.generic  'for (ii=0;ii<%d;ii++)' '\n',...
        o.indent.generic o.indent.generic 'z[ii] = x[ii]*y[ii];' '\n\n',...
        '}' '\n'],o.nc(o.K_nc{jj}(1)) )]; %#ok
    info.flops.mul = info.flops.mul+ o.nc(o.K_nc{jj}(1));
end
if ~isempty(o.K_nc)
    code = [code, sprintf([o.inline ' void minus_scale_nc(const ' o.real ...
        '* x, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj = 1:length(o.K_nc)
        code = [code, sprintf([o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_nc{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        if o.variable_stepSize.active
            code = [code, sprintf([o.indent.generic o.indent.generic 'scale_nc_%i( &res[0], &x[0], &alpha);' '\n'], jj)];%#ok
        else
            code = [code, sprintf([o.indent.generic o.indent.generic 'scale_nc_%i( &res[0], &x[0]);' '\n'], jj)]; %#ok
        end
        code = [code, sprintf([o.indent.generic '}' '\n'])]; %#ok
    end
    code = [code, sprintf(['}' '\n'])];
    
    code = [code, sprintf([o.inline ' void product_matlab_nc(const ' o.real ...
        '* x, const ' o.real '* y, const unsigned int k, ' o.real '* res){' '\n'])];
    for jj = 1:length(o.K_nc)
        code = [code, sprintf([o.indent.generic 'if(' falcopt.internal.generateRangeCheck(o.K_nc{jj}-1, 'k', 'precision', 'unsigned integer') ') {' '\n'])]; %#ok
        code = [code, sprintf([o.indent.generic o.indent.generic 'product_matlab_nc_%i(&x[0], &y[0], &res[0]);' '\n'], jj)]; %#ok
        code = [code, sprintf([o.indent.generic  '}' '\n'])]; %#ok
    end
    code = [code, sprintf(['}' '\n'])];
end

if (o.contractive || o.terminal)
        ter_struct = repmat({ones(nu,1)},1,o.N); 
else
        ter_struct = repmat({[]},1,o.N);
end

[c, in] = falcopt.generateConstraintInv(o.Jac_n_struct_hor, ter_struct, 'N', N, 'bounds', struct('lb', ~isinf(o.lb), 'ub', ~isinf(o.ub)),  'types', o.real,'precision', o.precision, ...
    'indent', o.indent, 'inline', o.inline, 'verbose', max(0,o.verbose-1), 'test', o.test);
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

% Gradient step function
code = [code, sprintf(['\n' '/* Gradient step function */' '\n'])];
code = [code, sprintf([o.inline ' void gradient_step(const ' o.real '* dot_J, const ' o.real '* u, const ' o.real '* sl,' '\n',...
    o.indent.generic  'const ' o.real '* sl_sqr, const ' o.real '* gps' c_psi_dot_dec ', ' o.real '* du, ' o.real '* dsl, ' o.real '* muG){' '\n\n'])];

if ~isempty(o.Jac_n_struct)
    
    % determine required size of Dn
    size_Dn = 0;
    for jj=1:length(o.K_n)
        current_size_Dn = max(max(o.Jac_n_struct{o.K_n{jj}(1)}))*length(o.K_n{jj});
        size_Dn = size_Dn + current_size_Dn;
    end
    
    code = [code, sprintf([o.indent.generic  o.real ' Dn[%d];' '\n'],size_Dn)];
end

code = [code, sprintf([o.indent.generic  o.real ' v_Nnc[%d], tmp_nu[%d], tmp_nc_m[%d], tmp_contr = 0.0;' '\n'],...
    sum(o.nc), nu, max(o.nc))];

if ~isempty(o.K_lb)
    code = [code, sprintf([o.indent.generic  o.real ' temp_lb[%i];' '\n'],o.nu)];
end
if ~isempty(o.K_ub)
    code = [code, sprintf([o.indent.generic  o.real ' temp_ub[%i];' '\n'],o.nu)];
end
if ~isempty(o.K_n)
    code = [code, sprintf([o.indent.generic  o.real ' temp_n[%i], temp_n2[%i];' '\n'], max(o.nn), o.nu)];
end


Dn_need = 0;
Dn_need_vec = zeros(1,o.N);
for ii = 1:o.N
    
    code = [code, sprintf(['\n' o.indent.generic  '/* loop unrolling: step %i of %i */' '\n'],ii-1,o.N-1)]; %#ok
    lb_size = 0;
    ub_size = 0;
    if ~isempty(o.K_lb)
        code = [code, sprintf([o.indent.generic  'build_vNnc_lb(&gps[%i], &dot_J[%i], %i, &v_Nnc[%i]);' '\n'],...
            sum(o.nc(1:ii-1)), (ii-1)*o.nu, ii-1, sum(o.nc(1:ii-1)))]; %#ok
        lb_size = sum(~isinf(o.lb(:,ii)));
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf([o.indent.generic  'build_vNnc_ub(&gps[%i], &dot_J[%i], %i, &v_Nnc[%i]);' '\n'],...
            sum(o.nc(1:ii-1)) + lb_size, (ii-1)*o.nu, ii-1, sum(o.nc(1:ii-1)) + lb_size )]; %#ok
        ub_size = sum(~isinf(o.ub(:,ii)));
    end
    if ~isempty(o.K_n)
        code = [code, sprintf([o.indent.generic  'build_Dn(&u[%d], %i, &Dn[%d]);' '\n'], (ii-1)*nu, ii-1, Dn_need)]; %#ok
        
        code = [code, sprintf([o.indent.generic  'Dntop_times_dotJ_n(&Dn[%i], &dot_J[%i], %i, &temp_n[0]);' '\n'],...
            Dn_need, (ii-1)*o.nu, ii-1)]; %#ok
        
        code = [code, sprintf([o.indent.generic  'build_vNnc_n(&gps[%i], &temp_n[0], %i, &v_Nnc[%i]);' '\n'],...
            sum(o.nc(1:ii-1)) + lb_size + ub_size, ii-1, sum(o.nc(1:ii-1)) + lb_size + ub_size )]; %#ok
        
        
        Dn_need = Dn_need + max(max(o.Jac_n_struct_hor{ii}));
        Dn_need_vec(ii) = Dn_need;
        
    end
end
    

if (o.contractive || o.terminal)
    if o.variable_stepSize.active
        code = [code, sprintf([o.indent.generic 'dot_product_Nnu(&tmp_contr,dot_psi_N,dot_J);' '\n',...
            o.indent.generic 'v_Nnc[%d] = alpha_inverse * gps[%d] - tmp_contr;' '\n'],sum(o.nc)-1, sum(o.nc)-1)];
    else
        code = [code, sprintf([o.indent.generic 'dot_product_Nnu(&tmp_contr,dot_psi_N,dot_J);' '\n',...
            o.indent.generic 'v_Nnc[%d] = ' falcopt.internal.num2str(1, o.precision) '/' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) '*gps[%d] - tmp_contr;' '\n'],sum(o.nc)-1, sum(o.nc)-1)];
    end
    
end
   
code = [code, sprintf(['\n' o.indent.generic 'solveConstraintSystem(&muG[0], '])];
if ~isempty(o.K_n)
    code = [code, sprintf('&Dn[0], ')];
end
if (o.contractive || o.terminal)
    code = [code, sprintf('&dot_psi_N[0], ')];
end
code = [code, sprintf(['&v_Nnc[0], &sl_sqr[0]);' '\n'])];

for ii = 1:o.N
    code = [code, sprintf(['\n' o.indent.generic  '/* loop unrolling: step %i of %i */' '\n'],ii-1,o.N-1)]; %#ok
    lb_size = 0;
    ub_size = 0;
    if ~isempty(o.K_lb)
        code = [code, sprintf([o.indent.generic  'minus_Ina_muG(&muG[%i], %i, &temp_lb[0]);' '\n'],...
            sum(o.nc(1:ii-1)), ii-1)]; %#ok
        lb_size = sum(~isinf(o.lb(:,ii)));
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf([o.indent.generic  'Inb_muG(&muG[%i], %i, &temp_ub[0]);' '\n'],...
            sum(o.nc(1:ii-1)) + lb_size, ii-1)]; %#ok
        ub_size = sum(~isinf(o.ub(:,ii)));
    end
    if ~isempty(o.K_n)
        code = [code, sprintf([o.indent.generic  'Dn_times_muG_n(&Dn[%i], &muG[%i], %i, &temp_n2[0]);' '\n'],...
            Dn_need_vec(ii) - Dn_need_vec(1), sum(o.nc(1:ii-1)) + lb_size + ub_size, ii-1)]; %#ok
    end
    
    code = [code, sprintf([o.indent.generic  'sum_nr_constr(&tmp_nu[0]'])]; %#ok
    
    if ~isempty(o.K_lb)
        code = [code, ', &temp_lb[0]']; %#ok
    end
    if ~isempty(o.K_ub)
        code = [code, ', &temp_ub[0]']; %#ok
    end
    if ~isempty(o.K_n)
        code = [code, ', &temp_n2[0]']; %#ok
    end
    code = [code, sprintf([', &dot_J[%i]);' '\n'], nu*(ii-1))]; %#ok
    
    if (o.contractive || o.terminal)
        code = [code, sprintf([o.indent.generic  'sum_terminal(&tmp_nu[0], &dot_psi_N[%i], &muG[%i]);'...
            '\n'], nu*(ii-1), sum(o.nc) - 1)]; %#ok
    end
    if o.variable_stepSize.active
         code = [code, sprintf([o.indent.generic 'scale_nu(&du[%i], &tmp_nu[0], &alpha);'...
            '\n'],(ii-1)*nu)];%#ok
    else
        code = [code, sprintf([o.indent.generic 'scale_nu(&du[%i], &tmp_nu[0]);'...
            '\n'],(ii-1)*nu)];%#ok
    end   
    code = [code, sprintf([o.indent.generic 'product_matlab_nc(&sl[%d], &muG[%d], %i, &tmp_nc_m[0]);'...
        '\n'],sum(o.nc(1:ii-1)), sum(o.nc(1:ii-1)), ii-1 )]; %#ok
    code = [code, sprintf([o.indent.generic  'minus_scale_nc(&tmp_nc_m[0], %i, &dsl[%d]);'...
        '\n'], ii-1, sum(o.nc(1:ii-1)) )]; %#ok
end
if (o.contractive || o.terminal)
    if o.variable_stepSize.active
        code = [code, sprintf([o.indent.generic 'dsl[%d] = -alpha * sl[%d] * muG[%d];' '\n'], sum(o.nc) - 1, sum(o.nc) - 1, sum(o.nc) - 1)];
    else
        code = [code, sprintf([o.indent.generic 'dsl[%d] = -' falcopt.internal.toDefineName(o.names.lineSearchAlphaMax) '*sl[%d]*muG[%d];' '\n'], sum(o.nc) - 1, sum(o.nc) - 1, sum(o.nc) - 1)];
    end
end
code = [code, sprintf(['}' '\n\n'])];   
info.flops = falcopt.internal.multFlops( info.flops, N);

end

function [code, data, info] = generate_build_gps(o)
code = [];
data = [];

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
    o.indent.generic  'unsigned int jj;' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code,  sprintf([o.indent.generic  'for (jj=0;jj< na;jj++){' '\n',...
        o.indent.generic o.indent.generic 'gps[jj]= amu[jj] + 0.5*sl_sqr[jj];' '\n',...
        o.indent.generic  '}' '\n'])];
end
if ~isempty(o.K_ub)
    code = [code,  sprintf([o.indent.generic  'for (jj=na;jj< nb;jj++){' '\n'...
        o.indent.generic o.indent.generic 'gps[jj]= umb[jj-na] + 0.5*sl_sqr[jj];' '\n',...
        o.indent.generic  '}' '\n'])];
end
if ~isempty(o.K_n)
    code = [code,  sprintf([o.indent.generic  'for (jj= nb;jj< nc;jj++) {' '\n'...
        o.indent.generic o.indent.generic 'gps[jj]= n[jj - nb] + 0.5*sl_sqr[jj];' '\n',...
        o.indent.generic  '}' '\n'])];
end
code = [code,  sprintf(['}' '\n\n'])];

code = [code, sprintf(['\n' '/* It computes gps = g + 0.5 * sl_sqr */' '\n'])];
code = [code, sprintf([o.inline ' void build_gpsl(const ' o.real '* u' c_psi_dec ...
    c_contr_dec ', const ' o.real '* sl_sqr, ' o.real '* gps){' '\n\n'])];

if ~isempty(o.K_lb)
    code = [code, sprintf([o.indent.generic   o.real ' amu[%d];' '\n'], max(sum(~isinf( o.lb))) )];
end
if ~isempty(o.K_ub)
    code = [code, sprintf([o.indent.generic   o.real ' umb[%d];' '\n'], max(sum(~isinf( o.ub))) )];
end
if ~isempty(o.K_n)
    code = [code, sprintf([o.indent.generic   o.real ' n[%d];' '\n'], max(o.nn) )];
end

if (o.contractive || o.terminal)
    code = [code, sprintf([o.indent.generic  o.real ' g_contr = 0.0;' '\n'])];
end


info.flops.mul = info.flops.mul+ N*3; % mul inside only first cycle % NOT CLEAR, ToDo


lb = sum(~isinf(o.lb),1);
ub = sum(~isinf(o.ub),1);

for k = 1:o.N
    code = [code, sprintf(['\n',...
        o.indent.generic  '/* Unrolling the for loop: iteration %i of %i */' '\n'], k-1, o.N -1)]; %#ok
    if ~isempty(o.K_lb)
        code = [code, sprintf([o.indent.generic  'build_amu(&u[%i],%i,&amu[0]);' '\n'], (k-1)*o.nu, k-1)]; %#ok
    end
    if ~isempty(o.K_ub)
        code = [code, sprintf([o.indent.generic  'build_umb(&u[%i],%i,&umb[0]);' '\n'], (k-1)*o.nu, k-1)]; %#ok
    end
    if ~isempty(o.K_n)
        code = [code, sprintf([o.indent.generic  'build_n(&u[%i],%i,&n[0]);' '\n'], (k-1)*o.nu, k-1)]; %#ok
    end
    
    code = [code, sprintf([o.indent.generic  'build_gpsl_lowLevel('])]; %#ok
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

if (o.contractive || o.terminal)
    code = [code, sprintf([o.indent.generic  'g_contr = *psi_N - c_contr;' '\n'...
        o.indent.generic  'gps[%d] = g_contr + 0.5*sl_sqr[%d];' '\n'],...
        sum(o.nc)-1, sum(o.nc)-1)];
    info.flops.add = info.flops.add + 1;
end


code = [code, sprintf(['}' '\n'])];


info.flops.mul = info.flops.mul+ 3*N;
% The flops of "build_g" are counted outside this function (in generateConverter)

%% build_sqr

code = [code, sprintf(['\n' '/* It computes x_sqr = x.*x */' '\n'])];
code = [code, sprintf([o.inline ' void build_sqr_Nnc( const ' o.real ' *x, ' o.real ' *x_sqr){' '\n\n'])];

code = [code, sprintf([o.indent.generic  'unsigned int ii = 0;' '\n\n'])];

code = [code, sprintf([o.indent.generic  'for (ii=0;ii<%d;ii++)' '\n',...
    o.indent.generic o.indent.generic 'x_sqr[ii] = x[ii]*x[ii];' '\n',...
    '}' '\n\n'],sum(o.nc))];

info.flops.mul = info.flops.mul+ sum(o.nc);
end

function [code, data, info] = generate_dot_product_Nnu(o)
% Damian please use your function to generate products

code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

code = [code, sprintf(['\n' '/* dot_product of dim N*nu */' '\n'])];
[d, c, in] = falcopt.generateMVMult(ones(1,o.N*o.nu), ...
    'names', struct('fun', 'dot_product_Nnu', 'M', 'u',...
    'v', 'du'), 'static', false, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end

function [code, data, info] = generate_dot_product_Nnc(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

code = [code, sprintf(['\n' '/* dot_product of dim N*nc */' '\n'])];
[d, c, in] = falcopt.generateMVMult(ones(1,sum(o.nc)), ...
    'names', struct('fun', 'dot_product_Nnc', 'M', 'u',...
    'v', 'du'), 'static', false, 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);
end

function [code, data] = generate_copy(o)
% copy_nx

nx = o.nx;
nu = o.nu;
N = o.N;

code = [];
data = [];

code = [code, sprintf(['\n' '/* it copies the content of a variable */' '\n'])];
[d, c] = falcopt.generateMVMult({eye(nx)}, ...
    'names', struct('fun', 'copy_nx', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

code = [code, sprintf(['\n' '/* it copies the content of a variable */' '\n'])];
[d, c] = falcopt.generateMVMult({eye(N*nx)}, ...
    'names', struct('fun', 'copy_Nnx', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

code = [code, sprintf(['\n' '/* it copies the content of a variable */' '\n'])];
[d, c] = falcopt.generateMVMult({eye(sum(o.nc))}, ...
    'names', struct('fun', 'copy_Nnc', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

code = [code, sprintf(['\n' '/* it copies the content of a variable */' '\n'])];
[d, c] = falcopt.generateMVMult({eye(N*nu)}, ...
    'names', struct('fun', 'copy_Nnu', 'M', {{'I'}}, 'v', {{'x'}}, 'r', 'res'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];

end

function [code, data, info, optCode] = generate_difference(o)
% generate  diff_Nnc

data = [];
code = [];

info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

code = [code, sprintf(['\n' '/* It computes dx = x - xref */' '\n'])];
[d, c, in] = falcopt.generateMVMult({eye(sum(o.nc)), -eye(sum(o.nc))}, ...
    'names', struct('fun', 'diff_Nnc', 'M', {{'I', 'mI'}},...
    'v', {{'x', 'xref'}}, 'r', 'dx'), 'types', o.real, 'precision', o.precision, 'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);

if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
info.flops = falcopt.internal.addFlops(info.flops, in.flops);

optCode = [];

end


function [code, data, info] = generate_det_phi(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

code = [code, sprintf(['\n' '/* It computes the merit function phi */' '\n'])];
code = [code, sprintf([o.inline ' void det_phi (const ' o.real ' J, const ' o.real '* gps, const ',...
    o.real '* mu, const ' o.real ' rho, ' o.real '* phi){' '\n\n'])];

code = [code, sprintf([o.indent.generic  o.real ' tmp[%d], pr = 0.0;' '\n'],sum(o.nc))];
code = [code, sprintf([o.indent.generic  'unsigned int ii = 0;' '\n\n'])];

code = [code, sprintf([o.indent.generic  'for (ii=%d; ii--; )' '\n'], sum(o.nc))];
code = [code, sprintf([o.indent.generic o.indent.generic 'tmp[ii] = mu[ii] + 0.50*rho*gps[ii];' '\n'], sum(o.nc))];
info.flops.mul = info.flops.mul+ 2*sum(o.nc);
info.flops.add = info.flops.add+ sum(o.nc);

code = [code, sprintf([o.indent.generic  'dot_product_Nnc(&pr, gps,tmp);' '\n\n'])];
info.flops.add = info.flops.add+ (sum(o.nc)-1); % addition of dot_product_Nnc
info.flops.mul = info.flops.mul+ sum(o.nc); % multiplication of dot_product_Nnc

code = [code, sprintf([o.indent.generic  '(*phi) = J + pr;' '\n'])];
info.flops.add = info.flops.add+ sum(o.nc);
code = [code, sprintf(['}' '\n\n'])];

end

function [code, data, info] = generate_det_dot_phi(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

code = [code, sprintf(['\n' '/* It computes the derivative of the merit function dot_phi */' '\n'])];
code = [code, sprintf([o.inline ' void det_dot_phi (const ' o.real '* du, const ' o.real '* DJ, const ' ...
    o.real ' rho, const '  o.real '* gps, ' '\n',...
    o.indent.generic ' const ' o.real '* mu, const ' o.real '* dm, ' o.real '* dot_phi){' '\n\n'])];

code = [code, sprintf([o.indent.generic  o.real ' tmp_prod[%d], prod_1 = 0.0, prod_2 = 0.0;' '\n'],sum(o.nc))];
code = [code, sprintf([o.indent.generic  'unsigned int ii = 0;' '\n\n'])];

code = [code, sprintf([o.indent.generic  'for (ii=%d; ii--; )' '\n'], sum(o.nc))];
code = [code, sprintf([o.indent.generic o.indent.generic 'tmp_prod[ii] = mu[ii] - dm[ii] + rho*gps[ii];' '\n'], sum(o.nc))];
info.flops.mul = info.flops.mul+ sum(o.nc);
info.flops.add = info.flops.add+ 2*sum(o.nc);

code = [code, sprintf([o.indent.generic  'dot_product_Nnu(&prod_1, du, DJ);' '\n\n'])];
code = [code, sprintf([o.indent.generic  'dot_product_Nnc(&prod_2, gps,tmp_prod);' '\n\n'])];

info.flops.add = info.flops.add+ o.N*(o.nu + sum(o.nc) - 1); % additions of dot_product_Nnc and _Nnu
info.flops.mul = info.flops.mul+ o.N*(o.nu + sum(o.nc)); % multiplications of dot_product_Nnc and _Nnu

code = [code, sprintf([o.indent.generic  '(*dot_phi) = prod_1 - prod_2;' '\n'])];

info.flops.add = info.flops.add+ 1;
code = [code, sprintf(['}' '\n\n'])];
end

function [ code, data, info] = generate_conditions_rho(o)
code = [];
data = [];
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

code = [code, sprintf(['\n' '/* It checks the decrease condition of the merit function */' '\n'])];
code = [code, sprintf([o.inline ' int conditions_rho_PM_simpler (const ' o.real ' dot_phi, const ' o.real ' du_sqr, const ' ...
    o.real ' dsl_sqr, const ' o.real ' alpha){' '\n\n'])];
code = [code, sprintf([o.indent.generic  'unsigned int res = 2;' '\n\n'])];

code = [code, sprintf([o.indent.generic 'if (dot_phi <= (-' falcopt.internal.num2str(0.5, o.precision) '/alpha)*(du_sqr + dsl_sqr))' '\n'])];
code = [code, sprintf([o.indent.generic o.indent.generic 'res = 1;' '\n',...
    o.indent.generic 'else' '\n',...
    o.indent.generic o.indent.generic 'res = 0;' '\n\n',...
    o.indent.generic 'return res;' '\n'...
    '}' '\n\n' ])];
info.flops.add = info.flops.add+ 1;
info.flops.comp = info.flops.comp +1;


end

function [code, data, info] = generate_Lagrangian_oracles_lp(o)

code = [];
data = [];
info.flops.it = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
info.flops.ls = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

%% define one norm and inf norm of dimension N*nc

if o.merit_function ~= 2
    
    code = [code, sprintf(['\n' '/* one_norm of a vector  */' '\n'])];
    code = [code, sprintf([o.inline ' ' o.real ' one_norm (const ' o.real '* g){' '\n\n',...
        o.indent.generic  'int ii = 0;' '\n'...
        o.indent.generic  o.real ' norm = 0.0;' '\n\n',...
        o.indent.generic  'for (ii = %d; ii-- >0; ){' '\n'],sum(o.nc))];
    
    code = [code, sprintf([o.indent.generic o.indent.generic 'norm += ' o.auxFunctions.abs '(g[ii]);' '\n',...
        o.indent.generic  '}' '\n',...
        o.indent.generic  'return norm;' '\n',...
        '}' '\n\n'])];
    
    code = [code, sprintf(['\n' '/* inf_norm of a vector  */' '\n'])];
    code = [code, sprintf([o.inline ' ' o.real ' inf_norm (const ' o.real '* g){' '\n\n',...
        o.indent.generic  'int ii = 0;' '\n'...
        o.indent.generic  o.real ' norm = 0.0;' '\n\n',...
        o.indent.generic  'for (ii = %d; ii-- >0; ){' '\n'],sum(o.nc))];
    
    code = [code, sprintf([o.indent.generic o.indent.generic 'norm = ' o.auxFunctions.max '(norm,' o.auxFunctions.abs '(g[ii]));' '\n',...
        o.indent.generic  '}' '\n',...
        o.indent.generic  'return norm;' '\n',...
        '}' '\n\n'])];
    
else
    code = [code, sprintf(['\n' '/* two_norm of a vector  */' '\n'])];
    code = [code, sprintf([o.inline ' ' o.real ' two_norm (const ' o.real '* g){' '\n\n',...
        o.indent.generic  'int ii = 0;' '\n'...
        o.indent.generic  o.real ' norm = 0.0, norm2 = 0.0;' '\n\n',...
        o.indent.generic  'for (ii = %d; ii-- >0; ){' '\n'],sum(o.nc))];
    
    code = [code, sprintf([o.indent.generic o.indent.generic 'norm2 += g[ii]*g[ii];' '\n',...
        o.indent.generic  '}' '\n',...
        o.indent.generic  'norm = sqrt(norm2);' '\n',...
        o.indent.generic  'return norm;' '\n',...
        '}' '\n\n'])];
    
end

%% det_phi
code = [code, sprintf(['\n' '/* it computes the merit function phi */' '\n'])];
code = [code, sprintf([o.inline ' void det_phi (const ' o.real ' J, const ' o.real '* gps, const ',...
    o.real ' rho, ' o.real '* phi){' '\n\n'])];

if o.merit_function == 1
    code = [code, sprintf([o.indent.generic  '(*phi) = J + rho * one_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
elseif o.merit_function == Inf
    code = [code, sprintf([o.indent.generic  '(*phi) = J + rho * inf_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
elseif o.merit_function == 2
    code = [code, sprintf([o.indent.generic  '(*phi) = J + rho * two_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
end

code = [code, sprintf(['}' '\n\n'])];
info.flops.ls = info.flops.it; %det_phi called in ls too.

%% det_dot_phi

code = [code, sprintf(['\n' '/* it computes the derivative of the merit function dot_phi */' '\n'])];
code = [code, sprintf([o.inline ' void det_dot_phi (const ' o.real '* du, const ' o.real '* DJ, const ' ...
    o.real ' rho, const '  o.real '* gps, ' o.real '* dot_phi){' '\n\n'])];

code = [code, sprintf([o.indent.generic  o.real ' pr = 0.0;' '\n'])];

code = [code, sprintf([o.indent.generic  'dot_product_Nnu(&pr, du, DJ);' '\n\n'])];
info.flops.it.add = info.flops.it.add+ o.N*(o.nu-1); % addition of dot_product_Nnu
info.flops.it.mul = info.flops.it.mul+ o.N * o.nu; % multiplication of dot_product_Nnu

if o.merit_function == 1
    code = [code, sprintf([o.indent.generic  '(*dot_phi) = pr - rho*one_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
elseif o.merit_function == Inf
    code = [code, sprintf([o.indent.generic  '(*dot_phi) = pr - rho*inf_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
elseif o.merit_function == 2
    code = [code, sprintf([o.indent.generic  '(*dot_phi) = pr - rho*two_norm(gps);' '\n'])];
    info.flops.it.add = info.flops.it.add+ 1;
    info.flops.it.mul = info.flops.it.mul+ 1;
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
end

code = [code, sprintf(['}' '\n\n'])];

%% conditions_rho

code = [code, sprintf(['\n' '/* it updates the penalty parameter rho */' '\n'])];
code = [code, sprintf([o.inline ' void update_rho (const ' o.real '* muG, ' o.real '* rho, ' o.real '* rho_hat){' '\n\n'])];

if o.merit_function == 1
    code = [code, sprintf([o.indent.generic  '(*rho_hat) = inf_norm(muG);' '\n'])];
    info.flops.it.comp = info.flops.it.comp+ 2*sum(o.nc); %flops of abs and max in inf_norm
elseif o.merit_function == Inf
    code = [code, sprintf([o.indent.generic  '(*rho_hat) = one_norm(muG);' '\n'])];
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
elseif o.merit_function == 2
    code = [code, sprintf([o.indent.generic  '(*rho_hat) = two_norm(muG);' '\n'])];
    info.flops.it.add = info.flops.it.add+ sum(o.nc); % flops of one_norm
    info.flops.it.comp = info.flops.it.comp+ sum(o.nc); %flops of abs in one_norm
end

code = [code, sprintf([o.indent.generic  '(*rho) = ' o.auxFunctions.max '( (*rho), (*rho_hat)); ' '\n'])];
info.flops.it.comp = info.flops.it.comp+ 1;

code = [code, sprintf(['}' '\n\n' ])];

end

function [code, data, info] = generate_weighted_sum_nc(o)
code = '';
data = '';
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);

code = [code, sprintf(['\n' '/*  dx = t*x + xref */' '\n'])];
[d, c, in] = falcopt.generateMVMult({eye(sum(o.nc)), eye(sum(o.nc))}, ...
    'names', struct('fun', 'weighted_sum_Nnc', 'M',...
    {'I'}, 'v', {{'x', 'xref'}}, 'r', 'dx'),...
    'static', struct('M',[false,true]),...
    'types', o.real, 'precision', o.precision, 'structure', 'unique', ...
    'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
info.flops = falcopt.internal.addFlops(info.flops, in.flops);
if ~isempty(d)
    data = [data, d, sprintf('\n')];
end
code = [code, c, sprintf('\n\n')];
end

function [code, data, info] = generate_weighted_sum_nu(o)
code = '';
data = '';
info.flops = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
code = [code, sprintf(['\n' '/*  du = t*u + uref */' '\n'])];
[d, c, in] = falcopt.generateMVMult({eye(o.N*o.nu), eye(o.N*o.nu)}, ...
    'names', struct('fun', 'weighted_sum_Nnu', 'M',...
    {'I'}, 'v', {{'x', 'xref'}}, 'r', 'dx'),...
    'static', struct('M',[false,true]),...
    'types', o.real, 'precision', o.precision, 'structure', 'unique', ...
    'verbose', o.verbose, 'test', o.test, 'inline', o.inline, 'indent', o.indent);
info.flops = falcopt.internal.addFlops(info.flops, in.flops);
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
code = [code, sprintf(['\n' '/* it generates a safeguarded quadratic interpolation */' '\n'])];
code = [code, sprintf([o.indent.code o.inline ' ' o.real ' quadratic_interp (const ' o.real ' f_l, const ', o.real ' g_l, const ' o.real ' t_u, const ' o.real ' f_u) {' '\n' '\n'])];
code = [code, sprintf([o.indent.code o.indent.generic o.real ' t_theo;' '\n' '\n'])];
code = [code, sprintf([o.indent.code o.indent.generic 't_theo = -0.5*(g_l*t_u*t_u)/(f_u - f_l - g_l*t_u);' '\n'])];
info.flops.mul = info.flops.mul+5;
info.flops.add = info.flops.add+2;
info.flops.div = info.flops.div +1;
% code = [code, sprintf([o.indent.code o.indent.generic 'a = 0.1*t_u;'
% '\n'])]; % TBD
% info.flops.mul = info.flops.mul+1;
code = [code, sprintf([o.indent.code o.indent.generic 'return ' o.auxFunctions.max '(' o.auxFunctions.min '(t_theo,' falcopt.internal.num2str(p.tau(2), o.precision) '*t_u), ' falcopt.internal.num2str(p.tau(1), o.precision) '*t_u);' '\n'])];
info.flops.mul = info.flops.mul+2;
info.flops.add = info.flops.add+1;
info.flops.comp = info.flops.comp+ 3;
code = [code, sprintf([o.indent.code '}'])];
end

function [code, data, info] = generate_compute_max(o)
code = '';
data = '';
code = [code, sprintf(['\n' '/*  compute maximum of a vector */' '\n'])];
code = [code, sprintf([o.indent.code o.inline ' ' o.real ' compute_max_Nnc( const ' o.real '* x) {' '\n' ...
    o.indent.code o.indent.generic 'unsigned int ii = 0;' '\n',...
    o.indent.code o.indent.generic o.real ' m = -100.0;' '\n' '\n'])];
code = [code, sprintf([o.indent.code o.indent.generic 'for (ii=%d; ii--; ) {' '\n'],sum(o.nc))];
code = [code, sprintf([o.indent.code o.indent.generic o.indent.generic 'm = ' o.auxFunctions.max '(m,x[ii]);' '\n' ...
    o.indent.code o.indent.generic '}' '\n' ...
    o.indent.code o.indent.generic 'return m;' '\n' ...
    o.indent.code '}'])];
info.flops.comp = sum(o.nc);
end

function [code, data, info] = general_objective_gradient_oracle(o)
nx = o.nx;
nu = o.nu;
nw = o.nw;
N = o.N;
trackRef = o.objective.trackReference;

code = [];
data = [];
info.flops.it = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
info.flops.ls = struct('add', 0, 'mul', 0, 'inv', 0, 'sqrt', 0, 'comp', 0);
info.src = {};
info.header = {};

%it generates: dot_product_nx_nx, dot_product_nu_nu,
%               product_and_sum_nu
[c, d, in] = generate_auxiliary_functions(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.internal.multFlops(in.flops, N);
info.flops.it = falcopt.internal.addFlops(info.flops.it,in.flops);
info.flops.ls = falcopt.internal.addFlops(info.flops.ls,in.flops);


[c, d, in] = generate_product_and_sum_nx(o);
code = [code, c];
data = [data, d];
in.flops = falcopt.internal.multFlops(in.flops, N-1);
info.flops.it = falcopt.internal.addFlops(info.flops.it,in.flops); % flops product_and_sum_nx in det_J_dot_J


if trackRef
    % generate diffX / diffU
    [c, d, in] = generate_diffXU(o);
    code = [code, c];
    data = [data, d];
    in.flops = falcopt.internal.multFlops(in.flops, N);
    info.flops.it = falcopt.internal.addFlops(info.flops.it,in.flops);
    info.flops.ls = falcopt.internal.addFlops(info.flops.ls,in.flops);
end

[c,d] = generate_copy(o);
% generate copy_Nnc copy_Nnu copy_Nnx
code = [code, c];
data = [data, d];

c_w_dec = argument_w(o,true);
c_tr_dec = argument_def(o,true);
c_psi_dec = argument_def_internal_psi_noconst(o,true);
c_psi_dot_dec = argument_def_internal_psi_dot_noconst(o,true);

% generate cost and cost_N functions + Jacobians
switch o.gradients
    case 'casadi'
        sxfcn = {}; %#ok
        
        % generate cost, cost_x, cost_u
        [d,c,i] = casadi_jacobians(o,o.objective.nonlinear,'cost','staticName','cost','jac',{'x','u'},...
                            'fileName','casadi_cost','jac_x',{'cost_x_data','in_cost_x','cost_x'},...
                            'jac_u',{'cost_u_data','in_cost_u','cost_u'},'generate_code',false, 'structure','dense');          
        code = [code, c];
        data = [data, d];
        cost_static = i.y.static; %#ok
        cost_x_static = i.in_cost_x.static;
        cost_u_static = i.in_cost_u.static;
        sxfcn = i.sxfcn;
        
        %flops
        info.flops.ls.add = info.flops.ls.add  + i.y.flops*(o.N-1);              %flops cost in det_J
        info.flops.it.add = info.flops.it.add  + i.y.flops*(o.N-1);              %flops cost in det_J_dot_J
        info.flops.it.add = info.flops.it.add  + i.in_cost_x.flops*(o.N-2);      %flops cost_x in det_J_dot_J
        info.flops.it.add = info.flops.it.add  + i.in_cost_u.flops*(o.N-1);      %flops cost_u in det_J_dot_J
        
        % generate cost_N and cost_N_x
        if ~isfield(o.objective,'nonlinearN')
            cost_N_static = 1;
            cost_N_x_static = 1;
            data = [ data, sprintf('static const cost_N_x_data[1] = {0.0};\n')];
        else
            [d,c,i] = casadi_jacobians(o,o.objective.nonlinearN,'cost_N','staticName','cost_N','jac',{'x'},...
                                'jac_x',{'cost_N_x_data','in_cost_N_x','cost_N_x'},...
                                'fileName','casadi_cost_N','generate_code',false,'structure','dense');   

            code = [code, c];
            data = [data, d];
            cost_N_static = i.y.static;
            cost_N_x_static = i.in_cost_N_x.static;
            sxfcn = [sxfcn, i.sxfcn];
            
            %flops
            info.flops.ls.add  = info.flops.ls.add  + i.y.flops;                  % flops cost_N in det_J
            info.flops.ls.add  = info.flops.ls.add  + i.in_cost_N_x.flops;        % flops cost_N_x in det_J
            info.flops.it.add  = info.flops.it.add  + i.y.flops;                  % flops cost_N in det_J_dot_J
            info.flops.it.add  = info.flops.it.add  + i.in_cost_N_x.flops;        % flops cost_N_x in det_J_dot_J
 
        end
        [info.src,info.header] = generate_casadi_c(o, 'casadi_cost', sxfcn);
        
    case {'matlab'}
        % generate cost, cost_x, cost_u
        [d,c,i] = matlab_jacobians(o,o.objective.nonlinear,'cost','staticName','cost','jac',{'x','u'},...
                            'fileName','casadi_cost','jac_x',{'cost_x_data','in_cost_x','cost_x'},...
                            'jac_u',{'cost_u_data','in_cost_u','cost_u'},'generate_code',false,'structure','dense');          
        code = [code, c];
        data = [data, d];
        cost_static = i.y.static; %#ok
        cost_x_static = i.in_cost_x.static;
        cost_u_static = i.in_cost_u.static;
        
        % flops
        info.flops.ls = falcopt.internal.addFlops(info.flops.ls,falcopt.internal.multFlops(i.y.flops,o.N-1));          % flops cost in det_J
        info.flops.it = falcopt.internal.addFlops(info.flops.it,falcopt.internal.multFlops(i.y.flops,o.N-1));          % flops cost in det_J_dot_J
        info.flops.it = falcopt.internal.addFlops(info.flops.it,falcopt.internal.multFlops(i.in_cost_x.flops,o.N-2));  % flops cost_x in det_J_dot_J 
        info.flops.it = falcopt.internal.addFlops(info.flops.it,falcopt.internal.multFlops(i.in_cost_u.flops,o.N-1));  % flops cost_u in det_J_dot_J 
        
        % generate cost_N and cost_N_x
        if ~isfield(o.objective,'nonlinearN')
            cost_N_static = 1;
            cost_N_x_static = 1;
            data = [ data, sprintf('static const cost_N_x_data[1] = {0.0};\n')];
        else
            [d,c,i] = matlab_jacobians(o,o.objective.nonlinearN,'cost_N','staticName','cost_N','jac',{'x'},...
                                'jac_x',{'cost_N_x_data','in_cost_N_x','cost_N_x'},...
                                'fileName','casadi_cost_N','generate_code',false,'structure','dense'); 
                            code = [code, c];
                            data = [data, d];
                            cost_N_static = i.y.static;
                            cost_N_x_static = i.in_cost_N_x.static;
                            
             % flops
             info.flops.ls = falcopt.internal.addFlops(info.flops.ls,i.y.flops);            % flops cost_N in det_J
             info.flops.it = falcopt.internal.addFlops(info.flops.it,i.y.flops);            % flops cost_N in det_J_dot_J
             info.flops.it = falcopt.internal.addFlops(info.flops.it,i.in_cost_N_x.flops);  % flops cost_N_x in det_J_dot_J
        end
    otherwise
        error(['nonlinear objective function with options ''manual'' or ''ccode'' not implemented yet. Use quadratic cost function with objective.Q '...
               'objective.R and objective.P instead.']);
end

%% generate det_J

code = [code, sprintf([o.inline ' void det_J(const ' o.real '* x0, const ' o.real ...
    '* u, const ' o.real '* x' c_tr_dec ', ' ...
    o.real '* J' c_psi_dec '){' '\n\n',...
    ])];
code = [code, sprintf([o.indent.generic o.real ' tmp_N,tmp_cost;' '\n'])];

if trackRef
    code = [code, sprintf([o.indent.generic o.real ' dx[%d], du[%d];' '\n'],...
        nx,nu)];
end
if o.contractive
    code = [code, sprintf([o.indent.generic o.real ' tmp_contr = 0.0;' '\n'], nx)];
    if trackRef
        code = [code, sprintf([o.indent.generic o.real ' dx_contr[%d];' '\n'], nx)];
    end
end
code = [code, sprintf([o.indent.generic 'unsigned int ii = 0;' '\n\n'])];

if trackRef
    code = [code, sprintf([o.indent.generic 'diffX(dx, x + %d, xref + %d);' '\n'],(N-1)*nx, (N-1)*nx)];
    if ~cost_N_static
        code = [code, sprintf([o.indent.generic 'cost_N( dx, &tmp_N);' '\n'])];
    end
    if o.contractive
        code = [code, sprintf([o.indent.generic 'diffX(dx_contr, x + (ind - 1)*%d, xref + (ind-1)*%d);' '\n'],nx, nx)];
        if ~cost_N_static
            code = [code, sprintf([o.indent.generic 'cost_N( dx_contr, &tmp_contr);' '\n'])];
        end
        code = [code, sprintf([o.indent.generic '(*psi_N) = tmp_contr;' '\n'])];
    end
else
    if ~cost_N_static
        code = [code, sprintf([o.indent.generic 'cost_N( x + %d, &tmp_N);' '\n'],(N-1)*nx)];
    end
    if o.contractive
        if ~cost_N_static
            code = [code, sprintf([o.indent.generic 'cost_N( x + (ind - 1)*%d, &tmp_contr);' '\n'],nx)];
        end
        code = [code, sprintf([o.indent.generic '(*psi_N) = tmp_contr;' '\n'])];
    end
end

code = [code, sprintf([o.indent.generic '(*J) = tmp_N;' '\n'])];
if o.terminal
    code = [code, sprintf([o.indent.generic '(*psi_N) = tmp_N;' '\n'])];
end

code = [code, sprintf(['\n' o.indent.generic 'for (ii=%d; ii-->0; ) {' '\n\n'],N-1)];


if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffX(dx, x + ii*%d, xref + ii*%d);' '\n'],nx, nx)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffU(du, u + (ii+1)*%d, uref + (ii+1)*%d);' '\n'],nu, nu)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'cost(dx, du, &tmp_cost);' '\n'])];
else
    code = [code, sprintf([o.indent.generic o.indent.generic 'cost(x + ii*%d, u + (ii+1)*%d, &tmp_cost);' '\n'],nx,nu)];
end

code = [code, sprintf([o.indent.generic o.indent.generic '(*J) += tmp_cost;' '\n'...
    o.indent.generic '}' '\n'])];
if trackRef
    code = [code, sprintf([o.indent.generic 'diffU(du, u , uref);' '\n'])];
    code = [code, sprintf([o.indent.generic 'cost(x0, du, &tmp_cost);' '\n'])];
    code = [code, sprintf([o.indent.generic '(*J) += tmp_cost;' '\n'])];
else
    code = [code, sprintf([o.indent.generic 'cost(x0, u, &tmp_cost);' '\n'])];
    code = [code, sprintf([o.indent.generic '(*J) += tmp_cost;' '\n'])];
end
code = [code, sprintf(['}' '\n\n'])];

%% generate det_J_and_dot_J

code = [code, sprintf([o.inline ' void det_J_and_dot_J(const ' o.real '* x0, const ' o.real ...
    '* u, const ' o.real '* x' c_w_dec c_tr_dec ', ' ...
    o.real '* J, ' o.real '* dot_J' c_psi_dec c_psi_dot_dec '){' '\n\n',...
    ])];

code = [code, sprintf([o.indent.generic o.real ' Px[%d], Ru[%d], Qx[%d], mem_tmp2[%d], ' '\n'...
    o.indent.generic o.indent.generic 'tmp_x = 0.0, tmp_u = 0.0;' '\n'],...
    nx,nu,nx,nx)];

if o.contractive
    code = [code, sprintf([o.indent.generic 'int index = 0;' '\n'])];
end

code = [code, sprintf([o.indent.generic o.real ' tmp_N, tmp_cost;' '\n'])];
if cost_N_x_static
    data_cost_N_x = 'cost_N_x_data';
else
    data_cost_N_x = 'Px_N';
    code = [code, sprintf([o.indent.generic o.real ' Px_N[%d];' '\n'],nx)];
end
if cost_x_static
    data_cost_x = 'cost_x_data';
else
    data_cost_x = 'Qx';
end
if cost_u_static
    data_cost_u = 'cost_u_data';
else
    data_cost_u = 'Ru';
end

if trackRef
    code = [code, sprintf([o.indent.generic o.real ' dx[%d], du[%d];' '\n'],...
        nx,nu)];
end
if o.contractive||o.terminal
    code = [code, sprintf([o.indent.generic o.real ' Px_contr[%d], mem_tmp_contr[%d], tmp_contr = 0.0;' '\n'], nx,nx)];
    if trackRef
        code = [code, sprintf([o.indent.generic o.real ' dx_contr[%d];' '\n'], nx)];
    end
end
code = [code, sprintf([o.indent.generic 'unsigned int ii = 0;' '\n\n'])];


if trackRef
    code = [code, sprintf([o.indent.generic 'diffX(dx, x + %d, xref + %d);' '\n'],(N-1)*nx, (N-1)*nx)];
    if ~cost_N_static
        code = [code, sprintf([o.indent.generic 'cost_N( dx, &tmp_N);' '\n'])];
    end
    if o.contractive
        code = [code, sprintf([o.indent.generic 'diffX(dx_contr, x + (ind - 1)*%d, xref + (ind-1)*%d);' '\n'],nx, nx)];
        if ~cost_N_static
            code = [code, sprintf([o.indent.generic 'cost_N( dx_contr, &tmp_contr);' '\n'])];
        end
        if ~cost_N_x_static
            code = [code, sprintf([o.indent.generic 'cost_N_x( dx_contr, Px_contr);' '\n'])];
        else
            code = [code, sprintf([o.indent.generic 'copy_nx( Px_contr, cost_N_x_data);' '\n'])];
        end
        code = [code, sprintf([o.indent.generic '(*psi_N) = tmp_contr;' '\n'])];
    end
else
    if ~cost_N_static
        code = [code, sprintf([o.indent.generic 'cost_N( x + %d, &tmp_N);' '\n'],(N-1)*nx)];
    end
    if o.contractive
        if ~cost_N_static
            code = [code, sprintf([o.indent.generic 'cost_N( x + (ind - 1)*%d, &tmp_contr);' '\n'],nx)];
        end
        if ~cost_N_x_static
            code = [code, sprintf([o.indent.generic 'cost_N_x( x + (ind - 1)*%d, Px_contr);' '\n'],nx)];
        else
            code = [code, sprintf([o.indent.generic 'copy_nx( Px_contr, cost_N_x_data);' '\n'])];
        end
        code = [code, sprintf([o.indent.generic '(*psi_N) = tmp_contr;' '\n'])];
    end
end

code = [code, sprintf([o.indent.generic '(*J) = tmp_N;' '\n'])];
if o.terminal
    code = [code, sprintf([o.indent.generic '(*psi_N) = tmp_N;' '\n'])];
end

if trackRef
    if ~cost_N_x_static
        code = [code, sprintf([o.indent.generic 'cost_N_x( dx, Px_N);' '\n'],(N-1)*nx)];
    end
else
    if ~cost_N_x_static
        code = [code, sprintf([o.indent.generic 'cost_N_x( x + %d, Px_N);' '\n'],(N-1)*nx)];
    end
end


if o.terminal
    code = [code, sprintf([o.indent.generic 'copy_nx(Px_contr,' data_cost_N_x ');' '\n'])]; 
end


if trackRef
    code = [code, sprintf([o.indent.generic 'diffU(du, u + %d, uref + %d);' '\n'],(N-1)*nu, (N-1)*nu)];
    code = [code, sprintf([o.indent.generic 'diffX(dx, x + %d, xref + %d);' '\n'],(N-2)*nx, (N-2)*nx)];
    if ~cost_u_static
        code = [code, sprintf([o.indent.generic 'cost_u(dx, du, Ru);' '\n'])];
    end
    code = [code, sprintf([o.indent.generic 'cost(dx, du, &tmp_cost);' '\n'])];
else
    if ~cost_u_static
        code = [code, sprintf([o.indent.generic 'cost_u(x + %d, u + %d, Ru);' '\n'],(N-2)*nx,(N-1)*nu)];
    end
    code = [code, sprintf([o.indent.generic 'cost(x + %d, u + %d, &tmp_cost);' '\n'],(N-2)*nx,(N-1)*nu)];
end
code = [code, sprintf([o.indent.generic '(*J) += tmp_cost;' '\n'])];

if o.nw > 0
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic 'Jacobian_u(x + %d,u + %d, w + %d, G);' '\n'],(N-2)*nx,(N-1)*nu,(N-1)*nw)];
    end
else
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic 'Jacobian_u(x + %d,u + %d, G);' '\n'],(N-2)*nx,(N-1)*nu)];
    end
end

code = [code, sprintf([o.indent.generic 'product_and_sum_nu(dot_J + %d,' data_cost_N_x ', ' data_cost_u ', G);' '\n\n'], (N-1)*nu)];

if o.terminal
    code = [code, sprintf([o.indent.generic 'product_contr_nu(&dot_psi_N[%d], ' data_cost_N_x ', G);' '\n\n'], (N-1)*nu)];
elseif o.contractive
    code = [code, sprintf([o.indent.generic 'if (ind == %d) ' '\n'...
        o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[%d], Px_contr, G);' '\n'], N, (N-1)*nu)];
    code = [code, sprintf([o.indent.generic 'else' '\n'...
        o.indent.generic o.indent.generic 'set_zero_nu(&dot_psi_N[%d]);' '\n\n'], (N-1)*nu)];
end



if o.contractive
    code = [code, sprintf([o.indent.generic 'if (ind == %d) ' '\n'...
        o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[%d], Px_contr, G);' '\n'], N, (N-1)*nu)];
    code = [code, sprintf([o.indent.generic 'else' '\n'...
        o.indent.generic o.indent.generic 'set_zero_nu(&dot_psi_N[%d]);' '\n\n'], (N-1)*nu)];
end

code = [code, sprintf([o.indent.generic 'copy_nx(Px,' data_cost_N_x ');' '\n'])]; 


code = [code, sprintf(['\n' o.indent.generic 'for (ii=%d; ii-->0; ) {' '\n\n'],N-1)];

if o.nw > 0
    if ~o.Jac_x_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'Jacobian_x(x + ii*%d,u + (ii+1)*%d, w + (ii+1)*%d, F);' '\n'],nx,nu,nw)];
        info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
        
    end
else
    if ~o.Jac_x_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'Jacobian_x(x + ii*%d,u + (ii+1)*%d, F);' '\n'],nx,nu)];
        info.flops.it.mul = info.flops.it.mul+ 2*(N-1);
    end
end

code = [code, sprintf([o.indent.generic o.indent.generic 'if (ii==0){' '\n'])];
if o.nw>0
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x0,u,w,G);' '\n'])];
    end
else
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x0,u,G);' '\n'])];
    end
end
if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'diffU(du, u + ii*%d, uref + ii*%d);' '\n'],nu, nu)]; 
    code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost(x0, du, &tmp_cost);' '\n'])];
    if ~cost_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost_u(x0, du, Ru);' '\n'])];
    end
else
    code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost(x0, u , &tmp_cost);' '\n'])];
    if ~cost_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost_u( x0, u , Ru);' '\n'])];
    end
end
code = [code, sprintf([o.indent.generic o.indent.generic '}\n'...
                o.indent.generic o.indent.generic 'else{' '\n'])];
if o.nw>0
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x + (ii-1)*%d, u + ii*%d, w + ii*%d, G);' '\n'],nx,nu,nw)];
    end
else
    if ~o.Jac_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'Jacobian_u(x + (ii-1)*%d, u + ii*%d, G);' '\n'],nx,nu)];
    end
end
if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'diffX(dx, x + (ii-1)*%d, xref + (ii-1)*%d);' '\n'],nx, nx)];
    code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'diffU(du, u + ii*%d, uref + ii*%d);' '\n'],nu, nu)]; 
    code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost(dx, du, &tmp_cost);' '\n'])];
    if ~cost_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost_u(dx, du, Ru);' '\n'])];
    end
else
    code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost(x + (ii-1)*%d, u + ii*%d, &tmp_cost);' '\n'],nx,nu)];
    if ~cost_u_static
        code = [code, sprintf([o.indent.generic o.indent.generic o.indent.generic 'cost_u( x + (ii-1)*%d, u + ii*%d, Ru);' '\n'],nx,nu)];
    end
end
code = [code, sprintf([o.indent.generic o.indent.generic '}\n\n'])];
code = [code, sprintf([o.indent.generic o.indent.generic '(*J) += tmp_cost;' '\n'])];

if o.terminal
    code = [code, sprintf([o.indent.generic o.indent.generic 'product_contr_nx(mem_tmp_contr, Px_contr, F);' '\n',...
        o.indent.generic o.indent.generic 'copy_nx(Px_contr,mem_tmp_contr);' '\n',...
        o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[ii*%d], Px_contr, G);' '\n\n'], nu)];
end



if trackRef
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffU(du, u + (ii+1)*%d, uref + (ii+1)*%d);' '\n'],nu, nu)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'diffX(dx, x + ii*%d, xref + ii*%d);' '\n'],nx, nx)];
    if ~cost_x_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'cost_x(dx, du, Qx);' '\n'])];
    end
else
    if ~cost_x_static
        code = [code, sprintf([o.indent.generic o.indent.generic 'cost_x( x + ii*%d, u + (ii+1)*%d, Qx);' '\n'],nx,nu)];
    end
end

if o.contractive
    code = [code, sprintf([o.indent.generic o.indent.generic 'index = ind - ii - 1;' '\n'])];
    
    code = [code, sprintf([o.indent.generic o.indent.generic 'if ( index >= 0){' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'if (index > 0){' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'product_contr_nx(mem_tmp_contr, Px_contr, F);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic o.indent.generic 'copy_nx(Px_contr,mem_tmp_contr);' '\n',...
        o.indent.generic o.indent.generic o.indent.generic '}' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'product_contr_nu(&dot_psi_N[ii*%d], Px_contr, G);' '\n'...
        o.indent.generic o.indent.generic '}' '\n'], nu)];
    code = [code, sprintf([o.indent.generic o.indent.generic 'else' '\n',...
        o.indent.generic o.indent.generic o.indent.generic 'set_zero_nu(&dot_psi_N[ii*%d]);' '\n\n'],nu)];
    
end

code = [code, sprintf([o.indent.generic o.indent.generic 'product_and_sum_nx(mem_tmp2, Px, ' data_cost_x ', F);' '\n',...
    o.indent.generic o.indent.generic 'copy_nx(Px,mem_tmp2);' '\n',...
    o.indent.generic o.indent.generic 'product_and_sum_nu(dot_J + ii*%d, Px, ' data_cost_u ', G);' '\n',...
    o.indent.generic '}' '\n',...
    '}' '\n\n'],nu)];
info.flops.it.mul = info.flops.it.mul+ 1*(N-1);



end
