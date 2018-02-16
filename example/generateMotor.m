function generateMotor(par)

    % options
    eps = 1e-3;                   % tolerance
    merit_function = 0;         % merit function
    debug = 3;                    % level of debug

    contractive = false;        % no constractive constraints
    terminal = false;             % no terminal constraints
    gradients = 'matlab';
    precision = 'double';

    %% Dynamics of the system
    dynamics = @(x,u) model_upd(x,u,par.Ts);

    % defining cost function
    J.Q = par.Q;
    J.R = par.R;
    J.P = par.P;

    J.trackReference = true;                     % to track a desired (possibly time-varying) reference

    %% code generation ( with 4 different ways to generate derivatives)
    warning('OFF', 'falcopt:MissingImplemementation:SymmetricMatrix');
    switch gradients
        case 'casadi' 
            % first option: Automatic differentiation via CasADi to generate derivatives

            % use of internally defined and variable step size alpha (default)

            variable_stepSize.active = true;

            info = falcopt.generateCode(dynamics,par.N,par.nx,par.nu, J,...
                'variable_stepSize',variable_stepSize,...
                'constraints', par.constraint,'nn',par.nn, 'gradients', gradients,...
                'lb',par.umin, 'ub', par.umax,...
                'contractive',contractive, 'terminal', terminal, ...
                'debug',debug,'merit_function', merit_function,...
                'eps',eps,'precision', precision, ...
                'name', 'Motor_example_FalcOpt', 'gendir', 'generatedCode');
        case 'matlab'
            % second option: Automatic differentiation via Matlab symbolic toolbox

            % use of user-defined step size alpha (may perform better but
            % requires tuning of "variable_stepSize.alpha_max")

            variable_stepSize.active = false;
            variable_stepSize.alpha_max = 5;


            info = falcopt.generateCode(dynamics,par.N,par.nx,par.nu, J,...
                'variable_stepSize',variable_stepSize,...
                'constraints', par.constraint,'nn',par.nn, 'gradients', gradients,...
                'lb',par.umin, 'ub', par.umax,...
                'contractive',contractive, 'terminal', terminal, ...
                'debug',debug,'merit_function', merit_function,...
                'verbose', 0, ...
                'eps',eps,'precision', precision,...
                'buildTypes', {'mex', 'production'}, ...
                'name', 'Motor_example_FalcOpt', 'gendir', 'generatedCode');

        case 'manual'

            % third option: specify derivatives by hand, with constant step size alpha
            external.manual.jacobian_x = @Jacobian_x;
            external.manual.jacobian_u = @Jacobian_u;
            external.manual.jacobian_n = @(u) u;

            variable_stepSize.active = false;
            variable_stepSize.alpha_max = 5;

            info = falcopt.generateCode(dynamics,par.N,par.nx,par.nu, J,...
                'constraints', par.constraint,'nn',par.nn, 'gradients', gradients,...
                'lb',par.umin, 'ub', par.umax,...
                'contractive',contractive, 'terminal', terminal, 'variable_stepSize', variable_stepSize, ...
                'debug',debug,'merit_function', merit_function,...
                'eps',eps,'precision', precision,...
                'name', 'Motor_example_FalcOpt', 'gendir', 'generatedCode',...
                'external',external);

        case 'ccode'

            % fourth option: write C code that evaluates the model and the
            % jacobians and their structure (only for experienced users)
            [jac_x_struct,jac_u_struct,jac_n_struct, K_n] = jacobian_structure(par);   % function returning structure of jacobians
            external.ccode.jacobian_x.structure = jac_x_struct;
            external.ccode.jacobian_u.structure = jac_u_struct;
            external.ccode.jacobian_n.structure = jac_n_struct;
            external.ccode.K_n = K_n;
            external.ccode.c = './external_functions.c';
            external.ccode.h = './external_functions.h';

            variable_stepSize.active = false;
            variable_stepSize.alpha_max = 5;

            info = falcopt.generateCode(dynamics,par.N,par.nx,par.nu, J,...
                'nn',par.nn,'contractive',contractive, 'terminal', terminal, 'gradients', gradients, ...
                'debug',debug,'merit_function', merit_function,...
                'eps',eps,'precision', precision,'variable_stepSize', variable_stepSize, ...
                'lb',par.umin, 'ub', par.umax,...
                'external', external, ...
                'name', 'Motor_example_FalcOpt', 'gendir', 'generatedCode');
    end
    warning('ON', 'falcopt:MissingImplemementation:SymmetricMatric');
end

