
plot_figures = 1;           % 1 or 0
generate_solver = 1;        % 1: generate c-file solver
                            % 0: run simulation without generating solver

%% Problem parameters
Rr = 1.1747;
Lrbar = 0.210;
Lhbar = 0.2011;
b = Rr*Lhbar/Lrbar;
a = -Rr/Lrbar;
J = 0.0053;
Tload = 0.3;
k = 1e-2;

%% MPC Problem setup

% reference state and input
problem.xref = [0.6; 80];
problem.uref = [ -a/b * problem.xref(1); ...
                 Lrbar/Lhbar/problem.xref(1)*(Tload+k*problem.xref(2))];

% number states, inputs and nonlinear constraints
problem.nx = 2;
problem.nu = 2;
problem.nn = 1;

% MPC options
problem.N = 10;
problem.Q = 1e-1*diag([0,1e0]);
problem.R = 1.1747e-2*eye(problem.nu);
problem.P = problem.Q;
problem.Ts = 0.01;

% nonlinear constraints
problem.imax = 7.5690;
problem.constraint = @(u) 0.5*u'*u - problem.imax;

% inputs box constraints
problem.umin = [0.1;0.1];
problem.umax = sqrt(2*problem.imax)*[1;1];


%% generate solver
if generate_solver
    generateMotor(problem);
end

%% simulation

t_min = 1;

tf = 6;
Ts = 0.01;
n_steps = tf/Ts;
time = zeros(1,n_steps+1);
state = zeros(problem.nx,n_steps+1);
input = zeros(problem.nu,n_steps);


% variables to store results
aggregate_time = 0;
worst_time = 0;
best_time = 1e5;
agg_iter = 0;
worst_iter = 0;
best_iter = 1000;
agg_ls = 0;
worst_ls = 0;
best_ls = 1000;


% initialization
state(:,1) = [0, 0]';
U_pred = zeros(problem.N*problem.nu,1);
xref = repmat(problem.xref,problem.N,1);
uref = zeros(problem.N*problem.nu,1);

for ii = 1:n_steps
    
    if (mod(ii,floor(n_steps/10)) == 0)
        fprintf('simulation progress: %d of 10 \n',ii/floor(n_steps/10));
    end
    
    % proposed algorithm
    
    [~, flag, info] = Motor_example_FalcOpt(state(:,ii),xref,uref,U_pred);
    
    
    if flag < 0
        disp('Some problem in the solver')
        keyboard
    end
    
    % get input
    U_pred = reshape(info.u,problem.nu*problem.N,1);
    input(:,ii) = U_pred(1:problem.nu);
    
    % shift input sequence
    U_pred = [U_pred(problem.nu+1:end);U_pred(end-problem.nu+1:end)];
    
    % get number iterations
    iter = double(info.iterations);
    agg_iter = agg_iter + iter;
    worst_iter = max(iter,worst_iter);
    best_iter = min(iter,best_iter);
    lsperiter = sum(info.lineSearch.iterations(1:info.iterations))/double(info.iterations);
    agg_ls = agg_ls + lsperiter;
    worst_ls = max(lsperiter,worst_ls);
    best_ls = min(lsperiter,best_ls);
    
    % get computational time
    delta_t = info.time;
    aggregate_time= aggregate_time + delta_t;
    worst_time = max(delta_t,worst_time);
    best_time = min(delta_t,best_time);
    
    % System simulation
    [~,X_nonlin] = ode45(@model_upd_TC,[0 Ts],state(:,ii),odeset('RelTol',1e-9),input(:,ii));
    state(:,ii+1) = X_nonlin(end,:)';
    
    time(:,ii+1) = Ts*(ii);
end

%% Figure and results
% display algorithms results
fprintf('Proposed algorithm, average time = %d s, best %d s, worst case %d s, closed loop cost %d \n',aggregate_time/n_steps,best_time, worst_time,CostComputer(state,input,problem));
fprintf('Proposed algorithm, average iter = %d , best %d, worst case %d \n',agg_iter/n_steps,best_iter, worst_iter);
fprintf('Proposed algorithm, average line search iters per iter = %d , best %d, worst case %d \n\n',agg_ls/n_steps,best_ls, worst_ls);
% plots
if plot_figures
    figure()
    hold on
    plot(time(1,1:size(state,2)),state(1,:),'b');
    title('flux x_1')

    
    figure
    hold on
    plot(time(1,1:size(state,2)),state(2,:),'b');
    plot(time,problem.xref(2)*ones(size(time)),'b:');
    title('speed x_2')

    figure()
    hold on
    plot(time(1,1:size(input,2)),input(1,:),'g');
    plot(time(1,1:size(input,2)),input(2,:),'b');
    title('Input')
end

