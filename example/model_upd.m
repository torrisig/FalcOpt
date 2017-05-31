function xp = model_upd (x,u,Ts)

Rr = 1.1747;                    % rotor resistance
Lrbar = 0.210;                  % rotor inductance
Lhbar = 0.2011;                 % magnetization inductance
J = 0.0053;                     % motor inertia
Tload = 0.3;                    % load torque
k = 1e-2;                       % friction coefficient
a = -Rr/Lrbar;   
b1 = Rr*Lhbar/Lrbar;
b2 = Lhbar/Lrbar;

xp = [ x(1) + ( a*x(1) + b1*u(1))*Ts;...
       x(2) + 1/J*( b2*x(1)*u(2) - Tload - k*x(2))*Ts];

end