function [dx] = Jacobian_x(x,u)

% parameters
Rr = 1.1747;
Lrbar = 0.210;
Lhbar = 0.2011;
b = Rr*Lhbar/Lrbar;
a = -Rr/Lrbar;
J = 0.0053;
Tload = 0.3;
k = 1e-2;
Ts = 0.01;

% Jacobian matrix
dx = [(1+Ts*a),                 0;...
      (Ts/J*Lhbar/Lrbar*u(2)),  (1-Ts*k/J)]; 

end