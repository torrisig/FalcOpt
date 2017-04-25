function xdot = model_upd_TC (t,x,u)

Rr = 1.1747;
Lrbar = 0.210;
Lhbar = 0.2011;
b = Rr*Lhbar/Lrbar;
a = -Rr/Lrbar;
J = 0.0053;
Tload = 0.3;
k = 1e-2;


xdot = [a*x(1) + b*u(1);...
    -k/J*x(2) + 1/J*(Lhbar/Lrbar*x(1)*u(2) -  Tload)];

end