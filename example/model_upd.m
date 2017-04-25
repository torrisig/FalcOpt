function xp = model_upd (x,u,Ts)

Rr = 1.1747;
Lrbar = 0.210;
Lhbar = 0.2011;
b = Rr*Lhbar/Lrbar;
a = -Rr/Lrbar;
J = 0.0053;
Tload = 0.3;
k = 1e-2;


xp = [(1+Ts*a)*x(1) + Ts*b*u(1);...
    (1-Ts*k/J)*x(2) + Ts/J*(Lhbar/Lrbar*x(1)*u(2) -  Tload)];

end