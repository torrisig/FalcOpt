
/* paramters for functions */
const double L_ratio = 0.957619047619048, b = 1.124915095238095, a = -5.593809523809524,
		J_inv = 1.886792452830189e+02, Tload = 0.3, k = 1e-2, Ts = 0.01;

/* dynamics of the system */
void model_mpc(const double* x,const double* u, double* xp){
	  xp[0] = (1 + Ts*a)*x[0] + Ts*b*u[0];
	  xp[1] = (1 - Ts*k *J_inv)*x[1] + Ts *J_inv *(L_ratio *x[0]*u[1] - Tload);
}

/* jacobian of dynamics w.r.t. x */
void Jacobian_x(const double* x, const double* u, double* F){

	  F[0] = (1 + Ts*a);
	  F[1] = Ts * J_inv * L_ratio * u[1];
	  F[2] = 1 - Ts * k * J_inv;
}

/* jacobian of dynamics w.r.t. u */
void Jacobian_u(const double* x, const double* u, double* G){

	  G[0] = Ts * b;
	  G[1] = Ts * J_inv * L_ratio * x[0] ;
}

/* nonlinear constraint. It can depend on stage k. For example:
  if( k<N-2){ 
    n[0] = ...;
  }
  else ...;                                                     */
void build_n(const double* u, const unsigned int k, double* n){
    n[0] = u[0]*u[0]*0.5 + u[1]*u[1]*0.5 - 7.569;
}

/* jacobian_u of nonlinear constraint */
void build_Dn(const double* u, const unsigned int k, double* Dn){
    Dn[0] = u[0];
    Dn[1] = u[1];
}
