
/* dynamics of the system*/
void model_mpc(const double* x,const double* u, double* xp);

/* jacobian of dynamics w.r.t. x */
void Jacobian_x(const double* x, const double* u, double* F);

/* jacobian of dynamics w.r.t. u */
void Jacobian_u(const double* x, const double* u, double* G);

/* nonlinear constraint. It can depend on stage k. For example:
  if( k<N-2){ 
    n[0] = ...;
  }
  else ...;                                                     */
void build_n(const double* u, const unsigned int k, double* n);

/* jacobian_u of nonlinear constraint */
void build_Dn(const double* u, const unsigned int k, double* Dn);
