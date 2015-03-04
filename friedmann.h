
double faexp_tilde(double a, double om, double ov);
double faexp(double a, double om, double ov);
double integ_da_dt_tilde(double a, double b, double om, double ov,double tol);
double integ_da_dt(double a, double b, double om, double ov,double tol);
double interp_aexp(double ttilde,double *tab_aexp,double *tab_ttilde);
void compute_friedmann(double amin, double amax, int ntab, double omegam, double omegav, double *tab_aexp, double *tab_ttilde, double *tab_t);
double dladt(double a, double omegam, double omegav);
double ddplus(double a,double omegam, double omegav);
double integ_ddplus(double a, double b, double om, double ov,double tol);
double dplus(double a, double omegam, double omegav);
double fomega(double a, double omegam, double omegav);

