/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2353261711785034300);
void inv_err_fun(double *nom_x, double *true_x, double *out_1289726089137086488);
void H_mod_fun(double *state, double *out_865462163854239133);
void f_fun(double *state, double dt, double *out_2529703530971468786);
void F_fun(double *state, double dt, double *out_7037179984384406310);
void h_25(double *state, double *unused, double *out_3284274841600165049);
void H_25(double *state, double *unused, double *out_9169984535917636411);
void h_24(double *state, double *unused, double *out_5173473705231232781);
void H_24(double *state, double *unused, double *out_3053204908593227049);
void h_30(double *state, double *unused, double *out_4862210440789913083);
void H_30(double *state, double *unused, double *out_1694529923267701359);
void h_26(double *state, double *unused, double *out_6528428809595139310);
void H_26(double *state, double *unused, double *out_6761496803203334989);
void h_27(double *state, double *unused, double *out_1850378850566814470);
void H_27(double *state, double *unused, double *out_470425778667975665);
void h_29(double *state, double *unused, double *out_4827364935328157768);
void H_29(double *state, double *unused, double *out_4380796137653531637);
void h_28(double *state, double *unused, double *out_7863686467677602205);
void H_28(double *state, double *unused, double *out_2868671003221819373);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
