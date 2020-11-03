/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2979432275415314763);
void inv_err_fun(double *nom_x, double *true_x, double *out_8210458818486868805);
void H_mod_fun(double *state, double *out_2095463338726842165);
void f_fun(double *state, double dt, double *out_431075374196780526);
void F_fun(double *state, double dt, double *out_8284872676528927080);
void h_3(double *state, double *unused, double *out_2102870914289254584);
void H_3(double *state, double *unused, double *out_7584385583185505203);
void h_4(double *state, double *unused, double *out_295645676549977109);
void H_4(double *state, double *unused, double *out_6504910247896389024);
void h_9(double *state, double *unused, double *out_3833117803990374915);
void H_9(double *state, double *unused, double *out_9134740444023733999);
void h_10(double *state, double *unused, double *out_1451576847657428427);
void H_10(double *state, double *unused, double *out_7646617256047821503);
void h_12(double *state, double *unused, double *out_1387734146316736601);
void H_12(double *state, double *unused, double *out_6745049292577392449);
void h_31(double *state, double *unused, double *out_8404574603050392738);
void H_31(double *state, double *unused, double *out_2375026812758704951);
void h_32(double *state, double *unused, double *out_3389019376852184025);
void H_32(double *state, double *unused, double *out_2087431086203915649);
void h_13(double *state, double *unused, double *out_8929451567799493985);
void H_13(double *state, double *unused, double *out_1493934397849623260);
void h_14(double *state, double *unused, double *out_3833117803990374915);
void H_14(double *state, double *unused, double *out_9134740444023733999);
void h_19(double *state, double *unused, double *out_3451257868663390707);
void H_19(double *state, double *unused, double *out_4882504419065250251);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);