
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2353261711785034300) {
   out_2353261711785034300[0] = delta_x[0] + nom_x[0];
   out_2353261711785034300[1] = delta_x[1] + nom_x[1];
   out_2353261711785034300[2] = delta_x[2] + nom_x[2];
   out_2353261711785034300[3] = delta_x[3] + nom_x[3];
   out_2353261711785034300[4] = delta_x[4] + nom_x[4];
   out_2353261711785034300[5] = delta_x[5] + nom_x[5];
   out_2353261711785034300[6] = delta_x[6] + nom_x[6];
   out_2353261711785034300[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1289726089137086488) {
   out_1289726089137086488[0] = -nom_x[0] + true_x[0];
   out_1289726089137086488[1] = -nom_x[1] + true_x[1];
   out_1289726089137086488[2] = -nom_x[2] + true_x[2];
   out_1289726089137086488[3] = -nom_x[3] + true_x[3];
   out_1289726089137086488[4] = -nom_x[4] + true_x[4];
   out_1289726089137086488[5] = -nom_x[5] + true_x[5];
   out_1289726089137086488[6] = -nom_x[6] + true_x[6];
   out_1289726089137086488[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_865462163854239133) {
   out_865462163854239133[0] = 1.0;
   out_865462163854239133[1] = 0.0;
   out_865462163854239133[2] = 0.0;
   out_865462163854239133[3] = 0.0;
   out_865462163854239133[4] = 0.0;
   out_865462163854239133[5] = 0.0;
   out_865462163854239133[6] = 0.0;
   out_865462163854239133[7] = 0.0;
   out_865462163854239133[8] = 0.0;
   out_865462163854239133[9] = 1.0;
   out_865462163854239133[10] = 0.0;
   out_865462163854239133[11] = 0.0;
   out_865462163854239133[12] = 0.0;
   out_865462163854239133[13] = 0.0;
   out_865462163854239133[14] = 0.0;
   out_865462163854239133[15] = 0.0;
   out_865462163854239133[16] = 0.0;
   out_865462163854239133[17] = 0.0;
   out_865462163854239133[18] = 1.0;
   out_865462163854239133[19] = 0.0;
   out_865462163854239133[20] = 0.0;
   out_865462163854239133[21] = 0.0;
   out_865462163854239133[22] = 0.0;
   out_865462163854239133[23] = 0.0;
   out_865462163854239133[24] = 0.0;
   out_865462163854239133[25] = 0.0;
   out_865462163854239133[26] = 0.0;
   out_865462163854239133[27] = 1.0;
   out_865462163854239133[28] = 0.0;
   out_865462163854239133[29] = 0.0;
   out_865462163854239133[30] = 0.0;
   out_865462163854239133[31] = 0.0;
   out_865462163854239133[32] = 0.0;
   out_865462163854239133[33] = 0.0;
   out_865462163854239133[34] = 0.0;
   out_865462163854239133[35] = 0.0;
   out_865462163854239133[36] = 1.0;
   out_865462163854239133[37] = 0.0;
   out_865462163854239133[38] = 0.0;
   out_865462163854239133[39] = 0.0;
   out_865462163854239133[40] = 0.0;
   out_865462163854239133[41] = 0.0;
   out_865462163854239133[42] = 0.0;
   out_865462163854239133[43] = 0.0;
   out_865462163854239133[44] = 0.0;
   out_865462163854239133[45] = 1.0;
   out_865462163854239133[46] = 0.0;
   out_865462163854239133[47] = 0.0;
   out_865462163854239133[48] = 0.0;
   out_865462163854239133[49] = 0.0;
   out_865462163854239133[50] = 0.0;
   out_865462163854239133[51] = 0.0;
   out_865462163854239133[52] = 0.0;
   out_865462163854239133[53] = 0.0;
   out_865462163854239133[54] = 1.0;
   out_865462163854239133[55] = 0.0;
   out_865462163854239133[56] = 0.0;
   out_865462163854239133[57] = 0.0;
   out_865462163854239133[58] = 0.0;
   out_865462163854239133[59] = 0.0;
   out_865462163854239133[60] = 0.0;
   out_865462163854239133[61] = 0.0;
   out_865462163854239133[62] = 0.0;
   out_865462163854239133[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_2529703530971468786) {
   out_2529703530971468786[0] = state[0];
   out_2529703530971468786[1] = state[1];
   out_2529703530971468786[2] = state[2];
   out_2529703530971468786[3] = state[3];
   out_2529703530971468786[4] = state[4];
   out_2529703530971468786[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2529703530971468786[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2529703530971468786[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7037179984384406310) {
   out_7037179984384406310[0] = 1;
   out_7037179984384406310[1] = 0;
   out_7037179984384406310[2] = 0;
   out_7037179984384406310[3] = 0;
   out_7037179984384406310[4] = 0;
   out_7037179984384406310[5] = 0;
   out_7037179984384406310[6] = 0;
   out_7037179984384406310[7] = 0;
   out_7037179984384406310[8] = 0;
   out_7037179984384406310[9] = 1;
   out_7037179984384406310[10] = 0;
   out_7037179984384406310[11] = 0;
   out_7037179984384406310[12] = 0;
   out_7037179984384406310[13] = 0;
   out_7037179984384406310[14] = 0;
   out_7037179984384406310[15] = 0;
   out_7037179984384406310[16] = 0;
   out_7037179984384406310[17] = 0;
   out_7037179984384406310[18] = 1;
   out_7037179984384406310[19] = 0;
   out_7037179984384406310[20] = 0;
   out_7037179984384406310[21] = 0;
   out_7037179984384406310[22] = 0;
   out_7037179984384406310[23] = 0;
   out_7037179984384406310[24] = 0;
   out_7037179984384406310[25] = 0;
   out_7037179984384406310[26] = 0;
   out_7037179984384406310[27] = 1;
   out_7037179984384406310[28] = 0;
   out_7037179984384406310[29] = 0;
   out_7037179984384406310[30] = 0;
   out_7037179984384406310[31] = 0;
   out_7037179984384406310[32] = 0;
   out_7037179984384406310[33] = 0;
   out_7037179984384406310[34] = 0;
   out_7037179984384406310[35] = 0;
   out_7037179984384406310[36] = 1;
   out_7037179984384406310[37] = 0;
   out_7037179984384406310[38] = 0;
   out_7037179984384406310[39] = 0;
   out_7037179984384406310[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7037179984384406310[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7037179984384406310[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7037179984384406310[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7037179984384406310[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7037179984384406310[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7037179984384406310[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7037179984384406310[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7037179984384406310[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7037179984384406310[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7037179984384406310[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7037179984384406310[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7037179984384406310[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7037179984384406310[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7037179984384406310[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7037179984384406310[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7037179984384406310[56] = 0;
   out_7037179984384406310[57] = 0;
   out_7037179984384406310[58] = 0;
   out_7037179984384406310[59] = 0;
   out_7037179984384406310[60] = 0;
   out_7037179984384406310[61] = 0;
   out_7037179984384406310[62] = 0;
   out_7037179984384406310[63] = 1;
}
void h_25(double *state, double *unused, double *out_3284274841600165049) {
   out_3284274841600165049[0] = state[6];
}
void H_25(double *state, double *unused, double *out_9169984535917636411) {
   out_9169984535917636411[0] = 0;
   out_9169984535917636411[1] = 0;
   out_9169984535917636411[2] = 0;
   out_9169984535917636411[3] = 0;
   out_9169984535917636411[4] = 0;
   out_9169984535917636411[5] = 0;
   out_9169984535917636411[6] = 1;
   out_9169984535917636411[7] = 0;
}
void h_24(double *state, double *unused, double *out_5173473705231232781) {
   out_5173473705231232781[0] = state[4];
   out_5173473705231232781[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3053204908593227049) {
   out_3053204908593227049[0] = 0;
   out_3053204908593227049[1] = 0;
   out_3053204908593227049[2] = 0;
   out_3053204908593227049[3] = 0;
   out_3053204908593227049[4] = 1;
   out_3053204908593227049[5] = 0;
   out_3053204908593227049[6] = 0;
   out_3053204908593227049[7] = 0;
   out_3053204908593227049[8] = 0;
   out_3053204908593227049[9] = 0;
   out_3053204908593227049[10] = 0;
   out_3053204908593227049[11] = 0;
   out_3053204908593227049[12] = 0;
   out_3053204908593227049[13] = 1;
   out_3053204908593227049[14] = 0;
   out_3053204908593227049[15] = 0;
}
void h_30(double *state, double *unused, double *out_4862210440789913083) {
   out_4862210440789913083[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1694529923267701359) {
   out_1694529923267701359[0] = 0;
   out_1694529923267701359[1] = 0;
   out_1694529923267701359[2] = 0;
   out_1694529923267701359[3] = 0;
   out_1694529923267701359[4] = 1;
   out_1694529923267701359[5] = 0;
   out_1694529923267701359[6] = 0;
   out_1694529923267701359[7] = 0;
}
void h_26(double *state, double *unused, double *out_6528428809595139310) {
   out_6528428809595139310[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6761496803203334989) {
   out_6761496803203334989[0] = 0;
   out_6761496803203334989[1] = 0;
   out_6761496803203334989[2] = 0;
   out_6761496803203334989[3] = 0;
   out_6761496803203334989[4] = 0;
   out_6761496803203334989[5] = 0;
   out_6761496803203334989[6] = 0;
   out_6761496803203334989[7] = 1;
}
void h_27(double *state, double *unused, double *out_1850378850566814470) {
   out_1850378850566814470[0] = state[3];
}
void H_27(double *state, double *unused, double *out_470425778667975665) {
   out_470425778667975665[0] = 0;
   out_470425778667975665[1] = 0;
   out_470425778667975665[2] = 0;
   out_470425778667975665[3] = 1;
   out_470425778667975665[4] = 0;
   out_470425778667975665[5] = 0;
   out_470425778667975665[6] = 0;
   out_470425778667975665[7] = 0;
}
void h_29(double *state, double *unused, double *out_4827364935328157768) {
   out_4827364935328157768[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4380796137653531637) {
   out_4380796137653531637[0] = 0;
   out_4380796137653531637[1] = 1;
   out_4380796137653531637[2] = 0;
   out_4380796137653531637[3] = 0;
   out_4380796137653531637[4] = 0;
   out_4380796137653531637[5] = 0;
   out_4380796137653531637[6] = 0;
   out_4380796137653531637[7] = 0;
}
void h_28(double *state, double *unused, double *out_7863686467677602205) {
   out_7863686467677602205[0] = state[5];
   out_7863686467677602205[1] = state[6];
}
void H_28(double *state, double *unused, double *out_2868671003221819373) {
   out_2868671003221819373[0] = 0;
   out_2868671003221819373[1] = 0;
   out_2868671003221819373[2] = 0;
   out_2868671003221819373[3] = 0;
   out_2868671003221819373[4] = 0;
   out_2868671003221819373[5] = 1;
   out_2868671003221819373[6] = 0;
   out_2868671003221819373[7] = 0;
   out_2868671003221819373[8] = 0;
   out_2868671003221819373[9] = 0;
   out_2868671003221819373[10] = 0;
   out_2868671003221819373[11] = 0;
   out_2868671003221819373[12] = 0;
   out_2868671003221819373[13] = 0;
   out_2868671003221819373[14] = 1;
   out_2868671003221819373[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
