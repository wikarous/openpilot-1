
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
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_455506380669311009) {
   out_455506380669311009[0] = delta_x[0] + nom_x[0];
   out_455506380669311009[1] = delta_x[1] + nom_x[1];
   out_455506380669311009[2] = delta_x[2] + nom_x[2];
   out_455506380669311009[3] = delta_x[3] + nom_x[3];
   out_455506380669311009[4] = delta_x[4] + nom_x[4];
   out_455506380669311009[5] = delta_x[5] + nom_x[5];
   out_455506380669311009[6] = delta_x[6] + nom_x[6];
   out_455506380669311009[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8624225469945293820) {
   out_8624225469945293820[0] = -nom_x[0] + true_x[0];
   out_8624225469945293820[1] = -nom_x[1] + true_x[1];
   out_8624225469945293820[2] = -nom_x[2] + true_x[2];
   out_8624225469945293820[3] = -nom_x[3] + true_x[3];
   out_8624225469945293820[4] = -nom_x[4] + true_x[4];
   out_8624225469945293820[5] = -nom_x[5] + true_x[5];
   out_8624225469945293820[6] = -nom_x[6] + true_x[6];
   out_8624225469945293820[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_4453914907973678378) {
   out_4453914907973678378[0] = 1.0;
   out_4453914907973678378[1] = 0.0;
   out_4453914907973678378[2] = 0.0;
   out_4453914907973678378[3] = 0.0;
   out_4453914907973678378[4] = 0.0;
   out_4453914907973678378[5] = 0.0;
   out_4453914907973678378[6] = 0.0;
   out_4453914907973678378[7] = 0.0;
   out_4453914907973678378[8] = 0.0;
   out_4453914907973678378[9] = 1.0;
   out_4453914907973678378[10] = 0.0;
   out_4453914907973678378[11] = 0.0;
   out_4453914907973678378[12] = 0.0;
   out_4453914907973678378[13] = 0.0;
   out_4453914907973678378[14] = 0.0;
   out_4453914907973678378[15] = 0.0;
   out_4453914907973678378[16] = 0.0;
   out_4453914907973678378[17] = 0.0;
   out_4453914907973678378[18] = 1.0;
   out_4453914907973678378[19] = 0.0;
   out_4453914907973678378[20] = 0.0;
   out_4453914907973678378[21] = 0.0;
   out_4453914907973678378[22] = 0.0;
   out_4453914907973678378[23] = 0.0;
   out_4453914907973678378[24] = 0.0;
   out_4453914907973678378[25] = 0.0;
   out_4453914907973678378[26] = 0.0;
   out_4453914907973678378[27] = 1.0;
   out_4453914907973678378[28] = 0.0;
   out_4453914907973678378[29] = 0.0;
   out_4453914907973678378[30] = 0.0;
   out_4453914907973678378[31] = 0.0;
   out_4453914907973678378[32] = 0.0;
   out_4453914907973678378[33] = 0.0;
   out_4453914907973678378[34] = 0.0;
   out_4453914907973678378[35] = 0.0;
   out_4453914907973678378[36] = 1.0;
   out_4453914907973678378[37] = 0.0;
   out_4453914907973678378[38] = 0.0;
   out_4453914907973678378[39] = 0.0;
   out_4453914907973678378[40] = 0.0;
   out_4453914907973678378[41] = 0.0;
   out_4453914907973678378[42] = 0.0;
   out_4453914907973678378[43] = 0.0;
   out_4453914907973678378[44] = 0.0;
   out_4453914907973678378[45] = 1.0;
   out_4453914907973678378[46] = 0.0;
   out_4453914907973678378[47] = 0.0;
   out_4453914907973678378[48] = 0.0;
   out_4453914907973678378[49] = 0.0;
   out_4453914907973678378[50] = 0.0;
   out_4453914907973678378[51] = 0.0;
   out_4453914907973678378[52] = 0.0;
   out_4453914907973678378[53] = 0.0;
   out_4453914907973678378[54] = 1.0;
   out_4453914907973678378[55] = 0.0;
   out_4453914907973678378[56] = 0.0;
   out_4453914907973678378[57] = 0.0;
   out_4453914907973678378[58] = 0.0;
   out_4453914907973678378[59] = 0.0;
   out_4453914907973678378[60] = 0.0;
   out_4453914907973678378[61] = 0.0;
   out_4453914907973678378[62] = 0.0;
   out_4453914907973678378[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6311072751643462934) {
   out_6311072751643462934[0] = state[0];
   out_6311072751643462934[1] = state[1];
   out_6311072751643462934[2] = state[2];
   out_6311072751643462934[3] = state[3];
   out_6311072751643462934[4] = state[4];
   out_6311072751643462934[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6311072751643462934[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6311072751643462934[7] = state[7];
}
void F_fun(double *state, double dt, double *out_8033668173644051700) {
   out_8033668173644051700[0] = 1;
   out_8033668173644051700[1] = 0;
   out_8033668173644051700[2] = 0;
   out_8033668173644051700[3] = 0;
   out_8033668173644051700[4] = 0;
   out_8033668173644051700[5] = 0;
   out_8033668173644051700[6] = 0;
   out_8033668173644051700[7] = 0;
   out_8033668173644051700[8] = 0;
   out_8033668173644051700[9] = 1;
   out_8033668173644051700[10] = 0;
   out_8033668173644051700[11] = 0;
   out_8033668173644051700[12] = 0;
   out_8033668173644051700[13] = 0;
   out_8033668173644051700[14] = 0;
   out_8033668173644051700[15] = 0;
   out_8033668173644051700[16] = 0;
   out_8033668173644051700[17] = 0;
   out_8033668173644051700[18] = 1;
   out_8033668173644051700[19] = 0;
   out_8033668173644051700[20] = 0;
   out_8033668173644051700[21] = 0;
   out_8033668173644051700[22] = 0;
   out_8033668173644051700[23] = 0;
   out_8033668173644051700[24] = 0;
   out_8033668173644051700[25] = 0;
   out_8033668173644051700[26] = 0;
   out_8033668173644051700[27] = 1;
   out_8033668173644051700[28] = 0;
   out_8033668173644051700[29] = 0;
   out_8033668173644051700[30] = 0;
   out_8033668173644051700[31] = 0;
   out_8033668173644051700[32] = 0;
   out_8033668173644051700[33] = 0;
   out_8033668173644051700[34] = 0;
   out_8033668173644051700[35] = 0;
   out_8033668173644051700[36] = 1;
   out_8033668173644051700[37] = 0;
   out_8033668173644051700[38] = 0;
   out_8033668173644051700[39] = 0;
   out_8033668173644051700[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8033668173644051700[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8033668173644051700[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8033668173644051700[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8033668173644051700[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8033668173644051700[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8033668173644051700[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8033668173644051700[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8033668173644051700[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8033668173644051700[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8033668173644051700[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8033668173644051700[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8033668173644051700[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8033668173644051700[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8033668173644051700[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8033668173644051700[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8033668173644051700[56] = 0;
   out_8033668173644051700[57] = 0;
   out_8033668173644051700[58] = 0;
   out_8033668173644051700[59] = 0;
   out_8033668173644051700[60] = 0;
   out_8033668173644051700[61] = 0;
   out_8033668173644051700[62] = 0;
   out_8033668173644051700[63] = 1;
}
void h_25(double *state, double *unused, double *out_6845862044070992032) {
   out_6845862044070992032[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7634934532946663080) {
   out_7634934532946663080[0] = 0;
   out_7634934532946663080[1] = 0;
   out_7634934532946663080[2] = 0;
   out_7634934532946663080[3] = 0;
   out_7634934532946663080[4] = 0;
   out_7634934532946663080[5] = 0;
   out_7634934532946663080[6] = 1;
   out_7634934532946663080[7] = 0;
}
void h_24(double *state, double *unused, double *out_1806551022244134745) {
   out_1806551022244134745[0] = state[4];
   out_1806551022244134745[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2103276594648636244) {
   out_2103276594648636244[0] = 0;
   out_2103276594648636244[1] = 0;
   out_2103276594648636244[2] = 0;
   out_2103276594648636244[3] = 0;
   out_2103276594648636244[4] = 1;
   out_2103276594648636244[5] = 0;
   out_2103276594648636244[6] = 0;
   out_2103276594648636244[7] = 0;
   out_2103276594648636244[8] = 0;
   out_2103276594648636244[9] = 0;
   out_2103276594648636244[10] = 0;
   out_2103276594648636244[11] = 0;
   out_2103276594648636244[12] = 0;
   out_2103276594648636244[13] = 1;
   out_2103276594648636244[14] = 0;
   out_2103276594648636244[15] = 0;
}
void h_30(double *state, double *unused, double *out_1047077450813531327) {
   out_1047077450813531327[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1978964378002520200) {
   out_1978964378002520200[0] = 0;
   out_1978964378002520200[1] = 0;
   out_1978964378002520200[2] = 0;
   out_1978964378002520200[3] = 0;
   out_1978964378002520200[4] = 1;
   out_1978964378002520200[5] = 0;
   out_1978964378002520200[6] = 0;
   out_1978964378002520200[7] = 0;
}
void h_26(double *state, double *unused, double *out_6048750718474376972) {
   out_6048750718474376972[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6293498399251015528) {
   out_6293498399251015528[0] = 0;
   out_6293498399251015528[1] = 0;
   out_6293498399251015528[2] = 0;
   out_6293498399251015528[3] = 0;
   out_6293498399251015528[4] = 0;
   out_6293498399251015528[5] = 0;
   out_6293498399251015528[6] = 0;
   out_6293498399251015528[7] = 1;
}
void h_27(double *state, double *unused, double *out_2136107127935848733) {
   out_2136107127935848733[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3266546365839145512) {
   out_3266546365839145512[0] = 0;
   out_3266546365839145512[1] = 0;
   out_3266546365839145512[2] = 0;
   out_3266546365839145512[3] = 1;
   out_3266546365839145512[4] = 0;
   out_3266546365839145512[5] = 0;
   out_3266546365839145512[6] = 0;
   out_3266546365839145512[7] = 0;
}
void h_29(double *state, double *unused, double *out_1990587303715193762) {
   out_1990587303715193762[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1404405764146984025) {
   out_1404405764146984025[0] = 0;
   out_1404405764146984025[1] = 1;
   out_1404405764146984025[2] = 0;
   out_1404405764146984025[3] = 0;
   out_1404405764146984025[4] = 0;
   out_1404405764146984025[5] = 0;
   out_1404405764146984025[6] = 0;
   out_1404405764146984025[7] = 0;
}
void h_28(double *state, double *unused, double *out_788741996809673489) {
   out_788741996809673489[0] = state[5];
   out_788741996809673489[1] = state[6];
}
void H_28(double *state, double *unused, double *out_9133830807137324866) {
   out_9133830807137324866[0] = 0;
   out_9133830807137324866[1] = 0;
   out_9133830807137324866[2] = 0;
   out_9133830807137324866[3] = 0;
   out_9133830807137324866[4] = 0;
   out_9133830807137324866[5] = 1;
   out_9133830807137324866[6] = 0;
   out_9133830807137324866[7] = 0;
   out_9133830807137324866[8] = 0;
   out_9133830807137324866[9] = 0;
   out_9133830807137324866[10] = 0;
   out_9133830807137324866[11] = 0;
   out_9133830807137324866[12] = 0;
   out_9133830807137324866[13] = 0;
   out_9133830807137324866[14] = 1;
   out_9133830807137324866[15] = 0;
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
