/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8563710128349221921);
void inv_err_fun(double *nom_x, double *true_x, double *out_5578910876619024124);
void H_mod_fun(double *state, double *out_8546593595067059713);
void f_fun(double *state, double dt, double *out_8277233601273549051);
void F_fun(double *state, double dt, double *out_574409882676075772);
void h_3(double *state, double *unused, double *out_9116543606686938938);
void H_3(double *state, double *unused, double *out_7922626234726440817);
void h_4(double *state, double *unused, double *out_5090413725090204401);
void H_4(double *state, double *unused, double *out_7647411873935716739);
void h_9(double *state, double *unused, double *out_1572386837559577218);
void H_9(double *state, double *unused, double *out_1721122958569996056);
void h_10(double *state, double *unused, double *out_642041116040745982);
void H_10(double *state, double *unused, double *out_1989108395848182367);
void h_12(double *state, double *unused, double *out_9091083175107654761);
void H_12(double *state, double *unused, double *out_5540856362196507170);
void h_31(double *state, double *unused, double *out_3464164092081030860);
void H_31(double *state, double *unused, double *out_8215946137339926991);
void h_32(double *state, double *unused, double *out_1528470715396010716);
void H_32(double *state, double *unused, double *out_3930432192109223439);
void h_13(double *state, double *unused, double *out_1979309428128721026);
void H_13(double *state, double *unused, double *out_3304654588672431724);
void h_14(double *state, double *unused, double *out_1572386837559577218);
void H_14(double *state, double *unused, double *out_1721122958569996056);
void h_19(double *state, double *unused, double *out_2559465969196434591);
void H_19(double *state, double *unused, double *out_6153802809672427237);
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