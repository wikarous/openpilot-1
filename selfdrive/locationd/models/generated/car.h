/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_455506380669311009);
void inv_err_fun(double *nom_x, double *true_x, double *out_8624225469945293820);
void H_mod_fun(double *state, double *out_4453914907973678378);
void f_fun(double *state, double dt, double *out_6311072751643462934);
void F_fun(double *state, double dt, double *out_8033668173644051700);
void h_25(double *state, double *unused, double *out_6845862044070992032);
void H_25(double *state, double *unused, double *out_7634934532946663080);
void h_24(double *state, double *unused, double *out_1806551022244134745);
void H_24(double *state, double *unused, double *out_2103276594648636244);
void h_30(double *state, double *unused, double *out_1047077450813531327);
void H_30(double *state, double *unused, double *out_1978964378002520200);
void h_26(double *state, double *unused, double *out_6048750718474376972);
void H_26(double *state, double *unused, double *out_6293498399251015528);
void h_27(double *state, double *unused, double *out_2136107127935848733);
void H_27(double *state, double *unused, double *out_3266546365839145512);
void h_29(double *state, double *unused, double *out_1990587303715193762);
void H_29(double *state, double *unused, double *out_1404405764146984025);
void h_28(double *state, double *unused, double *out_788741996809673489);
void H_28(double *state, double *unused, double *out_9133830807137324866);
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
