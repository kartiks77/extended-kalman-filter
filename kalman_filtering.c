# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>

# define N 4
# define M 1
# define Z 1

double T = 0.05 ;
double T22 = 0.05*0.05/2 ;
double CDEG2RAD = 22/7/180.0 ;
double kga = -21.265126512883949 ;
double as = 5 ;

double THETA_SERVO_MAX = 50 ;
double PWM_AMPLTD = 500 ;

# define kMax 300

double offset = 50 ;

# include "kf_functions.h"

int main(){

    /*Extracting Values from .txt files*/

    double t_abs[kMax], xPosMeas_v[kMax], xThetPWMcmd_v[kMax], xThetServSim_v[kMax], u[kMax] , y[kMax];

    char *filename1 = "AbsTime_v.txt" ;
    fileopener(filename1, t_abs) ;

    char *filename2 = "xPosMeas_v.txt" ;
    fileopener(filename2, xPosMeas_v) ;

    char *filename3 = "xThetPWMcmd_v.txt" ;
    fileopener_offset(filename3, xThetPWMcmd_v) ;

    calc_input_u(xThetPWMcmd_v, u) ;

    memcpy(y, xPosMeas_v, sizeof(y));

    /* Define System Matrices */

    double A4k[N][N] = {
        {0, 1, 0, 0},
        {0, 0, kga, kga},
        {0, 0, -as, 0},
        {0, 0, 0, 0}
    } ;

    double b4[N][M] = {
        {0},
        {0},
        {as},
        {0}
    } ;

    double H4[Z][N] = {
        {1, 0, 0, 0}
    } ;

    double D4 = 0 ;

    double I4[N][N] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} } ;

    double F4k[N][N] = {0} ;

    double g4k[N][M] = {0} ;
    
    double thetaSk_rad = 0 ;

    double deltathetaSk_rad = 0 ;

    double a23 = 0 ;

    double a24 = 0 ;

    double sigmaXposPxl = sqrt(0.5) ;

    double R[Z][Z] = { {sigmaXposPxl * sigmaXposPxl} } ;

    double qVarXpos = 5.0 ;
    double qVarXdot = 200 ;
    double qVarTheta = 5.0 ;
    double qVarThetOffs = 0.1 ;

    double Q4k[N][N] = {
        {qVarXpos, 0, 0, 0},
        {0, qVarXdot, 0, 0},
        {0, 0, qVarTheta, 0},
        {0, 0, 0, qVarThetOffs}
    } ;

    double P4k[N][N] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} } ;

    double P4k1[N][N] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} } ;

    double K4[N][Z] = { {0},{0},{0},{0} } ;

    double xcap4_post[N][1] = { {xPosMeas_v[0]},{0},{0},{0} } ;

    double xcap4_prio[N][1] = { {xPosMeas_v[0]},{0},{0},{0} } ;

    double xcap4_post_rec[N][kMax] = {0} ;
    double xcap4_prio_rec[N][kMax] = {0} ;
    double K4_rec[N][kMax] = {0} ;
    double P4k_rec[N][kMax] = {0} ;
    double P4k1_rec[N][kMax] = {0} ;

    printf("EKF4 Setup Complete\n\n") ;

    /* EKF Activation */

    printf("EKF4 Activated\n\n") ;

    for( int k=0 ; k < kMax ; k++ )
    {
        /* System Update */

        thetaSk_rad = xcap4_post[2][0]*CDEG2RAD ;
        deltathetaSk_rad = xcap4_post[3][0]*CDEG2RAD ;

        a23 = kga*cos(thetaSk_rad + deltathetaSk_rad) ;
        a24 = a23 ;

        A4k[1][2] = a23 ;
        A4k[1][3] = a24 ;

        /*F4k = eye(4) + A4k*T + A4k^2*T^2/2*/
        calc_F4k(A4k, I4, F4k) ;

        g4k[1][0] = a23*T22*as ;
        g4k[2][0] = as*T - as*T22 ;

        /* Time Update */

        /* P4k1 = F4k*P4k*F4k' + Q4k ; */
        calc_P4k1(F4k, P4k, Q4k, P4k1) ;

        /* xcap4_prio = F4k*xcap4_post + g4k*u(k) */
        calc_xcap4_prio(F4k, xcap4_post, g4k, u[k], xcap4_prio) ;

        /* Measurement Update */

        /* K4 = P4k1*H4'*(H4*P4k1*H4' + R)^-1 ; */
        calc_K4(P4k1, H4, R, K4) ;

        /* xcap4_post = xcap4_prio + K4*(y(k) - H4*xcap4_prio - D*u(k)) ; */
        calc_xcap4_post(xcap4_prio, K4, y[k], H4, D4, u[k], xcap4_post) ;

        /* P4k = (eye(4) - K4*H4)*P4k1 ; */
        calc_P4k(I4, K4, H4, P4k1, P4k) ;

        /* Record Values */

        record_Nx1_vector(xcap4_post, k, xcap4_post_rec) ;
        record_Nx1_vector(xcap4_prio, k, xcap4_prio_rec) ;
        record_Nx1_vector(K4, k, K4_rec) ;

        record_NxN_diagonals(P4k, k, P4k_rec) ;
        record_NxN_diagonals(P4k1, k, P4k1_rec) ;
    }

    printf("EKF4 Values Recorded\n\n") ;

    printf("Parameter Estimation delta offset = %0.2lf Degrees\n\n", xcap4_post_rec[3][kMax-1]) ;

    double ploty[kMax] = {0} ;

    extract_row(xcap4_post_rec, 3, ploty) ;
    
    return 0;
}