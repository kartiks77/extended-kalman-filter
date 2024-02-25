# include <stdio.h>
# include <string.h>
# include <stdlib.h>

# define matrows(mat) sizeof(mat)/sizeof(mat[0]);
# define matcols(mat) sizeof(mat[0]) / sizeof(mat[0][0]);

#define MAX_LINES 300
#define MAX_LEN 100

/*kMax is called from the Main C code*/

double fileopener(char *fileinput, double variable[kMax])
{
    char data[MAX_LINES][MAX_LEN] ;

	FILE *file ;
	
	file = fopen(fileinput, "r") ;

    if ( file == NULL )
    {
        printf("Error opening file yo\n") ;
    }

    int line = 0 ;

    while ( !feof(file) && !ferror(file) )
        if (fgets(data[line], MAX_LEN, file) != NULL )
        {
            variable[line] = strtod(data[line], NULL) ;
            line++ ;
        }

    fclose(file) ;

    printf("%s values Stored\n\n", fileinput) ;

	return *variable;	
}

double fileopener_offset(char *fileinput, double variable[kMax])
{
    char data[MAX_LINES][MAX_LEN] ;

	FILE *file ;
	
	file = fopen(fileinput, "r") ;

    if ( file == NULL )
    {
        printf("Error opening file yo\n") ;
    }

    int line = 0 ;

    while ( !feof(file) && !ferror(file) )
        if (fgets(data[line], MAX_LEN, file) != NULL )
        {
            variable[line] = strtod(data[line], NULL) - 50 + offset ;
            line++ ;
        }

    fclose(file) ;

    printf("%s values Stored\n\n", fileinput) ;

	return *variable;	
}

double calc_input_u(double xThetPWMcmd[kMax], double variable[kMax])
{
    for( int i=0 ; i < kMax ; i++ )
    {
        variable[i] = - THETA_SERVO_MAX * xThetPWMcmd[i] / PWM_AMPLTD;
    }

    return *variable;
}

double calc_F4k(double A[N][N], double I[N][N], double F[N][N])
{

    for( int i=0; i<N ; i++ ){
        for( int j=0 ;j<N ; j++ ){

            F[i][j] = I[i][j] + A[i][j]*T ;                      /*Fill up empty F4k matrix with these values*/

            for( int k=0 ; k<N ; k++ ){                          /*Then add amt multiplication values on top*/
                F[i][j] += A[i][k]*A[k][j]*T22 ;
            }

        }
    }
}

double calc_P4k1(double F4k[N][N], double P4k[N][N], double Q4k[N][N], double P4k1[N][N])
{

    double F4kP4k[N][N] = {0} ;                             /*Intermediate var*/
    
    for( int i=0; i<N ; i++ ){                              /* First F4k*P4k */
        for( int j=0 ;j<N ; j++ ){

            for( int k=0 ; k<N ; k++ ){                     /*Then add mat multiplication values on top*/
                F4kP4k[i][j] += F4k[i][k]*P4k[k][j] ;
            }

        }
    }

    for( int i=0; i<N ; i++ ){                              /*Here F4k*P4k*F4k'*/
        for( int j=0 ;j<N ; j++ ){

            P4k1[i][j] = 0 ;                                /*Fill up empty P4k1 mat */

            for( int k=0 ; k<N ; k++ ){                     /*Then add mat multiplication values on top*/
                P4k1[i][j] += F4kP4k[i][k]*F4k[j][k] ;      /*Here since F4k is transposed, F4k[j][k] is post multiplied instead of [j][k]*/
            }

        }
    }

    for( int i=0; i<N ; i++ ){                              /*Here F4k*P4k*F4k' + Q4k*/
        for( int j=0 ;j<N ; j++ ){

            P4k1[i][j] = P4k1[i][j] + Q4k[i][j] ;                                /*Add Q4k elements to F4k*P4k*F4k' */

        }
    }
}

double calc_xcap4_prio(double F4k[N][N], double xcap4_post[N][1], double g4k[N][M], double uk, double xcap4_prio[N][1])
{

    for( int i=0; i<N ; i++ ){                              /* First F4k*xcap4_post */

        xcap4_prio[i][0] = 0 ;                              /*Initialise*/
                                                            /*j loop is not needed as it is always zero*/
        for( int k=0 ; k<N ; k++ ){
            xcap4_prio[i][0] += F4k[i][k]*xcap4_post[k][0] ;
        }
    }

    for( int i=0; i<N ; i++ ){                              /* Now F4k*xcap4_post  + g4k*u(k) */
        
        xcap4_prio[i][0] += g4k[i][0]*uk ;
    }
}

double calc_K4(double P4k1[N][N], double H4[Z][N], double R[Z][Z], double K4[N][Z])
{

    double Pzz[Z][Z] = {0} ;                                /*Intermediate var*/

    double H4P4k1[Z][N] = {0} ;

    for( int i=0; i<Z ; i++ ){                              /* First H4*P4k1 */
        for( int j=0 ;j<N ; j++ ){

            for( int k=0 ; k<N ; k++ ){                     /* Then add mat multiplication values on top */
                H4P4k1[i][j] += H4[i][k]*P4k1[k][j] ;
            }

        }
    }

    for( int i=0; i<Z ; i++ ){                              /* Now Pzz = H4*P4k1*H4' + R */
        for( int j=0 ;j<Z ; j++ ){

            Pzz[i][j] = R[i][j] ;                           /* i j loops not needed as there is only one measurement */

            for( int k=0 ; k<N ; k++ ){                     /* Then add mat multiplication values on top */
                Pzz[i][j] += H4P4k1[i][k]*H4[j][k] ;        /* jk instead of kj because H4 transpose is post multiplied*/
            }

        }
    }

    for( int i=0; i<N ; i++ ){                              /* now P4k1*H4'/Pzz */
        for( int j=0 ;j<Z ; j++ ){

            K4[i][j] = 0 ;                                  /*j loop not needed since only one measurement*/

            for( int k=0 ; k<N ; k++ ){                     /*Then add mat multiplication values on top*/
                K4[i][j] += P4k1[i][k]*H4[j][k]/Pzz[0][0] ; /*hk instead of kj because H4 transpose is post multiplied*/
            }                                               /*Dividing by Pzz would be matrix multiplication of Pzz^-1 if multiply measurements and Pzz is a ZxZ matrix*/

        }
    }

}

double calc_xcap4_post(double xcap4_prio[N][1], double K4[N][Z], double y, double H4[Z][N], double D4, double uk, double xcap4_post[N][1])
{

    double yy[Z][1] = {0} ;                                 /* Intermediate var */

    for( int i=0; i<Z ; i++ ){                              /* First yy = y(k) - H4*xcap4_prio - D*u(k) */

        yy[i][0] = y - D4*uk ;                              /* D*uk would be matrix multiplication if multiple measurements or inputs */

        for( int k=0 ; k<N ; k++ ){                         /* Then add mat multiplication values on top */
            yy[i][0] -= H4[i][k]*xcap4_prio[k][0] ;
        }
    }

    for( int i=0; i<N ; i++ ){                                        /* Now xcap4_post = xcap4_prio + K4*yy */

        xcap4_post[i][0] = xcap4_prio[i][0] + K4[i][0]*yy[0][0] ;     /* Then add mat multiplication values on top */
    }                                                                 /* K4*yy would be matrix multiplication if more measurements */
}

double calc_P4k(double I4[N][N], double K4[N][Z], double H4[Z][N], double P4k1[N][N], double P4k[N][N])
{

    double I4_K4H4[N][N] = {0} ;                            /*Intermediate var*/
    
    for( int i=0; i<N ; i++ ){                              /* First I4_K4H4 = eye(4) - K4*H4 */
        for( int j=0 ;j<N ; j++ ){

            I4_K4H4[i][j] = I4[i][j] ;

            for( int k=0 ; k<Z ; k++ ){                     /*Then add mat multiplication values on top*/
                I4_K4H4[i][j] -= K4[i][k]*H4[k][j] ;        /* k loop not needed since only one measurement */
            }
        }
    }

    for( int i=0; i<N ; i++ ){                              /* Now P4k = I4_K4H4*P4k1 */
        for( int j=0 ;j<N ; j++ ){

            P4k[i][j] = 0 ;

            for( int k=0 ; k<N ; k++ ){                     /*Then add mat multiplication values on top*/
                P4k[i][j] += I4_K4H4[i][k]*P4k1[k][j] ;
            }
        }
    }
}

double record_Nx1_vector(double vector[N][1], int step, double record[N][kMax])
{

    for( int i=0 ; i<N ; i++ ){
        
        record[i][step] = vector[i][0] ;
    }
}

double record_NxN_diagonals(double P[N][N], int step, double record[N][kMax])
{

    for( int i=0 ; i<N ; i++ ){
        
        record[i][step] = P[i][i] ;
    }
}

double extract_row(double record[N][kMax], int row, double datout[kMax])
{
    for( int i=0 ; i<kMax ; i++ ){
        
        datout[i] = record[row][i] ;
    }
}