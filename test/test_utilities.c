#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#ifdef USE_SINGLE_PRECISION
typedef float _FLOAT_;
#else
typedef double _FLOAT_;
#endif
#include"test_utilities.h"

_FLOAT_ relErr(_FLOAT_ x,_FLOAT_ x_ex)
{
   return fabs((x-x_ex)/x_ex);
}

_FLOAT_ relErr2(_FLOAT_ x[2],_FLOAT_ x_ex[2])
{
   return sqrt(((x[0]-x_ex[0])*(x[0]-x_ex[0])+
   	        (x[1]-x_ex[1])*(x[1]-x_ex[1]))
		/(x_ex[0]*x_ex[0]+x_ex[1]*x_ex[1]));   
}

_FLOAT_ relL2Err2(int flag, _FLOAT_ x[2],_FLOAT_ x_ex[2])
{
	static int N;
	static _FLOAT_ nom, den;

	switch (flag) {
	case 0:	/*init*/
		N = 0;
		nom = 0; 
		den = 0;
		return 0.0;
	case 1: /*accumulate*/
		nom += (x[0]-x_ex[0])*(x[0]-x_ex[0])+
   	               (x[1]-x_ex[1])*(x[1]-x_ex[1]);
		den += (x_ex[0])*(x_ex[0])+
   	               (x_ex[1])*(x_ex[1]);
		return 0.0;
	}	
	/* finalize */
	return sqrt(nom/den);
}	

_FLOAT_ relErr3(_FLOAT_ x[3],_FLOAT_ x_ex[3])
{
   return sqrt(((x[0]-x_ex[0])*(x[0]-x_ex[0])+
   	        (x[1]-x_ex[1])*(x[1]-x_ex[1])+
   	        (x[2]-x_ex[2])*(x[2]-x_ex[2]))
		/(x_ex[0]*x_ex[0]+x_ex[1]*x_ex[1]+x_ex[2]*x_ex[2]));   
}

_FLOAT_ relL2Err3(int flag, _FLOAT_ x[3],_FLOAT_ x_ex[3])
{
	static _FLOAT_ nom, den;

	switch (flag) {
	case 0:	/*initialize*/
		nom = 0; 
		den = 0;
		return 0.0;
	case 1: /*accumulate*/
		nom += (x[0]-x_ex[0])*(x[0]-x_ex[0])+
   	               (x[1]-x_ex[1])*(x[1]-x_ex[1])+
   	               (x[2]-x_ex[2])*(x[2]-x_ex[2]);
		den += (x_ex[0])*(x_ex[0])+
   	               (x_ex[1])*(x_ex[1])+
   	               (x_ex[2])*(x_ex[2]);
		return 0.0;
	}	
	/* finalize */
	return sqrt(nom/den);
}	


_FLOAT_ relL2Err(int flag, _FLOAT_ x,_FLOAT_ x_ex)
{
	static _FLOAT_ nom, den;

	switch (flag) {
	case 0:	/*initialize*/
		nom = 0; 
		den = 0;
		return 0.0;
	case 1: /*accumulate*/
		nom += (x-x_ex)*(x-x_ex);
		den += (x_ex)*(x_ex);
		return 0.0;
	}	
	/* finalize */
	return sqrt(nom/den);
}	
		

unsigned int IDUM;

void my_srand(unsigned int SEED)
{
	IDUM = SEED;
}	
	
unsigned int my_rand(void)
{
	IDUM = 1664525L*IDUM + 1013904223L;
	return IDUM;
}	

#define MY_RAND_MAX 4294967295UL

void set_from_command_line_int(int argc, char*argv[], char* varname, int *var)
{
    int k;
    char *lhs;
    int var_old;

    for (k=1; k<argc; k++) {
        lhs = strchr(argv[k], '=');
        if (lhs) {
            *lhs ='\0';
            if (strcmp(argv[k], varname)==0) {
                  var_old = *var;
                  *var = atoi(lhs+1);
                  *lhs = '=';
                  printf("changed %s from %i to %i\n", varname, var_old, *var);
                  return;
            }
            *lhs = '=';
        }   
    }
}

void set_from_command_line_float(int argc, char*argv[], char* varname, float *var)
{
    int k;
    char *lhs;
    float var_old;

    for (k=1; k<argc; k++) {
        lhs = strchr(argv[k], '=');
        if (lhs) {
            *lhs ='\0';
            if (strcmp(argv[k], varname)==0) {
                  var_old = *var;
                  *var = atof(lhs+1);
                  *lhs = '=';
                  printf("changed %s from %g to %g\n", varname, var_old, *var);
                  return;
            }
            *lhs = '=';
        }   
    }
}

void set_from_command_line_double(int argc, char*argv[], char* varname, double *var)
{
    int k;
    char *lhs;
    double var_old;

    for (k=1; k<argc; k++) {
        lhs = strchr(argv[k], '=');
        if (lhs) {
            *lhs ='\0';
            if (strcmp(argv[k], varname)==0) {
                  var_old = *var;
                  *var = atof(lhs+1);
                  *lhs = '=';
                  printf("changed %s from %g to %g\n", varname, var_old, *var);
                  return;
            }
            *lhs = '=';
        }   
    }
}

void set_from_command_line_bool(int argc, char*argv[], char* varname, int *var)
{
    int k;

    for (k=1; k<argc; k++) {
        if (strcmp(argv[k], varname)==0) {
             *var = 1;
             printf("option '%s'\n", varname);
             return;
        }   
    }
}



/*
int main(int argc, char *argv[]) 
{
   int test_int = 123;
   float test_float = 1.23;
   double test_double = 12.3;
   int test_bool = 0;
   
   set_from_command_line_int(argc, argv, "test_int", &test_int);
   set_from_command_line_float(argc, argv, "test_float", &test_float);
   set_from_command_line_double(argc, argv, "test_double", &test_double);
   set_from_command_line_bool(argc, argv, "test_bool", &test_bool);
   printf("test_int    = %i\n", test_int);
   printf("test_float  = %g\n", test_float);
   printf("test_double = %g\n", test_double);
   printf("test_bool = %i\n", test_bool);
   return 0;
}
*/ 
