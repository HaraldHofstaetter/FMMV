// #ifdef USE_SINGLE_PRECISION
// typedef float _FLOAT_;
// #else		
// typedef double _FLOAT_;
// #endif	

_FLOAT_ relErr(_FLOAT_ x,_FLOAT_ x_ex);
_FLOAT_ relErr2(_FLOAT_ x[2],_FLOAT_ x_ex[2]);
_FLOAT_ relL2Err2(int flag, _FLOAT_ x[2],_FLOAT_ x_ex[2]);
_FLOAT_ relErr3(_FLOAT_ x[3],_FLOAT_ x_ex[3]);
_FLOAT_ relL2Err3(int flag, _FLOAT_ x[3],_FLOAT_ x_ex[3]);
_FLOAT_ relL2Err(int flag, _FLOAT_ x,_FLOAT_ x_ex);
void my_srand(unsigned int SEED);
unsigned int my_rand(void);
#define MY_RAND_MAX 4294967295UL

void set_from_command_line_int(int argc, char*argv[], char* varname, int *var);
void set_from_command_line_float(int argc, char*argv[], char* varname, float *var);
void set_from_command_line_double(int argc, char*argv[], char* varname, double *var);
void set_from_command_line_bool(int argc, char*argv[], char* varname, int *var);

