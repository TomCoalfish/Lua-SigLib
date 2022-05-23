%module siglib
%{
#include "siglib.h"
%}

%include "stdint.i"
%include "std_vector.i"
%include "std_math.i"

#define SIGLIB_FUNC_DECL
#define SIGLIB_OUTPUT_PTR_DECL
#define SIGLIB_INPUT_PTR_DECL
#define SIGLIB_HUGE_DECL
%ignore SUF_Debugvfprintf;
%ignore SUF_WavReadInt;
%ignore SUF_WavWriteInt;
%include "siglib.h"