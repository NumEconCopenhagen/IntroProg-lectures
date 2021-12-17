//////////////////////////
// 1. external includes //
//////////////////////////

// standard C++ libraries
#include <windows.h>
#include <stdio.h>
#include <cmath>
#include <omp.h>

///////////////
// 2. macros //
///////////////

#define EXPORT extern "C" __declspec(dllexport)

/////////////
// 3. main //
/////////////

EXPORT void setup_omp(){ // required when using visual studio
    SetEnvironmentVariable("OMP_WAIT_POLICY", "passive"); 
}

EXPORT void myfun_cpp(double *x1, double *x2, double *y, int N1, int N2, int threads){

    #pragma omp parallel num_threads(threads)
    {

    #pragma omp for      
    for(int i = 0; i < N1; i++){
        y[i] = 0.0;
        if(x1[i] < 0.5){
            for(int j = 0; j < N2; j++){
                y[i] += exp(x2[j]*x1[i]);
            }
        } else {
            for(int j = 0; j < N2; j++){
                y[i] += log(x2[j]*x1[i]);
            }
        }
    }

    } // omp parallel
    
}