#include <Closed_form.h>
#include <iostream>
#include <math.h>
#include <armadillo>
using namespace arma;
using namespace std;

Closed_form::Closed_form()
{
}

void Closed_form::Closed_form_solution(int n, int m, double t, double dx)
{
    vec U(m);
    int k=n;


    for (int i=0; i<m; i++){
        double sum=0.0;
        for (int j=0; j<k-1; j++) {sum=sum+(2/((j+1)*M_PI))*sin((j+1)*M_PI*i*dx)*exp(-M_PI*M_PI*(j+1)*(j+1)*t);
        }
        U(i)=1.0-i*dx-sum;
        }

    U(m-1)=0.0;     //for the sace of printing
    U.print("Closed_form=");

   return;
}
