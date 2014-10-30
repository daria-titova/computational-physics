#include <iostream>
#include <math.h>
#include <Cra-Nic.h>
#include <tridiag.h>
#include <armadillo>
using namespace std;
using namespace arma;

Crank_Nicolson::Crank_Nicolson(){}

void Crank_Nicolson::Crank_Nicolson_Scheme(mat &U, double alpha, int n, int m)
{ vec V(m);
    for (int i=0; i<m; i++)
    V(i)=U(0,i);    //initial condition, t=0;

  mat I(m,m);
    I.eye();

  mat B(m,m);
    B.diag(-1).fill(-1);
    B.diag(1).fill(-1);
    B.diag().fill(2);

    mat A(m,m);
    A=2*I+alpha*B;

    vec V_new(m);
    vec V1=V;

    int j=1;
    while (j<n)

   /* { V=inv(2*I+alpha*B)*(2*I-alpha*B)*V;
      V(0)=1.0; //boundary condition for x=0;
      j++;}*/

    {   // V=(2*I-alpha*B)*V;
         //V.print();

        for (int i=0; i<m; i++)
        {if (i==0) V_new(i)=(2.0 - 2.0*alpha)*V(i) + alpha*V(i+1);
            else{ if (i<(m-1)) V_new(i)=alpha*V(i-1) + (2.0 - 2.0*alpha)*V(i) + alpha*V(i+1);
                 else V_new(i)=alpha*V(i-1) + (2.0 - 2.0*alpha)*V(i);}
            } V=V_new;
           //V.print();

     tridiag solve;
     solve.tridiag_solver(A, V, m);
        V(0)=1.0;   //boundary condition for x=0;
        V(m-1)=0.0; //boundary condition for x=L=1;
       // V.print();
    j++;}

    V.print("Crank-Nicolson=");

return;}

