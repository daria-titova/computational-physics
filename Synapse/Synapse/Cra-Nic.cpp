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

    int j=0;
    while (j<n)

   /* { V=inv(2*I+alpha*B)*(2*I-alpha*B)*V;
      V(0)=1.0; //initial condition
      j++;}*/

    {V=(2*I-alpha*B)*V;
     tridiag solve;
        solve.tridiag_solver(A, V, m);
        V(0)=1.0; //initial condition
    j++;}

    cout<<"Crank-Nicolson="<<endl;
    V.print();

return;}

