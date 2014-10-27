#include <Implicit.h>
#include <iostream>
#include <math.h>
#include <armadillo>
#include <tridiag.h>
using namespace arma;
using namespace std;

Implicit::Implicit()
{
}

void Implicit::Implicit_Scheme (mat &U, double alpha, int n, int m)
{
    mat A(m,m);
    A.zeros();
    A.diag().fill(1+2*alpha);
    A.diag(1).fill(-alpha);
    A.diag(-1).fill(-alpha);
   // A.print();
   // A=inv(A);

    vec V(m);
    for (int i=0; i<m; i++)
    V(i)=U(0,i);    //initial condition, t=0;
   // V.print();

    int j=0;
    while (j<n){
       // V=inv(A)*V;
        tridiag solve;
        solve.tridiag_solver(A, V, m);
        V(0)=1.0;      //boundary condition for x=0;
       // V(m-1)=0.0;  //boundary condition for x=L=1;
        j++;}
    cout<<"Implicit="<<endl;
    V.print();

   return;
}
