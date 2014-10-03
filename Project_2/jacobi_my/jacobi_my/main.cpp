#include <iostream>
#include <math.h>
#include <fstream>
#include "time.h"
#include <armadillo>
using namespace std;
using namespace arma;


void Jacobi(mat &a, mat &R, int n);

int main(){
    // declaration of variables
        //First we should insert the number of electron (e_number, 1 eller 2)
        //the task will be solved for; the program asks to insert e_number
        //until it is equal to 1 eller 2
        //
        int e_number;
        cout<< "Insert number of elecrons\n";
        cin>>e_number;
        while (e_number!=1 && e_number!=2)
        { cout<< "The task is solved for 1 or 2 electrons. Try again\n";
            cin>>e_number;}
        cout<<"number of electrons="<<e_number<<"\n\n";


        // define the number of steps, n,
        //the interval ro_min - ro_max will be divided into
        //
        int n;
        cout<< "Insert number os steps - n\n";
        cin>>n;
        cout<<"n="<<n<<"\n\n";

        //ro_min, ro_max - the begining and the end of interval,
        //ro_min=0 by default; ro_max we are able to choose
        float ro_min=0.0, ro_max;
        cout<< "Insert the end of the interval, ro_max\n";
        cin>>ro_max;
        cout<<"ro_max="<<ro_max<<"\n";
        float h=(ro_max-ro_min)/(n+1);    //define the step length of spacing, h
        cout<<"\nStep length h="<<h<<"\n"<<endl;

        //define oscillation frequency omega only for the case with two electrons
        //
        double omega;
        if(e_number==2) {
            cout<< "Insert oscillation frequency\n";
            cin>>omega;
            cout<<"oscillation frequerncy omega="<<omega<<"\n";}

        cout<<"Start clock\n";  //declare start and finish time
        clock_t start, finish;
        start=clock();          //clock starts

        //vector of values ro[i]
        //intermediate values of ro[i] of the interval ro_min - ro_max
        //
        vec ro(n);
                for (int i=0; i<n; i++)
                {ro(i)=ro_min+(i+1)*h;}
              //ro.print("\nro[i]:");


        //harmonic oscillator potential V(i)
        //is different for one and two electrons
        //
        vec V(n);
              for (int i=0; i<n; i++){
              if (e_number==1) {
                  V(i)=ro(i)*ro(i);}
              else{
                  V(i)=ro(i)*ro(i)*omega*omega+(1.0/ro(i));}}
             //V.print("\nHarmonic oscillator potential V[i]:");


        //declaration of the matrix of the eigenvectors R
        mat R(n,n);

        //creation of the array of a(ij) elements according to the task
        //
        mat a(n,n);
        for (int i=0; i<n; i++){
           for (int j=0; j<n; j++){
               if (i==j) a(i,j)=2/(h*h)+V(i);
               else if (j==i-1) a(i,j)=-1/(h*h);
               else if (j==i+1) a(i,j)=-1/(h*h);
               else a(i,j)=0;}}
        //a.print("\nArray a[n][n]: ");


        //call for the function 'jacobi'
        //which implements Jacobi's algorithm for searching
        //for eigenvectors and eigenvalues
        //
        Jacobi(a, R, n);
        //a.print("");


        //Print results:
        //the vector of eigenvalues of the matrix a is the main diagonal
        //of the transformed with Jacobi's method matrix a
        //lambda - vector of the eigenvalues of the matrix a
        //
        vec lambda(n);
        lambda=a.diag();
        //we should sort the eigenvalues calculated
        //in order to get the first three of them
        vec lambda_sort(n);
        lambda_sort=sort(lambda);
        //print only the first three values of the lambda
        cout<<"Lambda(0)="<<lambda_sort(0)<<endl<<"Lambda(1)="<<lambda_sort(1)<<endl<<"Lambda(2)="<<lambda_sort(2)<<endl;
        cout<<"Lambda(0)="<<lambda(0)<<endl<<"Lambda(1)="<<lambda(1)<<endl<<"Lambda(2)="<<lambda(2)<<endl;
        //a.print("\na:");
        //R.print("\nMatrix of the eigenvectors of the matrix a:");


        //find the index of the minimum eigenvalue in order to
        //choose the corresponding eigenvector
        //
        uword index_min;
        lambda.min(index_min);
        cout<<"\nindex_min="<<index_min<<endl;


        //In order to built wave function or probability function (wave*wave) we should
        //normalaze the eigenvectors
        //For that purpose we should extract a vector (column) from the matrix R,
        //which corresponds to chosen eigenvalue
        //find the sum of the elements^2 of that vector and devide all elements of the
        //vector by the root of this sum
        //
        vec U1=R.col(index_min);
        double sum;
        sum=0.0;
        for (int i=0; i<n; i++)
            sum=sum+U1(i)*U1(i);
        double norm = sqrt(sum);

        vec U(n);
        for (int i=0; i<n; i++)
            U(i)=U1(i)/norm;
        // U.print("");

        //The solution - wave function itself, normalised
        ofstream myfile_U; //write results into txt file in order to be able to build a plot
       myfile_U.open ("U(ro)_1.txt");
       for (int i=0; i<n; i++)
           myfile_U <<ro(i)<<" "<<U(i)<<endl;
       myfile_U.close();


        //probability = U*U
        vec W(n);
        for (int i=0; i<n; i++)
            W(i)=U(i)*U(i);
        //W.print("");

        ofstream myfile; //write results into txt file in order to be able to build a plot
       myfile.open ("W(ro)_1.txt");
       for (int i=0; i<n; i++)
           myfile <<ro(i)<<" "<<W(i)<<endl;
       myfile.close();


        //Use standart armadillo function eig_sym in order to find
        //the eigenvalues and eigenvectors of the matrix a
        //Thus we can check the results got with Jacobi method
        //and compare these two methods
 /*      mat eigvec;
         vec eigval;
         eig_sym(eigval, eigvec, a);
        //   cout<<"\neigenvec(A):\n"<<eigvec<<endl;
        //   cout<<"\neigenval:\n"<<eigval<<endl<<endl;
        cout<<eigval(0)<<endl<<eigval(1)<<endl<<eigval(2)<<endl;  */

    cout<<"\nFinish clock\n\n";
    finish=clock();                       //fix the finish moment
    ((finish - start)/CLOCKS_PER_SEC);    //calculate the execution time and print it
    cout<<(finish - start)/CLOCKS_PER_SEC<<"sec"<<endl;
    cout<<(finish - start)*1000/CLOCKS_PER_SEC<<"msec"<<endl;

return 0;}





// Function to find the maximum matrix element
//and to extract the indeces k and l in order to rotate around them
//
void Jacobi (mat &a, mat &R, int n ){              //send parametres to the function by address &

    double epsilon = 1.0e-8;                //the chosen tolerance
    double max_number_iterations = n*n*n;
    int iterations = 0;
    double max=1;        //should be more than epsilon by default in order to start the while loop

    //R is the matrix for the eigenvectors,
    //we start from setting the main diagonal og the R equal to 1
    //while the non-diagonal elements are equal to 0
    R.eye();
    //R.print("\nR:");

    while ( fabs(max) > epsilon && (double) iterations < max_number_iterations )

   { max = 0.0;
    // k, l - the indeces of the element for rotation
    int k, l;
    //s - sin, c - cos
    double s, c;

      for(int i = 0; i < n; i++ )
         {for(int j = i + 1; j < n; j++ )
          {if ( fabs(a(i,j)) > max )
              {max = fabs(a(i,j));
                 l = i;
                 k = j;}}}

      if ( a(k,l) != 0.0 ){
           double t, tau;
           tau = (a(l,l) - a(k,k))/(2*a(k,l));
              if( tau > 0 )
                   t = 1.0/(tau + sqrt(1.0 + tau*tau));
              else t = -1.0/(-tau + sqrt(1.0 + tau*tau));
           c = 1/sqrt(1+t*t);
           s = c*t;}
      else {c = 1.0; s = 0.0;}


      double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
      a_kk = a(k,k);
      a_ll = a(l,l);

 // changing the matrix elements with indices k and l
      a(k,k) = c*c*a_kk - 2.0*c*s*a(k,l) + s*s*a_ll;
      a(l,l) = s*s*a_kk + 2.0*c*s*a(k,l) + c*c*a_ll;
      a(k,l) = 0.0;    // hard-coding of the zeros
      a(l,k) = 0.0;

 //changing the remaining elements
      for (int i = 0; i < n; i++ )
          {if ( i != k && i != l )
          {a_ik = a(i,k);
              a_il = a(i,l);
              a(i,k) = c*a_ik - s*a_il;
              a(k,i) = a(i,k);
              a(i,l) = c*a_il + s*a_ik;
              a(l,i) = a(i,l);}

          // Computing the new eigenvectors
          r_ik = R(i,k);
          r_il = R(i,l);
          R(i,k) = c*r_ik - s*r_il;
          R(i,l) = c*r_il + s*r_ik;}

    iterations++;}

    cout<<"Number of iterations="<<iterations<<endl;

return;}





