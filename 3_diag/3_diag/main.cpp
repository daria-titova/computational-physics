#include <iostream>
#include <math.h>
#include <fstream>
#include "time.h"
using namespace std;

int main(){
// declaration of variables
    int n;                    // n-size of grid
    cout<< "Insert n\n";      //insert n, read it
    cin>>n;
    cout<<"n="<<n<<"\n";

    float x1=0, xn=1;         //x1, xn - the begining and the end of interval
    float h=(xn-x1)/(n+1);    //declare and identify a step length of spacing - h
    cout<<"h="<<h<<"\n"<<endl;

    cout<<"Start clock\n\n";  //declare start and finish time
    clock_t start, finish;
    start=clock();            //clock starts

    float **a=new float *[n];{   //pointer on pointer, create space for array of pointers on rows
    for (int i=0; i<n; i++)      //loop for creating space for every row of the array
        a[i]=new float [n]; }    // each element of the array of pointers on rows
                                    //gets address of the beginning of the space created for row, n - columns

    cout<<"\nArray a[n][n]: "<<endl;
    for (int i=0; i<n; i++){       //creation of array of a(ij) elements according to the task
       for (int j=0; j<n; j++){
           if (i==j) a[i][j]=2/(h*h);
           else if (j==i-1) a[i][j]=-1/(h*h);
           else if (j==i+1) a[i][j]=-1/(h*h);
           else a[i][j]=0;
//          cout<<"array[" <<i<< "][" <<j<< "]="<<a[i][j]<<" ";
       }   //printing the array of a(ij) element by element for easy checking
          cout<<"\n";}

    float *f=new float[n];{        //creation of array of f(i) elements according to the task
        cout<<"\nArray of elements f(i)"<<endl;
    for (int i=0; i<n; i++){
         f[i]=100*exp(-10*(x1+(i+1)*h));
    cout<<f[i]<<"\n"; }}          //print the array f(i)


    float *a1=new float[n];    //defining of array a1(i) from array a(ij)
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++)
            a1[i]=a[i][i-1];}
    for (int i=0; i<n; i++){
    // cout<<"Array of a1(i)"<<a1[i]<<endl;  //print a1(i)
    }

    float *b=new float[n];    //defining of array b(i) from array a(ij)
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++)
            b[i]=a[i][i];}
    for (int i=0; i<n; i++){
    // cout<<"Array of b(i)"<<b[i]<<endl; //print the array b(i)
    }

    float *c=new float[n];    //defining of array c(i) from array a(ij)
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++)
            c[i]=a[i][i+1];}
    for (int i=0; i<n; i++){
    //cout<<"Array of c(i)"<<c[i]<<endl;  //print the array c(i)
    }

//forward substitution
    float *c1=new float[n];      //change the array c(i) into c1(i) according to the method for tridiagonal matrix
    for (int i=0; i<n-1; i++){
            if (i==0) c1[i]=c[i]/b[i];
            else c1[i]=c[i]/(b[i]-a1[i]*c1[i-1]);}
    for (int i=0; i<n; i++){
    //cout<<"Array of c1(i)"<<c1[i]<<endl;  //print new c1(i)
    }

    float *f1=new float[n];    //change the array f(i) into f1(i) according to the method for tridiagonal matrix
    for (int i=0; i<n; i++){
            if (i==0) f1[i]=f[i]/b[i];
            else f1[i]=(f[i]-a1[i]*f1[i-1])/(b[i]-a1[i]*c1[i-1]);}
    for (int i=0; i<n; i++){
    // cout<<"Array of f1(i)"<<f1[i]<<endl;   //print new f1(i)
    }

//backward substitution
    float *res=new float[n];         //declare array of results - res
    for (int i=n-1; i>=0; i--){      //solve the set, find the solution
            if (i==n-1) res[i]=f1[i];
            else res[i]=f1[i]-c1[i]*res[i+1];}

    cout<<"\nResult of LU"<<endl;
    for (int i=0; i<n; i++)        //print array of the results
    cout<<res[i]<<endl;

    float *u=new float [n];                 //closed-form function
    //cout<<"Closed-form solution u(i)\n"<<endl;
    for (int i=0; i<n; i++){
         u[i]=1-(1-exp(-10*(x1+(i+1)*h)))*(x1+(i+1)*h)-exp(-10*(x1+(i+1)*h));
    //cout<<u[i]<<"\n";                //print array of results for closed-form function
    }

    float *e=new float [n];            //computing relative error
    //cout<<"Relative error\n"<<endl;
    for (int i=0; i<n; i++){
        e[i]=log10(fabs((res[i]-u[i])/u[i]));
    //cout<<e[i]<<"\n";
    }

    cout<<"\nMax relative error for n="<<n<<endl;   //calculation of max relative error
    float maxe=e[0];
    for (int i=0; i<n; i++)
        if (e[i]>maxe)  maxe=e[i];     //compare med maxe
        cout<<maxe<<"\n";             //print max relative error


    cout<<"\nFinish clock\n\n";
    finish=clock();                       //fix the finish moment
    ((finish - start)/CLOCKS_PER_SEC);    //calculate the execution time and print it
    cout<<(finish - start)/CLOCKS_PER_SEC<<"sec"<<endl;
    cout<<(finish - start)*1000/CLOCKS_PER_SEC<<"msec"<<endl;

     ofstream myfile;              //write results into txt file in order to be able to build a plot
     myfile.open ("3_diag_10000.txt");
     for (int i=0; i<n; i++)
         myfile <<x1+(i+1)*h<<" "<<res[i]<<endl;
         myfile.close();

    for (int i=0; i<n; i++)   //free space
        delete [] a[i];
    delete [] a;
    delete [] f;
    delete [] res;
    delete [] f1;
    delete [] a1;
    delete [] b;
    delete [] c;
    delete [] c1;

    return 0;
}
