#include <iostream>
#include <math.h>
#include <fstream>
#include "time.h"
using namespace std;

int main(){
    int n;               //the number of grid
    cout<< "Insert n\n"; //insert n
    cin>>n;
    cout<<"n="<<n<<"\n";
    int num=n*n*n*2/3;
    cout<<"num="<<num<<"\n";   //number of operations

    float x1=0, xn=1;    //the boundaries of the interval (x1;xn)
    float h=(xn-x1)/(n+1);
    cout<<"h="<<h<<"\n";       //grid step

    cout<<"Start clock\n\n";
    clock_t start, finish;   //declaration of the clock
    start=clock();           //start point for the clock

    float **a=new float *[n];{ //declaration of an array of a(ij) elements
    for (int i=0; i<n; i++)
        a[i]=new float [n]; }

    cout<<"Array of a(ij)\n";
    for (int i=0; i<n; i++){       //fill out a(ij)
       for (int j=0; j<n; j++){
           if (i==j) a[i][j]=2/(h*h);
           else if (j==i-1) a[i][j]=-1/(h*h);
           else if (j==i+1) a[i][j]=-1/(h*h);
           else a[i][j]=0;
           cout<<a[i][j]<<" ";}
           cout<<"\n";}

    float *f=new float[n];{        //create array f(i) of right-side part of the set of equations
    //cout<<"\nMassive of elements f(i)\n";
    for (int i=0; i<n; i++){
         f[i]=100*exp(-10*(x1+(i+1)*h));
    //cout<<f[i]<<"\n";
    }}   //cout<<"\n";


//forward substitution
   for (int m=0; m<n-1; m++)          //change elements of a(ij), f(i);  m - number of a condition
     {for (int l=m+1; l<n; l++){       //rows
             //cout<<"m="<<m<<"\n";      //print all the temporary values in order to compare with the control matrix
             // cout<<"l="<<l<<"\n";
             // cout<<"f["<<l<<"]="<<f[l]<<"  ";
             // cout<<"f["<<m<<"]="<<f[m]<<"  ";
           double koef=a[l][m]/a[m][m];
            f[l]=f[l]-f[m]*koef;
            // cout<<"koef="<<koef<<"  ";         //print all the temporary values in order to compare with the control matrix
            // cout<<"koef*f["<<m<<"]="<<f[m]*koef<<"  ";
            // cout<<"f["<<l<<"]="<<f[l]<<"\n";
        for (int k=0; k<n; k++){
            a[l][k]=a[l][k]-a[m][k]*koef;
           // cout<<"["<<l<<"]["<<k<<"]="<<a[l][k]<<"  ";     //print all the temporary values in order to compare with the control matrix
           } //cout<<"\n ";
          }/*cout<<"\n ";*/
         }

   cout<<"\nChanged array a(ij)"<<endl;
    for (int i=0; i<n; i++){             //print changed a(ij) in order to check
        for (int j=0; j<n; j++){
             cout<<a[i][j]<<" ";
        }    cout<<"\n ";
    }

/*    cout<<"f(i) new\n"<<"";
    for (int i=0; i<n; i++){         //print changed f(i) in order to check
    cout<<f[i]<<"\n";} cout<<"\n";*/

//backward substitution
   float *result;             //declare an array of results
     result= new float [n];
   { for (int m=n-1; m>=0; m--) {
        double sum=0;
        //cout<<"\nm="<<m<<endl;
        {  for (int k=m+1; k<n; k++)
           {sum=sum+a[m][k]*result[k];
            //  cout<<"k="<<k<<endl;                          //print all the temporary values in order to compare with the control matrix
            //  cout<<"a["<<m<<"]["<<k<<"]"<<a[m][k]<<endl;
            //  cout<<"result["<<k<<"]="<<result[k]<<endl;
            }
            //  cout<<"sum="<<sum<<endl;

            result[m]=(f[m]-sum)/a[m][m];
               //cout<<"f["<<m<<"]"<<f[m]<<endl;             //print all the temporary values in order to compare with the control matrix
              // cout<<"a["<<m<<"]["<<m<<"]"<<a[m][m]<<endl;
              // cout<<"result["<<m<<"]="<<result[m]<<"\n";
        }

         }}
//     cout<<"\n";


     cout<<"\nResult"<<endl;
    for (int k=0; k<n; k++)                //print the result matrix in order to check
         cout<<"res["<<k<<"]="<<result[k]<< endl;


     cout<<"\nFinish clock\n\n";        //stop clock, calculate elapsed time
     finish=clock();
     ((finish - start)/CLOCKS_PER_SEC);
     cout<<(finish - start)/CLOCKS_PER_SEC<<"sec"<<endl;
     cout<<(finish - start)*1000/CLOCKS_PER_SEC<<"msec"<<endl;


     ofstream myfile;                    //write the results into a txt file in order to build a plot
     myfile.open ("Gauss_1000.txt");
     for (int i=0; i<n; i++)
         myfile <<x1+(i+1)*h<<" "<<result[i]<<endl;
         myfile.close();

    for (int i=0; i<n; i++) //free space
        delete [] a[i];
    delete [] a;
    delete [] f;
    delete [] result;

    return 0;
}
