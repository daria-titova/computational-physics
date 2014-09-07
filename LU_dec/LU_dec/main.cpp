#include <iostream>
#include <armadillo>
#include <math.h>
using namespace std;
using namespace arma;

int main () {
    int n;            //declare n - grid size, interval (x1, xn)
    float x1=0, xn=1;

    cout<< "Insert n\n";     //read n from terminal
    cin>>n;
    cout<<"n="<<n<<"\n";

    cout<<"Start clock\n\n";  //declare clock, fix start time of calculation
    clock_t start, finish;
    start=clock();

    float h=(xn-x1)/(n+1);    //calculate step h
    cout<<"h="<<h<<"\n\n";

    fmat A(n,n);
    for (int i=0; i<n; i++){       //creates array of a(ij) elements
       for (int j=0; j<n; j++){
           if (i==j) A(i,j)=2/(h*h);
           else if (j==i-1) A(i,j)=-1/(h*h);
           else if (j==i+1) A(i,j)=-1/(h*h);
           else A(i,j)=0;
           //cout<<"array[" <<i<< "][" <<j<< "]="<<A(i,j)<<" ";
           cout<<A(i,j)<<" ";
       }   cout<<"\n";
    }

/*    for (int i=0; i<n; i++){       //prints array of a(ij) elements
       for (int j=0; j<n; j++)
           cout<<A(i,j)<<" ";
           cout<<"\n";}*/

    fvec F(n);                                //creates array/vector of f(i) elements (0 - n)
//    cout<<"\nMassive of elements fi\n";
    for (int i=0; i<n; i++){
         F(i)=100*exp(-10*(x1+(i+1)*h));
//    cout<<F(i)<<"\n";
    }
//    cout<<"\n";                          //print f(i)

//LU decomposition
    fmat L, U, P;
    lu(L, U, P, A);
      //    cout<<"A:"<<" "<<A<<endl;cout<<"\n";      //print matrices got from LU decomposition in order to check with the control ones
      //    cout<<"L:"<<" "<<L<<endl;cout<<"\n";
      //    cout<<"U:"<<" "<<U<<endl;cout<<"\n";
      //    cout<<"P:"<<" "<<P<<endl;cout<<"\n";

    fmat Y=inv(L)*F;                  //Create auxiliary matrix Y
    //cout<<"Y:\n"<<""<<Y<<endl;
    fmat Res = inv(U)*Y;             //Find solution/result
    cout<<"\nResult:\n"<<""<<Res<<endl;

//Calculate results with solve function
            /* vec Result = solve(A, F);
             cout<<"Result:\n"<<""<<Result<<endl;cout<<"\n"; */  //print results

    cout<<"Finish clock\n\n";             //finish time, calculate elapsed time
    finish=clock();
    ((finish - start)/CLOCKS_PER_SEC);
    cout<<(finish - start)/CLOCKS_PER_SEC<<"sec"<<endl;   //in seconds
    cout<<(finish - start)*1000/CLOCKS_PER_SEC<<"msec"<<endl;  //in milliseconds

        ofstream myfile;             //Write the result into txt file
        myfile.open ("LU_100.txt");
        for (int i=0; i<n; i++)
            myfile <<x1+(i+1)*h<<" "<<Res(i)<<endl;
            myfile.close();

    return 0;
}
