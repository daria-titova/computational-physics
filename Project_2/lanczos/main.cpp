#include <iostream>
#include <math.h>
#include <fstream>
#include "time.h"
#include <armadillo>
using namespace std;
using namespace arma;


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
    //
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
    start=clock();            //clock starts

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
    // V.print("\nHarmonic oscillator potential V[i]:");


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


    //In order to find minimum eigenvalues for the matrix a
    //we will find the maximum ones for the inversed a
    //and then just inverse them
    //thus we need inversed a
    //
    mat a_inv=inv(a);
    //a_inv.print();


    //the matrix Q of the orthonornal vectors q_i
    //using to get the eigenvectors for a from the ones for T
    //q_0=0, q_1=1.
    //
    mat q(n,n+1);
    vec q0=q.col(0);
    q0.fill(0.0);
    q.col(0)=q0;
    vec q1=q.col(1);
    q1.fill(sqrt(1.0/n));
    q.col(1)=q1;

    mat r(n,n+1);
    r.col(0)=q1;

    //elements of the matrix T
    vec alfa(n+1);
    vec beta(n+1);
    beta(0)=1;

    mat I(n,n);
    I.eye();

    int iterations=0;

    int k=0;


    //the tolerance which is used to choose m
    //if the tolerance is equal to zero, m=>n,
    //which doesn't help much to increase the speed of the calculation
    //
    double epsilon=0.001;
    //cout<<"Starting while\n"<<endl;
        while (k!=n && fabs(beta(k))>epsilon){
        vec q_k_plus_1=q.col(k+1);
        vec r_k=r.col(k);
        q_k_plus_1=r_k/beta(k);
        q.col(k+1)=q_k_plus_1;

        k=k+1;

        vec q_k=q.col(k);
        mat q_k_t(1,n);
        q_k_t=trans(q_k);

        mat part1(1,n);
        part1= (q_k_t)*(a_inv);
        vec part2 = part1*q_k;
        //although the multiplication (q_k_t)*(a_inv)part1*q_k gives a number in paper,
        //here it gives a vec and we should just chose its first element
        alfa(k)=part2(0,0);

        vec q_k_minus_1=q.col(k-1);
        r_k=(a_inv-alfa(k)*I)*q_k-beta(k-1)*q_k_minus_1;

        double sum=0.0;   //normalizing
        vec r_star=r_k;
        for (int i=0; i<n; i++)
        {sum=sum+r_star(i)*r_star(i);}
        beta(k)=sqrt(sum);

        r.col(k)=r_k;
        iterations+=1;
    }

    //cout<<"\n\nEnd while"<<endl;

    int m;    //the size of the matrix T defined by beta>epsilon
    if (iterations==n) m=n;
    else m=iterations+1;
    cout<<"m="<<m<<endl;

    //elements of -/+1 diagonals of T
    vec beta_1(m-1);
    for (int i=0; i<m-1; i++)
    {beta_1(i)=beta(i+1);}

    //elements of the main diagonal of T
    vec alfa_1(m);
    for (int i=0; i<m; i++)
    {alfa_1(i)=alfa(i+1);}

    mat T(m,m);
    T.eye();
    //cout << "T dimension: " << T.n_rows << ", " << T.n_cols << endl;

    T.diag(0)=alfa_1;
    T.diag(-1)=beta_1;
    T.diag(1)=beta_1;
    //cout<<"\nT\n"<<T<<endl;

    cout<<"iterations="<<iterations<<endl;


    //to solve the eigenproblem for T we use
    //standard function from armadillo
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, T);

    //the found eigenvalues are the eigenvalues of the a inversed
    //now we should sort them and find the maximum three
    //and also get their indeces to be able to find
    //eigenvectors
    double max_1=0.0, max_2=0.0, max_3=0.0;
    int i_max_1=0, i_max_2=0, i_max_3=0;
    double eps=1.0e-5;   //use it to compare double

    //cout<<"\nStartimg_max"<<endl;
    for (int i=0; i<m; i++){
        if (eigval(i)>max_1 && fabs(eigval(i)-max_1)>eps){
        max_3=max_2;
        max_2=max_1;
        max_1=eigval(i);
        i_max_3=i_max_2;
        i_max_2=i_max_1;
        i_max_1=i;
        }
        else {
            if ((eigval(i)>max_2) && fabs(eigval(i)-max_1)>eps && fabs(eigval(i)-max_2)>eps) {
                max_3=max_2;
                max_2=eigval(i);
                i_max_3=i_max_2;
                i_max_2=i;}
            else {if ((eigval(i)>max_3) && fabs(eigval(i)-max_1)>eps && fabs(eigval(i)-max_2)>eps && fabs(eigval(i)-max_3)>eps)
                    max_3=eigval(i);
            }
        }}

   //cout<<"\nEnding_max"<<endl;

   // cout<<"max_1="<<max_1<<endl;
   // cout<<"max_2="<<max_2<<endl;
   // cout<<"max_3="<<max_3<<endl;

    //the vector of the first three minimum eigenvalues
    //of the matrix a
    vec lambda(3);
    lambda(0)=1/max_1;
    lambda(1)=1/max_2;
    lambda(2)=1/max_3;
    lambda.print("Lambda:");


    mat Q(n,m);
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            Q(i,j)=q(i, j+1);}}
    //Q.print();

    //the first three eigenvectors of the matrix a
    //corresponding to the three minimum eigenvalues
    //of a
    mat eigvec_1(n,3);
    eigvec_1.col(0)=Q*eigvec.col(i_max_1);
    eigvec_1.col(1)=Q*eigvec.col(i_max_2);
    eigvec_1.col(2)=Q*eigvec.col(i_max_3);
    //cout<<"i_1="<<i_max_1<<endl<<"i_2="<<i_max_2<<endl<<"i_3="<<i_max_3<<endl;
    // eigvec_1.print("eigvec_1\n");

    ofstream myfile; //write results into txt file in order to be able to build a plot
    myfile.open ("U(ro)_1.txt");  //the wave function
    for (int i=0; i<n; i++)
       myfile <<ro(i)<<" "<<eigvec_1(i,0)<<endl;
       myfile.close();

    cout<<"\nFinish clock\n\n";
    finish=clock();                       //fix the finish moment
    ((finish - start)/CLOCKS_PER_SEC);    //calculate the execution time and print it
    cout<<(finish - start)/CLOCKS_PER_SEC<<"sec"<<endl;
    cout<<(finish - start)*1000/CLOCKS_PER_SEC<<"msec"<<endl;

    return 0;}
