#include <Verlet.h>
#include <celestial_body.h>
#include <iostream>
#include <math.h>
#include <solar_system.h>
#include <fstream>

using namespace std;

Verlet::Verlet()
{
}

void Verlet::verlet_integration (Solar_System &system, double dt, int n)
{   //first of all we should get vx1, vy1, x1, y1
    //uisng euler method

    double AU=1.5e11;
    double year= 31557600;
    double M_sun=2e30;
    double G=6.674e-11;
    double G_star=G/((AU*AU*AU)/(M_sun*year*year));

    vec3 v1, x1;

    vec3 x1_1=system.objects[1].position;
    vec3 x1_2=system.objects[1].velocity*dt;
    x1=x1_1 + x1_2;
    for (int i=0; i<3; i++) cout<<"vec x1="<<x1[i]<<endl;

    vec3 v1_1=system.objects[1].velocity;
    cout<<v1_1[0]<<"  "<<v1_1[1]<<"  "<<v1_1[2]<<endl;
    vec3 v1_2= (x1 + system.objects[1].position)/pow(system.distance_between_objects(system.objects[0], system.objects[1]), 3);
    cout<<v1_2[0]<<"  "<<v1_2[1]<<"  "<<v1_2[2]<<endl;
    v1_2=v1_2*G_star;
    cout<<v1_2[0]<<"  "<<v1_2[1]<<"  "<<v1_2[2]<<endl;
    v1_2=(v1_2/2)*dt;
    cout<<v1_2[0]<<"  "<<v1_2[1]<<"  "<<v1_2[2]<<endl;
    v1=v1_1 - v1_2;
    cout<<v1[0]<<"  "<<v1[1]<<"  "<<v1[2]<<endl;

    double** velocity= new double*[n];
    for (int i=0; i<n; i++){
        velocity[i]= new double [3];}

    double** position= new double*[n];
    for (int i=0; i<n; i++){
        position[i]= new double [3];}


    //force


    for (int i=0; i<n; i++) {
        if (i==0) { for (int j=0; j<3; j++)
            {velocity[i][j]=system.objects[1].velocity[j];
            position[i][j]=system.objects[1].position[j];
            }
            cout<<"pos ["<<i<<"]"<<"["<<1<<"]="<<position[i][1]<<"   ";
            cout<<"vel ["<<i<<"]"<<"["<<1<<"]="<<velocity[i][1]<<"";
            cout<<endl;
            }
        else {
            if (i==1)  { for (int j=0; j<3; j++)
                {velocity[i][j]=v1[j];
                    position[i][j]=x1[j];}
                cout<<"pos ["<<i<<"]"<<"["<<1<<"]="<<position[i][1]<<"   ";
                cout<<"vel ["<<i<<"]"<<"["<<1<<"]="<<velocity[i][1]<<"";
                cout<<endl;
                }

            else { for (int j=0; j<3; j++)
                    {
                    vec3 pos;
                    for (int k=0; k<3; k++) {pos[k]=position[i-1][k];}
                    position[i][0]=position[i-1][0] + dt*velocity[i-1][0] - (1/2)*dt*dt*(position[i-1][0]*G_star)/1;
                    position[i][1]=position[i-1][1] + dt*velocity[i-1][1] + (1/2)*dt*dt*(position[i-1][1]*G_star)/system.distance_between_objects(system.objects[0], system.objects[1]);
                    position[i][2]=0;
                    velocity[i][j]=velocity[i-1][j] - (dt/2)*((position[i][j] + position[i-1][j])*G_star/system.distance_between_objects(system.objects[0], system.objects[1]));


                }
               // cout<<"pos ["<<i<<"]"<<"["<<0<<"]="<<position[i][0]<<"   ";
               // cout<<"vel ["<<i<<"]"<<"["<<1<<"]="<<velocity[i][1]<<"   ";
                    cout<<endl;

                 }

              }
     // cout<<"coordinates ["<<i<<"]"<<"["<<j<<"]"<<position[i][j]<<endl;
      }

  //  cout<<"new G="<<G_star<<endl;

    ofstream myfile; //write results into txt file in order to be able to build a plot
    myfile.open ("vel_x_1000.txt");  //the wave function
    for (int i=0; i<n; i++)
       myfile <<0+i*dt<<" "<<velocity[i][0]<<endl;
       myfile.close();


   return;
}
