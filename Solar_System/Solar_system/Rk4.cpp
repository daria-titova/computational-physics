#include <Rk4.h>
#include <celestial_body.h>
#include <iostream>
#include <math.h>
#include <solar_system.h>
#include <fstream>
using namespace std;


RK4::RK4()
{
}


void RK4::integrate_RK4(Solar_System &system, double dt, int n)
{// Do RK4 integration here

    double** velocity= new double*[n];
    for (int i=0; i<n; i++){
        velocity[i]= new double [3];}

    double** position= new double*[n];
    for (int i=0; i<n; i++){
        position[i]= new double [3];}

    double *kinetik_energy=new double[n];
    double *potential_energy=new double[n];
    double *momentum=new double[n];

    vec3 v;

    vec3  k1_vel, k2_vel, k3_vel, k4_vel;
    vec3  k1_pos, k2_pos, k3_pos, k4_pos;

    for (int i=0; i<n-1; i++) {
        if (i==0)
        {
            kinetik_energy[i]=system.E_kinetik(system.objects[1], system.objects[1].velocity.length());
            cout<<"Kinetic_energy="<<kinetik_energy[i]<<endl;

                for (int j=0; j<3; j++){

                velocity[i][j]=system.objects[1].velocity[j];
                position[i][j]=system.objects[1].position[j];

                k1_pos[j]=(system.objects[1].velocity[j])*dt;
                k1_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j])/system.objects[1].mass;

                k2_pos[j]=(system.objects[1].velocity[j]+k1_vel[j]*dt/2);
              //  k2_vel[j]=(-1)*(system.objects[1].position[j]*G_star+k1_pos[j]*dt/2);
                k2_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j]+k1_pos[j]*dt/2)/system.objects[1].mass;

                k3_pos[j]=(system.objects[1].velocity[j]+k2_vel[j]*dt/2);
              //  k3_vel[j]=(-1)*(system.objects[1].position[j]*G_star+k2_pos[j]*dt/2);
                k3_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j]+k2_pos[j]*dt/2)/system.objects[1].mass;

                k4_pos[j]=(system.objects[1].velocity[j]+k3_vel[j]*dt);
              //  k4_vel[j]=(-1)*(system.objects[1].position[j]*G_star+k3_pos[j]*dt);
                k4_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j]+k3_pos[j]*dt)/system.objects[1].mass;

                position[i+1][j]=position[i][j] + (k1_pos[j] + k2_pos[j]*2 + k3_pos[j]*2 + k4_pos[j])*dt/6;
                velocity[i+1][j]=velocity[i][j] + (k1_vel[j] + k2_vel[j]*2 + k3_vel[j]*2 + k4_vel[j])*dt/6;
                v[j]=velocity[i+1][j];
            }

                kinetik_energy[i+1]=system.E_kinetik (system.objects[1], v.length());
                cout<<"Kinetic_enerdy="<<kinetik_energy[i+1]<<endl;

        }

        else {for (int j=0; j<3; j++){
                k1_pos[j]=(velocity[i][j]);
                //k1_vel[j]=position[i][j]*(-G_star);
                k1_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j])/system.objects[1].mass;

                k2_pos[j]=(velocity[i][j]+k1_vel[j]*dt/2);
              //  k2_vel[j]=(position[i][j]+k1_pos[j]*dt/2)*(-G_star);
                k2_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j]+k1_pos[j]*dt/2)/system.objects[1].mass;

                k3_pos[j]=(velocity[i][j]+k2_vel[j]*dt/2);
               // k3_vel[j]=(position[i][j]+k2_pos[j]*dt/2)*(-G_star);
                k3_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j]+k2_pos[j]*dt/2)/system.objects[1].mass;

                k4_pos[j]=(velocity[i][j]+k3_vel[j]*dt);
              //  k4_vel[j]=(position[i][j]+k3_pos[j]*dt)*(-G_star);
                k4_vel[j]=system.Force(system.objects[0], system.objects[1], position[i][j]+k1_pos[j]*dt)/system.objects[1].mass;

                position[i+1][j]=position[i][j] + (k1_pos[j] + k2_pos[j]*2 + k3_pos[j]*2 + k4_pos[j])*dt/6;
                velocity[i+1][j]=velocity[i][j] + (k1_vel[j] + k2_vel[j]*2 + k3_vel[j]*2 + k4_vel[j])*dt/6;

                v[j]=velocity[i+1][j];
                system.objects[1].position[j]=position[i+1][j];
            //    cout<<"new_position="<<system.objects[1].position[j]<<endl;

            }

            kinetik_energy[i+1]=system.E_kinetik (system.objects[1], v.length());
            cout<<"total_velocity="<<v.length()<<endl;
           // cout<<"Kinetic_enerdy="<<kinetik_energy[i+1]<<endl;


            cout<<"total_distance="<<system.distance_between_objects(system.objects[0], system.objects[1])<<endl;
            potential_energy[i+1]=system.E_potential(system.objects[0], system.objects[1]);
            //cout<<"Potential_energy="<<potential_energy[i+1]<<endl;

            momentum[i+1]=system.Momentum(system.objects[1], v.length());
            cout<<"Angular_momentum="<<momentum[i+1]<<endl;
          }
         }

  /* for (int i=0; i<n; i++)
    {for (int j=0; j<2; j++)
        {cout<<position[i][j]<<"  "<<velocity[i][j]<<"  ";}cout<<endl;}*/

    for (int i=0; i<n; i++)
          {
        ofstream myfile; //write results into txt file in order to be able to build a plot
        myfile.open ("en_1000.txt");  //the wave function
        for (int i=0; i<n; i++)
           myfile <<0+i*dt<<" "<<kinetik_energy[i]+potential_energy[i]<<endl;
           myfile.close();}

return;
}
