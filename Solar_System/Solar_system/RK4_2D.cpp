#include <RK4_2D.h>
#include <celestial_body.h>
#include <iostream>
#include <math.h>
#include <solar_system.h>
using namespace std;


RK4_2D::RK4_2D()
{
}


void RK4_2D::integrate_RK4_2D(Solar_System &system, double **data, int k, double dt)
{// Do RK4 integration here

    double **k1_vel = new double* [k];
    for (int i=0; i<k; i++)
    {k1_vel[i] = new double [3];}

    double **k2_vel = new double *[k];
    for (int i=0; i<k; i++){ k2_vel[i]=new double [3];}

    double **k3_vel = new double*[k];
    for (int i=0; i<k; i++){ k3_vel[i]=new double [3];}

    double **k4_vel=new double *[k];
    for (int i=0; i<k; i++){ k4_vel[i]=new double [3];}

    double **k1_pos=new double *[k];
    for (int i=0; i<k; i++){ k1_pos[i]=new double [3];}

    double **k2_pos=new double*[k];
    for (int i=0; i<k; i++){ k2_pos[i]=new double [3];}

    double **k3_pos=new double *[k];
    for (int i=0; i<k; i++){ k3_pos[i]=new double [3];}

    double **k4_pos=new double *[k];
    for (int i=0; i<k; i++){ k4_pos[i]=new double [3];}


   // cout<<"distance="<<system.objects[1].position.length()<<endl;
    //cout<<"speed="<<system.objects[1].velocity.length()<<endl;


    //increase it on two planets
    //calculate forces

   // cout<<"mass="<<system.objects[1].mass<<endl;
    //cout<<"mass="<<system.objects[0].mass<<endl;


    for (int a=0; a<k-1; a++)
    {   if (a==0) for (int m=0; m<6; m++){data[0][m]=0;
                                k1_pos[0][m]=0;}

        for (int q=0; q<3; q++){
            k1_pos[a+1][q]=(data[a+1][q+3]);}

       /* k1_pos[a+1][0]=(data[a+1][3]);   //0=x, 1=y, 2=z
        k1_pos[a+1][1]=(data[a+1][4]);      //a+1 - check al planets
        k1_pos[a+1][2]=(data[a+1][5]);
       //cout<<"k1_pos"<<k1_pos[a+1][0]<<"__"<<k1_pos[a+1][1]<<"__"<<k1_pos[a+1][2]<<endl;
       // cout<<"distance="<<system.distance_between_objects (system.objects[0], system.objects[1])<<endl;*/

        k1_vel[a+1][0]=system.Total_Force(system, k, a+1, 0, 0)/system.objects[a+1].mass;
        k1_vel[a+1][1]=system.Total_Force(system, k, a+1, 1, 0)/system.objects[a+1].mass;
        k1_vel[a+1][2]=system.Total_Force(system, k, a+1, 2, 0)/system.objects[a+1].mass;
        //cout<<"k1_vel"<<k1_vel[a+1][0]<<"__"<<k1_vel[a+1][1]<<"__"<<k1_vel[a+1][2]<<endl;
        //cout<<system.Total_Force(system, k, a+1, 0, 0)<<endl;

        k2_pos[a+1][0]=(data[a+1][3]+k1_vel[a+1][0]*dt/2);
        k2_pos[a+1][1]=(data[a+1][4]+k1_vel[a+1][1]*dt/2);
        k2_pos[a+1][2]=(data[a+1][5]+k1_vel[a+1][2]*dt/2);

        k2_vel[a+1][0]=system.Total_Force(system, k, a+1, 0, k1_pos[a+1][0]*dt/2)/system.objects[a+1].mass;
        k2_vel[a+1][1]=system.Total_Force(system, k, a+1, 1, k1_pos[a+1][1]*dt/2)/system.objects[a+1].mass;
        k2_vel[a+1][2]=system.Total_Force(system, k, a+1, 2, k1_pos[a+1][2]*dt/2)/system.objects[a+1].mass;

        k3_pos[a+1][0]=(data[a+1][3]+k2_vel[a+1][0]*dt/2);
        k3_pos[a+1][1]=(data[a+1][4]+k2_vel[a+1][1]*dt/2);
        k3_pos[a+1][2]=(data[a+1][5]+k2_vel[a+1][2]*dt/2);

        k3_vel[a+1][0]=system.Total_Force(system, k, a+1, 0, k2_pos[a+1][0]*dt/2)/system.objects[a+1].mass;
        k3_vel[a+1][1]=system.Total_Force(system, k, a+1, 1, k2_pos[a+1][1]*dt/2)/system.objects[a+1].mass;
        k3_vel[a+1][2]=system.Total_Force(system, k, a+1, 2, k2_pos[a+1][2]*dt/2)/system.objects[a+1].mass;

        k4_pos[a+1][0]=(data[a+1][3]+k3_vel[a+1][0]*dt);
        k4_pos[a+1][1]=(data[a+1][4]+k3_vel[a+1][1]*dt);
        k4_pos[a+1][2]=(data[a+1][5]+k3_vel[a+1][2]*dt);

        k4_vel[a+1][0]=system.Total_Force(system, k, a+1, 0, k3_pos[a+1][0]*dt)/system.objects[a+1].mass;
        k4_vel[a+1][1]=system.Total_Force(system, k, a+1, 1, k3_pos[a+1][1]*dt)/system.objects[a+1].mass;
        k4_vel[a+1][2]=system.Total_Force(system, k, a+1, 2, k3_pos[a+1][2]*dt)/system.objects[a+1].mass;

        data[a+1][0]+= (k1_pos[a+1][0] + k2_pos[a+1][0]*2 + k3_pos[a+1][0]*2 + k4_pos[a+1][0])*dt/6;
        data[a+1][1]+= (k1_pos[a+1][1] + k2_pos[a+1][1]*2 + k3_pos[a+1][1]*2 + k4_pos[a+1][1])*dt/6;
        data[a+1][2]+= (k1_pos[a+1][2] + k2_pos[a+1][2]*2 + k3_pos[a+1][2]*2 + k4_pos[a+1][2])*dt/6;

        data[a+1][3]+= (k1_vel[a+1][0] + k2_vel[a+1][0]*2 + k3_vel[a+1][0]*2 + k4_vel[a+1][0])*dt/6;
        data[a+1][4]+= (k1_vel[a+1][1] + k2_vel[a+1][1]*2 + k3_vel[a+1][1]*2 + k4_vel[a+1][1])*dt/6;
        data[a+1][5]+= (k1_vel[a+1][2] + k2_vel[a+1][2]*2 + k3_vel[a+1][2]*2 + k4_vel[a+1][2])*dt/6;}


return;
}

