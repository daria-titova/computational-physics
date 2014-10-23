#include <iostream>
#include <math.h>
#include <solar_system.h>
#include <Verlet.h>
#include <Rk4.h>
#include <RK4_2D.h>
#include <fstream>
using namespace std;

int main () {

    int n;
    double a=0.0, b;
    cout<<"Insert number of steps n"<<endl;
    cin>>n;
    cout<<"Insert the end of the interval b"<<endl;
    cin>>b;

    double h;
    h=(b-a)/n;

    cout<<"n="<<n<<endl<<"a="<<a<<endl<<"b="<<b<<endl<<"h="<<h<<endl;

    //declare the solar system
    Solar_System my_Solar_System;

    // set up object for solar system

   // Celestial_Body sun(vec3(2.0, 0, 0), vec3(0, -0.009, 0), 1.0, "Sun");
   // Celestial_Body earth(vec3(3.0, 0, 0), vec3(0, 2*M_PI, 0), 3e-6, "Earth");
   // Celestial_Body jupiter(vec3(7.20, 0, 0), vec3(0, 2.75, 0), 0.95e-3, "Jupiter");

    Celestial_Body sun(vec3(0, 0, 0), vec3(0, 0, 0), 1.0, "Sun");
    Celestial_Body earth(vec3(1.0, 0, 0), vec3(0, 2*M_PI, 0), 3e-6, "Earth");
    Celestial_Body jupiter(vec3(5.20, 0, 0), vec3(0, 2.75, 0), 0.95, "Jupiter");
    Celestial_Body mars(vec3(1.52, 0, 0), vec3(0, 5.07, 0), 3.30e-7, "mars");
    Celestial_Body venus(vec3(0.72, 0, 0), vec3(0, 7.38, 0), 2.45e-6, "venus");
    Celestial_Body saturn(vec3(9.54, 0, 0), vec3(0, 2.038, 0), 2.75e-4, "saturn");
    Celestial_Body mercury(vec3(0.39, 0, 0), vec3(0, 10.07, 0), 1.2e-7, "mercury");
    Celestial_Body uranus(vec3(19.19, 0, 0), vec3(0, 1.43, 0), 4.4e-5, "uranus");
    Celestial_Body neptun(vec3(30.06, 0, 0), vec3(0, 1.14, 0), 0.515e-4, "neptun");
    Celestial_Body pluto(vec3(39.53, 0, 0), vec3(0, 0.98, 0), 0.655e-8, "pluto");

    //add objects (planets/celestial bodies) into my solar system
    my_Solar_System.addCelestialBody(sun);
    my_Solar_System.addCelestialBody(earth);
    my_Solar_System.addCelestialBody(jupiter);
   // my_Solar_System.addCelestialBody(mars);
   // my_Solar_System.addCelestialBody(venus);
   // my_Solar_System.addCelestialBody(saturn);
   // my_Solar_System.addCelestialBody(mercury);
   // my_Solar_System.addCelestialBody(uranus);
   // my_Solar_System.addCelestialBody(neptun);
   // my_Solar_System.addCelestialBody(pluto);

    int k=my_Solar_System.objects.size();
    double **data=new double * [k];
    for (int p=0; p<k; p++){
        data[p]= new double [6];}

    for (int i=0; i<k; i++){
            data[i][0]=my_Solar_System.objects[i].position[0];
            data[i][1]=my_Solar_System.objects[i].position[1];
            data[i][2]=my_Solar_System.objects[i].position[2];
            data[i][3]=my_Solar_System.objects[i].velocity[0];
            data[i][4]=my_Solar_System.objects[i].velocity[1];
            data[i][5]=my_Solar_System.objects[i].velocity[2];
            for (int l=0; l<6; l++){
                cout<<"data="<<data[i][l]<<"___";}
            cout<<endl;
        } // assaign matrix of data to the initial values

    cout<<"k="<<k<<endl;

    ofstream myfile; //write results into txt file in order to be able to build a plot
    myfile.open ("x_1000.txt");  //the wave function


    for (int j=0; j<n; j++){  // j - the number(order) of the RK's iteration
      // cout<<"starting rk4"<<endl;
        RK4_2D method;
        method.integrate_RK4_2D(my_Solar_System, data, k, h);
        //we should reassign positions and velocities
        for (int i=0; i<k-1; i++){
                my_Solar_System.objects[i+1].position[0]=data[i+1][0];
                my_Solar_System.objects[i+1].position[1]=data[i+1][1];
                my_Solar_System.objects[i+1].position[2]=data[i+1][2];
                my_Solar_System.objects[i+1].velocity[0]=data[i+1][3];
                my_Solar_System.objects[i+1].velocity[1]=data[i+1][4];
                my_Solar_System.objects[i+1].velocity[2]=data[i+1][5];
                   }



        for (int l=0; l<6; l++){
        // cout<<data[1][l]<<"__";
        }

       // cout<<endl;

           myfile <<0+j*h<<" "<<data[1][0]<<endl;

    }
    cout<<"finish rk4"<<endl;
    myfile.close();


   /* Verlet method;
    method.verlet_integration(my_Solar_System, h, n);
    cout << "Position: " << my_Solar_System.objects[0].position[0] << endl;*/

  /*  RK4 method;
    method.integrate_RK4(my_Solar_System, h, n);
    cout << "Position: " << my_Solar_System.objects[0].position[0] << endl;*/

return 0;}
