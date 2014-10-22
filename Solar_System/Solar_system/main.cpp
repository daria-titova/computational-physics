#include <iostream>
#include <math.h>
#include <solar_system.h>
#include <Verlet.h>
#include <Rk4.h>
using namespace std;

int main () {

    int n;
    double a, b;
    cout<<"Insert number of steps n"<<endl;
    cin>>n;
    cout<<"Insert the beginning of the interval a"<<endl;
    cin>>a;
    cout<<"Insert the end of the interval b"<<endl;
    cin>>b;

    double h;
    h=(b-a)/n;

    cout<<"n="<<n<<endl<<"a="<<a<<endl<<"b="<<b<<endl<<"h="<<h<<endl;

    //declare the solar system
    Solar_System my_Solar_System;

    // set up object for solar system
    Celestial_Body sun(vec3(0, 0, 0), vec3(0, 0, 0), 1.0, "Sun");
    Celestial_Body earth(vec3(1.0, 0, 0), vec3(0, 2*M_PI, 0), 3e-6, "Earth");
    /*Celestial_Body jupiter(5.20, 5.20, 0, 0, 0, 0.95e-3, "Jupiter");
    Celestial_Body mars(1.52, 1.52, 0, 0, 0, 3.30-7);
    Celestial_Body venus(0.72, 0.72, 0, 0, 0, 2.45e-6);
    Celestial_Body saturn(9.54, 9.54, 0, 0, 0, 2.75e-4);
    Celestial_Body mercury(0.39, 0.39, 0, 0, 0, 1.2e-7);
    Celestial_Body uranus(19.19, 19.19, 0, 0, 0, 4.4e-5);
    Celestial_Body neptun(30.06, 30.06, 0, 0, 0, 0.515e-4);
    Celestial_Body pluto(39.53, 39.53, 0, 0, 0, 0.655e-8);*/

  //  earth.print_position();

    //add objects (planets/celestial bodies) into my solar system
    my_Solar_System.addCelestialBody(sun);
    my_Solar_System.addCelestialBody(earth);

    for(int i = 0; i<my_Solar_System.objects.size(); i++) {
        Celestial_Body this_Body = my_Solar_System.objects[i];
        //cout << "The name of this object is " << this_Body.name<< " and position is " << this_Body.position[0] << " and mass is " << this_Body.mass << endl;
    }

    my_Solar_System.distance_between_objects(sun, earth);
    cout<<"The distance ="<<my_Solar_System.distance_between_objects(sun, earth)<<endl;

   // cout<<"test distance="<<my_Solar_System.distance_between_objects(sun, earth)<<endl;
   // cout<<"test force="<<my_Solar_System.force_between_objects(sun, earth)<<endl;

   // Verlet method;
    //method.verlet_integration(my_Solar_System, h, n);
    //cout << "Position: " << my_Solar_System.objects[0].position[0] << endl;

    RK4 method;
    method.integrate_RK4(my_Solar_System, h, n);
    //cout << "Position: " << my_Solar_System.objects[0].position[0] << endl;

return 0;}
