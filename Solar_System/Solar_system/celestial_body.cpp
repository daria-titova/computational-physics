#include <celestial_body.h>
#include <iostream>
using namespace std;


Celestial_Body::Celestial_Body(vec3 pos, vec3 vel, double mass_, char* name_of_object) {
position = pos;
velocity = vel;
mass = mass_;
name = name_of_object;
}

Celestial_Body::Celestial_Body(double x, double y, double z, double vx, double vy, double vz, double mass_, char* name_of_object) {
position = vec3(x,y,z);
velocity = vec3(vx,vy,vz);
mass = mass_;
name = name_of_object;
}

void Celestial_Body::resetForce() {
force.setToZero();
}

void Celestial_Body::print_position()
{  cout<<"["<<position[0]<<","<<position[1]<<"]"<<endl;
    return;}


/*Celestial_Body::Celestial_Body( double x, double y, double vx, double vy, double mass_, char *name_of_object) {
    distance_to_the_sun = d;

    position = new double[2];
    position[0] = x;
    position[1] = y;

    velocity = new double[2];
    velocity[0] = vx;
    velocity[1] = vy;

    mass = mass_;

    name = name_of_object;

    resetForce();
}*/


/*void Celestial_Body::resetForce() {
    // force.zeros(); // Armadillo
    force[0] = 0;
    force[1] = 0;
}*/



