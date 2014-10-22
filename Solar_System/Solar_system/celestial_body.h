#ifndef CELESTIAL_BODY_H
#define CELESTIAL_BODY_H
#include <vec3.h>

class Celestial_Body
{
public:
   // double* position;           //x, y
    //double* velocity;           //vx, vy
    //double force[2];

    vec3 position;
    vec3 velocity;
    vec3 force;
    double mass;
    char* name;                    //name_of_object

    Celestial_Body (vec3 position, vec3 velocity, double mass_, char* name_of_object);
    Celestial_Body (double x, double y, double z, double vx, double vy, double vz, double mass_, char* name_of_object);

    void print_position();
    void resetForce();
};

#endif // CELESTIAL_BODY_H
