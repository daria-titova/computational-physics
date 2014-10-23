#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H

#include <celestial_body.h>

#include <vector>
using std::vector;

class Solar_System
{  public:

    //members of theh class
    vector<Celestial_Body> objects;

    //constructor
    Solar_System();

    //functions
    void addCelestialBody(Celestial_Body newObject);
    double distance_between_objects (Celestial_Body object1, Celestial_Body object2);
    double Force (Celestial_Body body1, Celestial_Body body2, double x);
    double Total_Force(Solar_System &system, int k, int p, int j, double x);
    double E_kinetik(Celestial_Body object, double v);
    double E_potential(Celestial_Body object1, Celestial_Body object2);
    double E(double k, double p);
    double Momentum(Celestial_Body object, double v);
};

#endif // SOLAR_SYSTEM_H
