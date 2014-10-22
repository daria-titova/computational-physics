#ifndef VERLET_H
#define VERLET_H

#include <solar_system.h>

class Verlet
{public:

    //member
    Solar_System system;

    //constructor
    Verlet();

    //functions
    void verlet_integration (Solar_System &system, double dt, int n);

};

#endif // VERLET_H

