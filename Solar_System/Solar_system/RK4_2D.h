#ifndef RK4_2D_H
#define RK4_2D_H

#include <solar_system.h>

class RK4_2D
{
public:

//constructor
RK4_2D();

//member
Solar_System system;

void integrate_RK4_2D(Solar_System &system, double **data, int k, double dt);
};


#endif // RK4_2D_H
