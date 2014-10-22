#ifndef RK4_H
#define RK4_H
#include <solar_system.h>

class RK4
{
public:

//constructor
RK4();

//member
Solar_System system;

void integrate_RK4(Solar_System &system, double dt, int n);
};


#endif // RK4_H
