#include <solar_system.h>
#include <iostream>
#include <math.h>
#include <vec3.h>
using namespace std;

Solar_System::Solar_System()
{
}

void Solar_System::addCelestialBody(Celestial_Body newObject) {
    objects.push_back(newObject);
}

double Solar_System::distance_between_objects (Celestial_Body object1, Celestial_Body object2)
{ vec3 distance;
    distance=object1.position-object2.position;
    double r=distance.length();
return r;}

double Solar_System::Force (Celestial_Body object1, Celestial_Body object2, double x)
{ double force;
  double G=6.674e-11;
  double AU=1.5e11;
  double year= 31557600;
  double M_sun=2e30;
  double G_star=G/((AU*AU*AU)/(M_sun*year*year));
  double r=Solar_System::distance_between_objects(object1, object2);
  force=x*(-G_star)*(object1.mass*object2.mass)/r*r*r;
  return force;}


double Solar_System::E_kinetik(Celestial_Body object, double v)
{double K;
    K=(object.mass*v*v)/2;
return K;}

double Solar_System::E_potential(Celestial_Body object1, Celestial_Body object2)
{double P;
    double G=6.674e-11;
    double AU=1.5e11;
    double year= 31557600;
    double M_sun=2e30;
    double G_star=G/((AU*AU*AU)/(M_sun*year*year));
    double r=Solar_System::distance_between_objects(object1, object2);
P=-G_star*object1.mass*object2.mass/r;
return P;}

double Solar_System::E(double k, double p)
{double E;
E=k+p;
return E;}

double Solar_System::Momentum(Celestial_Body object, double v)
{
double M;
double r=Solar_System::distance_between_objects(objects[0], object);
M=object.mass*r*v;
return M;}
