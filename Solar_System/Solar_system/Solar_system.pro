TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    celestial_body.cpp \
    solar_system.cpp \
    Verlet.cpp \
    vec3.cpp \
    Rk4.cpp \
    RK4_2D.cpp

HEADERS += \
    celestial_body.h \
    solar_system.h \
    Verlet.h \
    vec3.h \
    Rk4.h \
    RK4_2D.h

