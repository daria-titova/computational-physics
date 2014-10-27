TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    Explicit.cpp \
    Implicit.cpp \
    Closed_form.cpp \
    tridiag.cpp \
    Cra-Nic.cpp

LIBS += -llapack -lblas -larmadillo

HEADERS += \
    tridiag.h \
    Explicit.h \
    Implicit.h \
    Closed_form.h \
    Cra-Nic.h

