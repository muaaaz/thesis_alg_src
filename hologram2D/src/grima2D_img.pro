TEMPLATE = app
CONFIG += console
CONFIG += -std=c++0x
CONFIG -= app_bundl qt gui core

CONFIG(debug, release|debug):DEFINES += _DEBUG

QMAKE_CXXFLAGS += -std=c++0x
# remove possible other optimization flags
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE -= -O3
# add the desired -O3 if not present
QMAKE_CXXFLAGS_RELEASE += -O3

SOURCES += main.cpp \
    grima/greader.cpp \
    grima/gdatabase.cpp \
    grima/ggraph.cpp \
    grima/gsparseset.cpp \
    grima/grima.cpp \
    grima/gextensioncollect.cpp \
    grima/gpattern.cpp \
    grima/gsubgraphiso.cpp \
    grima/gvocab.cpp \
    grima/gcanonicaldfscomputer.cpp

HEADERS += \
    global.hpp \
    grima/gglobal.hpp \
    grima/greader.hpp \
    grima/gdatabase.hpp \
    grima/ggraph.hpp \
    grima/gsparseset.hpp \
    grima/grima.hpp \
    grima/gextensioncollect.hpp \
    grima/gpattern.hpp \
    grima/gsubgraphiso.hpp \
    grima/gvocab.hpp \
    grima/gcanonicaldfscomputer.hpp
