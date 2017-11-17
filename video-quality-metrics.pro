QT       += core
QT       -= gui

TARGET    = video-quality-metrics
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

SOURCES += main.cpp

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS_RELEASE += -O3 -mavx -mavx2
