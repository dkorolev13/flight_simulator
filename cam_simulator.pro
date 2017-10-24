#-------------------------------------------------
#
# Project created by QtCreator 2017-04-26T00:25:51
#
#-------------------------------------------------

QT       += core gui serialport


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = cam_simulator
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    qcustomplot.cpp \
    mygraphicsview.cpp \
    binssender.cpp

HEADERS  += mainwindow.h \
    qcustomplot.h \
    mygraphicsview.h \
    structs.h \
    binssender.h

FORMS    += \
    mainwindow.ui

#Подключение заголовочных файлов OpenCV
INCLUDEPATH+=C:/opencv/build/include
INCLUDEPATH+=C:/release-1800-x64-gdal-2-1-3-mapserver-7-0-4/include

LIBS+=C:/opencv/build/x64/vc14/lib/opencv_world320.lib
LIBS+=C:/release-1800-x64-gdal-2-1-3-mapserver-7-0-4/lib/gdal_i.lib
