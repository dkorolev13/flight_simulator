#ifndef STRUCTS_H
#define STRUCTS_H

#include <QMetaType>

struct ISNoutSTRUCT //Пакет данных навигации и ориентации 70h
{
    unsigned int SS;

    float Ax;
    float Ay;
    float Az;

    float Wx;
    float Wy;
    float Wz;

    float gamma;
    float psi;
    float tetta;

    double B_wgs84; // широта
    double L_wgs84; // долгота

    float H_wgs84;
};
Q_DECLARE_METATYPE(ISNoutSTRUCT)

struct ISNinSTRUCT //Ввод данных для управления режимами работы 45h
{
    int B_wgs84;
    int L_wgs84;
    int H_wgs84;
    int cmd;
};
Q_DECLARE_METATYPE(ISNinSTRUCT)

#endif // STRUCTS_H
