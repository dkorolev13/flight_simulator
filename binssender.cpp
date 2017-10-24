#include "binssender.h"
#include <QDebug>
#include <QDataStream>

BINSsender::BINSsender(QObject *parent) : QObject(parent)
{
    qDebug() << __FUNCTION__;

    serial = new QSerialPort(this);

    connect(serial, static_cast<void (QSerialPort::*)(QSerialPort::SerialPortError)>(&QSerialPort::error), this, &BINSsender::handleError);

    serial->setPortName("COM1");
    serial->setBaudRate(115200);
    serial->setDataBits(QSerialPort::Data8);
    serial->setParity(QSerialPort::NoParity);
    serial->setStopBits(QSerialPort::OneStop);
    serial->setFlowControl(QSerialPort::NoFlowControl);
    if(serial->isOpen())
        serial->close();

    qDebug() << "serial port ISN -> open" << serial->open(QIODevice::ReadWrite);
}

void BINSsender::sendNavData(const ISNoutSTRUCT &isnStruct)
{
    //из isnStruct сделаем QByteArray
    //
    QByteArray byteArray;
    QDataStream ds(&byteArray, QIODevice::WriteOnly);

    ds.setVersion(QDataStream::Qt_4_4); //Обязательно. Иначе float некорректно разбирается.
    ds.setByteOrder(QDataStream::LittleEndian); //Обязательно. По умолчанию - BigEndian.

    int B_integer = isnStruct.B_wgs84 * 100;
    int L_integer = isnStruct.L_wgs84 * 100;
    unsigned short crc = 0;

    //ds << (unsigned int)0xAAAA3770;
    ds << (unsigned int)0x7037AAAA;
    ds << (unsigned int&)isnStruct.SS;
    //ds << (unsigned short&)0;
    ds << (float &)isnStruct.Ax;
    ds << (float &)isnStruct.Ay;
    ds << (float &)isnStruct.Az;
    ds << (float &)isnStruct.Wx;
    ds << (float &)isnStruct.Wy;
    ds << (float &)isnStruct.Wz;
    ds << (float &)isnStruct.gamma;
    ds << (float &)isnStruct.psi;
    ds << (float &)isnStruct.tetta;
    ds << (int &)B_integer;
    ds << (int &)L_integer;
    ds << (float &)isnStruct.H_wgs84;
    ds << (unsigned short&)crc;

    writeData(byteArray);
}

void BINSsender::writeData(const QByteArray &data)
{
    serial->write(data);
}

bool BINSsender::checkCRC(const QByteArray &data)
{
//    unsigned short crc = 0;
//    for(int i=0; i<data.size(); ++i)
//        crc = ccitt_crc16_table[(crc >> 8 ^ data.at(i)) & 0xffU] ^ (crc << 8);

//    qDebug() << "crc=" << crc;
    return true;
}

void BINSsender::handleError(QSerialPort::SerialPortError error)
{
    if(error == QSerialPort::ResourceError)
    {
        qDebug() << serial->errorString();
        serial->close();
    }
}
