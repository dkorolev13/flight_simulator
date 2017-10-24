#ifndef BINSSENDER_H
#define BINSSENDER_H

#include <QObject>
#include <QtSerialPort/QSerialPort>
#include "structs.h"

class BINSsender : public QObject
{
    Q_OBJECT
public:
    explicit BINSsender(QObject *parent = 0);

private:
    QSerialPort *serial;

private slots:

    void writeData(const QByteArray &data);
    void handleError(QSerialPort::SerialPortError error);
    bool checkCRC(const QByteArray &data);

signals:

public slots:
    void sendNavData(const ISNoutSTRUCT &isnStruct); // принимающий слот
};

#endif // BINSSENDER_H
