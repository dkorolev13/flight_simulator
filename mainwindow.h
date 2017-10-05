#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <QTime>
#include <QTimer>
#include "qcustomplot.h"
#include "gdal_utils.h"
#include "mygraphicsview.h"
using namespace cv;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_openButton_clicked();
    void on_RUNButton_clicked();
    void on_Calculation();
    void on_Pause_pushButton_clicked();

private:
    Ui::MainWindow *ui;
    Mat img, dst;
    QTimer* timer;

    QPoint point;

    double pixSizeLon, pixSizeLat, topLeftAngleLon, topLeftAngleLat, phix, phiy, Sx, Sy,
           Pos0x, Pos0y, Dx, Dy, W, psi, tetta, gamma, fokus, pixel, delta, Xa, Ya, topleftx_px, toplefty_px;

    int K, k;

    QVector<double> W_X, W_Y;
    QVector<double> S_X, S_Y;
    QVector<double> Kurs_X, Kurs_Y;
    QVector<double> Kren_X, Kren_Y;
    QVector<double> Tangazh_X, Tangazh_Y;

};

#endif // MAINWINDOW_H
