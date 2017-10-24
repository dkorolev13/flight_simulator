#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/videoio.hpp>
#include <QTime>
#include <QTimer>
#include <QFile>

#include "qcustomplot.h"
#include "gdal_utils.h"
#include "mygraphicsview.h"
#include "binssender.h"

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
    void on_RUNButton_clicked(bool checked);
    void on_Calculation();
    void on_Pause_pushButton_clicked();
    void on_X_spinBox_valueChanged(int arg1);
    void on_Y_spinBox_valueChanged(int arg1);

    void route1_W();
    void route2_W();

    void on_Scale_spinBox_valueChanged(double arg1);

    void on_openButton_2_clicked();

    void on_adjustment_checkbox_clicked(bool checked);


private:
    Ui::MainWindow *ui;
    Mat imgFull, img, dst;
    QTimer* timer;

    QPoint point;

    VideoWriter videoWriter;

    QLabel* qdst_label;
    QImage adjustmentImage;

    BINSsender *binsSender;

    QFile file_H;

    double pixSizeLon, pixSizeLat, topLeftAngleLon, topLeftAngleLat, phix, phiy, Sx, Sy,
           Pos0x, Pos0y, Dx, Dy, W, psi, gamma, tetta, fokus, pixel, Xa, Ya, topleftx_px, toplefty_px,
           S, sigma, Tc, H, H_shum, psi_shum, gamma_shum, Y_H, Y_psi, Y_gamma, alpha, X_H, X_psi, X_gamma, mW, aW, as, WL, X_tetta, Y_tetta, tetta_shum;

    double s0, h0, s1, H_max, s2, s3;

    int K, k;

    QVector<double> W_X, W_Y;
    QVector<double> S_X, S_Y;
    QVector<double> Kurs_X, Kurs_Y;
    QVector<double> Kren_X, Kren_Y;
    QVector<double> Tangazh_X, Tangazh_Y;

};

#endif // MAINWINDOW_H
