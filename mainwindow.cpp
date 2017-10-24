#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QtCore/qmath.h>
#include <QMessageBox>
#include <opencv\cv.h>
#include <QDebug>
#include <QFileInfo>
#include <QTextStream>
#include <QTime>
#include "gdal_priv.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Создаем таймер
    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(on_Calculation()));

    binsSender = new BINSsender(this);

    // Обнуляем CV Mat'ы картинок
    img = imgFull = dst = NULL;
    qdst_label = new QLabel(this);
    qdst_label->setWindowFlags(Qt::Dialog);
    qdst_label->show();

    //ui->X_spinBox->setMaximum(QApplication::desktop()->size().width());
    //ui->Y_spinBox->setMaximum(QApplication::desktop()->size().height());

    ui->customplot_W->addGraph();
    ui->customplot_W->yAxis->setLabel("W, м/с");
    ui->customplot_W->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui->customplot_W->xAxis->setRange(0, 3000);
    ui->customplot_W->yAxis->setRange(0, 100);
    ui->customplot_W->setInteraction(QCP::iRangeZoom,true);
    ui->customplot_W->setInteraction(QCP::iRangeDrag, true);

    ui->customplot_S_H->addGraph();
    ui->customplot_S_H->setInteraction(QCP::iRangeZoom,true);
    ui->customplot_S_H->setInteraction(QCP::iRangeDrag, true);
    ui->customplot_S_H->xAxis->setLabel("S, м");
    ui->customplot_S_H->yAxis->setLabel("H, м");
    ui->customplot_S_H->graph(0)->setLineStyle(QCPGraph::lsLine);

    ui->customplot_Kurs->addGraph();
    ui->customplot_Kurs->yAxis->setLabel("Курс, °");
    ui->customplot_Kurs->yAxis->setRange(0, 370);

    ui->customplot_Kren->addGraph()->setName("Kren");
    ui->customplot_Kren->addGraph()->setName("Tangazh");
    ui->customplot_Kren->graph(1)->setPen(QPen(Qt::red));
    ui->customplot_Kren->legend->setVisible(true); //Включаем Легенду графика
    ui->customplot_Kren->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop); // Устанавливаем Легенду в левый верхний угол графика
    ui->customplot_Kren->yAxis->setLabel("Крен, Тангаж, ");
    ui->customplot_Kren->yAxis->setRange(-180, 180);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Кнопка выбора эталонного изображения
void MainWindow::on_openButton_clicked()
{
    QString tifFileName = QFileDialog::getOpenFileName(this, "Выбор файла", "D:/Projects/Strelnya", "Изображения (*.tif);;Изображения (*.tiff)");
    if(tifFileName.isEmpty())
        return;
    ui->leFileName->setText(tifFileName);

    // Открываем geotif
    GDALDataset  *poDataset;
    GDALAllRegister();
    poDataset = (GDALDataset *) GDALOpen(tifFileName.toStdString().data(), GA_ReadOnly);
    if( poDataset == NULL )
        {
            QMessageBox::information(this, "Информация", "Выберите снимок!");
            return;
        }

    double adfGeoTransform[6];
    if( poDataset->GetGeoTransform(adfGeoTransform) == CE_None)
    {
        pixSizeLon = adfGeoTransform[1]; // размер пикселя по долготе в метрах (по горизонтали)
        pixSizeLat = adfGeoTransform[5]; // размер пикселя по широте в метрах (по вертикали)

        // Координаты левого верхнего угла эталона по долготе и широте в метрах
        topLeftAngleLon = adfGeoTransform[0];
        topLeftAngleLat = adfGeoTransform[3];
    }

    GDALClose((GDALDatasetH) poDataset);

    //Функция перевода geotif в формат OpenCV
    imgFull = imread(tifFileName.toLocal8Bit().data(), cv::IMREAD_LOAD_GDAL | cv::IMREAD_GRAYSCALE );

//    //Координата правого верхнего угла эталдона по долготе в метрах
//    double topRightAngleLon = topLeftAngleLon + (img.cols) * pixSizeLon;

//    //Координата левого нижнего угла эталона по широте в метрах
//    double bottomLeftAngleLat = topLeftAngleLat + (img.rows) * pixSizeLat;

//    //Координаты центра эталона в метрах
//    ui->Lon_doubleSpinBox->setRange(topLeftAngleLon, topRightAngleLon);
//    ui->Lat_doubleSpinBox->setRange(bottomLeftAngleLat, topLeftAngleLat);

//    ui->Lon_doubleSpinBox->setValue(topLeftAngleLon + ((topRightAngleLon - topLeftAngleLon) / 2)); // центр снимка по долготе
//    ui->Lat_doubleSpinBox->setValue(topLeftAngleLat - ((topLeftAngleLat - bottomLeftAngleLat)/ 2)); // центр снимка по широте

//    // Позиция центра площади обхода
//    Pos0x = ui->Lon_doubleSpinBox->value();
//    Pos0y = ui->Lat_doubleSpinBox->value();
}

// Кнопка START
void MainWindow::on_RUNButton_clicked(bool checked)
{
    if(checked)
    {
//        qDebug() << videoWriter.open(
//                        "D:/Projects/flight.avi",
//                        //VideoWriter::fourcc('D','I','V','X'),
//                        //VideoWriter::fourcc('X','V','I','D'),
//                        VideoWriter::fourcc('P','I','M','1'),
//                        10.0,
//                        cvSize(ui->N2_doubleSpinBox->value(), ui->N1_doubleSpinBox->value()));

        if(imgFull.empty())
        {
            ui->RUNButton->setChecked(false);
            QMessageBox::information(this, "Информация", "Выберите картинку!");
            return;
        }

        //qdst_label->resize(ui->N2_doubleSpinBox->value(), ui->N1_doubleSpinBox->value());

        ui->X_spinBox->setValue((QApplication::desktop()->size().width() - qdst_label->size().width())/2);
        ui->Y_spinBox->setValue((QApplication::desktop()->size().height() - qdst_label->size().height())/2);

        file_H.setFileName("D:/Projects/H.txt");
        file_H.open(QFile::WriteOnly);

        // Запуск qsrand
        qsrand(QTime(0,0,0).secsTo(QTime::currentTime()));

        k = 0.0;

        // Считываем данные из интерфейса
        // Положение ЛА
        gamma = qDegreesToRadians(ui->Kren_doubleSpinBox->value());
        psi = qDegreesToRadians(ui->Kurs_doubleSpinBox->value());
        tetta = qDegreesToRadians(ui->Tangazh_doubleSpinBox->value());

        // Позиция центра площади обхода
        Pos0x = ui->Lon_doubleSpinBox->value();
        Pos0y = ui->Lat_doubleSpinBox->value();

        // Зависимость H от S (убрать лишние переменные)
        //H = ui->H_max_doubleSpinBox->value();

        s1 = ui->s1_doubleSpinBox->value();
        s2 = ui->s2_doubleSpinBox->value();
        s3 = ui->s3_doubleSpinBox->value();
        H_max = ui->H_max_doubleSpinBox->value();

        //Для зашумления
        sigma = ui->sigma_doubleSpinBox->value();
        Tc = ui->Tc_SpinBox->value();

        // Параметры ОСН
        fokus = ui->Fokus_doubleSpinBox->value();
        pixel = ui->Pixel_doubleSpinBox->value();

        // Кол-во итераций
        K = ui->K_spinBox->value();

        //Размеры площади обхода, м
        Dx = ui->Dx_doubleSpinBox->value();
        Dy = ui->Dy_doubleSpinBox->value();

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Начальная фаза
        phix = ui->fx_doubleSpinBox->value() * 2 * M_PI * ui->k_spinBox->value() / ui->K_spinBox->value();
        phiy = (ui->fy_doubleSpinBox->value() * 2 * M_PI * ui->k_spinBox->value() / ui->K_spinBox->value()) +  qDegreesToRadians(double(ui->F_spinBox->value()));

        phix = phix - qFloor(phix / (2 * M_PI)) * 2 * M_PI;
        phiy = phiy - qFloor(phiy / (2 * M_PI)) * 2 * M_PI;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        // Путевая скорость
        W = ui->w_doubleSpinBox->value();

        // Стартовая позиция
        Sx = Dx * qSin(phix) + Pos0x;
        Sy = Dy * qSin(phiy) + Pos0y;

        // Поправка на ветер
        mW = ui->mW_doubleSpinBox->value();
        aW = ui->aW_doubleSpinBox->value();

        // Пройденный путь
        S = 0.0;

        Xa = (K * phix / (2 * M_PI)) + K / 4;
        Ya = (K * phiy / (2 * M_PI)) + K / 4;

        // Координаты левого верхнего угла в пикселях с учетом рамки 1000 px
        topleftx_px = ((Pos0x - Dx) - topLeftAngleLon) / pixSizeLon - 1000;
        toplefty_px = ((Pos0y + Dy) - topLeftAngleLat) / pixSizeLat - 1000;

        double width_px = 2*Dx / pixSizeLon + 2000;
        double height_px = -(2*Dy / pixSizeLat) + 2000;

        // Вывод площади обхода
        Rect roi(topleftx_px, toplefty_px, width_px, height_px);
        img = imgFull(roi);
        QImage qimg(img.data, img.cols, img.rows, static_cast<int>(img.step), QImage::Format_Grayscale8);

        ui->mygraphicsview->scene->clear();
        ui->mygraphicsview->scene->setSceneRect(qimg.rect());
        ui->mygraphicsview->scene->addPixmap(QPixmap::fromImage(qimg));
        ui->mygraphicsview->polygon = ui->mygraphicsview->scene->addPolygon(QPolygon(), QPen(Qt::blue,3));

        // Создаем картинку результата
        dst = Mat(ui->N1_doubleSpinBox->value(), ui->N2_doubleSpinBox->value(), CV_8UC1);

        // Обнуляем вектора
        W_X.clear();
        W_Y.clear();

        S_X.clear();
        S_Y.clear();

        Kurs_X.clear();
        Kurs_Y.clear();

        Kren_X.clear();
        Kren_Y.clear();

        Tangazh_X.clear();
        Tangazh_Y.clear();


        // Создаем оси графика S_H
        ui->customplot_S_H->xAxis->setRange(0, s3 + 25);
        ui->customplot_S_H->yAxis->setRange(0, H_max + 25);

        // Создаем полотно и оси графика W
        // ПОМЕНЯТЬ НА ПЕРЕМЕННУЮ

       // ui->customplot_W->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::blue, 2));

        // Создаем полотно и оси графика psi (курса)
        ui->customplot_Kurs->xAxis->setRange(0, K);

        ui->customplot_Kurs->graph(0)->setLineStyle(QCPGraph::lsLine);
        //ui->customplot_Kurs->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::blue, 2));

        // Создаем полотно и оси графиков tetta(крена) и gamma(тангажа)
        ui->customplot_Kren->xAxis->setRange(0, K);
        ui->customplot_Kren->yAxis->setRange(-5, 5);
        ui->customplot_Kren->graph(0)->setLineStyle(QCPGraph::lsLine);
       // ui->customplot_Kren->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::blue, 2));
         //->setName("0");
        ui->customplot_Kren->graph(1)->setLineStyle(QCPGraph::lsLine);
        ui->customplot_Kren->graph(1)->setPen(QPen(Qt::red));
      //  ui->customplot_Kren->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::red, Qt::red, 2));

        // Запускаем таймер
        timer->start(ui->deltaT_doubleSpinBox->value() * 1000);
        ui->RUNButton->setText("STOP");
    }
    else
    {
        timer->stop();
        qdst_label->hide();
        qdst_label->clear();
        ui->RUNButton->setText("START");
    }
}

// Кнопка PAUSE
void MainWindow::on_Pause_pushButton_clicked()
{
    if(ui->Pause_pushButton->isChecked())
        timer->stop();
    else
        timer->start(ui->deltaT_doubleSpinBox->value() * 1000);
}

// Слот вычислений
void MainWindow::on_Calculation()
{
    // Зашумление
    alpha = qExp(-ui->deltaT_doubleSpinBox->value() / ui->Tc_SpinBox->value());

    double rnd_H = qrand() * 1.0 / RAND_MAX;
    double rnd_gamma = qrand() * 1.0 / RAND_MAX;
    double rnd_psi = qrand() * 1.0 / RAND_MAX;
    double rnd_tetta = qrand() * 1.0 / RAND_MAX;

    X_H = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_H)) * qCos(2 * M_PI * rnd_H);
    X_gamma = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_gamma)) * qCos(2 * M_PI * rnd_gamma);
    X_psi = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_psi)) * qCos(2 * M_PI * rnd_psi);
    X_tetta = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_tetta)) * qCos(2 * M_PI * rnd_tetta);

    if(ui->Route_1_radioButton->isChecked())
        route1_W();
    else
        route2_W();

    ISNoutSTRUCT struct_out;
    //struct_out.SS = 0x0000802C;
    struct_out.SS = 0x2C800000;
    struct_out.Ax = 0;
    struct_out.Ay = 0;
    struct_out.Az = 0;
    struct_out.Wx = 0;
    struct_out.Wy = 0;
    struct_out.Wz = 0;
//    struct_out.gamma = gamma_shum;
//    struct_out.psi = psi_shum;
//    struct_out.tetta = tetta_shum;
    struct_out.gamma = qRadiansToDegrees(gamma);
    struct_out.psi = qRadiansToDegrees(psi);
    struct_out.tetta = qRadiansToDegrees(tetta);
    struct_out.B_wgs84 = Sy;
    struct_out.L_wgs84 = Sx;
//    struct_out.H_wgs84 = H_shum;
    struct_out.H_wgs84 = H;

    binsSender->sendNavData(struct_out);

    QTime curTime(0, 0);
    curTime = curTime.addMSecs(k*0.2*1000);
    ui->time_label->setText(curTime.toString("mm:ss:zzz"));

    // Вывод значений на графики
    W_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    W_Y.push_back(W);

    S_X.push_back(S);
    S_Y.push_back(H);

    Kurs_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    Kurs_Y.push_back(psi * 180/M_PI);

    Kren_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    Kren_Y.push_back(gamma * 180/M_PI);

    Tangazh_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    Tangazh_Y.push_back(tetta * 180/M_PI);

    ui->customplot_W->graph(0)->setData(W_X, W_Y);
    ui->customplot_S_H->graph(0)->setData(S_X, S_Y);
    ui->customplot_Kurs->graph(0)->setData(Kurs_X, Kurs_Y);
    ui->customplot_Kren->graph(0)->setData(Kren_X, Kren_Y);
    ui->customplot_Kren->graph(1)->setData(Tangazh_X, Tangazh_Y);

    ui->customplot_W->replot();
    ui->customplot_Kurs->replot();
    ui->customplot_Kren->replot();
    ui->customplot_S_H->replot();

    double gamma_Rx = -gamma;
    double tetta_Rz = -tetta;
    double psi_Ry = psi - M_PI_2;

    //Инициируем массивы
    double Rx[3][3], Rz[3][3], Rxz[3][3], Ry[3][3], R[3][3], r[3][3];

    Rx[0][0] = 1; Rx[0][1] = 0;              Rx[0][2] = 0;
    Rx[1][0] = 0; Rx[1][1] = qCos(gamma_Rx); Rx[1][2] = -qSin(gamma_Rx);
    Rx[2][0] = 0; Rx[2][1] = qSin(gamma_Rx); Rx[2][2] = qCos(gamma_Rx);

    Rz[0][0] = qCos(tetta_Rz); Rz[0][1] = -qSin(tetta_Rz);  Rz[0][2] = 0;
    Rz[1][0] = qSin(tetta_Rz); Rz[1][1] = qCos(tetta_Rz);   Rz[1][2] = 0;
    Rz[2][0] = 0;              Rz[2][1] = 0;                Rz[2][2] = 1;

    Ry[0][0] = qCos(psi_Ry);  Ry[0][1] = 0; Ry[0][2] = qSin(psi_Ry);
    Ry[1][0] = 0;             Ry[1][1] = 1; Ry[1][2] = 0;
    Ry[2][0] = -qSin(psi_Ry); Ry[2][1] = 0; Ry[2][2] = qCos(psi_Ry);

    Rxz[0][0] = 0; Rxz[0][1] = 0; Rxz[0][2] = 0;
    Rxz[1][0] = 0; Rxz[1][1] = 0; Rxz[1][2] = 0;
    Rxz[2][0] = 0; Rxz[2][1] = 0; Rxz[2][2] = 0;

    R[0][0] = 0; R[0][1] = 0; R[0][2] = 0;
    R[1][0] = 0; R[1][1] = 0; R[1][2] = 0;
    R[2][0] = 0; R[2][1] = 0; R[2][2] = 0;

    r[0][0] = 1; r[0][1] = 0; r[0][2] = 0;
    r[1][0] = 0; r[1][1] = 0; r[1][2] = -1;
    r[2][0] = 0; r[2][1] = 1; r[2][2] = 0;

    Mat M_Rx = Mat(3, 3, CV_64FC1, Rx);
    Mat M_Rz = Mat(3, 3, CV_64FC1, Rz);
    Mat M_Rxz = Mat(3, 3, CV_64FC1, Rxz);
    Mat M_Ry = Mat(3, 3, CV_64FC1, Ry);
    Mat M_R = Mat(3, 3, CV_64FC1, R);
    Mat M_r = Mat(3, 3, CV_64FC1, r);

    M_Rxz = M_Rx * M_Rz;
    M_R = M_Rxz * M_Ry;

    Mat vx = M_R.col(0);
    Mat vy = M_R.col(1);
    Mat vz = M_R.col(2);

    double N1 = ui->N1_doubleSpinBox->value();
    double N2 = ui->N2_doubleSpinBox->value();

    // Координаты вершин углов проективной плоскости
    double p0[3] = { 0.0,  -fokus/pixel,  0.0};
    double p1[3] = {-N1/2, -fokus/pixel, -N2/2};
    double p2[3] = {-N1/2, -fokus/pixel,  N2/2};
    double p3[3] = { N1/2, -fokus/pixel,  N2/2};
    double p4[3] = { N1/2, -fokus/pixel, -N2/2};

    Mat V_p0 = Mat(1, 3, CV_64FC1, p0);
    Mat V_p1 = Mat(1, 3, CV_64FC1, p1);
    Mat V_p2 = Mat(1, 3, CV_64FC1, p2);
    Mat V_p3 = Mat(1, 3, CV_64FC1, p3);
    Mat V_p4 = Mat(1, 3, CV_64FC1, p4);

    Mat vx0 = vx.t() * V_p0.t();
    Mat vx1 = vx.t() * V_p1.t();
    Mat vx2 = vx.t() * V_p2.t();
    Mat vx3 = vx.t() * V_p3.t();
    Mat vx4 = vx.t() * V_p4.t();

    Mat vy0 = vy.t() * V_p0.t();
    Mat vy1 = vy.t() * V_p1.t();
    Mat vy2 = vy.t() * V_p2.t();
    Mat vy3 = vy.t() * V_p3.t();
    Mat vy4 = vy.t() * V_p4.t();

    Mat vz0 = vz.t() * V_p0.t();
    Mat vz1 = vz.t() * V_p1.t();
    Mat vz2 = vz.t() * V_p2.t();
    Mat vz3 = vz.t() * V_p3.t();
    Mat vz4 = vz.t() * V_p4.t();

    // Координаты центральной проекций в плане углов снимка
    double P0[3] = {-fokus/pixel * vx0.at<double>(0, 0) / vy0.at<double>(0, 0), -fokus/pixel, -fokus/pixel * vz0.at<double>(0, 0) / vy0.at<double>(0, 0)};
    double P1[3] = {-fokus/pixel * vx1.at<double>(0, 0) / vy1.at<double>(0, 0), -fokus/pixel, -fokus/pixel * vz1.at<double>(0, 0) / vy1.at<double>(0, 0)};
    double P2[3] = {-fokus/pixel * vx2.at<double>(0, 0) / vy2.at<double>(0, 0), -fokus/pixel, -fokus/pixel * vz2.at<double>(0, 0) / vy2.at<double>(0, 0)};
    double P3[3] = {-fokus/pixel * vx3.at<double>(0, 0) / vy3.at<double>(0, 0), -fokus/pixel, -fokus/pixel * vz3.at<double>(0, 0) / vy3.at<double>(0, 0)};
    double P4[3] = {-fokus/pixel * vx4.at<double>(0, 0) / vy4.at<double>(0, 0), -fokus/pixel, -fokus/pixel * vz4.at<double>(0, 0) / vy4.at<double>(0, 0)};

    Mat V_P0_tmp = Mat(1, 3, CV_64FC1, P0);
    Mat V_P1_tmp = Mat(1, 3, CV_64FC1, P1);
    Mat V_P2_tmp = Mat(1, 3, CV_64FC1, P2);
    Mat V_P3_tmp = Mat(1, 3, CV_64FC1, P3);
    Mat V_P4_tmp = Mat(1, 3, CV_64FC1, P4);

    Mat V_P0 = Mat(1, 3, CV_64FC1);
    Mat V_P1 = Mat(1, 3, CV_64FC1);
    Mat V_P2 = Mat(1, 3, CV_64FC1);
    Mat V_P3 = Mat(1, 3, CV_64FC1);
    Mat V_P4 = Mat(1, 3, CV_64FC1);

    // Переход от осей XZ (маткада) к осям XY (экрана)
    V_P0 = M_r * V_P0_tmp.t();
    V_P1 = M_r * V_P1_tmp.t();
    V_P2 = M_r * V_P2_tmp.t();
    V_P3 = M_r * V_P3_tmp.t();
    V_P4 = M_r * V_P4_tmp.t();

    // Переводим в координаты OpenCV
    double deltaX = (Sx - ((Pos0x - Dx)  - 1000 * pixSizeLon)) / pixSizeLon;
    double deltaY = (Sy - ((Pos0y + Dy) - 1000 * pixSizeLat)) / pixSizeLat;

    double delta = (H * pixel) / fokus;

    // Координаты вершин снимка в OpenCV
    double C0x = V_P0.at<double>(0,0) * delta / pixSizeLon + deltaX;
    double C0y = V_P0.at<double>(1,0) * delta / pixSizeLat + deltaY;

    double C1x = V_P1.at<double>(0,0) * delta / pixSizeLon + deltaX;
    double C1y = V_P1.at<double>(1,0) * delta / pixSizeLat + deltaY;

    double C2x = V_P2.at<double>(0,0) * delta / pixSizeLon + deltaX;
    double C2y = V_P2.at<double>(1,0) * delta / pixSizeLat + deltaY;

    double C3x = V_P3.at<double>(0,0) * delta / pixSizeLon + deltaX;
    double C3y = V_P3.at<double>(1,0) * delta / pixSizeLat + deltaY;

    double C4x = V_P4.at<double>(0,0) * delta / pixSizeLon + deltaX;
    double C4y = V_P4.at<double>(1,0) * delta / pixSizeLat + deltaY;

    QPolygon polygon = QPolygon()
            << QPoint(C1x, C1y)
            << QPoint(C2x, C2y)
            << QPoint(C3x, C3y)
            << QPoint(C4x, C4y)
            << QPoint(C1x, C1y)
            << QPoint(C0x, C0y)
            << QPoint(C2x, C2y);

    // Выводим на QGraphicScene
    if(k>0)
        ui->mygraphicsview->scene->addLine(point.x(), point.y(), (Sx - topLeftAngleLon)/ pixSizeLon - topleftx_px, (Sy - topLeftAngleLat)/ pixSizeLat - toplefty_px, QPen(Qt::red, 5));

    point = QPoint((Sx - topLeftAngleLon)/ pixSizeLon - topleftx_px, (Sy - topLeftAngleLat)/ pixSizeLat - toplefty_px);
    ui->mygraphicsview->polygon->setPolygon(polygon);

    // Начинаем съемку с 30 м
    if(H >= 30)
    {
        // OpenCV точки
        Point2f srcQuad[4], dstQuad[4];

        Mat lambda( 2, 4, CV_32FC1 );
        lambda = Mat::zeros( img.rows, img.cols, img.type() );

        ///////////////////////////////////////////////////

        srcQuad[0].x = C4x;  //src Top left
        srcQuad[0].y = C4y;

        srcQuad[1].x = C1x;  //src Bottom left
        srcQuad[1].y = C1y;

        srcQuad[2].x = C2x;  //src Bottom right
        srcQuad[2].y = C2y;

        srcQuad[3].x = C3x;  //src Top right
        srcQuad[3].y = C3y;

        ///////////////////////////////////////////////////
        dstQuad[0].x = 0; //Top left
        dstQuad[0].y = 0;

        dstQuad[1].x = 0; //Bottom left
        dstQuad[1].y = dst.rows - 1;

        dstQuad[2].x = dst.cols - 1; //Bottom right
        dstQuad[2].y = dst.rows - 1;

        dstQuad[3].x = dst.cols - 1; //Top right
        dstQuad[3].y = 0;
        ///////////////////////////////////////////////////

        float xMin = +FLT_MAX, yMin = +FLT_MAX;
        float xMax = -FLT_MAX, yMax = -FLT_MAX;

        for (int i = 0; i < 4; i++)
        {
            xMin = qMin(srcQuad[i].x, xMin);
            yMin = qMin(srcQuad[i].y, yMin);

            xMax = qMax(srcQuad[i].x, xMax);
            yMax = qMax(srcQuad[i].y, yMax);
        }

        // Вырезаем ИОР из большого снимка
        Rect roi(qFloor(xMin), qFloor(yMin), qCeil(xMax - xMin), qCeil(yMax - yMin));
        Mat img_small = img(roi);
        if(!img_small.empty())
        {
            //imshow(QString("ROI").toLocal8Bit().data(), img_small);

            srcQuad[0].x -= xMin;  //src Top left
            srcQuad[0].y -= yMin;

            srcQuad[1].x -= xMin;  //src Top right
            srcQuad[1].y -= yMin;

            srcQuad[2].x -= xMin;  //src Bottom left
            srcQuad[2].y -= yMin;

            srcQuad[3].x -= xMin;  //src Bot right
            srcQuad[3].y -= yMin;

        //------------------------------------------------//

            lambda = getPerspectiveTransform(srcQuad,  dstQuad);

            warpPerspective(img_small, dst, lambda, dst.size());

            //imshow("Output", dst);
            //moveWindow("Output", 100, 100);

//            Mat colorDst = dst.clone();
//            cvtColor(dst, colorDst, CV_COLOr),
//            videoWriter.write(dst);

            if(!dst.empty() && !ui->adjustment_checkbox->isChecked())
            {
                QImage qdst( dst.data, dst.cols, dst.rows, static_cast<int>(dst.step), QImage::Format_Grayscale8);
                QPixmap qdst_pixmap = QPixmap::fromImage(qdst);
                qdst_pixmap = qdst_pixmap.scaled(qdst_pixmap.size() *= ui->Scale_spinBox->value(), Qt::KeepAspectRatio);
                qdst_label->setPixmap(qdst_pixmap);
                qdst_label->resize(qdst_pixmap.size());
            }

        }
    }



    k++;

    if (k == K || H == 0.0)
    {
//        videoWriter.release();
        timer->stop();
        file_H.close();
        QMessageBox::information(this, "Информация", "Полет завершен!");
    }
}

// Маршрут 1
void MainWindow::route1_W()
{
    double dphix = ui->fx_doubleSpinBox->value() * 2 * M_PI / ui->K_spinBox->value();
    double dphiy = ui->fy_doubleSpinBox->value() * 2 * M_PI / ui->K_spinBox->value();

    phix += dphix;
    phiy += dphiy;

    double Mx = Sx;
    double My = Sy;

    // Координаты ЛА
    Sx = Dx * qSin(phix) + Pos0x;
    Sy = Dy * qSin(phiy) + Pos0y;

    S += qSqrt(qPow((Sx - Mx), 2) + qPow((Sy - My), 2));
    ui->S_doubleSpinBox->setValue(S);
    ui->S_label->setText(QString("S = %1 м").arg(S, 0, 'f', 2));

    double Y_H_tmp = Y_H;

    // Рассчитываем высоту на момент S
    if(S <= s1)
        H = (S * H_max) / s1;
    if(S > s1 && S < s2)
        H = H_max;
    if(S >= s2 && S < s3)
        H = H_max * (s3 - S) / (s3 - s2);
    if(S >= s3)
        H = 0.0;

    //ui->H_doubleSpinBox->setValue(H);
    ui->H_label->setText(QString("H = %1 м").arg(H, 0, 'f', 2));

    //ui->Sx_doubleSpinBox->setValue(Sx);
    ui->Lon_label->setText(QString("Lon = %1 м").arg(Sx, 0, 'f', 2));
    //ui->Sy_doubleSpinBox->setValue(Sy);
    ui->Lat_label->setText(QString("Lat = %1 м").arg(Sy, 0, 'f', 2));

    //Добавляем шум
    Y_H = qSqrt(1 - qPow(alpha, 2)) * X_H + alpha * Y_H_tmp;
    H_shum = H + Y_H * ui->sigma_H_doubleSpinBox->value();

    QTextStream stream_H(&file_H);
    stream_H << H_shum << "\t" <<H << "\n";

    // Путевая скорость
    W = (qSqrt(qPow((Sx - Mx), 2) + qPow((Sy - My), 2))) / ui->deltaT_doubleSpinBox->value();
    ui->w_doubleSpinBox->setValue(W);
    ui->W_label->setText(QString("W = %1 м/c").arg(W, 0, 'f', 2));


    // Курс
    double psim = psi;

    double Y_psi_tmp = Y_psi;

    psi =  5 * M_PI_2-qAtan2((Sy - My), (Sx - Mx));

    if(psi < 0)
        psi += 2 * M_PI;

    while(psi > 2*M_PI)
        psi -= 2 * M_PI;

    // Учитываем ветер
    WL = qSqrt(qPow(mW, 2) + qPow(W, 2) - 2 * mW * W * qCos(qDegreesToRadians(aW) - psi )); // продольная скорость ИЛИ СДЕЛАТЬ ЛОКАЛЬНОЙ ????
    as = qAsin(mW * qSin(qDegreesToRadians(aW) - psi) / WL); // угол сноса РАЗОБРАТЬСЯ С РАДИАНАМИ И ГРАДУСАМИ

    qDebug() << "WL = " << WL;
    qDebug() << "W = " << W;
    qDebug() << "угол сноса" << qRadiansToDegrees(as);

    //ui->Ugol_snosa_doubleSpinBox->setValue(qRadiansToDegrees(as)) ;
    ui->Ugol_snosa_label->setText(QString("Угол сноса = %1 °").arg(qRadiansToDegrees(as), 0, 'f', 2));

    psi = psi - qDegreesToRadians(as);

    //Добавляем шум
    Y_psi = qSqrt(1 - qPow(alpha, 2)) * X_psi + alpha * Y_psi_tmp;
    psi_shum = psi + Y_psi * ui->sigma_Kurs_doubleSpinBox->value();

    //ui->Kurs_doubleSpinBox->setValue(qRadiansToDegrees(psi));
    ui->Kurs_label->setText(QString("Курс = %1 °").arg(qRadiansToDegrees(psi), 0, 'f', 2));

    // Крен
    double Y_gamma_tmp = Y_gamma;
    if (k > 1)
        gamma = (psi - psim) * 3;
    //ui->Kren_doubleSpinBox->setValue(qRadiansToDegrees(gamma));
    ui->Kren_label->setText(QString("Крен = %1 °").arg(qRadiansToDegrees(gamma), 0, 'f', 2));

    //Добавляем шум
    Y_gamma = qSqrt(1 - qPow(alpha, 2)) * X_gamma + alpha * Y_gamma_tmp;
    gamma_shum = psi + Y_psi * ui->sigma_Kren_doubleSpinBox->value();

    // Тангаж

    // Рассчитываем тангаж на момент S
    if(S <= s1)
//        tetta = qDegreesToRadians(10.0);
        tetta = qAtan(H_max / s1);
        qDebug() << "H_max =" << H_max;
        qDebug() << "s1 =" << s1;
        qDebug() << "tetta = " << qRadiansToDegrees(tetta);
    if(S > s1 && S < s2)
        tetta = qDegreesToRadians(0.0);
    if(S >= s2 && S < s3)
//        tetta = qDegreesToRadians(-10.0);
        tetta = -qAtan(H_max / (s3 - s2));
    if(S >= s3)
        tetta = qDegreesToRadians(0.0);

    double Y_tetta_tmp = Y_tetta;
    ui->Tangazh_label->setText(QString("Тангаж = %1 °").arg(qRadiansToDegrees(tetta), 0, 'f', 2));

    //Добавляем шум
    Y_tetta = qSqrt(1 - qPow(alpha, 2)) * X_tetta + alpha * Y_tetta_tmp;
    tetta_shum = tetta + Y_tetta * ui->sigma_Tangazh_doubleSpinBox->value();
}

// Маршрут 2
void MainWindow::route2_W()
{
    double aX = 0.0;
    double aY = 0.0;

    Xa = qFloor(Xa + ui->fx_doubleSpinBox->value()) % K;
    Ya = qFloor(Ya + ui->fy_doubleSpinBox->value()) % K;

    if (Xa < (K / 2))
    {
        aX = Xa - K / 4;
    }
    else aX = (K - Xa) - K / 4;

    if (Ya < (K / 2))
    {
        aY = Ya - K / 4;
    }
    else aY = (K - Ya) - K / 4;

    aX = Pos0x + aX * (4 * Dx / K);
    aY = Pos0y + aY * (4 * Dy / K);

    double Mx = Sx;
    double My = Sy;

    double cS = ui->cS_doubleSpinBox->value();

    // Координаты ЛА
    Sx = Mx + cS * (aX - Mx);
    Sy = My + cS * (aY - My);

    S += qSqrt(qPow((Sx - Mx), 2) + qPow((Sy - My), 2));
    ui->S_doubleSpinBox->setValue(S);

    double Y_H_tmp = Y_H;

    // Рассчитываем высоту на момент S
    if(S <= s1)
        H = (S * H_max) / s1;
    if(S > s1 && S < s2)
        H = H_max;
    if(S >= s2 && S < s3)
        H = H_max * (s3 - S) / (s3 - s2);
    if(S >= s3)
        H = 0.0;

    //Добавляем шум
    Y_H = qSqrt(1 - qPow(alpha, 2)) * X_H + alpha * Y_H_tmp;
    H_shum = H + Y_H * ui->sigma_H_doubleSpinBox->value();

    QTextStream stream_H(&file_H);
    stream_H << H_shum << "\t" <<H << "\n";

    ui->H_doubleSpinBox->setValue(H);

    ui->Sx_doubleSpinBox->setValue(Sx);
    ui->Sy_doubleSpinBox->setValue(Sy);

    // Путевая скорость
    W = (qSqrt(qPow((Sx - Mx), 2) + qPow((Sy - My), 2))) / ui->deltaT_doubleSpinBox->value();
    ui->w_doubleSpinBox->setValue(W);

    // Курс
    double Y_psi_tmp = Y_psi;

    double psim = psi;

    psi =  5 * M_PI_2-qAtan2((Sy - My), (Sx - Mx));

    if(psi < 0)
        psi += 2 * M_PI;

    while(psi > 2*M_PI)
        psi -= 2 * M_PI;

    //Добавляем шум
    Y_psi = qSqrt(1 - qPow(alpha, 2)) * X_psi + alpha * Y_psi_tmp;
    psi_shum = psi + Y_psi * ui->sigma_Kurs_doubleSpinBox->value();

    ui->Kurs_doubleSpinBox->setValue(qRadiansToDegrees(psi));

    // Крен
    double Y_gamma_tmp = Y_gamma;
    if (k > 0) // Почему для первого маршрута k > 1 ????
        gamma = (psi - psim) * ui->c_spinBox->value();
    ui->Kren_doubleSpinBox->setValue(qRadiansToDegrees(gamma));

    //Добавляем шум
    Y_gamma = qSqrt(1 - qPow(alpha, 2)) * X_gamma + alpha * Y_gamma_tmp;
    gamma_shum = psi + Y_psi * ui->sigma_Kren_doubleSpinBox->value();
}

void MainWindow::on_X_spinBox_valueChanged(int arg1)
{
    qdst_label->move(arg1,ui->Y_spinBox->value());
}

void MainWindow::on_Y_spinBox_valueChanged(int arg1)
{
    qdst_label->move(ui->X_spinBox->value(),arg1);
}

void MainWindow::on_Scale_spinBox_valueChanged(double arg1)
{
    if(ui->adjustment_checkbox->isChecked())
    {
        QPixmap qdst_pixmap = QPixmap::fromImage(adjustmentImage);
        qdst_pixmap = qdst_pixmap.scaled(qdst_pixmap.size() *= ui->Scale_spinBox->value(), Qt::KeepAspectRatio);
        qdst_label->setPixmap(qdst_pixmap);
        qdst_label->resize(qdst_pixmap.size());
    }
    else if(!dst.empty() && ui->Pause_pushButton->isChecked())
    {
        QImage qdst( dst.data, dst.cols, dst.rows, static_cast<int>(dst.step), QImage::Format_Grayscale8);
        QPixmap qdst_pixmap = QPixmap::fromImage(qdst);
        qdst_pixmap = qdst_pixmap.scaled(qdst_pixmap.size() *= arg1, Qt::KeepAspectRatio);
        qdst_label->setPixmap(qdst_pixmap);
        qdst_label->resize(qdst_pixmap.size());
    }
}

void MainWindow::on_openButton_2_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, "Выбор эталонного изображения", QApplication::applicationDirPath(), "Изображения (*.bmp *.jpg *.jpeg *.png)");
    if(fileName.isEmpty())
        return;

    Mat smallFrame = imread(fileName.toStdString().data(), IMREAD_GRAYSCALE);
    if(smallFrame.empty())
        return;

    uint meanVal = mean(smallFrame)[0];

    adjustmentImage = QImage(640, 480, QImage::Format_Grayscale8);
    adjustmentImage.fill(meanVal);

    QImage qImgSmallFrame(smallFrame.data, smallFrame.cols, smallFrame.rows, static_cast<int>(smallFrame.step), QImage::Format_Grayscale8);
    QPainter p;
    p.begin(&adjustmentImage);
    p.drawImage((640-qImgSmallFrame.width())/2, (480-qImgSmallFrame.height())/2, qImgSmallFrame);
    p.end();

    if(ui->adjustment_checkbox->isChecked())
    {
        QPixmap qdst_pixmap = QPixmap::fromImage(adjustmentImage);
        qdst_pixmap = qdst_pixmap.scaled(qdst_pixmap.size() *= ui->Scale_spinBox->value(), Qt::KeepAspectRatio);
        qdst_label->setPixmap(qdst_pixmap);
        qdst_label->resize(qdst_pixmap.size());
    }
}

void MainWindow::on_adjustment_checkbox_clicked(bool checked)
{
    if(checked)
    {
        QPixmap qdst_pixmap = QPixmap::fromImage(adjustmentImage);
        qdst_pixmap = qdst_pixmap.scaled(qdst_pixmap.size() *= ui->Scale_spinBox->value(), Qt::KeepAspectRatio);
        qdst_label->setPixmap(qdst_pixmap);
        qdst_label->resize(qdst_pixmap.size());
    }
}
