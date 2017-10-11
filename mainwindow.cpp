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
//#include <QKeyEvent>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Создаем таймер
    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(on_Calculation()));

    // Обнуляем CV Mat'ы картинок
    img = NULL;
    dst = NULL;  
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
    img = imread(tifFileName.toLocal8Bit().data(), cv::IMREAD_LOAD_GDAL | cv::IMREAD_GRAYSCALE );

    //Координата правого верхнего угла эталдона по долготе в метрах
    double topRightAngleLon = topLeftAngleLon + (img.cols) * pixSizeLon;

    //Координата левого нижнего угла эталона по широте в метрах
    double bottomLeftAngleLat = topLeftAngleLat + (img.rows) * pixSizeLat;

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
void MainWindow::on_RUNButton_clicked()
{   
    file_H.setFileName("D:/Projects/H.txt");
    file_H.open(QFile::WriteOnly);

    // Запуск qsrand
    qsrand(QTime(0,0,0).secsTo(QTime::currentTime()));

    ui->RUNButton->setEnabled(false);
    if(img.empty())
        return;

    k = 0;

    // Считываем данные из интерфейса
    // Положение ЛА
    gamma = qDegreesToRadians(ui->Kren_doubleSpinBox->value());
    psi = qDegreesToRadians(ui->Kurs_doubleSpinBox->value());
    tetta = qDegreesToRadians(ui->Tangazh_doubleSpinBox->value());

    // Позиция центра площади обхода
    Pos0x = ui->Lon_doubleSpinBox->value();
    Pos0y = ui->Lat_doubleSpinBox->value();

//    // Постоянная высота полета
    //double H = ui->H_doubleSpinBox->value();

    // Зависимость H от S (убрать лишние переменные)
    s0 = ui->s0_doubleSpinBox->value();
    H = ui->h0_doubleSpinBox->value();

    s1 = ui->s1_doubleSpinBox->value();
    h1 = ui->h1_doubleSpinBox->value();

    s2 = ui->s2_doubleSpinBox->value();
    s3 = ui->s3_doubleSpinBox->value();



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

    // Начальная фаза
    phix = qDegreesToRadians(ui->phix0_doubleSpinBox->value());
    phiy = qDegreesToRadians(ui->phiy0_doubleSpinBox->value());

    // Путевая скорость
    W = ui->w_doubleSpinBox->value();

    // Стартовая позиция
    Sx = Dx * qSin(phix) + Pos0x;
    Sy = Dy * qSin(phiy) + Pos0y;

    // Пройденный путь
    S = ui->S_doubleSpinBox->value();

    Xa = (K * phix / (2 * M_PI)) + K / 4;
    Ya = (K * phiy / (2 * M_PI)) + K / 4;


    // Координаты левого верхнего угла в пикселях с учетом рамки 1000 px
    topleftx_px = ((Pos0x - Dx) - topLeftAngleLon) / pixSizeLon - 1000;
    toplefty_px = ((Pos0y + Dy) - topLeftAngleLat) / pixSizeLat - 1000;

    double width_px = 2*Dx / pixSizeLon + 2000;
    double height_px = -(2*Dy / pixSizeLat) + 2000;

    // Вывод площади обхода
    Rect roi(topleftx_px, toplefty_px, width_px, height_px);
    img = img(roi);
    QImage qimg( img.data, img.cols, img.rows, static_cast<int>(img.step), QImage::Format_Grayscale8);

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

    // Создаем полотно и оси графика W
    ui->customplot_W->addGraph();
    ui->customplot_W->xAxis->setRange(0, K);
    ui->customplot_W->yAxis->setRange(0, 50); // ПОМЕНЯТЬ НА ПЕРЕМЕННУЮ

    // Создаем полотно и оси графика S
    ui->customplot_S->addGraph();
    ui->customplot_S->xAxis->setRange((Pos0x - Dx) + 20, (Pos0x + Dx) + 20);
    ui->customplot_S->yAxis->setRange((Pos0y - Dy) + 20, (Pos0y + Dy) + 20);
    ui->customplot_S->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->customplot_S->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::blue, 2));

    // Создаем полотно и оси графика psi (курса)
    ui->customplot_Kurs->addGraph();
    ui->customplot_Kurs->xAxis->setRange(0, K);
    ui->customplot_Kurs->yAxis->setRange(0, 360);
    ui->customplot_Kurs->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->customplot_Kurs->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::blue, 2));

    // Создаем полотно и оси графиков tetta(крена) и gamma(тангажа)
    ui->customplot_Kren->addGraph();
    ui->customplot_Kren->xAxis->setRange(0, K);
    ui->customplot_Kren->yAxis->setRange(-5, 5);
    ui->customplot_Kren->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->customplot_Kren->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::blue, 2));
    ui->customplot_Kren->addGraph();
    ui->customplot_Kren->graph(1)->setLineStyle(QCPGraph::lsNone);
    ui->customplot_Kren->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::red, Qt::red, 2));

    // Запускаем таймер
    timer->start(ui->deltaT_doubleSpinBox->value() * 1000);
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
    double alpha = qExp(-ui->deltaT_doubleSpinBox->value() / ui->Tc_SpinBox->value());

    double rnd_H = qrand() * 1.0 / RAND_MAX;
    double rnd_gamma = qrand() * 1.0 / RAND_MAX;
    double rnd_psi = qrand() * 1.0 / RAND_MAX;
    double rnd_tetta = qrand() * 1.0 / RAND_MAX;

    double X_H = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_H)) * qCos(2 * M_PI * rnd_H);
    double X_gamma = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_gamma)) * qCos(2 * M_PI * rnd_gamma);
    double X_psi = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_psi)) * qCos(2 * M_PI * rnd_psi);
    double X_tetta = ui->sigma_doubleSpinBox->value() * qSqrt(-2.0 * qLn(rnd_tetta)) * qCos(2 * M_PI * rnd_tetta);



    qDebug() << "rnd_H = " << rnd_H;
    qDebug() << "X_H = " << X_H;

    // Первый маршрут
    if (ui->Route_1_radioButton->isChecked())
    {
        double dFx = (ui->dax_doubleSpinBox->value() * 2 * M_PI) / K;
        double dFy = (ui->day_doubleSpinBox->value() * 2 * M_PI) / K;

        phix += dFx;
        phiy += dFy;

        double Mx = Sx;
        double My = Sy;

        // Координаты ЛА
        Sx = Dx * qSin(phix) + Pos0x;
        Sy = Dy * qSin(phiy) + Pos0y;

        S += qSqrt(qPow((Sx - Mx), 2) + qPow((Sy - My), 2));
        ui->S_doubleSpinBox->setValue(S);


        ui->Sx_doubleSpinBox->setValue(Sx);
        ui->Sy_doubleSpinBox->setValue(Sy);

        // Путевая скорость
        W = (qSqrt(qPow((Sx - Mx), 2) + qPow((Sy - My), 2))) / ui->deltaT_doubleSpinBox->value();
        ui->w_doubleSpinBox->setValue(W);

        // Курс
        double psim = psi;
        psi =  5 * M_PI_2-qAtan2((Sy - My), (Sx - Mx));

        if(psi < 0)
            psi += 2 * M_PI;

        while(psi > 2*M_PI)
            psi -= 2 * M_PI;

        ui->Kurs_doubleSpinBox->setValue(qRadiansToDegrees(psi));

        // Крен
        if (k > 1)
            gamma = (psi - psim) * 3;
        ui->Kren_doubleSpinBox->setValue(qRadiansToDegrees(gamma));
    }

    // Второй маршрут
    else
    {
        double aX = 0.0;
        double aY = 0.0;

        Xa = qFloor(Xa + ui->dax_doubleSpinBox->value()) % K;
        Ya = qFloor(Ya + ui->day_doubleSpinBox->value()) % K;

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
            H = (S * h1) / s1;
        if(S > s1 && S < s2)
            H = h1;
        if(S >= s2 && S < s3)
            H = h1 * (s3 - S) / (s3 - s2);
        if(S >= s3)
            H = 0.0;

       // qDebug() << "S = " << S;
        qDebug() << "H = " << H;

        //Добавляем шум
        Y_H = qSqrt(1 - qPow(alpha, 2)) * X_H + alpha * Y_H_tmp;
        H_shum = H + Y_H * ui->sigma_H_doubleSpinBox->value();

        QTextStream stream_H(&file_H);
        stream_H << H_shum << "\t" <<H << "\n";

//        QFile file("C:/Dev/log.txt");
//        if (file.open(QIODevice::Append)) {
//           file.write(line);
//        }
//        file.close();

        qDebug() << "Шум = " << Y_H;
        qDebug() << "H зашумленный = " << H_shum;

        ui->H_doubleSpinBox_2->setValue(H);

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

    QTime curTime(0, 0);
    curTime = curTime.addMSecs(k*0.2*1000);
    ui->time_label->setText(curTime.toString("mm:ss:zzz"));

    // Вывод значений на графики
    W_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    W_Y.push_back(W);

    S_X.push_back(Sx);
    S_Y.push_back(Sy);

    Kurs_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    Kurs_Y.push_back(psi * 180/M_PI);

    Kren_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    Kren_Y.push_back(gamma * 180/M_PI);

    Tangazh_X.push_back(k * ui->deltaT_doubleSpinBox->value());
    Tangazh_Y.push_back(tetta * 180/M_PI);

    ui->customplot_W->graph(0)->setData(W_X, W_Y);
    ui->customplot_S->graph(0)->setData(S_X, S_Y);
    ui->customplot_Kurs->graph(0)->setData(Kurs_X, Kurs_Y);
    ui->customplot_Kren->graph(0)->setData(Kren_X, Kren_Y);
    ui->customplot_Kren->graph(1)->setData(Tangazh_X, Tangazh_Y);

    ui->customplot_W->replot();
    ui->customplot_S->replot();
    ui->customplot_Kurs->replot();
    ui->customplot_Kren->replot();


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
        Rect roi(xMin, yMin, xMax - xMin, yMax - yMin);
        Mat img_small = img(roi);
        imshow(QString("ROI").toLocal8Bit().data(), img_small);

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

        imshow("Output", dst);

        QImage qdst( dst.data, dst.cols, dst.rows, static_cast<int>(dst.step), QImage::Format_Grayscale8);
        ui->qdst_label->setPixmap(QPixmap::fromImage(qdst));

    }

    k++;

    if (k == K || H == 0.0)
    {
        timer->stop();
        file_H.close();
        QMessageBox::information(this, "Информация", "Полет завершен!");
    }
}
