#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <algorithms.h>
#include <QTime>
#include <QMessageBox>
#include <time.h>

database a;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    UpdateInfo();
    UpdateInfo2();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::UpdateInfo(){
    ui->tableWidget->clear();
    ui->tableWidget->setColumnCount(6);
    ui->tableWidget->setRowCount(a.data.size());
    QTableWidgetItem * cnt = new QTableWidgetItem(tr("%1").arg("Organization:"));
    QTableWidgetItem * reg = new QTableWidgetItem(tr("%1").arg("Project:"));
    QTableWidgetItem * cit = new QTableWidgetItem(tr("%1").arg("Component:"));
    QTableWidgetItem * str = new QTableWidgetItem(tr("%1").arg("Class:"));
    QTableWidgetItem * hs = new QTableWidgetItem(tr("%1").arg("Method №:"));
    QTableWidgetItem * rm = new QTableWidgetItem(tr("%1").arg("Creator №:"));
    ui->tableWidget->setHorizontalHeaderItem(0, cnt);
    ui->tableWidget->setHorizontalHeaderItem(1, reg);
    ui->tableWidget->setHorizontalHeaderItem(2, cit);
    ui->tableWidget->setHorizontalHeaderItem(3, str);
    ui->tableWidget->setHorizontalHeaderItem(4, hs);
    ui->tableWidget->setHorizontalHeaderItem(5, rm);

    int k = 0;
    for (auto i = a.data.begin(); i != a.data.end(); i++){
        QTableWidgetItem * cnt = new QTableWidgetItem(tr("%1").arg(i->organization));
        QTableWidgetItem * reg = new QTableWidgetItem(tr("%1").arg(i->project));
        QTableWidgetItem * cit = new QTableWidgetItem(tr("%1").arg(i->component));
        QTableWidgetItem * str = new QTableWidgetItem(tr("%1").arg(i->cllass));
        QTableWidgetItem * hs = new QTableWidgetItem(tr("%1").arg(i->method));
        QTableWidgetItem * rm = new QTableWidgetItem(tr("%1").arg(i->creator));
        ui->tableWidget->setItem(k, 0, cnt);
        ui->tableWidget->setItem(k, 1, reg);
        ui->tableWidget->setItem(k, 2, cit);
        ui->tableWidget->setItem(k, 3, str);
        ui->tableWidget->setItem(k, 4, hs);
        ui->tableWidget->setItem(k, 5, rm);
        k++;
    }
}


void MainWindow::on_pushButton_clicked()
{
    a.add_rand();
    UpdateInfo();
}

void MainWindow::on_pushButton_2_clicked()
{
    if (ui->lineEdit->text() != nullptr && ui->lineEdit_2->text() != nullptr && ui->lineEdit_3->text() != nullptr && ui->lineEdit_4->text() != nullptr && ui->lineEdit_5->text() != nullptr && ui->lineEdit_6->text() != nullptr){
        a.add(ui->lineEdit->text(), ui->lineEdit_2->text(), ui->lineEdit_3->text(), ui->lineEdit_4->text(), ui->lineEdit_5->text().toInt(), ui->lineEdit_6->text().toInt());
        ui->lineEdit->clear();
        ui->lineEdit_2->clear();
        ui->lineEdit_3->clear();
        ui->lineEdit_4->clear();
        ui->lineEdit_5->clear();
        ui->lineEdit_6->clear();
        UpdateInfo();
    }
}

bool Greater(const codepart& a, const codepart& b)
{
    return a.method < b.method;
}

void MainWindow::on_pushButton_4_clicked()
{
    if (a.data.size() > 1){
        QVector<int> tmp;
        for (auto i = a.data.begin(); i != a.data.end(); i++){
            tmp.push_back(i->method);
        }
        QVector<int> tmp1 = tmp;
        QVector<int> tmp2 = tmp;
        QVector<int> tmp3 = tmp;
        QVector<int> tmp4 = tmp;

        clock_t start = clock();
        insertionSort(tmp);
        clock_t end = clock();
        double t = (double)(end - start);
        ui->time_1->setText(tr("%1t").arg(t));

        clock_t start1 = clock();
        HeapSort(tmp1);
        clock_t end1 = clock();
        double t1 = (double)(end1 - start1);
        ui->time_2->setText(tr("%1t").arg(t1));

        clock_t start2 = clock();
        MergeSort(tmp2);
        clock_t end2 = clock();
        double t2 = (double)(end2 - start2);
        ui->time_3->setText(tr("%1t").arg(t2));

        clock_t start3 = clock();
        countSort(tmp3);
        clock_t end3 = clock();
        double t3 = (double)(end3 - start3);
        ui->time_4->setText(tr("%1t").arg(t3));

        clock_t start4 = clock();
        std::sort(tmp4.begin(), tmp4.end());
        clock_t end4 = clock();
        double t4 = (double)(end4 - start4);
        ui->time_5->setText(tr("%1t").arg(t4));

        QString tmps;
        for (auto i = tmp1.begin(); i != tmp1.end(); i++){
            tmps.push_back(tr("%1 ").arg(*i));
        }
        std::sort(a.data.begin(), a.data.end(), Greater);
        UpdateInfo();
        QMessageBox::information(nullptr, "Information", tmps);
    }
}

void MainWindow::on_pushButton_room_clicked()
{
    if (a.data.size() > 1){
        QVector<int> tmp;
        for (auto i = a.data.begin(); i != a.data.end(); i++){
            tmp.push_back(i->creator);
        }
        QVector<int> tmp1 = tmp;
        QVector<int> tmp2 = tmp;
        QVector<int> tmp3 = tmp;
        QVector<int> tmp4 = tmp;

        clock_t start = clock();
        insertionSort(tmp);
        clock_t end = clock();
        double t = (double)(end - start);
        ui->time_1->setText(tr("%1t").arg(t));

        clock_t start1 = clock();
        HeapSort(tmp1);
        clock_t end1 = clock();
        double t1 = (double)(end1 - start1);
        ui->time_2->setText(tr("%1t").arg(t1));

        clock_t start2 = clock();
        MergeSort(tmp2);
        clock_t end2 = clock();
        double t2 = (double)(end2 - start2);
        ui->time_3->setText(tr("%1t").arg(t2));

        clock_t start3 = clock();
        countSort(tmp3);
        clock_t end3 = clock();
        double t3 = (double)(end3 - start3);
        ui->time_4->setText(tr("%1t").arg(t3));

        clock_t start4 = clock();
        std::sort(tmp4.begin(), tmp4.end());
        clock_t end4 = clock();
        double t4 = (double)(end4 - start4);
        ui->time_5->setText(tr("%1t").arg(t4));

        QString tmps;
        for (auto i = tmp1.begin(); i != tmp1.end(); i++){
            tmps.push_back(tr("%1 ").arg(*i));
        }
        std::sort(a.data.begin(), a.data.end());
        UpdateInfo();
        QMessageBox::information(nullptr, "Information", tmps);
    }
}

MyList a1;
MyList a2;


void MainWindow::on_pushButton_3_clicked()
{
    if (ui->spinBox->text().toInt()-1 >= 0 && ui->spinBox->text().toInt()-1 < a.data.size()){
        if (a1.implement(&a.data[ui->spinBox->text().toInt()-1]) == false){
            a1.add(&a.data[ui->spinBox->text().toInt()-1]);
        }
    }
    ui->label_14->setText(tr("Elements in a list №1: %1").arg(a1.size()));
}

void MainWindow::on_pushButton_5_clicked()
{
    if (ui->spinBox->text().toInt()-1 >= 0 && ui->spinBox->text().toInt()-1 < a.data.size()){
        if (a1.implement(&a.data[ui->spinBox->text().toInt()-1]) == true){
            a1.del(&a.data[ui->spinBox->text().toInt()-1]);
        }
    }
    ui->label_14->setText(tr("Elements in a list №1: %1").arg(a1.size()));
}

void MainWindow::on_pushButton_6_clicked()
{
    if (ui->spinBox->text().toInt()-1 >= 0 && ui->spinBox->text().toInt()-1 < a.data.size()){
        if (a1.implement(&a.data[ui->spinBox->text().toInt()-1]) == true){
            QMessageBox::information(nullptr, "Information", "True!!!");
        } else {
            QMessageBox::information(nullptr, "Information", "False!!!");
        }
    }
}

void MainWindow::on_pushButton_7_clicked()
{
    if (ui->spinBox_2->text().toInt()-1 >= 0 && ui->spinBox_2->text().toInt()-1 < a.data.size()){
        if (a2.implement(&a.data[ui->spinBox_2->text().toInt()-1]) == false){
            a2.add(&a.data[ui->spinBox_2->text().toInt()-1]);
        }
    }
    ui->label_15->setText(tr("Elements in a list №1: %1").arg(a2.size()));
}

void MainWindow::on_pushButton_9_clicked()
{
    if (ui->spinBox_2->text().toInt()-1 >= 0 && ui->spinBox_2->text().toInt()-1 < a.data.size()){
        if (a2.implement(&a.data[ui->spinBox_2->text().toInt()-1]) == true){
            a2.del(&a.data[ui->spinBox_2->text().toInt()-1]);
        }
    }
    ui->label_15->setText(tr("Elements in a list №1: %1").arg(a2.size()));
}

void MainWindow::on_pushButton_8_clicked()
{
    if (ui->spinBox_2->text().toInt()-1 >= 0 && ui->spinBox_2->text().toInt()-1 < a.data.size()){
        if (a2.implement(&a.data[ui->spinBox_2->text().toInt()-1]) == true){
            QMessageBox::information(nullptr, "Information", "True!!!");
        } else {
            QMessageBox::information(nullptr, "Information", "False!!!");
        }
    }
}

void MainWindow::on_pushButton_10_clicked()
{
    int k = a2.size();
    for (auto i = a1.head; i != nullptr; i = i->next){
        if (a2.implement(i->a) == false){
            k++;
        }
    }
    QMessageBox::information(nullptr, "Information", tr("Number of items in association of lists: %1").arg(k));
}

void MainWindow::on_pushButton_11_clicked()
{
    int k = 0;
    for (auto i = a1.head; i != nullptr; i = i->next){
        if (a2.implement(i->a) == true){
            k++;
        }
    }
    QMessageBox::information(nullptr, "Information", tr("Number of items in intersection of lists: %1").arg(k));
}

void MainWindow::on_pushButton_12_clicked()
{
    int k = a2.size();
    for (auto i = a1.head; i != nullptr; i = i->next){
        if (a2.implement(i->a) == true){
            k--;
        }
    }
    QMessageBox::information(nullptr, "Information", tr("Number of items in distinction of lists (№2 - №1): %1").arg(k));
}

void MainWindow::on_pushButton_13_clicked()
{
    int k = a2.size() + a1.size();
    for (auto i = a1.head; i != nullptr; i = i->next){
        if (a2.implement(i->a) == true){
            k = k - 2;
        }
    }
    QMessageBox::information(nullptr, "Information", tr("Number of items in symmetric distinction of lists: %1").arg(k));
}

hashmap<int, codepart*> ah;
QVector<pair<int,int>> keys;

void MainWindow::UpdateInfo2(){
    ui->tableWidget_2->clear();
    ui->tableWidget_2->setColumnCount(2);
    ui->tableWidget_2->setRowCount(keys.size());
    QTableWidgetItem * cnt = new QTableWidgetItem(tr("%1").arg("№ of item"));
    QTableWidgetItem * reg = new QTableWidgetItem(tr("%1").arg("Key"));
    ui->tableWidget_2->setHorizontalHeaderItem(0, cnt);
    ui->tableWidget_2->setHorizontalHeaderItem(1, reg);

    int k = 0;
    for (auto i = keys.begin(); i != keys.end(); i++){
        QTableWidgetItem * cnt = new QTableWidgetItem(tr("%1").arg(i->first));
        QTableWidgetItem * reg = new QTableWidgetItem(tr("%1").arg(i->second));
        ui->tableWidget_2->setItem(k, 0, cnt);
        ui->tableWidget_2->setItem(k, 1, reg);
        k++;
    }
}

void MainWindow::on_pushButton_14_clicked()
{
    if (ui->spinBox_3->text().toInt()-1 >= 0 && ui->spinBox_3->text().toInt()-1 < a.data.size()){
        if (ah.get(ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method)) == NULL){
            ah.insertNode(ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method), &a.data[ui->spinBox_3->text().toInt()-1]);
            keys.push_back(make_pair(ui->spinBox_3->text().toInt(), ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method)));
            UpdateInfo2();
            QMessageBox::information(nullptr, "Information", tr("Item added with key: %1").arg(ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method)));
        }
    }
}

void MainWindow::on_pushButton_15_clicked()
{
    if (ui->spinBox_3->text().toInt()-1 >= 0 && ui->spinBox_3->text().toInt()-1 < a.data.size()){
        if (ah.get(ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method)) != NULL){
            ah.deleteNode(ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method));
            for (auto i = 0; i < keys.size(); i++){
                if (keys[i].second == ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method)){
                    keys.erase(keys.begin()+i);
                }
            }
            UpdateInfo2();
            QMessageBox::information(nullptr, "Information", tr("Item deleted with key: %1").arg(ah.hashCode(a.data[ui->spinBox_3->text().toInt()-1].creator+a.data[ui->spinBox_3->text().toInt()-1].method)));
        }
    }
}

node *rootAVL = NULL;
vector<int> inserted;

void MainWindow::on_pushButton_17_clicked()
{
    if (ui->spinBox_4->text().toInt()-1 >= 0 && ui->spinBox_4->text().toInt()-1 < a.data.size()){
        rootAVL = insert(rootAVL, &a.data[ui->spinBox_4->text().toInt()-1]);
        inserted.push_back(ui->spinBox_4->text().toInt());
        QString tmp;
        for (auto i : inserted){
            tmp.push_back(tr("№%1 ").arg(i));
        }
        ui->label_18->setText("Inserted items: " + tmp);
        QMessageBox::information(nullptr, "Information", "Successfully inserted");
    }
}

void MainWindow::on_pushButton_16_clicked()
{
    if (ui->spinBox_4->text().toInt()-1 >= 0 && ui->spinBox_4->text().toInt()-1 < a.data.size()){
        rootAVL = delete_key(rootAVL, &a.data[ui->spinBox_4->text().toInt()-1]);
        for (auto i = 0; i < inserted.size(); i++){
            if (inserted[i] == ui->spinBox_4->text().toInt()){
                inserted.erase(inserted.begin()+i);
            }
        }
        QString tmp;
        for (auto i : inserted){
            tmp.push_back(tr("№%1 ").arg(i));
        }
        ui->label_18->setText("Inserted items: " + tmp);
        QMessageBox::information(nullptr, "Information", "Successfully deleted");
    }
}
