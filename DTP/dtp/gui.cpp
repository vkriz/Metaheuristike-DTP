#include "gui.h"
#include "population.h"
#include "solution.h"
#include "graph.h"
#include "ui_gui.h"
#include <iostream>
#include <QMessageBox>
#include <windows.h>

std::string name;
unsigned int pop, n;
int check;
double p, p_c, p_m, pm;


long long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 random_number_generator (seed);

gui::gui(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::gui)
{
    ui->setupUi(this);
}

void gui::ispis_u_text_box()
{
    Graph problem;
    std::string path=".."+name;
    problem.load(path);
    std::cout<<name<<std::endl;
    Population p1(pop,problem,2);
    for(unsigned int i=0; i<n; ++i){
        unsigned int rand=random_num(0,100);
        if(rand<(p*100)) p1.crossover(p_c);
        else p1.mutation(p_m, pm);
    }

    Solution best=p1.best();
    ui->gui_text->insert(QString::number(best.get_fitness()));
}

gui::~gui()
{
    delete ui;
}

void gui::on_gui_ok_clicked()
{
    int sve_ok=1, ok=1;
    pop=200; n=5000; check=1; p=0.8; p_c=0.8; p_m=0.8; pm=0.4;
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();    //postavimo seed (SAMO JEDNOM!)

    if(ui->gui_pop_size->text().toInt()) pop= static_cast<unsigned int>(ui->gui_pop_size->text().toInt());
    if(ui->gui_n->text().toInt()) n=static_cast<unsigned int>(ui->gui_n->text().toInt());
    if(ui->gui_p->text().toDouble()>static_cast<double>(0) || ui->gui_p->text().toDouble()<static_cast<double>(0)) p=ui->gui_p->text().toDouble();
    if(ui->gui_c->text().toDouble()>static_cast<double>(0) || ui->gui_c->text().toDouble()<static_cast<double>(0)) p_c=ui->gui_c->text().toDouble();
    if(ui->gui_m->text().toDouble()>static_cast<double>(0) || ui->gui_m->text().toDouble()<static_cast<double>(0)) p_m=ui->gui_m->text().toDouble();
    if(ui->gui_pm->text().toDouble()>static_cast<double>(0) || ui->gui_pm->text().toDouble()<static_cast<double>(0)) pm=ui->gui_pm->text().toDouble();

    std::cout<<pop<<" "<<n<<" "<<p<<" "<<p_c<<" "<<p_m<<" "<<pm<<std::endl;

    sve_ok=1; ok=1;
    if((ui->range_100->selectedItems()).size() && sve_ok && ok){
        name="/Range_100/"+ui->range_100->currentItem()->text().toStdString()+".txt";
        ispis_u_text_box();
        return;
    }

    else if((ui->range_125->selectedItems()).size() && sve_ok && ok){
        name="/Range_125/"+ui->range_125->currentItem()->text().toStdString()+".txt";
        ispis_u_text_box();
        return;
    }

    else if((ui->range_150->selectedItems()).size() && sve_ok && ok){
        name="/Range_150/"+ui->range_150->currentItem()->text().toStdString()+".txt";
        ispis_u_text_box();
        return;
    }

}

void gui::on_gui_ok_2_clicked()
{
    ui->gui_text->insert("n");
    ui->gui_text->clear();
    ui->range_100->clearSelection();
    ui->range_125->clearSelection();
    ui->range_150->clearSelection();

}
