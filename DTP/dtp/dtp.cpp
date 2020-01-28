#include "gui.h"
#include <iostream>
#include <QApplication>
#include <QDir>
#include <QTime>
#include <fstream>

extern std::string name;
extern unsigned int pop, n;
extern int check;
extern double p, p_c, p_m, pm;

int main(int argc, char* argv[]){

    QApplication app(argc, argv);

    gui w;
    w.show();
    app.exec();

    return 0;
}
