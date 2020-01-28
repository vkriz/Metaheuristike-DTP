TEMPLATE = app

QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += console c++11

SOURCES += \
    dtp.cpp \
    gui.cpp

HEADERS += \
    graph.h \
    population.h \
    solution.h \
    gui.h

DISTFILES += \
    Range_100/ins_050_1.txt \
    Range_100/ins_050_2.txt \
    Range_100/ins_050_3.txt \
    Range_100/ins_100_1.txt \
    Range_100/ins_100_2.txt \
    Range_100/ins_100_3.txt \
    Range_100/ins_200_1.txt \
    Range_100/ins_200_2.txt \
    Range_100/ins_200_3.txt \
    Range_100/ins_300_1.txt \
    Range_100/ins_300_2.txt \
    Range_100/ins_300_3.txt \
    Range_100/ins_400_1.txt \
    Range_100/ins_400_2.txt \
    Range_100/ins_400_3.txt \
    Range_100/ins_500_1.txt \
    Range_100/ins_500_2.txt \
    Range_100/ins_500_3.txt \
    Range_125/ins_100_1.txt \
    Range_125/ins_100_2.txt \
    Range_125/ins_100_3.txt \
    Range_125/ins_200_1.txt \
    Range_125/ins_200_2.txt \
    Range_125/ins_200_3.txt \
    Range_125/ins_300_1.txt \
    Range_125/ins_300_2.txt \
    Range_125/ins_300_3.txt \
    Range_125/ins_400_1.txt \
    Range_125/ins_400_2.txt \
    Range_125/ins_400_3.txt \
    Range_125/ins_500_1.txt \
    Range_125/ins_500_2.txt \
    Range_125/ins_500_3.txt \
    Range_125/ins_50_1.txt \
    Range_125/ins_50_2.txt \
    Range_125/ins_50_3.txt \
    Range_150/ins_100_1.txt \
    Range_150/ins_100_2.txt \
    Range_150/ins_100_3.txt \
    Range_150/ins_200_1.txt \
    Range_150/ins_200_2.txt \
    Range_150/ins_200_3.txt \
    Range_150/ins_300_1.txt \
    Range_150/ins_300_2.txt \
    Range_150/ins_300_3.txt \
    Range_150/ins_400_1.txt \
    Range_150/ins_400_2.txt \
    Range_150/ins_400_3.txt \
    Range_150/ins_500_1.txt \
    Range_150/ins_500_2.txt \
    Range_150/ins_500_3.txt \
    Range_150/ins_50_1.txt \
    Range_150/ins_50_2.txt \
    Range_150/ins_50_3.txt

FORMS += \
    gui.ui


INCLUDEPATH+= "C:\Program Files\boost\include\boost-1_69"
