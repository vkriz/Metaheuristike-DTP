#ifndef GUI_H
#define GUI_H

#include <QWidget>


namespace Ui {
class gui;
}

class gui : public QWidget
{
    Q_OBJECT

public:
    explicit gui(QWidget *parent = nullptr);
    void ispis_u_text_box();
    ~gui();

private slots:
    void on_gui_ok_clicked();

    void on_gui_ok_2_clicked();

private:
    Ui::gui *ui;
};

#endif // GUI_H
