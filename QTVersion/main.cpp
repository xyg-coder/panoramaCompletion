#include "imagecompletion.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ImageCompletion w;
	w.show();
	return a.exec();
}
