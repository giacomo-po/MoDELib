/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <QtWidgets/QApplication>
#include <QVTKOpenGLNativeWidget.h>
#include <DDqtMainWindow.h>


int main(int argc, char *argv[])
{
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat()); // needed to ensure appropriate OpenGL context is created for VTK rendering.
    QApplication a(argc, argv);
    model::DDqtMainWindow window;
    window.show();
    return a.exec();
}

