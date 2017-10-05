#include "mygraphicsview.h"

myGraphicsView::myGraphicsView(QWidget *parent)
	: QGraphicsView(parent)
{
	setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
	setResizeAnchor(QGraphicsView::AnchorUnderMouse);
	setDragMode(QGraphicsView::ScrollHandDrag);

    scene = new QGraphicsScene(this);
    setScene(scene);
}

myGraphicsView::~myGraphicsView()
{

}

void myGraphicsView::wheelEvent(QWheelEvent *event)
{
	if (event->delta() > 0)
		scale(1.1, 1.1);
	else
		scale(1 / 1.1, 1 / 1.1);
}
