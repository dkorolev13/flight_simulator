#ifndef MYGRAPHICSVIEW_H
#define MYGRAPHICSVIEW_H

#include <QGraphicsView>
#include <QWheelEvent>
#include <QGraphicsPolygonItem>

class myGraphicsView : public QGraphicsView
{
	Q_OBJECT

public:
	myGraphicsView(QWidget *parent);
    QGraphicsScene *scene;
    QGraphicsPolygonItem *polygon;
	~myGraphicsView();

protected:
	/*virtual*/ void wheelEvent(QWheelEvent *event);

private:
	
};

#endif // MYGRAPHICSVIEW_H
