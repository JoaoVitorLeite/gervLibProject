#ifndef QUEUEITEM_H
#define QUEUEITEM_H

//#include <BasicArrayObject.h>

//template <class DType>
//class QueueItem{
//    public:
//        QueueItem(double dist, BasicArrayObject<DType> featureVector){
//            this->dist = dist;
//            this->featureVector = featureVector;
//        }

//        double dist;
//        BasicArrayObject<DType> featureVector;

//        bool operator <(const QueueItem &item) const{
//            return dist < item.dist;
//        }
//        bool operator >(const QueueItem &item) const{
//            return dist > item.dist;
//        }
//};

template <class DType>
class QueueItem{
    public:
        QueueItem(double dist, DType featureVector){
            this->dist = dist;
            this->featureVector = featureVector;
        }

        double dist;
        DType featureVector;

        bool operator <(const QueueItem &item) const{
            return dist < item.dist;
        }
        bool operator >(const QueueItem &item) const{
            return dist > item.dist;
        }
};




#endif // QUEUEITEM_H
