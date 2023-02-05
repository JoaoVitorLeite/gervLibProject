#ifndef PAIR_H
#define PAIR_H

#include <vector>
#include <BasicArrayObject.h>
#include <cstdlib>

//template <class DType>
//class Pair{
//    private:
//        std::pair<BasicArrayObject<DType>, std::vector<double>> pair;

//    public:
//        Pair(BasicArrayObject<DType> obj, std::vector<double> distanceVector);

//        BasicArrayObject<DType> first();
//        std::vector<double> second();

//};


//template <class DType>
//Pair<DType>::Pair(BasicArrayObject<DType> obj, std::vector<double> distanceVector){
//    pair.first = obj;
//    pair.second = distanceVector;
//}

//template <class DType>
//BasicArrayObject<DType> Pair<DType>::first(){
//    return pair.first;
//}

//template <class DType>
//std::vector<double> Pair<DType>::second(){
//    return pair.second;
//}

template <class DType>
class Pair{
    private:
        std::pair<DType, std::vector<double>> pair;

    public:
        Pair(DType obj, std::vector<double> distanceVector);

        DType first();
        std::vector<double> second();
};

template <class DType>
Pair<DType>::Pair(DType obj, std::vector<double> distanceVector){
    pair.first = obj;
    pair.second = distanceVector;
}

template <class DType>
DType Pair<DType>::first(){
    return pair.first;
}

template <class DType>
std::vector<double> Pair<DType>::second(){
    return pair.second;
}


#endif // PAIR_H
