#ifndef BUCKET_H
#define BUCKET_H

#include <cstdint>
#include <vector>
#include <cstdlib>
#include <VpTree/Node/Bucket/Pair.h>
//#include <BasicArrayObject.h>

//template <class DType>
//class Bucket {
//    private:
//        std::vector<Pair<BasicArrayObject<DType>>> bucket;
//        uint_fast32_t currentPosition;

//    public:
//        Bucket();

//        void reset(){ bucket.clear(); }

//        void setPair(const Pair<BasicArrayObject<DType>> &pair, const uint_fast32_t idx);
//        Pair<BasicArrayObject<DType>> getPair(const uint_fast32_t idx);

//        void push_back(const Pair<BasicArrayObject<DType>> &pair);

//        Pair<BasicArrayObject<DType>> &operator[](const uint_fast32_t idx);
//        const Pair<BasicArrayObject<DType>> &operator[](const uint_fast32_t idx) const;

//        uint_fast32_t numberOfElements() const;

//        void setCurrentPosition(uint_fast32_t paramCurrentPosition){
//            currentPosition = paramCurrentPosition;
//        }

//};

//template <class DType>
//Bucket<DType>::Bucket(){
//    currentPosition = 0;
//}

//template <class DType>
//uint_fast32_t Bucket<DType>::numberOfElements() const{
//    return currentPosition;
//}

//template <class DType>
//void Bucket<DType>::setPair(const Pair<BasicArrayObject<DType>> &pair, const uint_fast32_t idx){
//    bucket[idx] = pair;
//}

//template <class DType>
//Pair<BasicArrayObject<DType>> Bucket<DType>::getPair(const uint_fast32_t idx){
//    return bucket[idx];
//}

//template <class DType>
//void Bucket<DType>::push_back(const Pair<BasicArrayObject<DType>> &pair){
//    bucket.push_back(pair);
//    currentPosition++;
//}

//template <class DType>
//Pair<BasicArrayObject<DType>> &Bucket<DType>::operator[](const uint_fast32_t idx){
//    return bucket[idx];
//}

//template <class DType>
//const Pair<BasicArrayObject<DType>> &Bucket<DType>::operator[](const uint_fast32_t idx) const{
//    return bucket[idx];
//}

template <class DType>
class Bucket {
    private:
        std::vector<Pair<DType>> bucket;
        uint_fast32_t currentPosition;

    public:
        Bucket();

        void reset(){ bucket.clear(); }

        void setPair(const Pair<DType> &pair, const uint_fast32_t idx);
        Pair<DType> getPair(const uint_fast32_t idx);

        void push_back(const Pair<DType> &pair);

        Pair<DType> &operator[](const uint_fast32_t idx);
        const Pair<DType> &operator[](const uint_fast32_t idx) const;

        uint_fast32_t numberOfElements() const;

        void setCurrentPosition(uint_fast32_t paramCurrentPosition){
            currentPosition = paramCurrentPosition;
        }

};

template <class DType>
Bucket<DType>::Bucket(){
    currentPosition = 0;
}

template <class DType>
uint_fast32_t Bucket<DType>::numberOfElements() const{
    return currentPosition;
}

template <class DType>
void Bucket<DType>::setPair(const Pair<DType> &pair, const uint_fast32_t idx){
    bucket[idx] = pair;
}

template <class DType>
Pair<DType> Bucket<DType>::getPair(const uint_fast32_t idx){
    return bucket[idx];
}

template <class DType>
void Bucket<DType>::push_back(const Pair<DType> &pair){
    bucket.push_back(pair);
    currentPosition++;
}

template <class DType>
Pair<DType> &Bucket<DType>::operator[](const uint_fast32_t idx){
    return bucket[idx];
}

template <class DType>
const Pair<DType> &Bucket<DType>::operator[](const uint_fast32_t idx) const{
    return bucket[idx];
}



#endif // BUCKET_H
