#ifndef LeafNodeVPTreeVPTreeVPTREE_H
#define LeafNodeVPTreeVPTreeVPTREE_H

#include <VpTree/Node/Node.h>
#include <VpTree/Node/Bucket/Bucket.h>

//template <class DType>
//class LeafNodeVPTreeVPTree: public Node<DType>{
//    private:
//        size_t maxElements;
//        bool isCircunscribed;
//        Bucket<DType> bucket;
//        std::vector<Node<DType>*> exceededNodes;
//        std::vector<BasicArrayObject<DType>> previousPivots;
//        std::vector<std::vector<double>> distanceVector;

//    public:
//        LeafNodeVPTreeVPTree();
//        LeafNodeVPTreeVPTree(size_t paramMaxElements);

//        ~LeafNodeVPTreeVPTree();

//        void setPreviousPivots(const std::vector<BasicArrayObject<DType>> &paramPreviousPivots);
//        const std::vector<BasicArrayObject<DType>> &getPreviousPivots() const;

//        void setDistanceVector(const std::vector<std::vector<double>> &paramDistanceVector);
//        const std::vector<std::vector<double>> &getDistanceVector() const;

//        void setBucket(const Bucket<DType> &paramBucket, const size_t numberOfElements);
//        Bucket<DType> &getBucket();

//        void resetBucket();

//        void setPair(Pair<DType> pair, const size_t idx);
//        Pair<DType> getPair(const size_t idx);

//        void push_back(const Pair<DType> &pair);

//        Pair<DType> &operator[](const size_t idx);
//        const Pair<DType> &operator[](const size_t idx) const;

//        void push_exceeded_node(Node<DType> *leaf);

//        std::vector<Node<DType>*> getExceededNodes();

//        void setCircunscribed(const bool paramIsCircunscribed);
//        bool circunscribed() const;
//        size_t numberOfElements() const;
//};


//template <class DType>
//LeafNodeVPTreeVPTree<DType>::LeafNodeVPTreeVPTree()
//{
//    maxElements = 0;

//}

//template <class DType>
//LeafNodeVPTreeVPTree<DType>::LeafNodeVPTreeVPTree(size_t paramMaxElements)
//{
//    maxElements = paramMaxElements;

//}


//template <class DType>
//LeafNodeVPTreeVPTree<DType>::~LeafNodeVPTreeVPTree()
//{



//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::setPreviousPivots(const std::vector<BasicArrayObject<DType>> &paramPreviousPivots)
//{
//    previousPivots = paramPreviousPivots;
//}

//template <class DType>
//const std::vector<BasicArrayObject<DType>>& LeafNodeVPTreeVPTree<DType>::getPreviousPivots() const
//{
//    return previousPivots;
//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::setDistanceVector(const std::vector<std::vector<double>> &paramDistanceVector)
//{
//    distanceVector = paramDistanceVector;
//}

//template <class DType>
//const std::vector<std::vector<double>>& LeafNodeVPTreeVPTree<DType>::getDistanceVector() const
//{
//    return distanceVector;
//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::setBucket(const Bucket<DType> &paramBucket, const size_t numberOfElements)
//{
//    bucket = paramBucket;
//    bucket.setCurrentPosition(numberOfElements);
//}

//template <class DType>
//Bucket<DType>& LeafNodeVPTreeVPTree<DType>::getBucket()
//{
//    return bucket;
//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::resetBucket()
//{
//    bucket.reset();
//    bucket.setCurrentPosition(0);
//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::setPair(Pair<DType> pair, const size_t idx)
//{
//    bucket[idx] = pair;
//}

//template <class DType>
//Pair<DType> LeafNodeVPTreeVPTree<DType>::getPair(const size_t idx)
//{
//    return bucket[idx];
//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::push_back(const Pair<DType> &pair)
//{
//    bucket.push_back(pair);
//}

//template <class DType>
//Pair<DType>& LeafNodeVPTreeVPTree<DType>::operator[](const size_t idx)
//{
//    return bucket[idx];
//}

//template <class DType>
//const Pair<DType>& LeafNodeVPTreeVPTree<DType>::operator[](const size_t idx) const
//{
//    return bucket[idx];
//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::push_exceeded_node(Node<DType> *leaf)
//{
//    exceededNodes.push_back(leaf);
//}

//template <class DType>
//std::vector<Node<DType>*> LeafNodeVPTreeVPTree<DType>::getExceededNodes()
//{
//    return exceededNodes;
//}

//template <class DType>
//void LeafNodeVPTreeVPTree<DType>::setCircunscribed(const bool paramIsCircunscribed)
//{
//    isCircunscribed = paramIsCircunscribed;
//}

//template <class DType>
//bool LeafNodeVPTreeVPTree<DType>::circunscribed() const
//{
//    return isCircunscribed;
//}

//template <class DType>
//size_t LeafNodeVPTreeVPTree<DType>::numberOfElements() const
//{
//    return bucket.numberOfElements();
//}

template <class DType>
class LeafNodeVPTree: public Node<DType>{
    private:
        uint_fast32_t maxElements;
        bool isCircunscribed;
        Bucket<DType> bucket;
        std::vector<Node<DType>*> exceededNodes;
        std::vector<DType> previousPivots;
        std::vector<std::vector<double>> distanceVector;

    public:
        LeafNodeVPTree();
        LeafNodeVPTree(size_t paramMaxElements);

        ~LeafNodeVPTree();

        void setPreviousPivots(const std::vector<DType> &paramPreviousPivots);
        const std::vector<DType> &getPreviousPivots() const;

        void setDistanceVector(const std::vector<std::vector<double>> &paramDistanceVector);
        const std::vector<std::vector<double>> &getDistanceVector() const;

        void setBucket(const Bucket<DType> &paramBucket,
                       const uint_fast32_t numberOfElements);
        Bucket<DType> &getBucket();

        void resetBucket();

        void setPair(Pair<DType> pair, const uint_fast32_t idx);
        Pair<DType> getPair(const uint_fast32_t idx);

        void push_back(const Pair<DType> &pair);

        Pair<DType> &operator[](const uint_fast32_t idx);
        const Pair<DType> &operator[](const uint_fast32_t idx) const;

        void push_exceeded_node(Node<DType> *leaf);

        std::vector<Node<DType>*> getExceededNodes();

        void setCircunscribed(const bool paramIsCircunscribed);
        bool circunscribed() const;
        uint_fast32_t numberOfElements() const;
};

#define IMPL_TEMPL template<class DType>

IMPL_TEMPL LeafNodeVPTree<DType>::LeafNodeVPTree(){
    maxElements = 0;
}

IMPL_TEMPL LeafNodeVPTree<DType>::LeafNodeVPTree(size_t paramMaxElements){
    maxElements = paramMaxElements;
}

IMPL_TEMPL LeafNodeVPTree<DType>::~LeafNodeVPTree(){}

IMPL_TEMPL void LeafNodeVPTree<DType>::setBucket(const Bucket<DType> &paramBucket,
                                           const uint_fast32_t numberOfElements){
    bucket = paramBucket;
    bucket.setCurrentPosition(numberOfElements);
}

IMPL_TEMPL Bucket<DType> &LeafNodeVPTree<DType>::getBucket() { return bucket; }

IMPL_TEMPL void LeafNodeVPTree<DType>::resetBucket(){
    bucket.reset();
    bucket.setCurrentPosition(0);
}

IMPL_TEMPL void LeafNodeVPTree<DType>::setPreviousPivots(const std::vector<DType> &paramPreviousPivots){
    previousPivots = paramPreviousPivots;
}

IMPL_TEMPL const std::vector<DType> &LeafNodeVPTree<DType>::getPreviousPivots() const{ return previousPivots; }

IMPL_TEMPL void LeafNodeVPTree<DType>::setDistanceVector(const std::vector<std::vector<double>> &paramDistanceVector){
    distanceVector = paramDistanceVector;
}

IMPL_TEMPL const std::vector<std::vector<double> > &LeafNodeVPTree<DType>::getDistanceVector() const{ return distanceVector; }

IMPL_TEMPL void LeafNodeVPTree<DType>::setCircunscribed(const bool paramIsCircunscribed){
    isCircunscribed = paramIsCircunscribed;
}

IMPL_TEMPL bool LeafNodeVPTree<DType>::circunscribed() const{
    return isCircunscribed;
}

IMPL_TEMPL uint_fast32_t LeafNodeVPTree<DType>::numberOfElements() const{ return bucket.numberOfElements(); }

IMPL_TEMPL void LeafNodeVPTree<DType>::setPair(Pair<DType> pair, const uint_fast32_t idx){
    bucket[idx] = pair;
}

IMPL_TEMPL Pair<DType> LeafNodeVPTree<DType>::getPair(const uint_fast32_t idx){
    return bucket[idx];
}

IMPL_TEMPL void LeafNodeVPTree<DType>::push_back(const Pair<DType> &pair){
    bucket.push_back(pair);
}

IMPL_TEMPL Pair<DType> &LeafNodeVPTree<DType>::operator[](const uint_fast32_t idx){
    return bucket[idx];
}

IMPL_TEMPL const Pair<DType> &LeafNodeVPTree<DType>::operator[](const uint_fast32_t idx) const{
    return bucket[idx];
}

IMPL_TEMPL void LeafNodeVPTree<DType>::push_exceeded_node(Node<DType> *leaf){
    exceededNodes.push_back(leaf);
}

IMPL_TEMPL std::vector<Node<DType>*> LeafNodeVPTree<DType>::getExceededNodes(){
    return exceededNodes;
}



#endif // LeafNodeVPTreeVPTreeVPTREE_H
