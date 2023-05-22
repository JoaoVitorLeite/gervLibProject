#ifndef NODE_H
#define NODE_H

//#include <BasicArrayObject.h>

//template <class DType>
//class Node{
//    private:
//        BasicArrayObject<DType> pivot;
//        double mu;
//        double dMin;
//        double dMax;
//        double coverage;

//    public:
//        virtual ~Node() = default;

//        bool wasVisited = false;

//        void setPivot(const BasicArrayObject<DType> &paramPivot);
//        const BasicArrayObject<DType> &getPivot() const;

//        void setMu(const double paramMu);
//        double getMu() const;

//        void setDMin(const double paramDMin);
//        double getDMin() const;

//        void setDMax(const double paramDMax);
//        double getDMax() const;

//        void setCoverage(double maxDistFromPivot);
//        double getCoverage();

//        bool operator <(const Node &node) const{
//            return dMin < node.dMin;
//        }

//        bool operator >(const Node &node) const{
//            return dMin > node.dMin;
//        }
//};


//template <class DType>
//const BasicArrayObject<DType> &Node<DType>::getPivot() const{ return pivot; }

//template <class DType>
//void Node<DType>::setPivot(const BasicArrayObject<DType> &paramPivot){ pivot = paramPivot; }

//template <class DType>
//void Node<DType>::setMu(const double paramMu){ mu = paramMu; }

//template <class DType>
//double Node<DType>::getMu() const{ return mu; }

//template <class DType>
//void Node<DType>::setDMin(const double paramDMin){ dMin = paramDMin; }

//template <class DType>
//double Node<DType>::getDMin() const{ return dMin; }

//template <class DType>
//void Node<DType>::setDMax(const double paramDMax){ dMax = paramDMax; }

//template <class DType>
//double Node<DType>::getDMax() const{ return dMax; }

//template <class DType>
//void Node<DType>::setCoverage(double maxDistFromPivot) { coverage = maxDistFromPivot; }

//template <class DType>
//double Node<DType>::getCoverage() { return coverage; }

#include <limits.h>

template <class DType>
class Node{
    private:
        DType pivot;
        double mu;
        double dMin;
        double dMax;
        double coverage;
        size_t pageID = std::numeric_limits<size_t>::max();

    public:
        virtual ~Node() = default;

        bool wasVisited = false;

        void setPivot(const DType &paramPivot);
        const DType &getPivot() const;

        void setMu(const double paramMu);
        double getMu() const;

        void setDMin(const double paramDMin);
        double getDMin() const;

        void setDMax(const double paramDMax);
        double getDMax() const;

        void setCoverage(double maxDistFromPivot);
        double getCoverage();

        void setPageID(size_t page)
        {

            pageID = page;

        }

        size_t getPageID()
        {

            return pageID;

        }

        bool operator <(const Node &node) const{
            //return dMin < node.dMin;
            if(dMin != node.dMin)
                return dMin < node.dMin;
            else
                return dMax < node.dMax;
        }

        bool operator >(const Node &node) const{
            //return dMin > node.dMin;
            if(dMin != node.dMin)
                return dMin > node.dMin;
            else
                return dMax > node.dMax;
        }
};

#undef IMPL_TEMPL
#define IMPL_TEMPL template<class DType>

IMPL_TEMPL const DType &Node<DType>::getPivot() const{ return pivot; }
IMPL_TEMPL void Node<DType>::setPivot(const DType &paramPivot){ pivot = paramPivot; }

IMPL_TEMPL void Node<DType>::setMu(const double paramMu){ mu = paramMu; }
IMPL_TEMPL double Node<DType>::getMu() const{ return mu; }

IMPL_TEMPL void Node<DType>::setDMin(const double paramDMin){ dMin = paramDMin; }
IMPL_TEMPL double Node<DType>::getDMin() const{ return dMin; }

IMPL_TEMPL void Node<DType>::setDMax(const double paramDMax){ dMax = paramDMax; }
IMPL_TEMPL double Node<DType>::getDMax() const{ return dMax; }

IMPL_TEMPL void Node<DType>::setCoverage(double maxDistFromPivot) { coverage = maxDistFromPivot; }
IMPL_TEMPL double Node<DType>::getCoverage() { return coverage; }



#endif // NODE_H
