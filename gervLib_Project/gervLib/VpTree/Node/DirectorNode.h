#ifndef DIRECTORNODE_H
#define DIRECTORNODE_H

#include <VpTree/Node/Node.h>
//#include <Histogram/PivotHistogram.h>

//template <class DType>
//class DirectorNode: public Node<DType>{
//    private:
//        Node<DType> *rightNode;
//        Node<DType> *leftNode;
////        PivotHistogram histogram;

//    public:
//        DirectorNode();
//        ~DirectorNode(){ delete(rightNode); rightNode = nullptr; delete(leftNode); leftNode = nullptr; }

////        void setHistogram(const PivotHistogram &paramHistogram);
////        const PivotHistogram &getHistogram() const;

//        void setRightNode(Node<DType> *paramRightNode);
//        Node<DType> *getRightNode() const;

//        void setLeftNode(Node<DType> *paramLeftNode);
//        Node<DType> *getLeftNode() const;
//};

//template<class DType>
//DirectorNode<DType>::DirectorNode(){
//    rightNode = nullptr;
//    leftNode = nullptr;
//}

//template<class DType>
//void DirectorNode<DType>::setHistogram(const PivotHistogram &paramHistogram){ histogram = paramHistogram; }

//template<class DType>
//const PivotHistogram &DirectorNode<DType>::getHistogram() const{ return histogram; }

//template<class DType>
//void DirectorNode<DType>::setRightNode(Node<DType> *paramRightNode){ rightNode = paramRightNode; }

//template<class DType>
//Node<DType> *DirectorNode<DType>::getRightNode() const{ return rightNode; }

//template<class DType>
//void DirectorNode<DType>::setLeftNode(Node<DType> *paramLeftNode){ leftNode = paramLeftNode; }

//template<class DType>
//Node<DType> *DirectorNode<DType>::getLeftNode() const{ return leftNode; }

template <class DType>
class DirectorNode: public Node<DType>{
    private:
        Node<DType> *rightNode;
        Node<DType> *leftNode;
//        PivotHistogram histogram;

    public:
        DirectorNode();
        ~DirectorNode(){ delete(rightNode); rightNode = nullptr; delete(leftNode); leftNode = nullptr; }

//        void setHistogram(const PivotHistogram &paramHistogram);
//        const PivotHistogram &getHistogram() const;

        void setRightNode(Node<DType> *paramRightNode);
        Node<DType> *getRightNode() const;

        void setLeftNode(Node<DType> *paramLeftNode);
        Node<DType> *getLeftNode() const;
};

#define IMPL_TEMPL template<class DType>

IMPL_TEMPL DirectorNode<DType>::DirectorNode(){
    rightNode = nullptr;
    leftNode = nullptr;
}

//IMPL_TEMPL void DirectorNode<DType>::setHistogram(const PivotHistogram &paramHistogram){ histogram = paramHistogram; }
//IMPL_TEMPL const PivotHistogram &DirectorNode<DType>::getHistogram() const{ return histogram; }

IMPL_TEMPL void DirectorNode<DType>::setRightNode(Node<DType> *paramRightNode){ rightNode = paramRightNode; }
IMPL_TEMPL Node<DType> *DirectorNode<DType>::getRightNode() const{ return rightNode; }

IMPL_TEMPL void DirectorNode<DType>::setLeftNode(Node<DType> *paramLeftNode){ leftNode = paramLeftNode; }
IMPL_TEMPL Node<DType> *DirectorNode<DType>::getLeftNode() const{ return leftNode; }

#endif // DIRECTORNODE_H
