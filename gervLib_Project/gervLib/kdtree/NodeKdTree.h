#ifndef NodeKdTreeKDTREE_H
#define NodeKdTreeKDTREE_H

#include <Dataset.h>
#include <Hermes.h>


//template <class Dtype>
//struct DynamicArray
//{

//    Dtype array[2];

//    DynamicArray(Dtype _array[2])
//    {

//        array[0] = _array[0];
//        array[1] = _array[1];

//    }

//    DynamicArray()
//    {



//    }

//};



//template<class Dtype>
//class NodeKdTree
//{

//    //Private
//    private:
//        NodeKdTree<Dtype> *left, *right, *parent;
//        std::vector<DynamicArray<Dtype>> bounds;
//        size_t nodeID;

//    //Public
//    public:

//    //Constructors and destructors
//        NodeKdTree();
//        NodeKdTree(NodeKdTree* _left, NodeKdTree* _right);
//        NodeKdTree(NodeKdTree* _left, NodeKdTree* _right, NodeKdTree* _parent);
//        NodeKdTree(NodeKdTree* _left, NodeKdTree* _right, NodeKdTree* _parent, std::vector<DynamicArray<Dtype>> _bounds);
//        virtual ~NodeKdTree();

//    //Public methods
//        bool isRoot();
//        bool equalBounds(std::vector<DynamicArray<Dtype>> other);
//        std::string showBounds();

//    //Getters
//        std::vector<DynamicArray<Dtype>> getBoundary();
//        NodeKdTree* getLeft();
//        NodeKdTree* getRight();
//        NodeKdTree* getParent();
//        size_t getNodeID();

//    //Setters
//        void setLeftNode(NodeKdTree* _left);
//        void setRightNode(NodeKdTree* _right);
//        void setParentNode(NodeKdTree* _parent);
//        void setBoundary(std::vector<DynamicArray<Dtype>> _bounds);
//        void setNodeID(size_t _nodeID);
//        void setBounds(DynamicArray<Dtype> dim, size_t pos);

//    //Virtual
//        virtual bool isDirectoryNode() = 0;
//        virtual bool isLeafNode() = 0;
//        virtual bool equals(NodeKdTree<Dtype>* other) = 0;
//        virtual std::string toString() = 0;
//        virtual u_int32_t getSerializeSize() = 0;
//        virtual u_char* serialize() = 0;
//        virtual void unserialize(u_char* data) = 0;

//};

//template<class Dtype>
//NodeKdTree<Dtype>::NodeKdTree()
//{

//    setLeftNode(nullptr);
//    setRightNode(nullptr);
//    setParentNode(nullptr);
//    setBoundary(std::vector<DynamicArray<Dtype>>());
//    setNodeID(0);

//}



//template<class Dtype>
//NodeKdTree<Dtype>::NodeKdTree(NodeKdTree* _left, NodeKdTree* _right)
//{

//    setLeftNode(_left);
//    setRightNode(_right);
//    setParentNode(nullptr);
//    setBoundary(std::vector<DynamicArray<Dtype>>());
//    setNodeID(0);

//}



//template<class Dtype>
//NodeKdTree<Dtype>::NodeKdTree(NodeKdTree* _left, NodeKdTree* _right, NodeKdTree* _parent)
//{

//    setLeftNode(_left);
//    setRightNode(_right);
//    setParentNode(_parent);
//    setBoundary(std::vector<DynamicArray<Dtype>>());
//    setNodeID(0);

//}



//template<class Dtype>
//NodeKdTree<Dtype>::NodeKdTree(NodeKdTree* _left, NodeKdTree* _right, NodeKdTree* _parent, std::vector<DynamicArray<Dtype>> _bounds)
//{

//    setLeftNode(_left);
//    setRightNode(_right);
//    setParentNode(_parent);
//    setBoundary(_bounds);
//    setNodeID(0);

//}



//template<class Dtype>
//NodeKdTree<Dtype>::~NodeKdTree<Dtype>()
//{



//}



//template<class Dtype>
//bool NodeKdTree<Dtype>::isRoot()
//{

//    return getParent() == nullptr;

//}



//template<class Dtype>
//bool NodeKdTree<Dtype>::equalBounds(std::vector<DynamicArray<Dtype>> other)
//{

//    if(other.size() != bounds.size())
//        return false;
//    else
//    {

//        for(size_t x = 0; x < bounds.size(); x++)
//        {

//            if((bounds[x].array[0] != other[x].array[0]) || (bounds[x].array[1] != other[x].array[1]))
//                return false;

//        }

//    }

//    return true;

//}



//template<class Dtype>
//std::string NodeKdTree<Dtype>::showBounds()
//{

//    std::stringstream result;

//    result << "[";

//    for(size_t x = 0; x < bounds.size(); x++)
//    {

//        result << "<" << std::to_string(bounds[x].array[0]) << ", " << std::to_string(bounds[x].array[1]) << ">";

//        if(x != (bounds.size() - 1))
//            result << ", ";

//    }

//    result << "]";

//    return result.str();

//}


//template<class Dtype>
//std::vector<DynamicArray<Dtype>> NodeKdTree<Dtype>::getBoundary()
//{

//    return bounds;

//}



//template<class Dtype>
//NodeKdTree<Dtype>* NodeKdTree<Dtype>::getLeft()
//{

//    return left;

//}



//template<class Dtype>
//NodeKdTree<Dtype>* NodeKdTree<Dtype>::getRight()
//{

//    return right;

//}



//template<class Dtype>
//NodeKdTree<Dtype>* NodeKdTree<Dtype>::getParent()
//{

//    return parent;

//}



//template<class Dtype>
//size_t NodeKdTree<Dtype>::getNodeID()
//{

//    return nodeID;

//}



//template<class Dtype>
//void NodeKdTree<Dtype>::setLeftNode(NodeKdTree* _left)
//{

//    left = _left;

//}



//template<class Dtype>
//void NodeKdTree<Dtype>::setRightNode(NodeKdTree* _right)
//{

//    right = _right;

//}



//template<class Dtype>
//void NodeKdTree<Dtype>::setParentNode(NodeKdTree* _parent)
//{

//    parent = _parent;

//}



//template<class Dtype>
//void NodeKdTree<Dtype>::setBoundary(std::vector<DynamicArray<Dtype>> _bounds)
//{

//    bounds = _bounds;

//}



//template<class Dtype>
//void NodeKdTree<Dtype>::setNodeID(size_t _nodeID)
//{

//    nodeID = _nodeID;

//}

//template<class Dtype>
//void NodeKdTree<Dtype>::setBounds(DynamicArray<Dtype> dim, size_t pos)
//{

//    if((pos < 0) || (pos > getBoundary().size()))
//    {

//        throw std::out_of_range("Out of range");

//    }
//    else
//    {

//        bounds[pos] = dim;

//    }

//}


struct Bound
{

    double array[2];

    Bound(double arr[2])
    {

        array[0] = arr[0];
        array[1] = arr[1];

    }

    Bound()
    {



    }


};


class NodeKdTree
{

    private:
        NodeKdTree *left, *right, *parent;
        std::vector<Bound> bounds;
        size_t nodeID;

    public:
        NodeKdTree(){

            setLeftNode(nullptr);
            setRightNode(nullptr);
            setParentNode(nullptr);
            setBoundary(std::vector<Bound>());
            setNodeID(0);

        }

        NodeKdTree(NodeKdTree *left_, NodeKdTree *right_)
        {

            setLeftNode(left_);
            setRightNode(right_);
            setParentNode(nullptr);
            setBoundary(std::vector<Bound>());
            setNodeID(0);

        }

        NodeKdTree(NodeKdTree *left_, NodeKdTree *right_, NodeKdTree* parent_)
        {

            setLeftNode(left_);
            setRightNode(right_);
            setParentNode(parent_);
            setBoundary(std::vector<Bound>());
            setNodeID(0);

        }

        NodeKdTree(NodeKdTree *left_, NodeKdTree *right_, NodeKdTree* parent_, std::vector<Bound> bounds_)
        {

            setLeftNode(left_);
            setRightNode(right_);
            setParentNode(parent_);
            setBoundary(bounds_);
            setNodeID(0);

        }

        ~NodeKdTree()
        {



        }

        bool isRoot()
        {

            return getParent() == nullptr;

        }

        bool equalsBounds(std::vector<Bound> other)
        {

            if(bounds.size() != other.size())
            {

                return false;

            }


            for(size_t i = 0; i < other.size(); i++)
            {

                if(bounds[i].array[0] != other[i].array[0] || bounds[i].array[1] != other[i].array[1])
                {

                    return false;

                }

            }

            return true;

        }

        std::string showBounds()
        {

            std::stringstream st;
            st << "[";

            for(size_t i = 0; i < bounds.size(); i++)
            {

                st << "<" << std::to_string(bounds[i].array[0]) << ", " << std::to_string(bounds[i].array[1]) << ">";

                if(i != (bounds.size()-1))
                {

                    st << " ,";

                }

            }

            st << "]";

            return st.str();

        }

        std::vector<Bound> getBoundary()
        {

            return bounds;

        }

        NodeKdTree* getLeft()
        {

            return left;

        }

        NodeKdTree* getRight()
        {

            return right;

        }

        NodeKdTree* getParent()
        {

            return parent;

        }

        size_t getNodeID()
        {

            return nodeID;

        }

        void setLeftNode(NodeKdTree* left_)
        {

            left = left_;

        }

        void setRightNode(NodeKdTree* right_)
        {

            right = right_;

        }

        void setParentNode(NodeKdTree* parent_)
        {

            parent = parent_;

        }

        void setBoundary(std::vector<Bound> bounds_)
        {

            bounds = bounds_;

        }

        void setNodeID(size_t id_)
        {

            nodeID = id_;

        }

        void setBound(Bound bound, size_t pos)
        {


            if((pos >= 0) && (pos < bounds.size()))
            {


                bounds[pos] = bound;

            }
            else
            {

                throw std::invalid_argument("Out of bounds !_!");

            }

        }

        virtual bool isDirectoryNode() = 0;
        virtual bool isLeafNode() = 0;
        virtual bool equals(NodeKdTree* node_) = 0;
        virtual std::string toString() = 0;
        virtual size_t getSerializedSize() = 0;
        virtual unsigned char* serialize() = 0;
        virtual void unserialize(unsigned char* data) = 0;

};



#endif // NodeKdTreeKDTREE_H
