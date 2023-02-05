#ifndef DIRECTORYNODEKDTREE_H
#define DIRECTORYNODEKDTREE_H

#include <NodeKdTree.h>

//template<class Dtype>
//class DirectoryNodeKdTree : public NodeKdTree<Dtype>
//{

//    //Private
//    private:
//        size_t s_bounds_size;

//        //Getters
//        size_t getBoundsSize();

//        //Setters
//        void setSerializeVariables(size_t _bounds_size);


//    //Public
//    public:
//        DirectoryNodeKdTree();
//        DirectoryNodeKdTree(NodeKdTree<Dtype>* _parent, std::vector<DynamicArray<Dtype>> _bounds);
//        ~DirectoryNodeKdTree();

//        //Virtual
//        bool isDirectoryNode();
//        bool isLeafNode();
//        bool equals(NodeKdTree<Dtype>* other);
//        std::string toString();
//        u_int32_t getSerializeSize();
//        u_char* serialize();
//        void unserialize(u_char* data);


//};



//template<class Dtype>
//size_t DirectoryNodeKdTree<Dtype>::getBoundsSize()
//{

//    return s_bounds_size;

//}



//template<class Dtype>
//void DirectoryNodeKdTree<Dtype>::setSerializeVariables(size_t _s_bounds_size)
//{

//    s_bounds_size = _s_bounds_size;

//}



//template<class Dtype>
//DirectoryNodeKdTree<Dtype>::DirectoryNodeKdTree()
//{

//    this->setLeftNode(nullptr);
//    this->setRightNode(nullptr);
//    this->setParentNode(nullptr);
//    this->setBoundary(std::vector<DynamicArray<Dtype>>());
//    this->setNodeID(0);

//}



//template<class Dtype>
//DirectoryNodeKdTree<Dtype>::DirectoryNodeKdTree(NodeKdTree<Dtype>* _parent, std::vector<DynamicArray<Dtype>> _bounds)
//{

//    this->setLeftNode(nullptr);
//    this->setRightNode(nullptr);
//    this->setParentNode(_parent);
//    this->setBoundary(_bounds);
//    this->setNodeID(0);

//}



//template<class Dtype>
//DirectoryNodeKdTree<Dtype>::~DirectoryNodeKdTree()
//{

//    if(this->getBoundary().size() > 0)
//    {

//        this->getBoundary().clear();

//    }

//}



//template<class Dtype>
//bool DirectoryNodeKdTree<Dtype>::isDirectoryNode()
//{

//    return true;

//}



//template<class Dtype>
//bool DirectoryNodeKdTree<Dtype>::isLeafNode()
//{

//    return false;

//}



//template<class Dtype>
//bool DirectoryNodeKdTree<Dtype>::equals(NodeKdTree<Dtype>* other)
//{

//    bool ans = true;

//    if(other->isLeafNode())
//    {

//        ans = false;

//    }
//    else
//    {

//        ans = this->equalBounds(other->getBoundary());

//    }

//    return ans;

//}



//template<class Dtype>
//std::string DirectoryNodeKdTree<Dtype>::toString()
//{

//    std::stringstream st;
//    st << "Node ID = " << this->getNodeID() << "\n";
//    st << "Bounds : " << this->showBounds() << "\n";
//    return st.str();

//}



//template<class Dtype>
//u_int32_t DirectoryNodeKdTree<Dtype>::getSerializeSize()
//{

//    setSerializeVariables(this->getBoundary().size());
//    return sizeof(u_int32_t) * 2 + getBoundsSize() * (sizeof(Dtype) * 2);

//}



//template<class Dtype>
//u_char* DirectoryNodeKdTree<Dtype>::serialize()
//{

//    u_char* data = nullptr;

//    if(data == nullptr)
//    {

//        data = new u_char[getSerializeSize()];
//        u_int32_t nodeID = (u_int32_t)this->getNodeID(), boundsSize = (u_int32_t)s_bounds_size, total = 0;

//        memcpy(data + total, &nodeID, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);
//        memcpy(data + total, &boundsSize, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);

//        for(size_t x = 0; x < boundsSize; x++)
//        {

//            Dtype arr[2];
//            arr[0] = this->getBoundary()[x].array[0];
//            arr[1] = this->getBoundary()[x].array[1];

//            memcpy(data + total, &arr[0], sizeof(Dtype));
//            total += sizeof(Dtype);
//            memcpy(data + total, &arr[1], sizeof(Dtype));
//            total += sizeof(Dtype);

//        }

//    }

//    return data;

//}



//template<class Dtype>
//void DirectoryNodeKdTree<Dtype>::unserialize(u_char* data)
//{

//    u_int32_t nodeID = 0, boundsSize = 0, total = 0;

//    memcpy(&nodeID, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&boundsSize, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);

//    this->setNodeID((size_t)nodeID);
//    s_bounds_size = (size_t)boundsSize;
//    setSerializeVariables(getBoundsSize());

//    std::vector<DynamicArray<Dtype>> vecBounds;

//    for(size_t x = 0; x < getBoundsSize(); x++)
//    {

//        Dtype arr[2];
//        memcpy(&arr[0], data + total, sizeof(Dtype));
//        total += sizeof(Dtype);
//        memcpy(&arr[1], data + total, sizeof(Dtype));
//        total += sizeof(Dtype);

//        vecBounds.push_back(DynamicArray<Dtype>(arr));

//    }

//    this->setBoundary(vecBounds);

//}


class DirectoryNodeKdTree : public NodeKdTree
{

public:
    DirectoryNodeKdTree()
    {

        this->setLeftNode(nullptr);
        this->setRightNode(nullptr);
        this->setParentNode(nullptr);
        this->setBoundary(std::vector<Bound>());
        this->setNodeID(0);

    }

    DirectoryNodeKdTree(NodeKdTree* parent_, std::vector<Bound> bounds_)
    {

        this->setLeftNode(nullptr);
        this->setRightNode(nullptr);
        this->setParentNode(parent_);
        this->setBoundary(bounds_);
        this->setNodeID(0);

    }

    ~DirectoryNodeKdTree()
    {

        if(this->getBoundary().size() > 0)
        {

            this->getBoundary().clear();

        }

    }

    bool isDirectoryNode()
    {

        return true;

    }

    bool isLeafNode()
    {

        return false;

    }

    bool equals(NodeKdTree* node_)
    {

        if(node_->isLeafNode())
        {

            return false;

        }


        return this->equalsBounds(node_->getBoundary());

    }


    std::string toString()
    {

        std::stringstream st;
        st << "Node ID = " << this->getNodeID() << "\n";
        st << "Bounds = " << this->showBounds() << "\n";
        return st.str();

    }

    size_t getSerializedSize()
    {

        return 2 * sizeof(size_t) + this->getBoundary().size() * (sizeof(double) * 2);

    }

    unsigned char* serialize()
    {

        unsigned char* data = new unsigned char[getSerializedSize()];

        size_t total = 0, boundsSize = this->getBoundary().size(), node_id = this->getNodeID();

        memcpy(data + total, &node_id, sizeof(size_t));
        total += sizeof(size_t);
        memcpy(data + total, &boundsSize, sizeof(size_t));
        total += sizeof(size_t);

        for(size_t i = 0; i < boundsSize; i++)
        {

            double arr[2];
            arr[0] = this->getBoundary()[i].array[0];
            arr[1] = this->getBoundary()[i].array[1];

            memcpy(data + total, &arr[0], sizeof(double));
            total += sizeof(double);
            memcpy(data + total, &arr[1], sizeof(double));
            total += sizeof(double);

        }

        return data;

    }

    void unserialize(unsigned char* data)
    {

        size_t total = 0, boundsSize = 0, node_id = 0;

        memcpy(&node_id, data + total, sizeof(size_t));
        total += sizeof(size_t);
        memcpy(&boundsSize, data + total, sizeof(size_t));
        total += sizeof(size_t);

        this->setNodeID(node_id);
        this->setBoundary(std::vector<Bound>(boundsSize));

        for(size_t i = 0; i < boundsSize; i++)
        {

            double arr[2];

            memcpy(&arr[0], data + total, sizeof(double));
            total += sizeof(double);
            memcpy(&arr[1], data + total, sizeof(double));
            total += sizeof(double);

            this->setBound(Bound(arr), i);

        }

    }


};


#endif // DIRECTORYNODEKDTREE_H
