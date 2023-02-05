#ifndef LEAFNODEKDTREE_H
#define LEAFNODEKDTREE_H

#include <NodeKdTree.h>

//template<class Dtype>
//class LeafNodeKdTree : public NodeKdTree<Dtype>
//{

//    //Private
//    private:
//        Dataset<Dtype>* dataset;
//        size_t s_dataset_cardinality, s_dataset_dimensionality, s_dataset_seed, s_bounds_size;

//        //Getters
//        size_t getDatasetCardinality();
//        size_t getDatasetDimensionality();
//        size_t getDatasetSeed();
//        size_t getBoundsSize();

//        //Setters
//        void setSerializeVariables(size_t _s_dataset_cardinality, size_t _s_dataset_dimensionality, size_t _s_dataset_seed, size_t _s_bounds_size);


//    //Public
//    public:

//        //Constructors and destructors
//        LeafNodeKdTree();
//        LeafNodeKdTree(Dataset<Dtype>* _dataset);
//        LeafNodeKdTree(Dataset<Dtype>* _dataset, NodeKdTree<Dtype>* _parent, std::vector<DynamicArray<Dtype>> _bounds);
//        ~LeafNodeKdTree();

//       //Getters
//        Dataset<Dtype>* getDataset();
//        BasicArrayObject<Dtype>* getInstance(size_t pos);

//       //Setters
//        void setDataset(Dataset<Dtype>* _dataset);

//       //Virtual
//        bool isDirectoryNode();
//        bool isLeafNode();
//        bool equals(NodeKdTree<Dtype>* other);
//        std::string toString();
//        u_int32_t getSerializeSize();
//        u_char* serialize();
//        void unserialize(u_char* data);

//};



//template<class Dtype>
//size_t LeafNodeKdTree<Dtype>::getDatasetCardinality()
//{

//    return s_dataset_cardinality;

//}



//template<class Dtype>
//size_t LeafNodeKdTree<Dtype>::getDatasetDimensionality()
//{

//    return s_dataset_dimensionality;

//}



//template<class Dtype>
//size_t LeafNodeKdTree<Dtype>::getDatasetSeed()
//{

//    return s_dataset_seed;

//}



//template<class Dtype>
//size_t LeafNodeKdTree<Dtype>::getBoundsSize()
//{

//    return s_bounds_size;

//}



//template<class Dtype>
//void LeafNodeKdTree<Dtype>::setSerializeVariables(size_t _s_dataset_cardinality, size_t _s_dataset_dimensionality, size_t _s_dataset_seed, size_t _s_bounds_size)
//{

//    s_dataset_cardinality = _s_dataset_cardinality;
//    s_dataset_dimensionality = _s_dataset_dimensionality;
//    s_dataset_seed = _s_dataset_seed;
//    s_bounds_size = _s_bounds_size;

//}



//template<class Dtype>
//LeafNodeKdTree<Dtype>::LeafNodeKdTree()
//{

//    this->setLeftNode(nullptr);
//    this->setRightNode(nullptr);
//    this->setParentNode(nullptr);
//    this->setBoundary(std::vector<DynamicArray<Dtype>>());
//    this->setNodeID(0);
//    setDataset(nullptr);

//}



//template<class Dtype>
//LeafNodeKdTree<Dtype>::LeafNodeKdTree(Dataset<Dtype>* _dataset)
//{

//    this->setLeftNode(nullptr);
//    this->setRightNode(nullptr);
//    this->setParentNode(nullptr);
//    this->setBoundary(std::vector<DynamicArray<Dtype>>());
//    this->setNodeID(0);
//    setDataset(_dataset);

//}



//template<class Dtype>
//LeafNodeKdTree<Dtype>::LeafNodeKdTree(Dataset<Dtype>* _dataset, NodeKdTree<Dtype>* _parent, std::vector<DynamicArray<Dtype>> _bounds)
//{

//    this->setLeftNode(nullptr);
//    this->setRightNode(nullptr);
//    this->setParentNode(_parent);
//    this->setBoundary(_bounds);
//    this->setNodeID(0);
//    setDataset(_dataset);

//}



//template<class Dtype>
//LeafNodeKdTree<Dtype>::~LeafNodeKdTree()
//{

//    if(dataset != nullptr)
//    {

//        delete dataset;

//    }

//    if(this->getBoundary().size() > 0)
//    {

//        this->getBoundary().clear();

//    }


//}



//template<class Dtype>
//Dataset<Dtype>* LeafNodeKdTree<Dtype>::getDataset()
//{

//    return dataset;

//}



//template<class Dtype>
//BasicArrayObject<Dtype>* LeafNodeKdTree<Dtype>::getInstance(size_t pos)
//{

//    return dataset->getInstance(pos);

//}



//template<class Dtype>
//void LeafNodeKdTree<Dtype>::setDataset(Dataset<Dtype>* _dataset)
//{

//    dataset = _dataset;

//}



//template<class Dtype>
//bool LeafNodeKdTree<Dtype>::isDirectoryNode()
//{

//    return false;

//}



//template<class Dtype>
//bool LeafNodeKdTree<Dtype>::isLeafNode()
//{

//    return true;

//}



//template<class Dtype>
//bool LeafNodeKdTree<Dtype>::equals(NodeKdTree<Dtype>* other)
//{

//    bool ans = true;

//    if(other->isDirectoryNode())
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
//std::string LeafNodeKdTree<Dtype>::toString()
//{

//    std::stringstream st;
//    st << "Node ID = " << this->getNodeID() << "\n";
//    st << "Bounds : " << this->showBounds() << "\n";
//    st << "Dataset : \n";

//    for(size_t x = 0; x < getDataset()->getCardinality(); x++)
//    {

//        BasicArrayObject<Dtype>* inst = getDataset()->getInstance(x);
//        st << inst->getOID() << " -> ";

//        for(size_t y = 0; y < getDataset()->getDimensionality(); y++)
//        {

//            st << inst->get(y) << " ";

//        }

//        st << "\n";

//    }

//    return st.str();

//}



//template<class Dtype>
//u_int32_t LeafNodeKdTree<Dtype>::getSerializeSize()
//{

//    setSerializeVariables(getDataset()->getCardinality(), getDataset()->getDimensionality(), getDataset()->getSeed(), this->getBoundary().size());
//    u_int32_t datasetSize = (sizeof(u_int32_t) + sizeof(double) * (u_int32_t)getDatasetDimensionality()) * (u_int32_t)getDatasetCardinality();
//    u_int32_t boundsSize = (u_int32_t)getBoundsSize() * (sizeof(Dtype) * 2);
//    return datasetSize + boundsSize + sizeof(u_int32_t) * 5;

//}



//template<class Dtype>
//u_char* LeafNodeKdTree<Dtype>::serialize()
//{

//    u_char* data = nullptr;

//    if(data == nullptr)
//    {

//        data = new u_char[getSerializeSize()];

//        u_int32_t nodeID = (u_int32_t)this->getNodeID(), datasetCardinality = (u_int32_t)getDatasetCardinality(), datasetDimensionality = (u_int32_t)getDatasetDimensionality(), datasetSeed = (u_int32_t)getDatasetSeed(), boundsSize = (u_int32_t)getBoundsSize(), total = 0;

//        memcpy(data + total, &nodeID, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);
//        memcpy(data + total, &datasetCardinality, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);
//        memcpy(data + total, &datasetDimensionality, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);
//        memcpy(data + total, &datasetSeed, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);
//        memcpy(data + total, &boundsSize, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);

//        for(size_t x = 0; x < getBoundsSize(); x++)
//        {

//            memcpy(data + total, &this->getBoundary()[x].array[0], sizeof(Dtype));
//            total += sizeof(Dtype);
//            memcpy(data + total, &this->getBoundary()[x].array[1], sizeof(Dtype));
//            total += sizeof(Dtype);

//        }

//        for(size_t x = 0; x < getDatasetCardinality(); x++)
//        {


//            BasicArrayObject<Dtype>* inst = getDataset()->getInstance(x);
//            u_int32_t oid = inst->getOID();

//            memcpy(data + total, &oid, sizeof(u_int32_t));
//            total += sizeof(u_int32_t);

//            for(size_t y = 0; y < getDatasetDimensionality(); y++)
//            {

//                memcpy(data + total, &inst->getData()[y], sizeof(double));
//                total += sizeof(double);

//            }

//        }

//    }

//    return data;

//}



//template<class Dtype>
//void LeafNodeKdTree<Dtype>::unserialize(u_char* data)
//{

//    u_int32_t nodeID = 0, datasetCardinality = 0, datasetDimensionality = 0, datasetSeed = 0, boundsSize = 0, total = 0;
//    double dl;
//    Dtype arr[2];

//    memcpy(&nodeID, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&datasetCardinality, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&datasetDimensionality, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&datasetSeed, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&boundsSize, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);

//    this->setNodeID((size_t)nodeID);
//    setSerializeVariables((size_t)datasetCardinality, (size_t)datasetDimensionality, (size_t)datasetSeed, (size_t)boundsSize);
//    this->setBoundary(std::vector<DynamicArray<Dtype>>(getBoundsSize()));

//    BasicArrayObject<Dtype>** dataAux = new BasicArrayObject<Dtype>*[getDatasetCardinality()];

//    for(size_t x = 0; x < getDatasetCardinality(); x++)
//    {

//        dataAux[x] = new BasicArrayObject<Dtype>(x, getDatasetDimensionality());

//    }

//    for(size_t x = 0; x < getBoundsSize(); x++)
//    {

//        memcpy(&arr[0], data + total, sizeof(Dtype));
//        total += sizeof(Dtype);
//        memcpy(&arr[1], data + total, sizeof(Dtype));
//        total += sizeof(Dtype);

//        this->setBounds(DynamicArray<Dtype>(arr), x);

//    }

//    for(size_t x = 0; x < getDatasetCardinality(); x++)
//    {

//        u_int32_t oid = 0;
//        memcpy(&oid, data + total, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);

//        dataAux[x]->setOID(oid);

//        for(size_t y = 0; y < getDatasetDimensionality(); y++)
//        {

//            memcpy(&dl, data + total, sizeof(double));
//            total += sizeof(double);

//            dataAux[x]->set(y, dl);

//        }

//    }

//    setDataset(new Dataset(dataAux, getDatasetCardinality(), getDatasetDimensionality()));
//    dataset->setSeed(getDatasetSeed());


//}

class LeafNodeKdTree : public NodeKdTree
{

    private:
        Dataset<double>* dataset;


    public:
        LeafNodeKdTree()
        {

            this->setLeftNode(nullptr);
            this->setRightNode(nullptr);
            this->setParentNode(nullptr);
            this->setBoundary(std::vector<Bound>());
            this->setNodeID(0);
            setDataset(nullptr);


        }

        LeafNodeKdTree(Dataset<double>* dataset_)
        {

            this->setLeftNode(nullptr);
            this->setRightNode(nullptr);
            this->setParentNode(nullptr);
            this->setBoundary(std::vector<Bound>());
            this->setNodeID(0);
            setDataset(dataset_);

        }

        LeafNodeKdTree(Dataset<double>* dataset_, NodeKdTree* parent_, std::vector<Bound> bounds_)
        {

            this->setLeftNode(nullptr);
            this->setRightNode(nullptr);
            this->setParentNode(parent_);
            this->setBoundary(bounds_);
            this->setNodeID(0);
            setDataset(dataset_);

        }

        ~LeafNodeKdTree()
        {

            if(dataset != nullptr)
            {

                delete dataset;

            }


            if(this->getBoundary().size() > 0)
            {

                this->getBoundary().clear();

            }

        }

        Dataset<double>* getDataset()
        {

            return dataset;

        }

        BasicArrayObject<double>* getInstance(size_t pos)
        {

            return dataset->instance(pos);

        }

        void setDataset(Dataset<double>* dataset_)
        {

            dataset = dataset_;

        }

        bool isDirectoryNode()
        {

            return false;

        }

        bool isLeafNode()
        {

            return true;

        }

        bool equals(NodeKdTree* other)
        {

            if(other->isDirectoryNode())
            {

                return false;

            }

            return other->equalsBounds(this->getBoundary());


        }

        std::string toString()
        {

            std::stringstream st;
            st << "Node ID = " << this->getNodeID() << "\n";
            st << "Bounds = " << this->showBounds() << "\n";
            st << "Dataset = \n";

            for(size_t i = 0; i < dataset->getCardinality(); i++)
            {

                st << dataset->instance(i)->toStringWithOID(" ") << "\n";

            }

            return st.str();

        }

        size_t getSerializedSize()
        {

            return sizeof(size_t) + (sizeof(size_t) + dataset->getSerializedSize()) + (sizeof(size_t) + this->getBoundary().size() * (sizeof(double) * 2));

        }

        unsigned char* serialize()
        {

            size_t total = 0, boundsSize = this->getBoundary().size(), node_id = this->getNodeID(), datasetSize = dataset->getSerializedSize();
            unsigned char* data = new unsigned char[getSerializedSize()];

            memcpy(data + total, &node_id, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(data + total, &boundsSize, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(data + total, &datasetSize, sizeof(size_t));
            total += sizeof(size_t);

            for(size_t i = 0; i < boundsSize; i++)
            {

                memcpy(data + total, &this->getBoundary()[i].array[0], sizeof(double));
                total += sizeof(double);
                memcpy(data + total, &this->getBoundary()[i].array[1], sizeof(double));
                total += sizeof(double);

            }

            memcpy(data + total, dataset->serialize(), datasetSize);
            total += datasetSize;

            return data;

        }


        void unserialize(unsigned char* data)
        {

            size_t total = 0, boundsSize = 0, node_id = 0, datasetSize = 0;
            double arr[2];

            memcpy(&node_id, data + total, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(&boundsSize, data + total, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(&datasetSize, data + total, sizeof(size_t));
            total += sizeof(size_t);

            this->setNodeID(node_id);
            this->setBoundary(std::vector<Bound>(boundsSize));

            for(size_t i = 0; i < boundsSize; i++)
            {

                memcpy(&arr[0], data + total, sizeof(double));
                total += sizeof(double);
                memcpy(&arr[1], data + total, sizeof(double));
                total += sizeof(double);

                this->setBound(Bound(arr), i);

            }

            unsigned char* datasetChar = new unsigned char[datasetSize];
            memcpy(datasetChar, data + total, datasetSize);

            this->dataset = new Dataset<double>();
            this->dataset->unserialize(datasetChar);

            delete [] datasetChar;

        }


};

#endif // LEAFNODEKDTREE_H
