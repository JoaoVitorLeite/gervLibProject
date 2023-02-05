#ifndef KDTREE_H
#define KDTREE_H

//#include <Node.h>
//#include <DirectoryNodeKdTree.h>
//#include <LeafNodeKdTree.h>
//#include <Dataset.h>
////#include <Index.h>
//#include <Hermes.h>
////#include <DistanceFunction.h>
////#include <DistanceBrowsing.h>
//#include <Pivots.h>
//#include <vector>
//#include <stack>
//#include <queue>
//#include <Util.h>
//#include <sys/stat.h>



//template<class Dtype>
//struct Partition
//{

//    Node<Dtype>* node;
//    double min, max;

//    Partition(Node<Dtype>* _node, double _min, double _max)
//    {

//        node = _node;
//        min = _min;
//        max = _max;

//    }

//    Partition()
//    {



//    }

//};

//template<class Dtype>
//struct ComparePartition
//{

//    bool operator()(Partition<Dtype> const& a, Partition<Dtype> const& b)
//    {

//        bool ans;

//        if(a.min != b.min)
//        {

//            ans = a.min > b.min;

//        }
//        else
//        {

//            ans = a.max > b.max;

//        }

//        return ans;

//    }

//};

//struct TupleResult
//{

//    size_t index;
//    double value1, value2;
//    size_t orderBy;

//    TupleResult(size_t index_, double value1_, double value2_, size_t orderBy_ = 2)
//    {

//        index = index_;
//        value1 = value1_;
//        value2 = value2_;
//        orderBy = orderBy_;

//    }

//};

//struct CompareTupleResult
//{

//    bool operator()(TupleResult const& a, TupleResult const& b)
//    {

//        bool ans;

//        if((a.orderBy == 1) && (b.orderBy == 1))
//        {

//            ans = a.value1 > b.value1;

//        }
//        else if((a.orderBy == 2) && (b.orderBy == 2))
//        {

//            ans = a.value2 > b.value2;

//        }
//        else
//        {

//            throw std::invalid_argument("Can't sort");

//        }

//        return ans;

//    }

//};

//struct ComparePairResult
//{

//    bool operator()(PairResult const& a, PairResult const& b)
//    {

//        return a.distance > b.distance;

//    }

//};

//struct ComparePairResult2
//{

//    bool operator()(PairResult const& a, PairResult const& b)
//    {

//        return a.distance < b.distance;

//    }

//};

//template<class Type, class Comp>
//std::vector<Type> dequeueInOrder(std::priority_queue<Type, std::vector<Type>, Comp> pq)
//{

//    std::vector<Type> ans;
//    std::priority_queue<Type, std::vector<Type>, Comp> pqClone = pq;

//    while(!pqClone.empty())
//    {

//        ans.push_back(pqClone.top());
//        pqClone.pop();

//    }

//    return ans;

//}




//template<class Dtype>
//class KdTree
//{

//    //PRIVATE
//    private:
//        Dataset<Dtype>* laesaMatrix;
//        DistanceFunction<BasicArrayObject<Dtype>>* df;
//        Node<Dtype>* root;
//        Pivot<Dtype>* pivot;
//        //DistanceBrowsing* distanceBrowsing;
//        size_t numberOfNodes;
//        size_t nodesSize, datasetSize, pivotSize;

//    //Private methods
//        void setSerializeVariables();
//        void setSerializeVariables(size_t number_of_nodes, size_t nodes_size, size_t dataset_size, size_t pivot_size);
//        void setNumberOfNodes();
//        void setRoot(Node<Dtype>* _root);
//        void buildTree(size_t numPerLeaf);
//        void buildTree(std::vector<Node<Dtype>*> vecNodes);
//        std::pair<std::vector<DynamicArray<Dtype>>, std::vector<DynamicArray<Dtype>>> splitBounds(std::vector<DynamicArray<Dtype>> _bounds, Dtype _split, size_t _pos);
//        double medianBy(Dataset<Dtype>* dataset, size_t pos);
//        void setNodeID();
//        TupleResult calculateV1V2(Instance* sq_, BasicArrayObject<Dtype>* inst);
//        void insertPriorityQueueResult(TupleResult inst, std::priority_queue<TupleResult, std::vector<TupleResult>, CompareTupleResult> & pqResult, size_t k, double & radius);
//        Instance* pivotCoordinates(Instance* inst);


//    //PUBLIC
//    public:
//        //Constructors and destructors
//        KdTree();
//        KdTree(Dataset* dataset, DistanceFunction<Instance>* distanceFunction, Pivot* pivot, size_t numPerLeaf);
//        KdTree(std::string serializedPath, DistanceFunction<Instance>* distanceFunction);
//        ~KdTree();

//        //Public methods
//        std::string toString();
//        void knn(Dataset* train, Instance* query, size_t k, std::vector<PairResult> & ans);
//        u_int8_t* serialize();
//        void unserialize(u_int8_t* data);
//        u_int32_t getSerializedSize();
//        void saveToFile();
//        void loadFromFile(std::string fileName);
//        void resetDistanceCount();

//        //Getters
//        Node<Dtype>* getRoot();
//        size_t getDiskAccess();
//        size_t getDistanceCount();
//        Dataset* getLaesaMatrix();
//        DistanceFunction<Instance>* getDistanceFunction();
//        Pivot* getPivot();


//        //Setters
//        void setDistanceFunction(DistanceFunction<Instance>* distanceFunction);



//};



//template<class Dtype>
//void KdTree<Dtype>::setSerializeVariables()
//{

//    datasetSize = laesaMatrix->getSerializedSize();
//    pivotSize = pivot->getSerializedSize();

//}



//template<class Dtype>
//void KdTree<Dtype>::setSerializeVariables(size_t number_of_nodes, size_t nodes_size, size_t dataset_size, size_t pivot_size)
//{

//    numberOfNodes = number_of_nodes;
//    nodesSize = nodes_size;
//    datasetSize = dataset_size;
//    pivotSize = pivot_size;

//}



//template<class Dtype>
//void KdTree<Dtype>::setNumberOfNodes()
//{

//    if(getRoot() == nullptr)
//    {

//        numberOfNodes = 0;
//        return;

//    }

//    std::queue<Node<Dtype>*> queueNode;
//    queueNode.push(getRoot());
//    Node<Dtype>* curr = nullptr;
//    numberOfNodes = 0;
//    nodesSize = 0;

//    while(!queueNode.empty())
//    {

//        curr = queueNode.front();
//        queueNode.pop();
//        numberOfNodes++;
//        nodesSize += curr->getSerializeSize();

//        if(curr->getLeft() != nullptr)
//        {

//            queueNode.push(curr->getLeft());

//        }

//        if(curr->getRight() != nullptr)
//        {

//            queueNode.push(curr->getRight());

//        }

//    }

//}


//template<class Dtype>
//void KdTree<Dtype>::setRoot(Node<Dtype>* _root)
//{

//    root = _root;

//}



//template<class Dtype>
//void KdTree<Dtype>::buildTree(size_t numPerLeaf)
//{

//    size_t dPartition = 0;
//    std::stack<std::tuple<Node<Dtype>*, Dataset*, size_t>> stackNode;
//    std::vector<DynamicArray<Dtype>> boundaries(laesaMatrix->getDimensionality());

//    for(size_t x = 0; x < laesaMatrix->getDimensionality(); x++)
//    {

//        Dtype arr[2] = {0.0, std::numeric_limits<Dtype>::max()};
//        boundaries[x] = DynamicArray<Dtype>(arr);

//    }

//    Node<Dtype>* node = new DirectoryNode<Dtype>(nullptr, boundaries);
//    setRoot(node);
//    stackNode.push(std::make_tuple(node, laesaMatrix, dPartition));

//    while(!stackNode.empty())
//    {

//        numberOfNodes++;
//        std::tuple<Node<Dtype>*, Dataset*, size_t> hd = stackNode.top();
//        stackNode.pop();
//        Node<Dtype>* node = std::get<0>(hd);
//        Dataset* dataAux = std::get<1>(hd);
//        dPartition = std::get<2>(hd) % laesaMatrix->getDimensionality();
//        double median = medianBy(dataAux, dPartition);
//        std::vector<Instance*> l, r;

//        for(size_t x = 0; x < dataAux->getCardinality(); x++)
//        {

//            Instance* inst = dataAux->getInstance(x);

//            if(inst->get(dPartition) < median)
//            {

//                l.push_back(inst);

//            }
//            else
//            {

//                r.push_back(inst);

//            }

//        }

//        if((l.size() == 0) || (r.size() == 0))
//        {


//            if(node->getParent()->getLeft()->equals(node))
//            {

//                node->getParent()->setLeftNode(new LeafNode<Dtype>(dataAux, node->getParent(), node->getBoundary()));

//            }
//            else
//            {

//                node->getParent()->setRightNode(new LeafNode<Dtype>(dataAux, node->getParent(), node->getBoundary()));

//            }

//        }
//        else
//        {


//            std::pair<std::vector<DynamicArray<Dtype>>, std::vector<DynamicArray<Dtype>>> splitBoundaries = splitBounds(node->getBoundary(), median, dPartition);

//            Instance** auxLeft = new Instance*[l.size()];
//            Instance** auxRight = new Instance*[r.size()];

//            for(size_t x = 0; x < l.size(); x++)
//            {

//                auxLeft[x] = l[x];

//            }

//            for(size_t x = 0; x < r.size(); x++)
//            {

//                auxRight[x] = r[x];

//            }


//            if(l.size() <= numPerLeaf)
//            {

//                node->setLeftNode(new LeafNode<Dtype>(new Dataset(auxLeft, l.size(), l[0]->size()), node, splitBoundaries.first));

//            }
//            else
//            {

//                Node<Dtype>* nextL = new DirectoryNode<Dtype>(node, splitBoundaries.first);
//                node->setLeftNode(nextL);
//                stackNode.push(std::make_tuple(nextL, new Dataset(auxLeft, l.size(), l[0]->size()), dPartition + 1));

//            }

//            if(r.size() <= numPerLeaf)
//            {

//                node->setRightNode(new LeafNode<Dtype>(new Dataset(auxRight, r.size(), r[0]->size()), node, splitBoundaries.second));

//            }
//            else
//            {

//                Node<Dtype>* nextR = new DirectoryNode<Dtype>(node, splitBoundaries.second);
//                node->setRightNode(nextR);
//                stackNode.push(std::make_tuple(nextR, new Dataset(auxRight, r.size(), r[0]->size()), dPartition + 1));

//            }

//        }


//    }

//    setNumberOfNodes();
//    setNodeID();

//}



//template<class Dtype>
//void KdTree<Dtype>::buildTree(std::vector<Node<Dtype>*> vecNodes)
//{

//    struct CompareNode{

//        bool operator()(Node<Dtype>* a, Node<Dtype>* b)
//        {

//            return a->getNodeID() < b->getNodeID();

//        }

//    };

//    std::sort(vecNodes.begin(), vecNodes.end(), CompareNode());

//    if(vecNodes[0]->getNodeID() == 1)
//    {

//        setRoot(vecNodes[0]);

//    }
//    else
//    {

//        throw std::invalid_argument("Root node not found");

//    }

//    std::queue<Node<Dtype>*> queueNode;
//    queueNode.push(getRoot());
//    size_t pos = 0, index = 1, flag = numberOfNodes ;
//    Node<Dtype>* curr = nullptr;

//    while(flag)
//    {

//        curr = queueNode.front();
//        queueNode.pop();

//        if(curr == nullptr)
//        {

//            index += 2;
//            queueNode.push(nullptr);
//            queueNode.push(nullptr);

//        }
//        else
//        {

//            index++;
//            pos = -1;
//            for(size_t x = 0; x < vecNodes.size(); x++)
//            {

//                if(vecNodes[x]->getNodeID() == index)
//                {

//                    pos = x;
//                    break;

//                }

//            }
//            if(pos != (size_t)-1)
//            {

//                queueNode.push(vecNodes[pos]);
//                curr->setLeftNode(vecNodes[pos]);
//                vecNodes[pos]->setParentNode(curr);

//            }
//            else
//            {

//                queueNode.push(nullptr);

//            }

//            index++;
//            pos = -1;
//            for(size_t x = 0; x < vecNodes.size(); x++)
//            {

//                if(vecNodes[x]->getNodeID() == index)
//                {

//                    pos = x;
//                    break;

//                }

//            }
//            if(pos != (size_t)-1)
//            {

//                queueNode.push(vecNodes[pos]);
//                curr->setRightNode(vecNodes[pos]);
//                vecNodes[pos]->setParentNode(curr);

//            }
//            else
//            {

//                queueNode.push(nullptr);

//            }

//            flag--;

//        }

//    }

//}



//template<class Dtype>
//std::pair<std::vector<DynamicArray<Dtype>>, std::vector<DynamicArray<Dtype>>> KdTree<Dtype>::splitBounds(std::vector<DynamicArray<Dtype>> _bounds, Dtype _split, size_t _pos)
//{

//    std::vector<DynamicArray<Dtype>> cp1 = _bounds;
//    std::vector<DynamicArray<Dtype>> cp2 = _bounds;
//    Dtype arr_cp1[2] = {_bounds[_pos].array[0], -_split};
//    Dtype arr_cp2[2] = {_split, _bounds[_pos].array[1]};
//    cp1[_pos] = DynamicArray<Dtype>(arr_cp1);
//    cp2[_pos] = DynamicArray<Dtype>(arr_cp2);
//    return std::make_pair(cp1, cp2);

//}



//template<class Dtype>
//double KdTree<Dtype>::medianBy(Dataset* dataset, size_t pos)
//{

//    std::vector<double> vec(dataset->getCardinality());

//    for(size_t x = 0; x < dataset->getCardinality(); x++)
//    {

//        vec[x] = dataset->getInstance(x)->get(pos);

//    }

//    std::sort(vec.begin(), vec.end());

//    if(vec.size() % 2 == 0)
//    {

//        return (vec[vec.size()/2 - 1] + vec[vec.size()/2])/2.0;

//    }
//    else
//    {

//        return vec[vec.size()/2];

//    }


//}



//template<class Dtype>
//void KdTree<Dtype>::setNodeID()
//{

//    std::queue<Node<Dtype>*> queueNode;
//    getRoot()->setNodeID(1);
//    queueNode.push(getRoot());
//    size_t index = 1, flag = numberOfNodes;
//    Node<Dtype>* curr = nullptr;

//    while(flag)
//    {

//        curr = queueNode.front();
//        queueNode.pop();

//        if(curr == nullptr)
//        {

//            index += 2;
//            queueNode.push(nullptr);
//            queueNode.push(nullptr);

//        }
//        else
//        {

//            index++;

//            if(curr->getLeft() != nullptr)
//            {

//                queueNode.push(curr->getLeft());
//                curr->getLeft()->setNodeID(index);

//            }
//            else
//            {

//                queueNode.push(nullptr);

//            }

//            index++;

//            if(curr->getRight() != nullptr)
//            {

//                queueNode.push(curr->getRight());
//                curr->getRight()->setNodeID(index);

//            }
//            else
//            {

//                queueNode.push(nullptr);

//            }

//            flag--;

//        }

//    }


//}



//template<class Dtype>
//TupleResult KdTree<Dtype>::calculateV1V2(Instance* sq_, Instance* inst)
//{

//    double ansV1 = std::numeric_limits<double>::max(), ansV2 = std::numeric_limits<double>::max(), inf = 0, sup = 0;

//    for(size_t x = 0; x < inst->size(); x++)
//    {

//        inf = std::abs(sq_->get(x) - laesaMatrix->getInstance(inst->getOID())->get(x));
//        sup = sq_->get(x) + laesaMatrix->getInstance(inst->getOID())->get(x);

//        if(inf < ansV1 && sup < ansV2)
//        {

//            ansV1 = inf;
//            ansV2 = sup;

//        }
//        else if(inf < ansV1)
//        {

//            ansV1 = inf;

//        }
//        else if(sup < ansV2)
//        {

//            ansV2 = sup;

//        }

//    }

//    return TupleResult(inst->getOID(), ansV1, ansV2);

//}



//template<class Dtype>
//void KdTree<Dtype>::insertPriorityQueueResult(TupleResult inst, std::priority_queue<TupleResult, std::vector<TupleResult>, CompareTupleResult> & pqResult, size_t k, double & radius)
//{

//    if(pqResult.size() == k)
//    {

//        if(inst.value2 <= radius)
//        {

//            pqResult.push(inst);
//            std::vector<TupleResult> vec = dequeueInOrder<TupleResult, CompareTupleResult>(pqResult);
//            radius = vec[k-1].value2;
//            vec.clear();

//        }
//        else
//        {

//            if(inst.value1 <= radius)
//            {

//                pqResult.push(inst);

//            }

//        }

//    }
//    else
//    {

//        pqResult.push(inst);

//        if(pqResult.size() == k)
//        {

//            std::vector<TupleResult> vec = dequeueInOrder<TupleResult, CompareTupleResult>(pqResult);
//            radius = vec[k-1].value2;
//            vec.clear();

//        }

//    }

//}



//template<class Dtype>
//Instance* KdTree<Dtype>::pivotCoordinates(Instance* inst)
//{

//    Instance* ans = new Instance(0, pivot->getNumberOfPivots());

//    for(size_t x = 0; x < pivot->getNumberOfPivots(); x++)
//    {

//        ans->set(x, df->getDistance(*inst, *pivot->getPivot(x)));

//    }

//    return ans;

//}



//template<class Dtype>
//KdTree<Dtype>::KdTree()
//{

//   laesaMatrix = nullptr;
//   df = nullptr;
//   root = nullptr;
//   pivot = nullptr;
//   distanceBrowsing = nullptr;
//   numberOfNodes = 0;

//}



//template<class Dtype>
//KdTree<Dtype>::KdTree(Dataset* dataset, DistanceFunction<Instance>* distanceFunction, Pivot* pivot, size_t numPerLeaf)
//{

//    df = distanceFunction;
//    this->pivot = pivot;
//    distanceBrowsing = new DistanceBrowsing(df);
//    numberOfNodes = 0;
//    root = nullptr;

//    Instance** laesa = new Instance*[dataset->getCardinality()];

//    for(size_t x = 0; x < dataset->getCardinality(); x++)
//    {

//        laesa[x] = new Instance(x, pivot->getNumberOfPivots());

//        for(size_t y = 0; y < pivot->getNumberOfPivots(); y++)
//        {


//            laesa[x]->set(y, distanceFunction->getDistance(*dataset->instance(x), *pivot->getPivot(y)));

//        }

//    }

//    laesaMatrix = new Dataset(laesa, dataset->getCardinality(), pivot->getNumberOfPivots());

//    buildTree(numPerLeaf);

//}



//template<class Dtype>
//KdTree<Dtype>::KdTree(std::string serializedPath, DistanceFunction<Instance>* distanceFunction)
//{

//    loadFromFile(serializedPath);
//    setDistanceFunction(distanceFunction);
//    distanceBrowsing = new DistanceBrowsing(distanceFunction);

//}



//template<class Dtype>
//KdTree<Dtype>::~KdTree()
//{

//    if(root != nullptr)
//    {

//        std::queue<Node<Dtype>*> queueNode;
//        queueNode.push(root);

//        while(!queueNode.empty())
//        {

//            Node<Dtype>* curr = queueNode.front();
//            queueNode.pop();

//            if(curr->getLeft() != nullptr)
//            {

//                queueNode.push(curr->getLeft());

//            }

//            if(curr->getRight() != nullptr)
//            {

//                queueNode.push(curr->getRight());

//            }

//            delete curr;

//        }

//    }

//}



//template<class Dtype>
//std::string KdTree<Dtype>::toString()
//{

//    std::queue<std::pair<Node<Dtype>*, size_t>> queue;
//    queue.push(std::make_pair(getRoot(), 0));
//    std::stringstream ss;

//    while(!queue.empty())
//    {

//        std::pair<Node<Dtype>*, size_t> element = queue.front();
//        queue.pop();

//        if(element.first->getLeft() != nullptr)
//        {

//            queue.push(std::make_pair(element.first->getLeft(), element.second + 1));

//        }

//        if(element.first->getRight() != nullptr)
//        {

//            queue.push(std::make_pair(element.first->getRight(), element.second + 1));

//        }

//        ss << "Level = " << element.second << "\n";
//        ss << "Node ID = " << element.first->getNodeID() << "\n";

//        if(element.first->isLeafNode())
//        {

//            ss << "Bounds :" << element.first->showBounds() << "\n";
//            ss << "Dataset :\n";
//            LeafNode<Dtype>* leaf = (LeafNode<Dtype>*)element.first;
//            Dataset* data = leaf->getDataset();
//            for(size_t x = 0; x < data->getCardinality(); x++)
//            {

//                ss << data->getInstance(x)->getOID() << " -> ";

//                for(size_t z = 0; z < data->getDimensionality(); z++)
//                {

//                    ss << data->getInstance(x)->get(z) << " ";

//                }

//                ss << "\n";

//            }

//        }
//        else
//        {

//            ss << "Bounds :" << element.first->showBounds() << "\n";

//        }

//        ss << "\n\n\n";

//    }

//    return ss.str();

//}



//template<class Dtype>
//void KdTree<Dtype>::knn(Dataset* train, Instance* query, size_t k, std::vector<PairResult> & ans)
//{

//    df->resetStatistics();
//    std::priority_queue<Partition<Dtype>, std::vector<Partition<Dtype>>, ComparePartition<Dtype>> pqPartition;
//    std::priority_queue<TupleResult, std::vector<TupleResult>, CompareTupleResult> pqResult;
//    Instance* sq_ = pivotCoordinates(query);
//    double radius = std::numeric_limits<double>::max();
//    pqPartition.push(Partition<Dtype>(getRoot(), 0.0, 0.0));
//    LeafNode<Dtype>* lf = nullptr;

//    while(!pqPartition.empty())
//    {

//        Node<Dtype>* node = pqPartition.top().node;
//        pqPartition.pop();

//        if(node->isDirectoryNode())
//        {

//            double minL = distanceBrowsing->minDist(sq_, node->getLeft()->getBoundary());
//            double maxL = distanceBrowsing->maxDist(sq_, node->getLeft()->getBoundary());
//            double minR = distanceBrowsing->minDist(sq_, node->getRight()->getBoundary());
//            double maxR = distanceBrowsing->maxDist(sq_, node->getRight()->getBoundary());

//            pqPartition.push(Partition<Dtype>(node->getLeft(), minL, maxL));
//            pqPartition.push(Partition<Dtype>(node->getRight(), minR, maxR));

//        }
//        else
//        {

//            lf = (LeafNode<Dtype>*)node;

//            for(size_t x = 0; x < lf->getDataset()->getCardinality(); x++)
//            {

//                insertPriorityQueueResult(calculateV1V2(sq_, lf->getDataset()->getInstance(x)), pqResult, k, radius);

//            }

//            lf = nullptr;

//        }

//    }

//    std::vector<TupleResult> orderByV1 = dequeueInOrder<TupleResult, CompareTupleResult>(pqResult);

//    std::sort(orderByV1.begin(), orderByV1.end(), [](TupleResult const& a, TupleResult const& b){

//        return a.value1 > b.value1;

//    });

//    std::priority_queue<PairResult, std::vector<PairResult>, ComparePairResult> pqAns;

//    size_t i = 0;

//    while(i != orderByV1.size())
//    {

//        TupleResult pop = orderByV1[i];

//        if(pop.value1 <= radius)
//        {

//            double dist = df->getDistance(*query, *train->getInstance(pop.index));

//            if(dist <= radius)
//            {

//                if(pqAns.size() == k)
//                {

//                    pqAns.push(PairResult(train->getInstance(pop.index)->getOID(), dist));
//                    std::vector<PairResult> aux = dequeueInOrder<PairResult, ComparePairResult>(pqAns);
//                    radius = aux[k-1].distance;

//                    pqAns = std::priority_queue<PairResult, std::vector<PairResult>, ComparePairResult>();

//                    for(size_t x = 0; x < k; x++)
//                    {

//                        pqAns.push(aux[x]);

//                    }

//                }
//                else
//                {

//                    pqAns.push(PairResult(train->getInstance(pop.index)->getOID(), dist));

//                    if(pqAns.size() == k)
//                    {

//                        std::vector<PairResult> aux = dequeueInOrder<PairResult, ComparePairResult>(pqAns);
//                        radius = aux[k-1].distance;

//                        pqAns = std::priority_queue<PairResult, std::vector<PairResult>, ComparePairResult>();

//                        for(size_t x = 0; x < k; x++)
//                        {

//                            pqAns.push(aux[x]);

//                        }

//                    }

//                }

//            }

//        }

//        i++;

//    }

//    ans = dequeueInOrder<PairResult, ComparePairResult>(pqAns);
//    ans.resize(k);



//    orderByV1.clear();

//    while(!pqResult.empty())
//    {

//        pqResult.pop();

//    }

//}



//template<class Dtype>
//u_int8_t* KdTree<Dtype>::serialize()
//{

//    setSerializeVariables();
//    u_int8_t* data = new u_int8_t[getSerializedSize()];
//    u_int32_t nodesSize_ = (u_int32_t)nodesSize, datasetSize_ = (u_int32_t)datasetSize, pivotSize_ = (u_int32_t)pivotSize, numberOfNodes_ = (u_int32_t)numberOfNodes, total = 0;
//    PIVOT_TYPE pvt_type = pivot->getPivotType();

//    memcpy(data + total, &numberOfNodes_, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(data + total, &nodesSize_, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(data + total, &datasetSize_, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(data + total, &pivotSize_, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(data + total, &pvt_type, sizeof(PIVOT_TYPE));
//    total += sizeof(PIVOT_TYPE);

//    std::queue<Node<Dtype>*> queueNode;
//    queueNode.push(getRoot());
//    Node<Dtype>* curr = nullptr;
//    u_int32_t nodeType, nodeSize;
//    u_int8_t* nodeData = nullptr;

//    while(!queueNode.empty())
//    {

//        curr = queueNode.front();
//        queueNode.pop();

//        if(curr->getLeft() != nullptr)
//        {

//            queueNode.push(curr->getLeft());

//        }

//        if(curr->getRight() != nullptr)
//        {

//            queueNode.push(curr->getRight());

//        }

//        nodeType = curr->isDirectoryNode() ? 1 : 0;
//        nodeSize = curr->getSerializeSize();

//        memcpy(data + total, &nodeType, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);
//        memcpy(data + total, &nodeSize, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);

//        nodeData = curr->serialize();
//        memcpy(data + total, nodeData, nodeSize);
//        total += nodeSize;

//    }

//    memcpy(data + total, (u_int8_t*)laesaMatrix->serialize(), datasetSize);
//    total += datasetSize;
//    memcpy(data + total, (u_int8_t*)pivot->serialize(), pivotSize);
//    total += pivotSize;


//    return data;

//}



//template<class Dtype>
//void KdTree<Dtype>::unserialize(u_int8_t* data)
//{

//    u_int32_t nodesSize_ = 0, datasetSize_ = 0, pivotSize_ = 0, numberOfNodes_ = 0, total = 0;
//    PIVOT_TYPE pvt_type;

//    memcpy(&numberOfNodes_, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&nodesSize_, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&datasetSize_, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&pivotSize_, data + total, sizeof(u_int32_t));
//    total += sizeof(u_int32_t);
//    memcpy(&pvt_type, data + total, sizeof(PIVOT_TYPE));
//    total += sizeof(PIVOT_TYPE);

//    setSerializeVariables(numberOfNodes_, nodesSize_, datasetSize_, pivotSize_);

//    u_int32_t nodeType, nodeSize;
//    std::vector<Node<Dtype>*> vecNodes(numberOfNodes);
//    Node<Dtype>* curr = nullptr;
//    u_int8_t* nodeData = nullptr;

//    for(size_t x = 0; x < numberOfNodes; x++)
//    {

//        memcpy(&nodeType, data + total, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);
//        memcpy(&nodeSize, data + total, sizeof(u_int32_t));
//        total += sizeof(u_int32_t);

//        if(nodeType == 1)
//        {

//            curr = new DirectoryNode<Dtype>();

//        }
//        else
//        {

//            curr = new LeafNode<Dtype>();

//        }

//        nodeData = new u_int8_t[nodeSize];
//        memcpy(nodeData, data + total, nodeSize);
//        total += nodeSize;

//        curr->unserialize(nodeData);

//        vecNodes[x] = curr;

//    }

//    u_int8_t* dataset_char = new u_int8_t[datasetSize],
//            *pivot_char = new u_int8_t[pivotSize];

//    memcpy(dataset_char, data + total, datasetSize);
//    total += datasetSize;
//    memcpy(pivot_char, data + total, pivotSize);
//    total += pivotSize;

//    laesaMatrix = new Dataset();
//    laesaMatrix->unserialize((char*)dataset_char);

//    std::string names[] = {"RANDOM", "GNAT", "CONVEX", "KMEDOIDS", "MAXSEPARETED", "MAXVARIANCE", "SELECTION", "PCA", "SSS"};

//    if(names[pvt_type] == "RANDOM")
//        pivot = new RandomPivots();
//    else if(names[pvt_type] == "GNAT")
//        pivot = new GnatPivots();
//    else if(names[pvt_type] == "CONVEX")
//        pivot = new ConvexPivots();
//    else if(names[pvt_type] == "KMEDOIDS")
//        pivot = new KmedoidsPivots();
//    else if(names[pvt_type] == "MAXSEPARETED")
//        pivot = new MaxSeparetedPivots();
//    else if(names[pvt_type] == "MAXVARIANCE")
//        pivot = new MaxVariancePivots();
//    else if(names[pvt_type] == "SELECTION")
//        pivot = new SelectionPivots();
//    else if(names[pvt_type] == "PCA")
//        pivot = new PCAPivots();
//    else if(names[pvt_type] == "SSS")
//        pivot = new SSSPivots();
//    else
//        throw std::invalid_argument("Can't read pivot type !_!");

//    pivot->unserialize((char*)pivot_char);

//    buildTree(vecNodes);

//}



//template<class Dtype>
//u_int32_t KdTree<Dtype>::getSerializedSize()
//{

//    setSerializeVariables();
//    return 4 * sizeof(u_int32_t) + nodesSize + datasetSize + pivotSize + numberOfNodes * 2 * sizeof(u_int32_t) + sizeof(PIVOT_TYPE);

//}



//template<class Dtype>
//void KdTree<Dtype>::saveToFile()
//{

//    setSerializeVariables();

//    int check = -1;

//    std::string fileName;

//    while(check == -1)
//    {

//        std::string directory = "../kdtree/kdtree_files/",
//                name = selectFolderName(directory, "kdtree"),
//                folderName = directory + name;

//        fileName = folderName + "/" + name + ".dat";

//        check = mkdir(folderName.c_str(), 0777);

//    }

//    std::ofstream file(fileName, std::ios::out | std::ios::binary);

//    u_int8_t* seri = new u_int8_t[sizeof(u_int32_t) + getSerializedSize()];
//    u_int32_t size = getSerializedSize();
//    memcpy(seri, &size, sizeof(u_int32_t));
//    memcpy(seri + sizeof(u_int32_t), serialize(), size);

//    file.write((char*)seri, size + sizeof(u_int32_t));
//    file.close();


//}



//template<class Dtype>
//void KdTree<Dtype>::loadFromFile(std::string fileName)
//{

//    std::ifstream file(fileName, std::ios::in | std::ios::binary);

//    u_int8_t* seri = new u_int8_t[sizeof(u_int32_t)];
//    u_int32_t size;

//    file.read((char*)seri, sizeof(u_int32_t));

//    memcpy(&size, seri, sizeof(u_int32_t));

//    u_int8_t* data = new u_int8_t[size + sizeof(u_int32_t)];

//    file.seekg(0);
//    file.read((char *)data, size + sizeof(u_int32_t));

//    unserialize(data + sizeof(u_int32_t));

//    file.close();

//}



//template<class Dtype>
//void KdTree<Dtype>::resetDistanceCount()
//{

//    df->resetStatistics();

//}



//template<class Dtype>
//Node<Dtype>* KdTree<Dtype>::getRoot()
//{

//    return root;

//}



//template<class Dtype>
//size_t KdTree<Dtype>::getDiskAccess()
//{

//    return 0;

//}



//template<class Dtype>
//size_t KdTree<Dtype>::getDistanceCount()
//{

//    return (size_t)df->getDistanceCount();

//}



//template<class Dtype>
//Dataset* KdTree<Dtype>::getLaesaMatrix()
//{

//    return laesaMatrix;

//}



//template<class Dtype>
//DistanceFunction<Instance>* KdTree<Dtype>::getDistanceFunction()
//{

//    return df;

//}



//template<class Dtype>
//Pivot* KdTree<Dtype>::getPivot()
//{

//    return pivot;

//}



//template<class Dtype>
//void KdTree<Dtype>::setDistanceFunction(DistanceFunction<Instance>* distanceFunction)
//{

//    df = distanceFunction;

//}


#include <NodeKdTree.h>
#include <LeafNodeKdTree.h>
#include <DirectoryNodeKdTree.h>
#include <Pivots.h>
#include <stack>
#include <tuple>
#include <limits.h>
#include <cstdlib>
#include <sys/stat.h>


struct PairResult
{

    size_t index;
    double distance;

    PairResult(size_t index_, double distance_)
    {

        index = index_;
        distance = distance_;

    }

    PairResult()
    {

        index = 0;
        distance = -1;

    }

};

struct Partition
{

    NodeKdTree* node;
    double min, max;

    Partition(NodeKdTree* _node, double _min, double _max)
    {

        node = _node;
        min = _min;
        max = _max;

    }

    Partition()
    {



    }

};

struct ComparePartition
{

    bool operator()(Partition const& a, Partition const& b)
    {

        bool ans;

        if(a.min != b.min)
        {

            ans = a.min > b.min;

        }
        else
        {

            ans = a.max > b.max;

        }

        return ans;

    }

};

struct TupleResult
{

    size_t index;
    double value1, value2;
    size_t orderBy;

    TupleResult(size_t index_, double value1_, double value2_, size_t orderBy_ = 2)
    {

        index = index_;
        value1 = value1_;
        value2 = value2_;
        orderBy = orderBy_;

    }

};

struct CompareTupleResult
{

    bool operator()(TupleResult const& a, TupleResult const& b)
    {

        bool ans;

        if((a.orderBy == 1) && (b.orderBy == 1))
        {

            ans = a.value1 > b.value1;

        }
        else if((a.orderBy == 2) && (b.orderBy == 2))
        {

            ans = a.value2 > b.value2;

        }
        else
        {

            throw std::invalid_argument("Can't sort");

        }

        return ans;

    }

};

struct ComparePairResult
{

    bool operator()(PairResult const& a, PairResult const& b)
    {

        return a.distance > b.distance;

    }

};

struct ComparePairResult2
{

    bool operator()(PairResult const& a, PairResult const& b)
    {

        return a.distance < b.distance;

    }

};

template<class Type, class Comp>
std::vector<Type> dequeueInOrder(std::priority_queue<Type, std::vector<Type>, Comp> pq)
{

    std::vector<Type> ans;
    std::priority_queue<Type, std::vector<Type>, Comp> pqClone = pq;

    while(!pqClone.empty())
    {

        ans.push_back(pqClone.top());
        pqClone.pop();

    }

    return ans;

}

template <class DType>
class KdTree
{

    private:
        Dataset<double>* laesaMatrix;
        NodeKdTree* root;
        Pivot<DType>* pvt;
        DistanceFunction<BasicArrayObject<DType>>* df;
        size_t numberOfNodes, nodesSize;

    private:

        void setRoot(NodeKdTree* root_)
        {

            root = root_;

        }

        double medianBy(Dataset<double>* data, size_t dim)
        {

            std::vector<double> vec(data->getCardinality());

            for(size_t i = 0; i < data->getCardinality(); i++)
            {

                vec[i] = *data->getInstance(i)->get(dim);

            }

            std::sort(vec.begin(), vec.end());

            if(vec.size() % 2 == 0)
            {

                return (vec[vec.size()/2 - 1] + vec[vec.size()/2])/2.0;

            }
            else
            {

                return vec[vec.size()/2];

            }

        }

        std::pair<std::vector<Bound>, std::vector<Bound>> splitBounds(std::vector<Bound> bounds, double split, size_t pos)
        {

            std::vector<Bound> cp1 = bounds;
            std::vector<Bound> cp2 = bounds;
            double arr_cp1[2] = {bounds[pos].array[0], -split};
            double arr_cp2[2] = {split, bounds[pos].array[1]};
            cp1[pos] = Bound(arr_cp1);
            cp2[pos] = Bound(arr_cp2);
            return std::make_pair(cp1, cp2);

        }

        void setNodeID()
        {

            std::queue<NodeKdTree*> queueNode;
            getRoot()->setNodeID(1);
            queueNode.push(getRoot());
            size_t index = 1, flag = numberOfNodes;
            NodeKdTree* curr = nullptr;

            while(flag)
            {

                curr = queueNode.front();
                queueNode.pop();

                if(curr == nullptr)
                {

                    index += 2;
                    queueNode.push(nullptr);
                    queueNode.push(nullptr);

                }
                else
                {

                    index++;

                    if(curr->getLeft() != nullptr)
                    {

                        queueNode.push(curr->getLeft());
                        curr->getLeft()->setNodeID(index);

                    }
                    else
                    {

                        queueNode.push(nullptr);

                    }

                    index++;

                    if(curr->getRight() != nullptr)
                    {

                        queueNode.push(curr->getRight());
                        curr->getRight()->setNodeID(index);

                    }
                    else
                    {

                        queueNode.push(nullptr);

                    }

                    flag--;


                }

            }


        }

        void setNumberOfNodes()
        {

            if(getRoot() == nullptr)
            {

                numberOfNodes = 0;
                return;

            }

            std::queue<NodeKdTree*> queueNode;
            queueNode.push(getRoot());
            NodeKdTree* curr = nullptr;
            numberOfNodes = 0;
            nodesSize = 0;

            while(!queueNode.empty())
            {

                curr = queueNode.front();
                queueNode.pop();
                numberOfNodes++;
                nodesSize += curr->getSerializedSize();

                if(curr->getLeft() != nullptr)
                {

                    queueNode.push(curr->getLeft());

                }

                if(curr->getRight() != nullptr)
                {

                    queueNode.push(curr->getRight());

                }

            }


        }




        void buildTree(size_t numPerLeaf)
        {

            size_t dPartition = 0;
            std::stack<std::tuple<NodeKdTree*, Dataset<double>*, size_t>> stackNode;
            std::vector<Bound> boundaries = std::vector<Bound>(laesaMatrix->getDimensionality());

            for(size_t i = 0; i < laesaMatrix->getDimensionality(); i++)
            {

                double arr[2] = {0.0, std::numeric_limits<double>::max()};
                boundaries[i] = Bound(arr);

            }

            NodeKdTree* node = new DirectoryNodeKdTree(nullptr, boundaries);
            setRoot(node);
            stackNode.push(std::make_tuple(node, laesaMatrix, dPartition));

            while(!stackNode.empty())
            {

                numberOfNodes++;
                std::tuple<NodeKdTree*, Dataset<double>*, size_t> hd = stackNode.top();
                stackNode.pop();
                NodeKdTree* node = std::get<0>(hd);
                Dataset<double>* dataAux = std::get<1>(hd);
                dPartition = std::get<2>(hd) % laesaMatrix->getDimensionality();
                double median = medianBy(dataAux, dPartition);

                std::vector<BasicArrayObject<double>> left, right;

                for(size_t i = 0; i < dataAux->getCardinality(); i++)
                {

                    BasicArrayObject<double> inst = dataAux->getFeatureVector(i);

                    if(*inst.get(dPartition) < median)
                    {

                        left.push_back(inst);

                    }
                    else
                    {

                        right.push_back(inst);

                    }

                }


                if(left.size() == 0 || right.size() == 0)
                {

                    if(node->getParent()->getLeft()->equals(node))
                    {

                        node->getParent()->setLeftNode(new LeafNodeKdTree(dataAux, node->getParent(), node->getBoundary()));

                    }
                    else
                    {

                        node->getParent()->setRightNode(new LeafNodeKdTree(dataAux, node->getParent(), node->getBoundary()));

                    }

                }
                else
                {

                    std::pair<std::vector<Bound>, std::vector<Bound>> splitBoundaries = splitBounds(node->getBoundary(), median, dPartition);

                    if(left.size() <= numPerLeaf)
                    {

                        node->setLeftNode(new LeafNodeKdTree(new Dataset<double>(left, left.size(), left[0].size()), node, splitBoundaries.first));

                    }
                    else
                    {

                        NodeKdTree* nextL = new DirectoryNodeKdTree(node, splitBoundaries.first);
                        node->setLeftNode(nextL);
                        stackNode.push(std::make_tuple(nextL, new Dataset<double>(left, left.size(), left[0].size()), dPartition+1));

                    }

                    if(right.size() <= numPerLeaf)
                    {

                        node->setRightNode(new LeafNodeKdTree(new Dataset<double>(right, right.size(), right[0].size()), node, splitBoundaries.second));

                    }
                    else
                    {

                        NodeKdTree* nextR = new DirectoryNodeKdTree(node, splitBoundaries.second);
                        node->setRightNode(nextR);
                        stackNode.push(std::make_tuple(nextR, new Dataset<double>(right, right.size(), right[0].size()), dPartition+1));

                    }


                }


            }
            setNumberOfNodes();
            setNodeID();

        }

        void buildTree(std::vector<NodeKdTree*> vecNodes)
        {

            struct CompareNode{

                bool operator()(NodeKdTree* a, NodeKdTree* b)
                {

                    return a->getNodeID() < b->getNodeID();

                }

            };

            std::sort(vecNodes.begin(), vecNodes.end(), CompareNode());

            if(vecNodes[0]->getNodeID() == 1)
            {

                setRoot(vecNodes[0]);

            }
            else
            {

                throw std::invalid_argument("Root node not found");

            }

            std::queue<NodeKdTree*> queueNode;
            queueNode.push(getRoot());
            size_t pos = 0, index = 1, flag = numberOfNodes;
            NodeKdTree* curr = nullptr;

            while(flag)
            {

                curr = queueNode.front();
                queueNode.pop();

                if(curr == nullptr)
                {

                    index += 2;
                    queueNode.push(nullptr);
                    queueNode.push(nullptr);

                }
                else
                {

                    index++;
                    pos = -1;
                    for(size_t x = 0; x < vecNodes.size(); x++)
                    {

                        if(vecNodes[x]->getNodeID() == index)
                        {

                            pos = x;
                            break;

                        }

                    }

                    if(pos != -1)
                    {

                        queueNode.push(vecNodes[pos]);
                        curr->setLeftNode(vecNodes[pos]);
                        vecNodes[pos]->setParentNode(curr);

                    }
                    else
                    {

                        queueNode.push(nullptr);

                    }

                    index++;
                    pos = -1;

                    for(size_t x = 0; x < vecNodes.size(); x++)
                    {

                        if(vecNodes[x]->getNodeID() == index)
                        {

                            pos = x;
                            break;

                        }

                    }

                    if(pos != -1)
                    {

                        queueNode.push(vecNodes[pos]);
                        curr->setRightNode(vecNodes[pos]);
                        vecNodes[pos]->setParentNode(curr);

                    }
                    else
                    {

                        queueNode.push(nullptr);

                    }


                    flag--;
                }

            }

        }




    public:
        KdTree()
        {

            laesaMatrix = nullptr;
            root = nullptr;
            pvt = nullptr;
            df = nullptr;
            numberOfNodes = 0;

        }

        KdTree(Dataset<DType>* dataset, DistanceFunction<BasicArrayObject<DType>>* df_, Pivot<DType>* pvt_, size_t numPerLeaf)
        {

            numberOfNodes = 0;
            df = df_;
            pvt = pvt_;
            root = nullptr;
            std::vector<BasicArrayObject<double>> elements(dataset->getCardinality());

            for(size_t i = 0; i < dataset->getCardinality(); i++)
            {

                BasicArrayObject<double> basic = BasicArrayObject<double>(i, pvt->getNumberOfPivots());

                for(size_t j = 0; j < pvt->getNumberOfPivots(); j++)
                {

                    basic.set(j, df->getDistance(*dataset->instance(i), *pvt->getPivot(j)));

                }

                elements[i] = basic;

            }

            laesaMatrix = new Dataset<double>(elements, dataset->getCardinality(), pvt->getNumberOfPivots());

            buildTree(numPerLeaf);

        }

        KdTree(std::string serializedPath, DistanceFunction<BasicArrayObject<DType>>* df_)
        {

            loadFromFile(serializedPath);
            setDistanceFunction(df_);

        }

        void setDistanceFunction(DistanceFunction<BasicArrayObject<DType>>* df_)
        {

            df = df_;

        }

        NodeKdTree* getRoot()
        {

            return root;

        }

        size_t getDistanceCount()
        {

            return (size_t)df->distCount;

        }

        void resetDistanceCount()
        {

            df->resetStatistics();

        }

        Dataset<double>* getLaesaMatrix()
        {

            return laesaMatrix;

        }

        Pivot<DType>* getPivot()
        {

            return pvt;

        }

        std::string toString()
        {

            std::queue<std::pair<NodeKdTree*, size_t>> queue;
            queue.push(std::make_pair(getRoot(), 0));
            std::stringstream st;

            while(!queue.empty())
            {

                std::pair<NodeKdTree*, size_t> element = queue.front();
                queue.pop();

                if(element.first->getLeft() != nullptr)
                {

                    queue.push(std::make_pair(element.first->getLeft(), element.second + 1));

                }

                if(element.first->getRight() != nullptr)
                {

                    queue.push(std::make_pair(element.first->getRight(), element.second + 1));

                }

                st << "Level = " << element.second << "\n";
                //st << "Node ID = " << element.first->getNodeID() << "\n";

                if(element.first->isLeafNode())
                {

                    LeafNodeKdTree* leaf = (LeafNodeKdTree*)element.first;
                    st << leaf->toString() << "\n\n\n";

                }
                else
                {

                    DirectoryNodeKdTree* dir = (DirectoryNodeKdTree*)element.first;
                    st << dir->toString() << "\n\n\n";

                }


            }

            return st.str();

        }

        size_t getSerializedSize()
        {

            size_t datasetSize = sizeof(size_t) + laesaMatrix->getSerializedSize();
            size_t pivotSize = sizeof(size_t) + pvt->getSerializedSize() + sizeof(PIVOT_TYPE);
            size_t nodesSize_ = nodesSize + sizeof(size_t)*2 + numberOfNodes * 2 * sizeof(size_t);
            return datasetSize + pivotSize + nodesSize_;

        }

        unsigned char* serialize()
        {

            unsigned char* data = new unsigned char[getSerializedSize()];
            size_t total = 0, datasetSize = laesaMatrix->getSerializedSize(), pivotSize = pvt->getSerializedSize();
            PIVOT_TYPE pvtType = pvt->getPivotType();

            memcpy(data + total, &numberOfNodes, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(data + total, &nodesSize, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(data + total, &datasetSize, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(data + total, &pivotSize, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(data + total, &pvtType, sizeof(PIVOT_TYPE));
            total += sizeof(PIVOT_TYPE);

            std::queue<NodeKdTree*> queueNode;
            queueNode.push(getRoot());
            NodeKdTree* curr = nullptr;

            size_t nodeType, nodeSize;
            unsigned char* nodeData = nullptr;

            while(!queueNode.empty())
            {

                curr = queueNode.front();
                queueNode.pop();

                if(curr->getLeft() != nullptr)
                {

                    queueNode.push(curr->getLeft());

                }

                if(curr->getRight() != nullptr)
                {

                    queueNode.push(curr->getRight());

                }

                nodeType = curr->isDirectoryNode() ? 1 : 0;
                nodeSize = curr->getSerializedSize();

                memcpy(data + total, &nodeType, sizeof(size_t));
                total += sizeof(size_t);
                memcpy(data + total, &nodeSize, sizeof(size_t));
                total += sizeof(size_t);

                nodeData = curr->serialize();

                memcpy(data + total, nodeData, nodeSize);
                total += nodeSize;

            }

            memcpy(data + total, laesaMatrix->serialize(), datasetSize);
            total += datasetSize;
            memcpy(data + total, pvt->serialize(), pivotSize);
            total += pivotSize;

            return data;

        }


        void unserialize(unsigned char* data)
        {

            size_t total = 0, datasetSize = 0, pivotSize = 0;
            PIVOT_TYPE pvtType;

            memcpy(&numberOfNodes, data + total, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(&nodesSize, data + total, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(&datasetSize, data + total, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(&pivotSize, data + total, sizeof(size_t));
            total += sizeof(size_t);
            memcpy(&pvtType, data + total, sizeof(PIVOT_TYPE));
            total += sizeof(PIVOT_TYPE);

            size_t nodeType, nodeSize;
            std::vector<NodeKdTree*> vecNodes(numberOfNodes);
            NodeKdTree* curr = nullptr;
            unsigned char* nodeData = nullptr;

            for(size_t i = 0; i < numberOfNodes; i++)
            {

                memcpy(&nodeType, data + total, sizeof(size_t));
                total += sizeof(size_t);
                memcpy(&nodeSize, data + total, sizeof(size_t));
                total += sizeof(size_t);

                if(nodeType == 1)
                {

                    curr = new DirectoryNodeKdTree();

                }
                else
                {

                    curr = new LeafNodeKdTree();

                }

                nodeData = new unsigned char[nodeSize];
                memcpy(nodeData, data + total, nodeSize);
                total += nodeSize;

                curr->unserialize(nodeData);

                vecNodes[i] = curr;

            }

            unsigned char* dataset_char = new unsigned char[datasetSize];
            unsigned char* pivot_char = new unsigned char[pivotSize];

            memcpy(dataset_char, data + total, datasetSize);
            total += datasetSize;
            memcpy(pivot_char, data + total, pivotSize);
            total += pivotSize;

            laesaMatrix = new Dataset<double>();
            laesaMatrix->unserialize(dataset_char);

            std::string names[] = {"RANDOM", "GNAT", "CONVEX", "KMEDOIDS", "MAXSEPARETED", "MAXVARIANCE", "SELECTION", "PCA", "SSS", "FFT", "IS", "HFI", "WDR"};

            if(names[pvtType] == "RANDOM")
                pvt = new RandomPivots<DType>();
            else if(names[pvtType] == "GNAT")
                pvt = new GnatPivots<DType>();
            else if(names[pvtType] == "CONVEX")
                pvt = new ConvexPivots<DType>();
            else if(names[pvtType] == "KMEDOIDS")
                pvt = new KmedoidsPivots<DType>();
            else if(names[pvtType] == "MAXSEPARETED")
                pvt = new MaxSeparetedPivots<DType>();
            else if(names[pvtType] == "MAXVARIANCE")
                pvt = new MaxVariancePivots<DType>();
            else if(names[pvtType] == "SELECTION")
                pvt = new SelectionPivots<DType>();
            else if(names[pvtType] == "PCA")
                pvt = new PCAPivots<DType>();
            else if(names[pvtType] == "SSS")
                pvt = new SSSPivots<DType>();
            else if(names[pvtType] == "FFT")
                pvt = new FFTPivots<DType>();
            else if(names[pvtType] == "IS")
                pvt = new ISPivots<DType>();
            else if(names[pvtType] == "HFI")
                pvt = new HFIPivots<DType>();
            else if(names[pvtType] == "WDR")
                pvt = new WDRPivots<DType>();
            else
                throw std::invalid_argument("Can't read pivot type !_!");


            pvt->unserialize(pivot_char);

            buildTree(vecNodes);

        }

        void saveToFile()
        {

            int check = -1;

            std::string fileName;

            while(check == -1)
            {

                std::string directory = "../../gervLib/kdtree/kdtree_files/",
                        name = selectFolderName(directory, "kdtree"),
                        folderName = directory + name;

                fileName = folderName + "/" + name + ".dat";

                check = mkdir(folderName.c_str(), 0777);

            }

            std::ofstream file(fileName, std::ios::out | std::ios::binary);

            unsigned char* data = new unsigned char[sizeof(size_t) + getSerializedSize()];
            size_t sizeData = getSerializedSize();

            //std::cout << "TAM SAVE: " << sizeof(size_t) + sizeData << "\n";

            memcpy(data, &sizeData, sizeof(size_t));
            memcpy(data + sizeof(size_t), serialize(), sizeData);

            file.write((char*)data, sizeof(size_t) + sizeData);
            file.close();

//            int check = -1;

//                std::string fileName;

//                while(check == -1)
//                {

//                    std::string directory = "../kdtree/kdtree_files/",
//                            name = selectFolderName(directory, "kdtree"),
//                            folderName = directory + name;

//                    fileName = folderName + "/" + name + ".dat";

//                    check = mkdir(folderName.c_str(), 0777);

//                }

//                std::ofstream file(fileName, std::ios::out | std::ios::binary);

//                u_int8_t* seri = new u_int8_t[sizeof(u_int32_t) + getSerializedSize()];
//                u_int32_t size = getSerializedSize();
//                memcpy(seri, &size, sizeof(u_int32_t));
//                memcpy(seri + sizeof(u_int32_t), serialize(), size);

//                file.write((char*)seri, size + sizeof(u_int32_t));
//                file.close();

        }

        void loadFromFile(std::string fileName)
        {

            std::ifstream file(fileName, std::ios::in | std::ios::binary);

            unsigned char* dataSize_char = new unsigned char[sizeof(size_t)];
            size_t dataSize;

            file.read((char*)dataSize_char, sizeof(size_t));
            memcpy(&dataSize, dataSize_char, sizeof(size_t));

            unsigned char* data = new unsigned char[sizeof(size_t) + dataSize];

            //std::cout << "TAM LOAD: " << sizeof(size_t) + dataSize << "\n";

            file.seekg(0);
            file.read((char*)data, dataSize + sizeof(size_t));

            unserialize(data + sizeof(size_t));

            file.close();

//            std::ifstream file(fileName, std::ios::in | std::ios::binary);

//                u_int8_t* seri = new u_int8_t[sizeof(u_int32_t)];
//                u_int32_t size;

//                file.read((char*)seri, sizeof(u_int32_t));

//                memcpy(&size, seri, sizeof(u_int32_t));

//                u_int8_t* data = new u_int8_t[size + sizeof(u_int32_t)];

//                file.seekg(0);
//                file.read((char *)data, size + sizeof(u_int32_t));

//                unserialize(data + sizeof(u_int32_t));

//                file.close();

        }

        bool isInterval(double infBound, double supBound, double test)
        {

            return ((test >= infBound) && (test <= supBound));

        }

        double minDist(BasicArrayObject<double> sq_, std::vector<Bound> bounds)
        {

            double limInfCase3 = -1.0;
            double limInfCase2 = -1.0;
            double answer = -1.0;
            bool within = true;

            for(size_t x = 0; x < bounds.size(); x++)
            {

                if(!isInterval(std::abs(bounds[x].array[0]), std::abs(bounds[x].array[1]), *sq_.get(x)))
                {

                    within = false;

                    limInfCase3 = std::max(limInfCase3,

                                                       std::min(

                                                           std::abs(*sq_.get(x) - std::abs(bounds[x].array[0])),

                                                       std::abs(*sq_.get(x) - std::abs(bounds[x].array[1]))

                                        )

                                        );
                }
                else
                {

                    limInfCase2 = std::min(limInfCase2,

                                                       std::min(

                                                           std::abs(*sq_.get(x) - std::abs(bounds[x].array[0])),

                                                       std::abs(*sq_.get(x) - std::abs(bounds[x].array[1]))

                                        )

                                        );

                }

            }

            if(within)
            {

                answer = 0.0;

            }
            else
            {

                if(limInfCase2 != -1.0)
                {

                    answer = limInfCase2;

                }
                else
                {

                    answer = limInfCase3;

                }

            }

            return answer;

        }

        double maxDist(BasicArrayObject<double> sq_, std::vector<Bound> bounds)
        {

            double answer = -1.0;

            answer = std::numeric_limits<double>::max();

            for(size_t x = 0; x < bounds.size(); x++)
            {

                if (std::numeric_limits<double>::max() -  std::abs(bounds[x].array[0]) >= *sq_.get(x))
                {

                    answer = std::min(answer, *sq_.get(x) + std::abs(bounds[x].array[0]));

                }

                if (std::numeric_limits<double>::max() -  std::abs(bounds[x].array[1]) >= *sq_.get(x))
                {

                    answer = std::min(answer, *sq_.get(x) + std::abs(bounds[x].array[1]));

                }

            }

            return answer;

        }

        BasicArrayObject<double> pivotCoordinates(BasicArrayObject<DType>* inst)
        {

            BasicArrayObject<double> ans = BasicArrayObject<double>(0, pvt->getNumberOfPivots());

            for(size_t x = 0; x < pvt->getNumberOfPivots(); x++)
            {

                ans.set(x, df->getDistance(*inst, *pvt->getPivot(x)));

            }

            return ans;

        }

        void kNN(Dataset<DType>* dataset, BasicArrayObject<DType>* query, size_t k, std::vector<PairResult>& ans)
        {

            //resetDiskAccess();
            df->resetStatistics();
            Dataset<double>* dataLeafNode;
            LeafNodeKdTree* leaf;
            std::priority_queue<Partition, std::vector<Partition>, ComparePartition> pqPartition;
            std::priority_queue<PairResult, std::vector<PairResult>, ComparePairResult> pqCandidates;
            std::priority_queue<PairResult, std::vector<PairResult>, ComparePairResult2> pqAns;
            BasicArrayObject<double> sq_ = pivotCoordinates(query);
            pqPartition.push(Partition(getRoot(), 0.0, 0.0));
            NodeKdTree* node = nullptr;
            Partition partition;

            while(!pqPartition.empty() || pqCandidates.size() > 0)
            {

                if(pqCandidates.size() == 0)
                {

                    partition = pqPartition.top();
                    node = partition.node;
                    pqPartition.pop();

                    if(node->isLeafNode())
                    {

                        leaf = (LeafNodeKdTree*)node;
                        dataLeafNode = leaf->getDataset();

                        for(size_t x = 0; x < dataLeafNode->getCardinality(); x++)
                        {

                            pqCandidates.push(PairResult(dataLeafNode->getInstance(x)->getOID(), df->getDistance(*query, *dataset->getInstance(dataLeafNode->getInstance(x)->getOID()))));

                        }

                        //delete dataLeafNode;

                    }
                    else
                    {

                        double minL = minDist(sq_, node->getLeft()->getBoundary());
                        double maxL = maxDist(sq_, node->getLeft()->getBoundary());
                        double minR = minDist(sq_, node->getRight()->getBoundary());
                        double maxR = maxDist(sq_, node->getRight()->getBoundary());

                        pqPartition.push(Partition(node->getLeft(), minL, maxL));
                        pqPartition.push(Partition(node->getRight(), minR, maxR));

                    }

                }
                else if((pqPartition.size() > 0) && pqPartition.top().min < df->getDistance(*query, *dataset->getInstance(pqCandidates.top().index)))
                {

                    node = pqPartition.top().node;
                    pqPartition.pop();

                    if(node->isLeafNode())
                    {

                        leaf = (LeafNodeKdTree*)node;
                        dataLeafNode = leaf->getDataset();

                        for(size_t x = 0; x < dataLeafNode->getCardinality(); x++)
                        {

                            pqCandidates.push(PairResult(dataLeafNode->getInstance(x)->getOID(), df->getDistance(*query, *dataset->getInstance(dataLeafNode->getInstance(x)->getOID()))));

                        }

//                        delete dataLeafNode;

                    }
                    else
                    {

                        double minL = minDist(sq_, node->getLeft()->getBoundary());
                        double maxL = maxDist(sq_, node->getLeft()->getBoundary());
                        double minR = minDist(sq_, node->getRight()->getBoundary());
                        double maxR = maxDist(sq_, node->getRight()->getBoundary());

                        pqPartition.push(Partition(node->getLeft(), minL, maxL));
                        pqPartition.push(Partition(node->getRight(), minR, maxR));

                    }

                }
                else
                {

                    if(!pqAns.empty() && !pqCandidates.empty() && pqAns.size() >= k && pqCandidates.top().distance > pqAns.top().distance)
                    {

                        break;

                    }

                    pqAns.push(pqCandidates.top());
                    pqCandidates.pop();

                    while(pqAns.size() > k)
                    {

                        pqAns.pop();

                    }

                }

            }


            ans = dequeueInOrder(pqAns);
            std::reverse(ans.begin(), ans.end());

            while(!pqCandidates.empty())
            {

                pqCandidates.pop();

            }

            while(!pqAns.empty())
            {

                pqAns.pop();

            }

            while(!pqPartition.empty())
            {

                pqPartition.pop();

            }

        }

};
















#endif // KDTREE_H
