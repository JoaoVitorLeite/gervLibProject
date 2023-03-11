#ifndef OMNIKDTREE_H
#define OMNIKDTREE_H

#include <KdTree.h>


template <class DType>
class OmniKdTree
{

    private:
        KdTree<DType>* tree;
        size_t diskAccess;
        std::string path;
        Pivot<DType>* pivot;
        DistanceFunction<BasicArrayObject<DType>>* df;
        Dataset<double>* laesaMatrix;
        //Dataset<DType>* dataset;

    private:

        void incrementDiskAccess()
        {

            diskAccess++;

        }

        void setPath()
        {

            int check = -1;

            while(check == -1)
            {

                std::string directory = "../omni/omni_files/",
                        name = selectFolderName(directory, "omni");

                path = directory + name + "/";

                check = mkdir(path.c_str(), 0777);

            }

        }

        bool checkIfFileExists(std::string path)
        {

            std::ifstream file;
            file.open(path);

            if(file)
            {

                return true;

            }
            else
            {

                return false;

            }

        }

        void saveLeafNode(LeafNodeKdTree* node)
        {

            std::string fileName = path + "leafnode_" + std::to_string(node->getNodeID()) + ".dat";

            if(!checkIfFileExists(fileName))
            {

                size_t size = node->getDataset()->getSerializedSize();
                unsigned char* data = new unsigned char[size + sizeof(size_t)];

                memcpy(data, &size, sizeof(size_t));
                memcpy(data + sizeof(size_t), node->getDataset()->serialize(), size);

                node->getDataset()->clear();
                node->getDataset()->setCardinality(0);
                node->getDataset()->setDimensionality(0);
                node->getDataset()->setSeed(0);
//                node->setDataset(nullptr);

                std::ofstream file(fileName, std::ios::out | std::ios::binary);

                file.write((char*)data, size + sizeof(size_t));
                file.close();

                delete [] data;

            }
//            else
//            {

//                node->setDataset(nullptr);

//            }

        }

        void saveAllLeafNodes()
        {

            std::queue<NodeKdTree*> queueNode;
            queueNode.push(tree->getRoot());
            NodeKdTree* curr = nullptr;

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

                if(curr->isLeafNode())
                {

                    saveLeafNode((LeafNodeKdTree*)curr);

                }

            }

        }

//        Dataset<double>* readLeafNode(LeafNodeKdTree* node)
//        {

//            Dataset<double>* dataset = new Dataset<double>();
//            std::string fileName = path + "leafnode_" + std::to_string(node->getNodeID()) + ".dat";

//            if(checkIfFileExists(fileName))
//            {

//                std::ifstream file(fileName, std::ios::in | std::ios::binary);

//                size_t size;
//                unsigned char* dataSize = new unsigned char[sizeof(size_t)];

//                file.read((char*)dataSize, sizeof(size_t));

//                memcpy(&size, dataSize, sizeof(size_t));

//                unsigned char* data = new unsigned char[size + sizeof(size_t)];

//                file.seekg(0);
//                file.read((char*)data, size + sizeof(size_t));

//                dataset->unserialize(data + sizeof(size_t));

//                file.close();

//                delete [] data;
//                delete [] dataSize;

//            }
//            else
//            {

//                throw std::invalid_argument("File does not exist");

//            }

//            return dataset;

//        }

//        void readLeafNode(LeafNodeKdTree* node, Dataset<double>* dataset)
//        {

//            //Dataset<double>* dataset = new Dataset<double>();
//            std::string fileName = path + "leafnode_" + std::to_string(node->getNodeID()) + ".dat";

//            if(checkIfFileExists(fileName))
//            {

//                std::ifstream file(fileName, std::ios::in | std::ios::binary);

//                size_t size;
//                unsigned char* dataSize = new unsigned char[sizeof(size_t)];

//                file.read((char*)dataSize, sizeof(size_t));

//                memcpy(&size, dataSize, sizeof(size_t));

//                unsigned char* data = new unsigned char[size + sizeof(size_t)];

//                file.seekg(0);
//                file.read((char*)data, size + sizeof(size_t));

//                dataset->unserialize(data + sizeof(size_t));

//                file.close();

//                delete [] data;
//                delete [] dataSize;

//            }
//            else
//            {

//                throw std::invalid_argument("File does not exist");

//            }

////            return dataset;

//        }

        void readLeafNode(LeafNodeKdTree* node)
        {

//            Dataset<double>* dataset = new Dataset<double>();
            std::string fileName = path + "leafnode_" + std::to_string(node->getNodeID()) + ".dat";

            if(checkIfFileExists(fileName))
            {

                std::ifstream file(fileName, std::ios::in | std::ios::binary);

                size_t size;
                u_char* dataSize = new u_char[sizeof(size_t)];

                file.read((char*)dataSize, sizeof(size_t));

                memcpy(&size, dataSize, sizeof(size_t));

                u_char* data = new u_char[size + sizeof(size_t)];

                file.seekg(0);
                file.read((char*)data, size + sizeof(size_t));

                node->getDataset()->unserialize(data + sizeof(size_t));

                file.close();

                delete [] data;
                delete [] dataSize;

            }
            else
            {

                throw std::invalid_argument("File does not exist");

            }

//            return dataset;

        }


    public:
        OmniKdTree()
        {

            tree = nullptr;
            diskAccess = 0;
            pivot = nullptr;
            df = nullptr;
            laesaMatrix = nullptr;
            setPath();

        }

        OmniKdTree(Dataset<DType>* dataset_, DistanceFunction<BasicArrayObject<DType>>* df_, Pivot<DType>* pivot_, size_t numPerLeaf)
        {

            tree = new KdTree<DType>(dataset_, df_, pivot_, numPerLeaf);
            //std::cout << tree->toString() << std::endl;
            diskAccess = 0;
            pivot = pivot_;
            df = df_;
            laesaMatrix = tree->getLaesaMatrix();
            //std::cout << "\n\n";
            //laesaMatrix->printDataset();
            setPath();
            saveToFile();
            saveAllLeafNodes();

        }

        OmniKdTree(std::string serializedPath, DistanceFunction<BasicArrayObject<DType>>* df_)
        {

            loadFromFile(serializedPath);
            saveAllLeafNodes();
            setDistanceFunction(df_);
            diskAccess = 0;
            path = serializedPath.substr(0, serializedPath.find_last_of("/") + 1);

        }

        size_t getDiskAccess()
        {

            return diskAccess;

        }

        size_t getDistanceCount()
        {

            return (size_t)df->getDistanceCount();

        }

        KdTree<DType>* getTree()
        {

            return tree;

        }

        std::string toString()
        {

            return tree->toString();

        }

        void setDistanceFunction(DistanceFunction<BasicArrayObject<DType>>* df_)
        {

            df = df_;

        }

        void saveToFile()
        {

            std::string nameAux = path.substr(path.rfind("/", path.size()-2)+1, path.find_last_of("/")),
                    name = nameAux.substr(0, nameAux.size()-1),
                    fileName = path + name + ".dat";

            std::ofstream file(fileName, std::ios::out | std::ios::binary);

            size_t size = tree->getSerializedSize();
            unsigned char* data = new unsigned char[size + sizeof(size_t)];

            memcpy(data, &size, sizeof(size_t));
            memcpy(data + sizeof(size_t), tree->serialize(), size);

            file.write((char*)data, size + sizeof(size_t));
            file.close();

        }

        void loadFromFile(std::string fileName)
        {

            std::ifstream file(fileName, std::ios::in | std::ios::binary);

            unsigned char* seri = new unsigned char[sizeof(size_t)];
            size_t size;

            file.read((char*)seri, sizeof(size_t));

            memcpy(&size, seri, sizeof(size_t));

            unsigned char* data = new unsigned char[size + sizeof(size_t)];

            file.seekg(0);
            file.read((char*)data, size + sizeof(size_t));

            tree = new KdTree<DType>();
            tree->unserialize(data + sizeof(size_t));

            pivot = tree->getPivot();
            //df = tree->getDistanceFunction();
            laesaMatrix = tree->getLaesaMatrix();


        }

        void resetDiskAccess()
        {

            diskAccess = 0;

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

            BasicArrayObject<double> ans = BasicArrayObject<double>(0, pivot->getNumberOfPivots());

            for(size_t x = 0; x < pivot->getNumberOfPivots(); x++)
            {

                ans.set(x, df->getDistance(*inst, *pivot->getPivot(x)));

            }

            return ans;

        }

        void kNN(Dataset<DType>* dataset, BasicArrayObject<DType>* query, size_t k, std::vector<PairResult>& ans)
        {

            resetDiskAccess();
            df->resetStatistics();
            Dataset<double>* dataLeafNode;
            LeafNodeKdTree* leaf;
            std::priority_queue<Partition, std::vector<Partition>, ComparePartition> pqPartition;
            std::priority_queue<PairResult, std::vector<PairResult>, ComparePairResult> pqCandidates;
            std::priority_queue<PairResult, std::vector<PairResult>, ComparePairResult2> pqAns;
            BasicArrayObject<double> sq_ = pivotCoordinates(query);
            //std::cout << sq_.toString() << "\n";
            pqPartition.push(Partition(tree->getRoot(), 0.0, 0.0));
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

                        readLeafNode((LeafNodeKdTree*)node);
                        leaf = (LeafNodeKdTree*)node;
                        dataLeafNode = leaf->getDataset();

                        //Dataset<double>* dataLeafNode = readLeafNode((LeafNodeKdTree*)node);
                        //Dataset<double>* dataLeafNode = static_cast<LeafNodeKdTree*>(node)->getDataset();
                        //readLeafNode((LeafNodeKdTree*)node, dataLeafNode);
                        incrementDiskAccess();

                        for(size_t x = 0; x < dataLeafNode->getCardinality(); x++)
                        {

                            pqCandidates.push(PairResult(dataLeafNode->getInstance(x)->getOID(), df->getDistance(*query, *dataset->getInstance(dataLeafNode->getInstance(x)->getOID()))));

                        }

                        leaf->getDataset()->clear();
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

                        readLeafNode((LeafNodeKdTree*)node);
                        leaf = (LeafNodeKdTree*)node;
                        dataLeafNode = leaf->getDataset();
                        //Dataset<double>* dataLeafNode = readLeafNode((LeafNodeKdTree*)node);
                        incrementDiskAccess();

                        for(size_t x = 0; x < dataLeafNode->getCardinality(); x++)
                        {

                            pqCandidates.push(PairResult(dataLeafNode->getInstance(x)->getOID(), df->getDistance(*query, *dataset->getInstance(dataLeafNode->getInstance(x)->getOID()))));

                        }

                        leaf->getDataset()->clear();
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

#endif // OMNIKDTREE_H
