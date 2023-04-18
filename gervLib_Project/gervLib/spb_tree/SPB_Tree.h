#ifndef SPB_TREE_H
#define SPB_TREE_H

#include <Hilbert.h>
#include <Hermes.h>
#include <btree_multimap.h>
#include <Dataset.h>
#include <fstream>
#include <chrono>
#include <EquiDepth.h>
#include <set>
#include <Pivots.h>

//using std::less, std::ofstream, std::floor, std::pair, std::make_pair, std::sort, std::set, std::pair, std::numeric_limits;

typedef stx::btree_multimap<ull, size_t, std::less<ull> > btree_type;

struct SPBPartition
{

    btree_type::Node* node;
    double min, max;

    SPBPartition()
    {



    }

    //SPBPartition(btree_type::Node* _node, ll _min, ll _max)
    SPBPartition(btree_type::Node* _node, double _min, double _max)
    {

        node = _node;
        min = _min;
        max = _max;

    }

};

struct ComparePartitionSPB
{

    bool operator()(SPBPartition const& a, SPBPartition const& b)
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

template <class T>
struct KnnSPB
{

    BasicArrayObject<T> element;
    double distance;

    KnnSPB()
    {



    }

    KnnSPB(BasicArrayObject<T> _element, double _distance)
    {

        element = _element;
        distance = _distance;

    }

    bool operator<(const KnnSPB<T>& item) const
    {

        return distance < item.distance;

    }

    bool operator>(const KnnSPB<T>& item) const
    {

        return distance > item.distance;

    }

};

template<class Type, class Comp>
std::vector<Type> dequeueInOrderSPB_Results(std::priority_queue<Type, std::vector<Type>, Comp> pq)
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



template <class T>
class SPBTree
{

private:
    HilbertCurve hc;
    EquiDepth<double> ed;
    DistanceFunction<BasicArrayObject<T>>* df;
    Pivot<T>* pvt;
//    vector<size_t> pivotsIds;
//    vector<BasicArrayObject<T>> pivots;
    btree_type bt;
    ll buildTime, leafNodeAccess;

private:
    void bulk_load(Dataset<T>* dataset)
    {

        vector<vector<double>> pivot_mapping;
        pivot_mapping.resize(dataset->getCardinality(), vector<double>(PIVOT_NUM));

        vector<vector<ull>> disc;
        disc.resize(dataset->getCardinality(), vector<ull>(PIVOT_NUM));

        vector<ull> keys(dataset->getCardinality());

        std::ofstream file_pivot_mapping(baseFilePath + fs::path::preferred_separator + "bulk_load_pivot_mapping.txt");
        std::ofstream file_disc(baseFilePath + fs::path::preferred_separator + "bulk_load_disc.txt");
        std::ofstream file_sfc(baseFilePath + fs::path::preferred_separator + "bulk_load_sfc.txt");

        for(size_t i = 0; i < dataset->getCardinality(); i++)
        {

            for(size_t j = 0; j < PIVOT_NUM; j++)
            {

                pivot_mapping[i][j] = df->getDistance(dataset->getFeatureVector(i), *pvt->getPivot(j));
                file_pivot_mapping << pivot_mapping[i][j] << " ";

            }

            file_pivot_mapping << "\n";

        }

        ed.build(pivot_mapping);
        ed.saveToFile();
        file_pivot_mapping.close();

        for(size_t i = 0; i < dataset->getCardinality(); i++)
        {

            for(size_t j = 0; j < PIVOT_NUM; j++)
            {

                disc[i][j] = ed.getBin(j, pivot_mapping[i][j]);
                file_disc << disc[i][j] << " ";

            }

            file_disc << "\n";

            keys[i] = hc.distance_from_point(disc[i]);
            file_sfc << keys[i] << "\n";

        }

        pivot_mapping.clear();
        disc.clear();
        ed.clear();
        file_disc.close();
        file_sfc.close();

        std::vector<std::pair<ull, size_t>> insertValues(keys.size());

        for(size_t i = 0; i < keys.size(); i++)
        {

            insertValues[i] = std::make_pair(keys[i], i);

        }

        keys.clear();

        std::sort(insertValues.begin(), insertValues.end());
        bt.setHilbertCurve(hc);
        bt.bulk_load(insertValues.begin(), insertValues.end());

        insertValues.clear();

        bt.init_Key_Minmax();
        bt.initDisk(dataset);

        delete dataset;

//        vector<vector<double>> pointsDist;
//        vector<vector<ull>> pointsDisc;
//        vector<ull> keys;

//        //ofstream file(baseFilePath + std::filesystem::path::preferred_separator + "bulk_load_disc.txt");

//        for(size_t i = 0; i < dataset->getCardinality(); i++)
//        {

//            pointsDist.push_back(vector<double>(PIVOT_NUM));

//            for(size_t j = 0; j < (size_t)PIVOT_NUM; j++)
//            {

//                pointsDist[i][j] = df->getDistance(dataset->getFeatureVector(pivotsIds[j]), dataset->getFeatureVector(i));

//            }

//        }

//        ed.build(pointsDist);

//        for(size_t i = 0; i < dataset->getCardinality(); i++)
//        {

//            pointsDisc.push_back(vector<ull>(PIVOT_NUM));

//            for(size_t j = 0; j < (size_t)PIVOT_NUM; j++)
//            {

//                pointsDisc[i][j] = ed.getBin(j, pointsDist[i][j]);
//                //file << pointsDisc[i][j] << " ";

//            }

//            //file << endl;

//        }

//        //file.close();

//        keys = hc.distances_from_points(pointsDisc);

//        //pointsDist.clear();
//        //pointsDisc.clear();

//        //ofstream file2(baseFilePath + std::filesystem::path::preferred_separator + "bulk_load_sfc.txt");

//        vector<pair<ull, size_t>> insertValues;

//        for(size_t i = 0; i < keys.size(); i++)
//        {

//            insertValues.push_back(make_pair(keys[i], i));
//            //file2 << keys[i] << endl;

//        }

//        sort(insertValues.begin(), insertValues.end());
//        //bt.setHilbertCurve(hc);
//        //bt.bulk_load(insertValues.begin(), insertValues.end());
//        //bt.init_Key_Minmax();
//        //bt.initDisk(dataset);

//        ed.saveToFile();

        //keys.clear();
        //insertValues.clear();
        //file2.close();
        //delete dataset;

    }

public:
    SPBTree()
    {

        setBaseFilePath();
//        pivotsIds = {};
//        pivots = {};
        df = nullptr;
        hc = HilbertCurve();
        ed = EquiDepth<double>();

    }

    SPBTree(Dataset<T>* dataset, DistanceFunction<BasicArrayObject<T>>* _df, Pivot<T>* pvt_, ull pivot_num, long num_bins)
    {

        auto start = std::chrono::steady_clock::now();

        baseFilePath = "../spb_tree/spb_files";
        setBaseFilePath();
//        pivotsIds = {(size_t)2, (size_t)12};
//        pivots = {dataset->getFeatureVector(2), dataset->getFeatureVector(12)};
        df = _df;
        pvt = pvt_;
        pvt->generatePivots(dataset, df, pivot_num);
        PIVOT_NUM = pivot_num;
        ed = EquiDepth<double>(num_bins, pivot_num);
        p = (ull)log2(num_bins-1) + 1;
        hc = HilbertCurve(p, pivot_num);
        GRID_L = ((1ull << p) - (ull)1);

        bulk_load(dataset);

        auto end = std::chrono::steady_clock::now();
        buildTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    }

    void dump_key_min_max()
    {

        bt.dump_Key_Minmax();

    }

    void dumpHilbertCurve()
    {

        cout << hc << endl;

    }

    template <typename U>
    bool isInterval(U infBound, U supBound, U test)
    {

        return ((test >= infBound) && (test <= supBound));

    }

    ll getLeafNodeAccess()
    {

        return leafNodeAccess;

    }

    size_t getDistanceCount()
    {

        return df->getDistanceCount();

    }

    double minDist(vector<double> sq_, ull** mbr_)
    {

        double** mbr = new double*[PIVOT_NUM];

        for(size_t i = 0; i < PIVOT_NUM; i++)
        {

            mbr[i] = new double[2];
            //mbr[i][0] = numeric_limits<double>::max();
            //mbr[i][1] = numeric_limits<double>::lowest();

            std::pair<double, double> pairMin = ed.getInterval(i, mbr_[i][0]);
            std::pair<double, double> pairMax = ed.getInterval(i, mbr_[i][1]);

            //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
              //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMax.first});
            mbr[i][0] = pairMin.first;

            //mbr[i][1] = std::max({mbr[i][1], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
              //mbr[i][1] = std::max({mbr[i][1], pairMin.second, pairMax.second});
            mbr[i][1] = pairMax.second;


        }

        double limInfCase3 = -1.0;
        double limInfCase2 = -1.0;
        double answer = -1.0;
        bool within = true;

        for(size_t i = 0; i < PIVOT_NUM; i++)
        {

            if(!isInterval(mbr[i][0], mbr[i][1], sq_[i]))
            {

                within = false;

                limInfCase3 = std::max(limInfCase3,
                                       std::min(
                                           std::abs(sq_[i] - mbr[i][0]),
                                       std::abs(sq_[i] - mbr[i][1])
                        )
                        );

            }
            else
            {

                limInfCase2 = std::min(limInfCase2,
                                       std::min(
                                           std::abs(sq_[i] - mbr[i][0]),
                                       std::abs(sq_[i] - mbr[i][1])
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

            if(limInfCase2 != -1)
            {

                answer = limInfCase2;

            }
            else
            {

                answer = limInfCase3;

            }

        }

        for(size_t i = 0; i < PIVOT_NUM; i++)
        {

            delete [] mbr[i];

        }
        delete [] mbr;
        sq_.clear();

//        cout << "MIN D = " << answer << endl;

        return answer;

    }

    double maxDist(vector<double> sq_, ull** mbr_)
    {


        double** mbr = new double*[PIVOT_NUM];

        for(size_t i = 0; i < PIVOT_NUM; i++)
        {

            mbr[i] = new double[2];
            //mbr[i][0] = numeric_limits<double>::max();
            //mbr[i][1] = numeric_limits<double>::lowest();

            std::pair<double, double> pairMin = ed.getInterval(i, mbr_[i][0]);
            std::pair<double, double> pairMax = ed.getInterval(i, mbr_[i][1]);

            //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
              //mbr[i][0] = std::min({mbr[i][0], pairMin.first, pairMax.first});
            mbr[i][0] = pairMin.first;

            //mbr[i][1] = std::max({mbr[i][1], pairMin.first, pairMin.second, pairMax.first, pairMax.second});
              //mbr[i][1] = std::max({mbr[i][1], pairMin.second, pairMax.second});
            mbr[i][1] = pairMax.second;

        }

        double answer = std::numeric_limits<double>::max();

        for(size_t i = 0; i < PIVOT_NUM; i++)
        {

            if(std::numeric_limits<double>::max() - mbr[i][0] >= sq_[i])
            {

                answer = std::min(answer, sq_[i] + mbr[i][0]);

            }

            if(std::numeric_limits<double>::max() - mbr[i][1] >= sq_[i])
            {

                answer = std::min(answer, sq_[i] + mbr[i][1]);

            }


        }

        for(size_t i = 0; i < PIVOT_NUM; i++)
        {

            delete [] mbr[i];

        }
        delete [] mbr;
        sq_.clear();

        //cout << "MAX D = " << answer << endl;

        return answer;

    }

    size_t getPageID(btree_type::Leaf* node, size_t pos)
    {

        if(pos >= 0 && pos <= node->slotuse)
        {

            return (size_t)(node->memory_cost[pos] + sizeof(size_t))/PAGE_SIZE;

        }
        else
        {

            throw std::invalid_argument("Leaf node has fewer elements than the requested position");

        }

    }

    void knn(BasicArrayObject<T> query, size_t k, std::vector<KnnSPB<T>>& ans)
    {

        //ed.readFromFile();
        ed.load();
        df->resetStatistics();
        leafNodeAccess = 0;
        IOread = 0;
        ans.clear();

        std::priority_queue<SPBPartition, std::vector<SPBPartition>, ComparePartitionSPB> nodeQueue;
        std::priority_queue<KnnSPB<T>, std::vector<KnnSPB<T>>, std::greater<KnnSPB<T>>> candidatesQueue;
        std::priority_queue<KnnSPB<T>, std::vector<KnnSPB<T>>, std::less<KnnSPB<T>>> resultQueue;
        nodeQueue.push(SPBPartition(bt.getRoot(), 0.0, std::numeric_limits<double>::infinity()));

        vector<double> sq_ = vector<double>(PIVOT_NUM);

        for(size_t i = 0; i < PIVOT_NUM; i++)
        {

            sq_[i] = df->getDistance(query, *pvt->getPivot(i));
//            cout << "SQ_ = " << sq_[i] << endl;
//            cout << "BIN = " << ed.getBin(i, sq_[i]) << endl;

        }

        SPBPartition partition;
        btree_type::Node* node;
        btree_type::Leaf *leafnode;
        btree_type::Inner* innernode;
        std::vector<BasicArrayObject<T>*> dataLeaf;
        std::set<size_t> globalPagesID;

        while((!nodeQueue.empty() || candidatesQueue.size() > 0) && resultQueue.size() < k)
        {

            if(candidatesQueue.size() == 0)
            {

                partition = nodeQueue.top();
                node = partition.node;
                nodeQueue.pop();

                if(node->isleafnode())
                {

                    leafnode = static_cast<btree_type::Leaf*>(node);
                    leafNodeAccess++;

                    for(size_t i = 0; i < leafnode->slotuse; i++)
                    {

                        size_t pid = getPageID(leafnode, i);
                        std::pair<std::set<size_t>::iterator, bool> insertPID = globalPagesID.insert(pid);

                        if(insertPID.second)
                        {

                            read_from_disk(dataLeaf, pid);

                        }

                    }

                    for(size_t i = 0; i < dataLeaf.size(); i++)
                    {

                        candidatesQueue.push(KnnSPB<T>(*dataLeaf[i], df->getDistance(query, *dataLeaf[i])));
                        delete dataLeaf[i];

                    }

                    dataLeaf.clear();

                }
                else
                {

                    innernode = static_cast<btree_type::Inner*>(node);

                    for(size_t i = 0; i < (size_t)(innernode->slotuse + 1); i++)
                    {

                        nodeQueue.push(SPBPartition(innernode->childid[i], minDist(sq_, innernode->childid[i]->mbr), maxDist(sq_, innernode->childid[i]->mbr)));

                    }

                }

            }
            else if((nodeQueue.size() > 0) && nodeQueue.top().min < candidatesQueue.top().distance)
            {

                partition = nodeQueue.top();
                node = partition.node;
                nodeQueue.pop();

                if(node->isleafnode())
                {

                    leafnode = static_cast<btree_type::Leaf*>(node);
                    leafNodeAccess++;

                    for(size_t i = 0; i < leafnode->slotuse; i++)
                    {

                        size_t pid = getPageID(leafnode, i);
                        std::pair<std::set<size_t>::iterator, bool> insertPID = globalPagesID.insert(pid);

                        if(insertPID.second)
                        {

                            read_from_disk(dataLeaf, pid);

                        }

                    }

                    for(size_t i = 0; i < dataLeaf.size(); i++)
                    {

                        candidatesQueue.push(KnnSPB<T>(*dataLeaf[i], df->getDistance(query, *dataLeaf[i])));
                        delete dataLeaf[i];

                    }

                    dataLeaf.clear();

                }
                else
                {

                    innernode = static_cast<btree_type::Inner*>(node);

                    for(size_t i = 0; i < (size_t)(innernode->slotuse + 1); i++)
                    {

                        nodeQueue.push(SPBPartition(innernode->childid[i], minDist(sq_, innernode->childid[i]->mbr), maxDist(sq_, innernode->childid[i]->mbr)));

                    }

                }

            }
            else
            {

                if(!resultQueue.empty() && !candidatesQueue.empty() && resultQueue.size() >= k && candidatesQueue.top().distance > resultQueue.top().distance)
                {

                    break;

                }

                resultQueue.push(candidatesQueue.top());
                candidatesQueue.pop();

                while(resultQueue.size() > k)
                {

                    resultQueue.pop();

                }

            }

        }

        ans = dequeueInOrderSPB_Results(resultQueue);
        std::reverse(ans.begin(), ans.end());

        while(!candidatesQueue.empty())
        {

            candidatesQueue.pop();

        }

        while(!resultQueue.empty())
        {

            resultQueue.pop();

        }

        while(!nodeQueue.empty())
        {

            nodeQueue.pop();

        }

    }

    void dumpEquiDepth()
    {

        ed.load();
        ed.print();

    }

    void test()
    {

        bt.test();

    }

};

//////#include <Hilbert.h>
////#include <Hermes.h>
////#include <Dataset.h>
////#include <btree_multimap.h>
////#include <fstream>
////#include <cstdlib>
////#include <chrono>
////#include <set>

////using std::less, std::ofstream, std::floor, std::pair, std::make_pair, std::sort, std::set;

////typedef stx::btree_multimap<ull, size_t, less<ull> > btree_type;
//////typedef typename btree_type::Node NODE;
//////typedef typename btree_type::Leaf LEAF;
//////typedef typename btree_type::Inner INNER;

////struct SPBPartition
////{

////    btree_type::Node* node;
////    double min, max;

////    SPBPartition()
////    {



////    }

////    SPBPartition(btree_type::Node* _node, ll _min, ll _max)
////    {

////        node = _node;
////        min = _min;
////        max = _max;

////    }

////};

////struct ComparePartition
////{

////    bool operator()(SPBPartition const& a, SPBPartition const& b)
////    {

////        bool ans;

////        if(a.min != b.min)
////        {

////            ans = a.min > b.min;

////        }
////        else
////        {

////            ans = a.max > b.max;

////        }

////        return ans;

////    }

////};

////template <class T>
////struct KnnSPB
////{

////    BasicArrayObject<T> element;
////    double distance;

////    KnnSPB()
////    {



////    }

////    KnnSPB(BasicArrayObject<T> _element, double _distance)
////    {

////        element = _element;
////        distance = _distance;

////    }

////    bool operator<(const KnnSPB<T>& item) const
////    {

////        return distance < item.distance;

////    }

////    bool operator>(const KnnSPB<T>& item) const
////    {

////        return distance > item.distance;

////    }

////};

////template<class Type, class Comp>
////std::vector<Type> dequeueInOrderSPB_Results(std::priority_queue<Type, std::vector<Type>, Comp> pq)
////{

////    std::vector<Type> ans;
////    std::priority_queue<Type, std::vector<Type>, Comp> pqClone = pq;

////    while(!pqClone.empty())
////    {

////        ans.push_back(pqClone.top());
////        pqClone.pop();

////    }

////    return ans;

////}


////template <class type>
////class SPBTree
////{

////public:
////    SPBTree()
////    {

////        setBaseFilePath();
////        pivotsIds = {};
////        pivots = {};
////        df = nullptr;
////        hc = HilbertCurve();

////    }

////    SPBTree(Dataset<type>* dataset, DistanceFunction<BasicArrayObject<type>>* _df, ull pivot_num, double _MAXDIST = -1.0)
////    {

////        auto start = std::chrono::steady_clock::now();

////        setBaseFilePath();
////        pivotsIds = {(size_t)0, (size_t)6};
////        pivots = {dataset->getFeatureVector(0), dataset->getFeatureVector(6)};
////        df = _df;
////        PIVOT_NUM = pivot_num;

////        if(_MAXDIST == -1.0)
////        {

////            findMAXDIST = true;

////        }
////        else
////        {

////            findMAXDIST = false;
////            MAXDIST = _MAXDIST;

////        }

////        bulk_load(dataset);

////        auto end = std::chrono::steady_clock::now();
////        buildTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

////    }

////    void bulk_load(Dataset<type>* dataset)
////    {

////        vector<vector<double>> pointsDist;
////        vector<vector<ull>> pointsDisc;
////        vector<ull> keys;

////        ofstream file(baseFilePath + std::filesystem::path::preferred_separator + "bulk_load_disc.txt");

////        if(findMAXDIST)
////        {

////            MAXDIST = 0.0;

////            for(size_t i = 0; i < dataset->getCardinality(); i++)
////            {

////                pointsDist.push_back(vector<double>(PIVOT_NUM));

////                for(size_t j = 0; j < (size_t)PIVOT_NUM; j++)
////                {

////                    pointsDist[i][j] = df->getDistance(dataset->getFeatureVector(pivotsIds[j]), dataset->getFeatureVector(i));
////                    MAXDIST = std::max(MAXDIST, pointsDist[i][j]);

////                }

////            }

////            EPS = MAXDIST/PREC;
////            MAXINT = MAXDIST/EPS;
////            p = ((ull)log2(MAXINT) + (ull)1);
////            GRID_L = ((1ull << p) - (ull)1);

////            hc = HilbertCurve(p, PIVOT_NUM);

////            for(size_t i = 0; i < dataset->getCardinality(); i++)
////            {

////                pointsDisc.push_back(vector<ull>(PIVOT_NUM));

////                for(size_t j = 0; j < (size_t)PIVOT_NUM; j++)
////                {

////                    pointsDisc[i][j] = (ull)floor(pointsDist[i][j]/EPS);
////                    file << pointsDisc[i][j] << " ";

////                }

////                file << endl;

////            }

////            file.close();

////            keys = hc.distances_from_points(pointsDisc);

////        }
////        else
////        {

////            EPS = MAXDIST/PREC;
////            MAXINT = MAXDIST/EPS;
////            p = ((ull)log2(MAXINT) + (ull)1);
////            GRID_L = ((1ull << p) - (ull)1);

////            hc = HilbertCurve(p, PIVOT_NUM);

////            for(size_t i = 0; i < dataset->getCardinality(); i++)
////            {

////                pointsDisc.push_back(vector<ull>(PIVOT_NUM));

////                for(size_t j = 0; j < (size_t)PIVOT_NUM; j++)
////                {

////                    pointsDisc[i][j] = (ull)floor(df->getDistance(dataset->getFeatureVector(pivotsIds[j]), dataset->getFeatureVector(i))/EPS);
////                    file << pointsDisc[i][j] << " ";

////                }

////                file << endl;

////            }

////            file.close();

////            keys = hc.distances_from_points(pointsDisc);

////        }

////        pointsDist.clear();
////        pointsDisc.clear();

////        ofstream file2(baseFilePath + std::filesystem::path::preferred_separator + "bulk_load_sfc.txt");

////        vector<pair<ull, size_t>> insertValues;

////        for(size_t i = 0; i < keys.size(); i++)
////        {

////            insertValues.push_back(make_pair(keys[i], i));
////            file2 << keys[i] << endl;

////        }

////        sort(insertValues.begin(), insertValues.end());
////        bt.setHilbertCurve(hc);
////        bt.bulk_load(insertValues.begin(), insertValues.end());
////        bt.init_Key_Minmax();
////        bt.initDisk(dataset);

////        keys.clear();
////        insertValues.clear();
////        file2.close();
////        delete dataset;

////    }

////    void dump_key_min_max()
////    {

////        bt.dump_Key_Minmax();

////    }

////    void dumpHilbertCurve()
////    {

////        cout << hc << endl;

////    }

////    bool isInterval(double infBound, double supBound, double test)
////    {

////        return ((test >= infBound) && (test <= supBound));

////    }

////    double minDist(vector<double> sqDisc, double** mbr)
////    {

////        double limInfCase3 = -1.0;
////        double limInfCase2 = -1.0;
////        double answer = -1.0;
////        bool within = true;

////        for(size_t i = 0; i < PIVOT_NUM; i++)
////        {

////            if(!isInterval(mbr[i][0], mbr[i][1], sqDisc[i]))
////            {

////                within = false;

////                limInfCase3 = std::max(limInfCase3,
////                                       std::min(
////                                           std::abs(sqDisc[i] - mbr[i][0]),
////                                           std::abs(sqDisc[i] - mbr[i][1])
////                                           )
////                                       );

////            }
////            else
////            {

////                limInfCase2 = std::min(limInfCase2,
////                                       std::min(
////                                           std::abs(sqDisc[i] - mbr[i][0]),
////                                           std::abs(sqDisc[i] - mbr[i][1])
////                                           )
////                                       );

////            }

////        }

////        if(within)
////        {

////            answer = 0.0;

////        }
////        else
////        {

////            if(limInfCase2 != -1)
////            {

////                answer = limInfCase2;

////            }
////            else
////            {

////                answer = limInfCase3;

////            }

////        }

////        return answer;

////    }

////    double maxDist(vector<double> sqDisc, double** mbr)
////    {

////        double answer = std::numeric_limits<double>::max();

////        for(size_t i = 0; i < PIVOT_NUM; i++)
////        {

////            if(std::numeric_limits<double>::max() - mbr[i][0] >= sqDisc[i])
////            {

////                answer = std::min(answer, sqDisc[i] + mbr[i][0]);

////            }

////            if(std::numeric_limits<double>::max() - mbr[i][1] >= sqDisc[i])
////            {

////                answer = std::min(answer, sqDisc[i] + mbr[i][1]);

////            }


////        }

////        return answer;

////    }

////    size_t getPageID(btree_type::Leaf* node, size_t pos)
////    {

//////        cout << "ACUM = " << node->memory_cost[pos] << endl;

////        if(pos >= 0 && pos <= node->slotuse)
////        {

////            return (size_t)(node->memory_cost[pos] + sizeof(size_t))/PAGE_SIZE;

////        }
////        else
////        {

////            throw std::invalid_argument("Leaf node has fewer elements than the requested position");

////        }

////    }

////    ll getPageAccess()
////    {

////        return pageAccess;

////    }

////    size_t getDistanceCount()
////    {

////        return df->getDistanceCount();

////    }

////    void knn(BasicArrayObject<type> query, size_t k, std::vector<KnnSPB<type>>& ans)
////    {

//////        vector<double> sqDisc = vector<double>(PIVOT_NUM);

//////        for(size_t i = 0; i < PIVOT_NUM; i++)
//////        {

//////            //sqDisc[i] = (ull)std::floor(df->getDistance(query, pivots[i])/EPS);
//////            sqDisc[i] = df->getDistance(query, pivots[i]);
//////            cout << "SQ " << i << " / " << sqDisc[i] << "\n";

//////        }
//////        cout << endl;

//////        btree_type::Inner* innernode = static_cast<btree_type::Inner*>(bt.getRoot());
//////        for(size_t i = 0; i < innernode->slotuse + 1; i++)
//////        {

//////            btree_type::Leaf* leafnode = static_cast<btree_type::Leaf*>(innernode->childid[i]);
//////            cout << "Min = " << minDist(sqDisc, leafnode->mbr) << " / Max = " << maxDist(sqDisc, leafnode->mbr) << endl;

//////        }

//////------------------------------------------------------------------------------------------------------------------------------------

////        df->resetStatistics();
////        pageAccess = 0;
////        ans.clear();

////        std::priority_queue<SPBPartition, std::vector<SPBPartition>, ComparePartition> nodeQueue;
////        std::priority_queue<KnnSPB<type>, std::vector<KnnSPB<type>>, std::greater<KnnSPB<type>>> candidatesQueue;
////        std::priority_queue<KnnSPB<type>, std::vector<KnnSPB<type>>, std::less<KnnSPB<type>>> resultQueue;
////        nodeQueue.push(SPBPartition(bt.getRoot(), 0.0, std::numeric_limits<double>::infinity()));

////        vector<double> sq_ = vector<double>(PIVOT_NUM);

////        for(size_t i = 0; i < PIVOT_NUM; i++)
////        {

////            //sqDisc[i] = (ull)std::floor(df->getDistance(query, pivots[i])/EPS);
////            sq_[i] = df->getDistance(query, pivots[i]);
//////            cout << "SQ " << i << " / " << sqDisc[i] << " ";

////        }
//////        cout << endl;

////        SPBPartition partition;
////        btree_type::Node* node;
////        btree_type::Leaf *leafnode;
////        btree_type::Inner* innernode;
////        std::vector<BasicArrayObject<type>*> dataLeaf;
////        set<size_t> globalPagesID;

////        while(!nodeQueue.empty() || candidatesQueue.size() > 0)
////        {

////            if(candidatesQueue.size() == 0)
////            {

////                partition = nodeQueue.top();
////                node = partition.node;
////                nodeQueue.pop();

////                if(node->isleafnode())
////                {

////                    leafnode = static_cast<btree_type::Leaf*>(node);
////                    pageAccess++;

////                    for(size_t i = 0; i < leafnode->slotuse; i++)
////                    {

////                        size_t pid = getPageID(leafnode, i);
////                        pair<set<size_t>::iterator, bool> insertPID = globalPagesID.insert(pid);

////                        if(insertPID.second)
////                        {

//////                            pageAccess++;
////                            read_from_disk(dataLeaf, pid);

////                        }

////                    }

////                    for(size_t i = 0; i < dataLeaf.size(); i++)
////                    {

////                        candidatesQueue.push(KnnSPB<type>(*dataLeaf[i], df->getDistance(query, *dataLeaf[i])));
////                        delete dataLeaf[i];

////                    }

////                    dataLeaf.clear();

////                }
////                else
////                {

////                     innernode = static_cast<btree_type::Inner*>(node);

////                     for(size_t i = 0; i < innernode->slotuse + 1; i++)
////                     {

////                         nodeQueue.push(SPBPartition(innernode->childid[i], minDist(sq_, innernode->childid[i]->mbr), maxDist(sq_, innernode->childid[i]->mbr)));

////                     }

////                }

////            }
////            else if((nodeQueue.size() > 0) && nodeQueue.top().min < candidatesQueue.top().distance)
////            {

////                partition = nodeQueue.top();
////                node = partition.node;
////                nodeQueue.pop();

////                if(node->isleafnode())
////                {

////                    leafnode = static_cast<btree_type::Leaf*>(node);
////                    pageAccess++;

////                    for(size_t i = 0; i < leafnode->slotuse; i++)
////                    {

////                        size_t pid = getPageID(leafnode, i);
////                        pair<set<size_t>::iterator, bool> insertPID = globalPagesID.insert(pid);

////                        if(insertPID.second)
////                        {

//////                            pageAccess++;
////                            read_from_disk(dataLeaf, pid);

////                        }

////                    }

////                    for(size_t i = 0; i < dataLeaf.size(); i++)
////                    {

////                        candidatesQueue.push(KnnSPB<type>(*dataLeaf[i], df->getDistance(query, *dataLeaf[i])));
////                        delete dataLeaf[i];

////                    }

////                    dataLeaf.clear();

////                }
////                else
////                {

////                     innernode = static_cast<btree_type::Inner*>(node);

////                     for(size_t i = 0; i < innernode->slotuse + 1; i++)
////                     {

////                         nodeQueue.push(SPBPartition(innernode->childid[i], minDist(sq_, innernode->childid[i]->mbr), maxDist(sq_, innernode->childid[i]->mbr)));

////                     }

////                }

////            }
////            else
////            {

////                if(!resultQueue.empty() && !candidatesQueue.empty() && resultQueue.size() >= k && candidatesQueue.top().distance > resultQueue.top().distance)
////                {

////                    break;

////                }

////                resultQueue.push(candidatesQueue.top());
////                candidatesQueue.pop();

////                while(resultQueue.size() > k)
////                {

////                    resultQueue.pop();

////                }

////            }

////        }

////        ans = dequeueInOrderSPB_Results(resultQueue);
////        std::reverse(ans.begin(), ans.end());

////        while(!candidatesQueue.empty())
////        {

////            candidatesQueue.pop();

////        }

////        while(!resultQueue.empty())
////        {

////            resultQueue.pop();

////        }

////        while(!nodeQueue.empty())
////        {

////            nodeQueue.pop();

////        }

////    }

////    void test()
////    {

////        bt.test();

//////        cout << bt.getRoot()->isleafnode() << "\n";

//////        cout << "\n\n";
//////        cout << bt.minleafslots << "\n";
//////        cout << bt.begin()->first << endl;

//////        for(ull t : hc.point_from_distance(39594))
//////            cout << t << " ";
//////        for(ull t : hc.point_from_distance(42840))
//////            cout << t << " ";

////    }

////private:
////    //Pivot<type>* pvt;
////    HilbertCurve hc;
////    DistanceFunction<BasicArrayObject<type>>* df;
////    vector<size_t> pivotsIds;
////    vector<BasicArrayObject<type>> pivots;
////    btree_type bt;
////    bool findMAXDIST;
////    ll buildTime, pageAccess;


////};


#endif // SPB_TREE_H
