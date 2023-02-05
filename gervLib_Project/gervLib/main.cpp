#include "mvptree.h"
#include <iostream>
#include <datapoint.h>
#include <Hermes.h>
#include <Pivots.h>
#include <Dataset.h>
#include <OmniKdTree.h>
#include <VpTree.h>
#include <PM_Tree.h>
#include <chrono>
#include <unistd.h>
#include <filesystem>

using namespace std;
using namespace mvp;

const int BF = 2;   //branchfactor
const int PL = 8;   // pathlength
const int LC = 50; // leafcap
const int LPN = 2;  // levelspernode
const int FO = 4; //fanout bf^lpn
const int NS = 2; //numsplits (bf-1)^lpn

typedef MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, RandomPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> MVPTREE_DOUBLE_RANDOM;
typedef MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxSeparetedPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> MVPTREE_DOUBLE_MAXSEPARETED;
typedef MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, KmedoidsPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> MVPTREE_DOUBLE_KMEDOIDS;

int main(int argc, char *argv[])
{

//    Dataset<double>* train = new Dataset<double>();
//    Dataset<double>::loadNumericDataset(train, "/home/joaovictor/Documents/TCC/Code/Project/gervLib/datasets/cities_norm.csv", ",");
//    Dataset<double>* test = new Dataset<double>();
//    Dataset<double>::loadNumericDataset(test, "/home/joaovictor/Documents/TCC/Code/Project/gervLib/datasets/cities_norm.csv", ",");
//    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
//    RandomPivots<double>* pvt = new RandomPivots<double>();
//    pvt->generatePivots(train, df, 2);
//    size_t numPerLeaf = 55;
//    std::vector<PairResult> ans;
//    std::vector<KnnEntry<double>> ans2;
//    Dataset<double>* ans3 = new Dataset<double>();

////    OmniKdTree<double> index = OmniKdTree<double>(train, df, pvt, numPerLeaf);
//    //VpTree<double, DistanceFunction<BasicArrayObject<double>>> index = VpTree<double, DistanceFunction<BasicArrayObject<double>>>(false, 0.0, numPerLeaf, pvt, train, df);
//    PM_Tree<double> index = PM_Tree<double>(train, df, pvt, numPerLeaf, 2);

//    std::ofstream file("/home/joaovictor/Documents/TCC/Code/cities_pmtree.csv");
//    file << "k,Time,Count,Leaf\n";

//    for(size_t k = 5; k <= 100; k+=5)
//    {

//        for(size_t i = 0; i < test->getCardinality(); i++)
//        {

//            df->ResetStatistics();
//            //ans.erase(ans.begin(), ans.end());
//            //ans.clear();
//            ans2.erase(ans2.begin(), ans2.end());
//            auto start = std::chrono::steady_clock::now();
////            index.kNN(train, test->instance(i), k, ans);
////            index.kNNInc(test->getFeatureVector(i), k, index.getRoot(), ans3, df);
//            index.kNN(test->getFeatureVector(i), k, ans2);
//            auto end = std::chrono::steady_clock::now();
//            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

//            file << std::to_string(k) << "," << std::to_string(elapsed.count()) << "," << std::to_string(df->getDistanceCount()) << "," << index.getLeafNodeAccess() << "\n";

//        }

//    }

//    file.close();

    Dataset<double>* train = new Dataset<double>();
    Dataset<double>::loadNumericDataset(train, "../../gervLib/datasets/train_nasa.csv", ",");
    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSampleSize(0.2);
    pvt->generatePivots(train, df, 2);

    OmniKdTree<double> index = OmniKdTree<double>(train, df, pvt, 360);

//    Dataset<double>* sample = train->sampleDataset(std::ceil(train->getCardinality()*0.2), false, 0);

//    for(size_t i = 0; i < train->getCardinality(); i++)
//        vec.push_back({static_cast<long long>(i), train->getFeatureVector(i)});

//    tree.Add(vec);

    return 0;

}
