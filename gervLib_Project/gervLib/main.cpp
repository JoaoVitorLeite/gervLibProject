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
#include <typeinfo>
#include <config_spb.h>
#include <MemoryManagerUtils.h>
#include <SPB_Tree.h>
#include <PivotExperiments.h>
#include <experimental/filesystem>
#include <VPExperiments.h>
#include <SPB_Tree.h>

using namespace std;
using namespace mvp;

const int BF = 2;   //branchfactor
const int PL = 8;   // pathlength
const int LC = 55; // leafcap
const int LPN = 2;  // levelspernode
const int FO = 4; //fanout bf^lpn
const int NS = 2; //numsplits (bf-1)^lpn

typedef MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, RandomPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> MVPTREE_DOUBLE_RANDOM;
typedef MVPTree<BasicArrayObject<vector<char>>, EditDistance<BasicArrayObject<vector<char>>>, RandomPivots<vector<char>>, Dataset<vector<char>>, BF,PL,LC,LPN,FO,NS> MVPTREE_STRING_RANDOM;
typedef MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxSeparatedPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> MVPTREE_DOUBLE_MAXSEPARETED;
typedef MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, KmedoidsPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> MVPTREE_DOUBLE_KMEDOIDS;

typedef std::vector<char> str;


int main(int argc, char *argv[])
{

    Dataset<double>* train = new Dataset<double>();
    Dataset<double>::loadNumericDataset(train, "../../gervLib/datasets/Dataset1.csv", " ");
    Pivot<double>* pvt = new MaxVariancePivots<double>();
    //pvt->setSeed(200);
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PM_Tree<double> pm = PM_Tree<double>(train, df, pvt, 5, 2);
    pm.get_root();

//    cout << train->getCardinality() << "\t" << train->getDimensionality() << endl;
//    cout << train->getFeatureVector(0).getSerializedSize() + sizeof(size_t) << endl;

//    pvt->setSeed(200);
//    VpTree<double, DistanceFunction<BasicArrayObject<double>>>* index = new VpTree<double, DistanceFunction<BasicArrayObject<double>>>(false, 0.0, 5, pvt, train, df);
//    Node<BasicArrayObject<double>>* root = index->getRoot();
//    Dataset<double>* ans = new Dataset<double>();
//    BasicArrayObject<double> b = BasicArrayObject<double>(-1, {6.0, 3.0});
//    index->kNNInc(b, 3, root, ans, df);

//    for(size_t i = 0; i < 3; i++)
//        cout << ans->getFeatureVector(i).toStringWithOID() << "\t" << df->getDistance(ans->getFeatureVector(i), b) << endl;

    //    MVPTREE_STRING_RANDOM mvp = MVPTREE_STRING_RANDOM(df, train);
    //    mvp.test();

    //    train->printDataset();

    //    cout << endl << endl;

    //    for(size_t i = 0; i < train->getCardinality(); i++)
    //    {

    //        for(size_t j = 0; j < train->getCardinality(); j++)
    //        {

    //            //cout << "d(" << i << "," << j << ")" << " = " << df->getDistance(*train->getInstance(i), *train->getInstance(j)) << endl;
    //            cout << df->getDistance(*train->getInstance(i), *train->getInstance(j)) << "\t";

    //        }

    //        cout << endl;

    //    }


    //    cout << "Query: " << test->getFeatureVector(0).toStringWithOID() << endl;

    //    std::vector<double> dist;
    //    for(size_t i = 0; i < test->getCardinality(); i++)
    //        dist.push_back(df->getDistance(*test->instance(0), *test->instance(i)));

    //    sort(dist.begin(), dist.end());

    //    index.knn(test->getFeatureVector(0), k, ans);
    //    for(size_t i = 0; i < k; i++)
    //        cout << dist[i] << endl;

    //    for(size_t i = 0; i < k; i++)
    //        cout << ans[i].element.toStringWithOID() << " / " << ans[i].distance << endl;

    //    cout << endl << endl;

    //    index.test();
    //    cout << "CARD = " << train->getCardinality() << endl;


    //*************************************************************************

    //    Dataset<std::vector<char>>* train = new Dataset<std::vector<char>>();
    //    Dataset<std::vector<char>>::loadTextDataset(train, "../../gervLib/datasets/sgb-words.csv", " ");
    //    Dataset<std::vector<char>>* test = new Dataset<std::vector<char>>();
    //    Dataset<std::vector<char>>::loadTextDataset(test, "../../gervLib/datasets/sgb-words.csv", " ");
    //    EditDistance<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    //    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    //    srand(175978);
    //    size_t k = 5, numPerLeaf = 55;
    //    MVPTree<BasicArrayObject<std::vector<char>>, EditDistance<BasicArrayObject<std::vector<char>>>, RandomPivots<std::vector<char>>, Dataset<std::vector<char>>, BF,PL,LC,LPN,FO,NS> index
    //        = MVPTree<BasicArrayObject<std::vector<char>>, EditDistance<BasicArrayObject<std::vector<char>>>, RandomPivots<std::vector<char>>, Dataset<std::vector<char>>, BF,PL,LC,LPN,FO,NS>(df, train);
    //    std::vector<KnnEntryMVP<BasicArrayObject<std::vector<char>>>> ans;

    ////    for(size_t x = 0; x < test->getCardinality(); x++)
    ////    {

    ////        ans.clear();
    ////        index.knn(test->getFeatureVector(x), k, ans);
    ////        std::vector<double> dist;
    ////        for(size_t i = 0; i < test->getCardinality(); i++)
    ////            dist.push_back(df->getDistance(*test->instance(x), *test->instance(i)));

    ////        sort(dist.begin(), dist.end());

    ////        for(size_t z = 0; z < k; z++)
    ////            if(dist[z] != ans[z].distance)
    ////                cout << "ERRO EM: " << x << endl;

    ////    }

    //    cout << "Query: " << test->getFeatureVector(0).toStringWithOID() << endl;

    //    std::vector<double> dist;
    //    for(size_t i = 0; i < test->getCardinality(); i++)
    //        dist.push_back(df->getDistance(*test->instance(0), *test->instance(i)));

    //    sort(dist.begin(), dist.end());

    //    index.knn(test->getFeatureVector(0), k, ans);
    //    for(size_t i = 0; i < k; i++)
    //        cout << dist[i] << endl;

    //    for(size_t i = 0; i < k; i++)
    //        cout << ans[i].element.toStringWithOID() << " / " << ans[i].distance << endl;

    //    cout << endl << endl;

    //    index.test();
    //    cout << "CARD = " << train->getCardinality() << endl;



    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //    Dataset<double>* data = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(data, "../datasets/open2/Dataset3.csv", ",");
    //    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    Pivot<double>* pvt = new PCAPivots<double>();
    //    //pvt->setSampleSize(0.5);
    //    pvt->generatePivots(data, df, 7);
    //    for(size_t i = 0; i < 7; i++)
    //        std::cout << pvt->get(i).getOID() << std::endl;

    //    PivotExperiments<double> expt = PivotExperiments<double>();
    //    expt.setDistanceFunctionName("EUCLIDEAN");
    //    expt.setTrainDataset(data);
    //    expt.setOutputPath("../results/");
    //    expt.setDistanceFunction(df);
    //    expt.setPivotMethod(pvt);
    //    expt.setSeed(157);
    //    expt.setSampleSize(0.5);
    //    expt.modifySeed();
    //    expt.modifySampleSize();
    //    expt.setPivotNum({100,150,200});
    //    expt.setSavePivot(false);
    //    expt.runExperiment();
    //    expt.setSavePivot(true);
    //    expt.runExperimentWithRepetitions(1);

    //-----------------------------------------------------------------------------------------------------------------------------------------

    //    Dataset<double>* data = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(data, "../datasets/Dataset1.csv", " ");
    //    Dataset<double>* test = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(test, "../datasets/Dataset1.csv", " ");
    //    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    Pivot<double>* pvt = new BPPPivots<double>();
    //    pvt->setSeed(4554);
    //    pvt->generatePivots(data, df, 6);

    //    for(size_t i = 0; i < 6; i++)
    //    {

    //        cout << pvt->get(i).toStringWithOID() << endl;

    //    }

    //    SPBTree<double> spb = SPBTree<double>(data, df, pvt, 2, 4);
    ////    spb.dump_key_min_max();
    ////    spb.test();

    //    BasicArrayObject<double> query = BasicArrayObject<double>(-1,2);
    //    query.set(0, 6.0);
    //    query.set(1, 3.0);

    //    std::vector<KnnSPB<double>> ans;
    //    spb.knn(query, 3, ans);

    //    for(auto &e : ans)
    //        cout << e.element.getOID() << " / " << e.distance << endl;

    //    cout << ans.size() << endl;
    //    cout << spb.getLeafNodeAccess() << endl;
    //    cout << spb.getDistanceCount() << endl;
    //    cout << IOread << endl;
    //    cout << p << endl;

    //----------------------------------------------------------------------------------------------------------------------------------------

    //    Dataset<str>* data = new Dataset<str>();
    //    Dataset<str>::loadTextDataset(data, "../datasets/names.csv", " ");
    //    Dataset<str>* test = new Dataset<str>();
    //    Dataset<str>::loadTextDataset(test, "../datasets/names.csv", " ");
    //    DistanceFunction<BasicArrayObject<str>>* df = new EditDistance<BasicArrayObject<str>>();
    //    Pivot<str>* pvt = new MaxVariancePivots<str>();
    //    SPBTree<str> spb = SPBTree<str>(data, df, pvt, 3, 3);
    //    spb.dump_key_min_max();

    //    Dataset<double>* data = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(data, "../datasets/cities_norm.csv", ",");
    //    Dataset<double>* test = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(test, "../datasets/cities_norm.csv", ",");
    //    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    Pivot<double>* pvt = new MaxVariancePivots<double>();
    //    SPBTree<double> spb = SPBTree<double>(data, df, pvt, 3, 200);
    //    //spb.dump_key_min_max();

    //    size_t leafMean = 0, distCnt = 0;
    //    size_t num_queries = test->getCardinality();

    //    for(size_t j = 0; j < num_queries; j++)
    //    {

    //        BasicArrayObject<double> query = test->getFeatureVector(j);
    //        //BasicArrayObject<str> query = test->getFeatureVector(j);

    //        size_t k = 100;

    //        std::vector<KnnSPB<double>> ans;
    //        //std::vector<KnnSPB<str>> ans;

    //        spb.knn(query, k, ans);

    //        //        cout << ans.size() << endl;
    //        //        cout << spb.getLeafNodeAccess() << endl;
    //        //        cout << spb.getDistanceCount() << endl;
    //        leafMean += spb.getLeafNodeAccess();
    //        distCnt += spb.getDistanceCount();

    //        //    for(auto i : ans)
    //        //        cout << i.element.toStringWithOID() << " / " << i.distance << endl;

    //        vector<pair<size_t, double>> dist;results/

    //        for(size_t i = 0; i < test->getCardinality(); i++)
    //        {

    //            dist.push_back(std::make_pair(i, df->getDistance(query, test->getFeatureVector(i))));

    //        }

    //        std::sort(dist.begin(), dist.end(), [](const std::pair<size_t, double>& lhs, const std::pair<size_t, double>& rhs){
    //            return lhs.second < rhs.second;
    //        });

    //        //    cout << "\n\n";

    //        //    for(size_t i = 0; i < k; i++)
    //        //        cout << dist[i].first << " / " << dist[i].second << endl;

    //        for(size_t z = 0; z < k; z++)
    //        {

    //            if(ans[z].distance != dist[z].second)
    //                cout << "ERRO EM: " << j << endl;

    //        }

    //    }

    //    cout << "LEAF MEAN : " << (leafMean*1.0)/num_queries << endl;
    //    cout << "DIST MEAN : " << (distCnt*1.0)/num_queries << endl;
    //    cout << "LEAF 2 : " << IOread/num_queries << endl;

    //-------------------------------------------------------------------------------------------------------

    //    Dataset<double>* train = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(train, "../datasets/Dataset1.csv", " ");
    //    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    //    //pvt->generatePivots(train, df, 2);

    ////    MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxVariancePivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> index
    ////            = MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxVariancePivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS>((EuclideanDistance<BasicArrayObject<double>>*)df, train);

    ////    OmniKdTree<double> index = OmniKdTree<double>(train, df, pvt, 5);
    //    PM_Tree<double> index = PM_Tree<double>(train, df, pvt, 5, 2);

    //    BasicArrayObject<double> query = BasicArrayObject<double>(-1,2);
    //    query.set(0, 6.0);
    //    query.set(1, 3.0);

    //    std::vector<KnnEntry<double>> ans;
    //    index.kNN(query, 3, ans);

    ////    std::vector<PairResult> ans;
    ////    index.kNN(train, &query, 3, ans);

    ////    std::vector<KnnEntryMVP<BasicArrayObject<double>>> ans;

    ////    index.knn(query, 3, ans);

    //    for(auto &i : ans)
    //        cout << i.element.toStringWithOID() << " / " << i.distance << endl;

    //    cout << endl << index.getLeafNodeAccess() << endl;


    //-------------------------------------------------------------------------------------------------------

    //    std::string s = std::bitset< 10 >( 2 ).to_string();
    //    cout << s << endl;


    //    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    WDRPivots<double> pvt1 = WDRPivots<double>();
    //    pvt1.setSeed(94819);
    ////    pvt1.setAlpha(0.1);

    //    //pvt1.setSampleSize(0.5);
    //    pvt1.generatePivots(train, df, 6);

    //    for(size_t i = 0; i < 6; i++)
    //        cout << pvt1.getPivot(i)->getOID() << "\n";

    //    Dataset<double>* train = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(train, "../datasets/Dataset1.csv", " ");
    //    Dataset<double>* test = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(test, "../datasets/Dataset1.csv", " ");
    //    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    //    pvt->generatePivots(train, df, 3);

    //    for(size_t x = 0; x < 3; x++)
    //        cout << pvt->getPivot(x)->toString() << endl;

    //    pvt->setPath("../results/maxvar.pvt");
    //    pvt->writePivotsToFile();


    ////    VpTree<double, DistanceFunction<BasicArrayObject<double>>> index = VpTree<double, DistanceFunction<BasicArrayObject<double>>>(false, 0.0, 5, pvt, train, df);

    //    MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxVariancePivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> index
    //            = MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxVariancePivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS>((EuclideanDistance<BasicArrayObject<double>>*)df, train);

    //    BasicArrayObject<double> query = BasicArrayObject<double>(-1,2);
    //    query.set(0, 6.0);
    //    query.set(1, 3.0);

    ////    Dataset<double>* ans3 = new Dataset<double>();
    ////    index.kNNInc(query, 3, index.getRoot(), ans3, df);

    //    std::vector<KnnEntryMVP<BasicArrayObject<double>>> ans;
    //    index.knn(query, 3, ans);
    //    cout << index.getLeafNodeAccess() << endl;

    //    for(size_t x = 0; x < ans.size(); x++)
    //        cout << ans[x].element.toStringWithOID() << endl;

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
    ////    //VpTree<double, DistanceFunction<BasicArrayObject<double>>> index = VpTree<double, DistanceFunction<BasicArrayObject<double>>>(false, 0.0, numPerLeaf, pvt, train, df);
    //    PM_Tree<double> index = PM_Tree<double>(train, df, pvt, numPerLeaf, 2);

    //    std::vector<std::pair<double, size_t>> ansVec;
    //    for(size_t x = 0; x < test->getCardinality(); x++)
    //    {

    //        ans2.clear();
    //        index.kNN(test->getFeatureVector(x), 100, ans2);

    //        std::vector<double> dist;
    //        for(size_t i = 0; i < train->getCardinality(); i++)
    //            dist.push_back(df->getDistance(*test->instance(x), *train->instance(i)));

    //        sort(dist.begin(), dist.end());

    //        for(size_t z = 0; z < 100; z++)
    //            if(dist[z] != ans2[z].distance)
    //                cout << "ERRO EM : " << x << endl;

    //    }

    //    std::ofstream file("/home/joaovictor/Documents/TCC/Code/cities_pmtree_range.csv");
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
    //            //index.kNN(test->getFeatureVector(i), k, ans2);
    //            index.rang
    //            auto end = std::chrono::steady_clock::now();
    //            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    //            file << std::to_string(k) << "," << std::to_string(elapsed.count()) << "," << std::to_string(df->getDistanceCount()) << "," << index.getLeafNodeAccess() << "\n";

    //        }

    //    }

    //    file.close();

    //    Dataset<double>* train = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(train, "../datasets/train_nasa.csv", ",");
    //    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    WDRPivots<double>* pvt = new WDRPivots<double>();
    //    pvt->setSampleSize(0.2);
    //    pvt->generatePivots(train, df, 2);

    ////    OmniKdTree<double> index = OmniKdTree<double>(train, df, pvt, 360);
    //    VpTree<double, DistanceFunction<BasicArrayObject<double>>> index = VpTree<double, DistanceFunction<BasicArrayObject<double>>>(false, 0.0, 360, pvt, train, df);

    //    cout << "OI\n";

    ////    Dataset<double>* sample = train->sampleDataset(std::ceil(train->getCardinality()*0.2), false, 0);

    ////    for(size_t i = 0; i < train->getCardinality(); i++)
    ////        vec.push_back({static_cast<long long>(i), train->getFeatureVector(i)});

    ////    tree.Add(vec);
    ///




    //**********************************************************************************************************************************

    //    Dataset<double>* train = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(train, "../datasets/cities_norm.csv", ",");
    //    Dataset<double>* test = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(test, "../datasets/cities_norm.csv", ",");
    //    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    RandomPivots<double>* pvt = new RandomPivots<double>();
    //    pvt->generatePivots(train, df, 2);
    //    size_t numPerLeaf = 55;
    //    std::vector<PairResult> ans;
    //    std::vector<KnnEntry<double>> ans2;
    //    Dataset<double>* ans3 = new Dataset<double>();

    //    OmniKdTree<double> index = OmniKdTree<double>(train, df, pvt, numPerLeaf);

    //    std::vector<double> dist;

    //    for(size_t i = 0; i < test->getCardinality(); i++)
    //    {

    //        BasicArrayObject<double> query = test->getFeatureVector(i);

    //        ans.clear();
    //        ans2.clear();
    //        ans3->clear();

    //        index.kNN(train, test->instance(i), 100, ans);
    //        train->getDistanceVector(query, df, &dist);
    //        std::sort(dist.begin(), dist.end());

    //        for(size_t z = 0; z < 100; z++)
    //            if(ans[z].distance != dist[z])
    //                cout << "ERRO EM : " << i << " " << ans[z].distance << " " << dist[z] << endl;

    //    }

    //**********************************************************************************************************************************

    //    Dataset<vector<char>>* train = new Dataset<vector<char>>();
    //    Dataset<vector<char>>::loadTextDataset(train, "../datasets/names.csv", " ");
    //    Dataset<vector<char>>* train2 = new Dataset<vector<char>>();
    //    *train2 = *train;
    //    Dataset<vector<char>>* test = new Dataset<vector<char>>();
    //    Dataset<vector<char>>::loadTextDataset(test, "../datasets/names.csv", " ");
    //    EditDistance<BasicArrayObject<vector<char>>>* df = new EditDistance<BasicArrayObject<vector<char>>>();
    //    RandomPivots<vector<char>>* pvt = new RandomPivots<vector<char>>();
    //    pvt->generatePivots(train, df, 2);
    //    size_t numPerLeaf = 55;
    //    std::vector<PairResult> ans;
    //    std::vector<KnnEntry<vector<char>>> ans2;
    //    Dataset<vector<char>>* ans3 = new Dataset<vector<char>>();


    //    //OmniKdTree<vector<char>> index = OmniKdTree<vector<char>>(train, df, pvt, numPerLeaf);
    //    //VpTree<vector<char>, DistanceFunction<BasicArrayObject<vector<char>>>> index = VpTree<vector<char>, DistanceFunction<BasicArrayObject<vector<char>>>>(false, 0.0, numPerLeaf, pvt, train, df);
    //    PM_Tree<vector<char>> index = PM_Tree<vector<char>>(train, df, pvt, numPerLeaf, 2);


    //    std::vector<double> dist;

    //    for(size_t i = 0; i < test->getCardinality(); i++)
    //    {

    //        BasicArrayObject<vector<char>> query = test->getFeatureVector(i);

    //        ans.clear();
    //        ans2.clear();
    //        ans3->clear();

    //        dist.clear();

    ////        index.kNN(train, test->instance(i), 100, ans);
    ////        index.kNNInc(test->getFeatureVector(i), 100, index.getRoot(), ans3, df);
    //        index.kNN(test->getFeatureVector(i), 100, ans2);

    //        train2->getDistanceVector(query, df, &dist);
    //        std::sort(dist.begin(), dist.end());

    //        for(size_t z = 0; z < 100; z++)
    //            if(ans2[z].distance != dist[z])
    //                cout << "ERRO EM : " << i << endl;

    //    }

    /////////////////////////////////////////////MVP

    //    Dataset<double>* train = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(train, "../datasets/train_nasa.csv", ",");
    //    Dataset<double>* test = new Dataset<double>();
    //    Dataset<double>::loadNumericDataset(test, "../datasets/test_nasa.csv", ",");
    //    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    //    PCAPivots<double> pvt1 = PCAPivots<double>();
    //    pvt1.setSampleSize(0.01);
    //    pvt1.generatePivots(train, df, 7);

    //    for(size_t i = 0; i < 7; i++)
    //        cout << pvt1.getPivot(i)->toStringWithOID() << "\n";


    //    PCAPivots<double>* pvt = new PCAPivots<double>();
    //    pvt->setSampleSize(0.01);
    //    pvt->generatePivots(train, df, 2);

    //    MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, PCAPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> mvp
    //            = MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, PCAPivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS>((EuclideanDistance<BasicArrayObject<double>>*)df);
    //    std::vector<datapoint_t<BasicArrayObject<double>, PL>> addPoints;
    //    for(size_t x = 0; x < train->getCardinality(); x++)
    //    {

    //        addPoints.push_back({static_cast<long long>(x), train->getFeatureVector(x)});

    //    }
    //    mvp.Add(addPoints);
    //    size_t k = 100;

    //    vector<size_t> distcnt, leafcnt;

    //    for(size_t j = 0; j < test->getCardinality(); j++)
    //    {

    //        std::vector<std::pair<size_t, double>> dist;
    //        for(size_t i = 0; i < train->getCardinality(); i++)
    //            dist.push_back(std::make_pair(i, df->getDistance(*test->instance(j), *train->instance(i))));
    //        std::sort(dist.begin(), dist.end(), [](std::pair<size_t, double> &a, std::pair<size_t, double> &b){
    //            return a.second < b.second;
    //        });

    //        std::vector<KnnEntryMVP<BasicArrayObject<double>>> ans;
    //        df->resetStatistics();
    //        mvp.knn(test->getFeatureVector(j), k, ans);

    //        for(size_t z = 0; z < k; z++)
    //        {

    //            if(ans[z].distance != dist[z].second)
    //                cout << "ERRO EM : " << j << endl;

    //        }

    //        distcnt.push_back(df->getDistanceCount());
    //        leafcnt.push_back(mvp.getLeafNodeAccess());

    //    }

    //    if(distcnt.size() % 2 == 0)
    //    {

    //        nth_element(distcnt.begin(),
    //                            distcnt.begin() + distcnt.size() / 2,
    //                            distcnt.end());

    //        nth_element(distcnt.begin(),
    //                            distcnt.begin() + (distcnt.size() - 1) / 2,
    //                            distcnt.end());


    //        cout << "DIST MEDIAN : " << (double)(distcnt[(distcnt.size() - 1) / 2]+ distcnt[distcnt.size() / 2])/ 2.0 << endl;

    //    }
    //    else
    //    {

    //        nth_element(distcnt.begin(),
    //                            distcnt.begin() + distcnt.size() / 2,
    //                            distcnt.end());

    //        cout << "DIST MEDIAN : " << (double)distcnt[distcnt.size() / 2] << endl;

    //    }

    //    if(leafcnt.size() % 2 == 0)
    //    {

    //        nth_element(leafcnt.begin(),
    //                            leafcnt.begin() + leafcnt.size() / 2,
    //                            leafcnt.end());

    //        nth_element(leafcnt.begin(),
    //                            leafcnt.begin() + (leafcnt.size() - 1) / 2,
    //                            leafcnt.end());


    //        cout << "LEAF MEDIAN : " << (double)(leafcnt[(leafcnt.size() - 1) / 2]+ leafcnt[leafcnt.size() / 2])/ 2.0 << endl;

    //    }
    //    else
    //    {

    //        nth_element(leafcnt.begin(),
    //                            leafcnt.begin() + leafcnt.size() / 2,
    //                            leafcnt.end());

    //        cout << "LEAF MEDIAN : " << (double)leafcnt[leafcnt.size() / 2] << endl;

    //    }

    //    Dataset<vector<char>>* train = new Dataset<vector<char>>();
    //    Dataset<vector<char>>::loadTextDataset(train, "../datasets/names.csv", " ");
    //    Dataset<vector<char>>* test = new Dataset<vector<char>>();
    //    Dataset<vector<char>>::loadTextDataset(test, "../datasets/names.csv", " ");
    //    EditDistance<BasicArrayObject<vector<char>>>* df = new EditDistance<BasicArrayObject<vector<char>>>();
    //    RandomPivots<vector<char>>* pvt = new RandomPivots<vector<char>>();
    //    pvt->generatePivots(train, df, 2);

    //    MVPTREE_STRING_RANDOM mvp = MVPTREE_STRING_RANDOM(df, pvt);
    //    std::vector<datapoint_t<BasicArrayObject<vector<char>>, PL>> addPoints;
    //    for(size_t x = 0; x < train->getCardinality(); x++)
    //    {

    //        addPoints.push_back({static_cast<long long>(x), train->getFeatureVector(x)});

    //    }
    //    mvp.Add(addPoints);
    //    size_t k = 100;

    //    vector<size_t> distcnt, leafcnt;

    //    for(size_t j = 0; j < test->getCardinality(); j++)
    //    {

    //        std::vector<std::pair<size_t, double>> dist;
    //        for(size_t i = 0; i < train->getCardinality(); i++)
    //            dist.push_back(std::make_pair(i, df->getDistance(*test->instance(j), *train->instance(i))));
    //        std::sort(dist.begin(), dist.end(), [](std::pair<size_t, double> &a, std::pair<size_t, double> &b){
    //            return a.second < b.second;
    //        });

    //        std::vector<KnnEntryMVP<BasicArrayObject<vector<char>>>> ans;
    //        df->resetStatistics();
    //        mvp.knn(test->getFeatureVector(j), k, ans);

    //        for(size_t z = 0; z < k; z++)
    //        {

    //            if(ans[z].distance != dist[z].second)
    //                cout << "ERRO EM : " << j << endl;

    //        }

    //        distcnt.push_back(df->getDistanceCount());
    //        leafcnt.push_back(mvp.getLeafNodeAccess());

    //    }

    //    if(distcnt.size() % 2 == 0)
    //    {

    //        nth_element(distcnt.begin(),
    //                            distcnt.begin() + distcnt.size() / 2,
    //                            distcnt.end());

    //        nth_element(distcnt.begin(),
    //                            distcnt.begin() + (distcnt.size() - 1) / 2,
    //                            distcnt.end());


    //        cout << "DIST MEDIAN : " << (double)(distcnt[(distcnt.size() - 1) / 2]+ distcnt[distcnt.size() / 2])/ 2.0 << endl;

    //    }
    //    else
    //    {

    //        nth_element(distcnt.begin(),
    //                            distcnt.begin() + distcnt.size() / 2,
    //                            distcnt.end());

    //        cout << "DIST MEDIAN : " << (double)distcnt[distcnt.size() / 2] << endl;

    //    }

    //    if(leafcnt.size() % 2 == 0)
    //    {

    //        nth_element(leafcnt.begin(),
    //                            leafcnt.begin() + leafcnt.size() / 2,
    //                            leafcnt.end());

    //        nth_element(leafcnt.begin(),
    //                            leafcnt.begin() + (leafcnt.size() - 1) / 2,
    //                            leafcnt.end());


    //        cout << "LEAF MEDIAN : " << (double)(leafcnt[(leafcnt.size() - 1) / 2]+ leafcnt[leafcnt.size() / 2])/ 2.0 << endl;

    //    }
    //    else
    //    {

    //        nth_element(leafcnt.begin(),
    //                            leafcnt.begin() + leafcnt.size() / 2,
    //                            leafcnt.end());

    //        cout << "LEAF MEDIAN : " << (double)leafcnt[leafcnt.size() / 2] << endl;

    //    }


    return 0;

}
