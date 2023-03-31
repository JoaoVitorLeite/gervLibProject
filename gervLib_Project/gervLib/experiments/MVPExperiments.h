#ifndef MVPEXPERIMENTS_H
#define MVPEXPERIMENTS_H

#include <IndexExperiments.h>
#include <mvptree.h>

using namespace mvp;

template <class DType, typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
class MVPExperiments : public IndexExperiments<DType>
{

private:
    MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* index;
    std::vector<KnnEntryMVP<T>> ans;

public:
    MVPExperiments()
    {

    }

    ~MVPExperiments()
    {

    }

    void buildIndex()
    {

        auto start = std::chrono::steady_clock::now();
        index = new MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>((F*)this->df, this->train);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = this->df->getDistanceCount();

    }

    void runQuery(BasicArrayObject<T>* query, size_t k)
    {

        ans.clear();
        auto start = std::chrono::steady_clock::now();
        index->knn(*query, k, ans);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = this->df->getDistanceCount();
        this->diskAccess = index->getLeafNodeAccess();

    }

    std::string indexName()
    {

        return "MVPTREE";

    }

};

#endif // MVPEXPERIMENTS_H
