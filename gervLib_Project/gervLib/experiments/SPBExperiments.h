#ifndef SPBEXPERIMENTS_H
#define SPBEXPERIMENTS_H

#include <IndexExperiments.h>
#include <SPB_Tree.h>

template <class T>
class SPBExperiments : public IndexExperiments<T>
{

private:
    SPBTree<T>* index;
    unsigned long long num_bins;
    std::vector<KnnSPB<T>> ans;

public:
    SPBExperiments()
    {

    }

    ~SPBExperiments()
    {


    }

    void setNumBins(unsigned long long num)
    {

        num_bins = num;

    }

    void buildIndex()
    {

        auto start = std::chrono::steady_clock::now();
        index = new SPBTree<T>(this->train, this->df, this->pvt, this->pivot_num, this->num_bins);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = index->getDistanceCount();

    }

    void runQuery(BasicArrayObject<T>* query, size_t k)
    {

        ans.clear();
        auto start = std::chrono::steady_clock::now();
        index->knn(query, k, ans);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = index->getDistanceCount();
        this->diskAccess = IOread;

    }

    std::string indexName()
    {

        return "SPBTREE";

    }

};

#endif // SPBEXPERIMENTS_H
