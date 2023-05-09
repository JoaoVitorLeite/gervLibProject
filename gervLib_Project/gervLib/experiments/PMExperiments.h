#ifndef PMEXPERIMENTS_H
#define PMEXPERIMENTS_H

#include <IndexExperiments.h>
#include <PM_Tree.h>

template <class T>
class PMExperiments : public IndexExperiments<T>
{

private:
    PM_Tree<T>* index = nullptr;
    std::vector<KnnEntry<T>> ans;

public:
    PMExperiments()
    {

    }

    ~PMExperiments()
    {

    }

    void clean()
    {

        if(index != nullptr)
            delete index;
        buildIndex();
        this->saveBuildStats();

    }

    void buildIndex()
    {

        auto start = std::chrono::steady_clock::now();
        index = new PM_Tree<T>(this->train, this->df, this->pvt, this->numPerLeaf, this->pivot_num);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = this->df->getDistanceCount();

    }

    void runQuery(BasicArrayObject<T>* query, size_t k)
    {

        ans.clear();
        auto start = std::chrono::steady_clock::now();
        index->kNN(*query, k, ans);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = this->df->getDistanceCount();
        this->diskAccess = index->getLeafNodeAccess();

    }

    std::string indexName()
    {

        return "PMTREE";

    }

};

#endif // PMEXPERIMENTS_H
