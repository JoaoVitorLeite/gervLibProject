#ifndef OMNIEXPERIMENTS_H
#define OMNIEXPERIMENTS_H

#include <IndexExperiments.h>
#include <OmniKdTree.h>

template <class T>
class OmniExperiments : public IndexExperiments<T>
{

private:
    OmniKdTree<T>* index = nullptr;
    std::vector<PairResult> ans;

public:
    OmniExperiments()
    {

    }

    ~OmniExperiments()
    {

    }

    void deleteIndex()
    {

        if(index != nullptr)
            delete index;

    }

    void buildIndex()
    {

        auto start = std::chrono::steady_clock::now();
        index = new OmniKdTree<T>(this->train, this->df, this->pvt, this->numPerLeaf, this->pivot_num);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = index->getDistanceCount();

    }

    void runQuery(BasicArrayObject<T>* query, size_t k)
    {

        ans.clear();
        auto start = std::chrono::steady_clock::now();
        index->kNN(this->train, query, k, ans);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = index->getDistanceCount();
        this->diskAccess = index->getDiskAccess();

    }

    std::string indexName()
    {

        return "OMNIKDTREE";

    }



};

#endif // OMNIEXPERIMENTS_H
