#ifndef MVPEXPERIMENTS_H
#define MVPEXPERIMENTS_H

#include <IndexExperiments.h>
#include <MVPTree.h>

template <class T>
class MVPExperiments : public IndexExperiments<T>
{

private:
    MVPTree<T>* index = nullptr;
    std::vector<KnnEntryMVP<T>> ans;

public:
    MVPExperiments()
    {

    }

    ~MVPExperiments()
    {

    }

//    void deleteIndex()
//    {

//        if(index != nullptr)
//            delete index;

//    }

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
        index = new MVPTree<T>(this->train, this->df, this->pvt, 2, 8, this->numPerLeaf, 2, 4, 2);
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
