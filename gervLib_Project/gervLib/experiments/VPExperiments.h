#ifndef VPEXPERIMENTS_H
#define VPEXPERIMENTS_H

#include <IndexExperiments.h>
#include <VpTree.h>

template <class T>
class VPExperiments : public IndexExperiments<T>
{

private:
    VpTree<T, DistanceFunction<BasicArrayObject<T>>>* index = nullptr;
    Dataset<T>* ans = new Dataset<T>();

public:
    VPExperiments()
    {

    }

    ~VPExperiments()
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
        index = new VpTree<T, DistanceFunction<BasicArrayObject<T>>>(false, 0.0, this->numPerLeaf, this->pvt, new Dataset<T>(this->train->getElements(), this->train->getCardinality(), this->train->getDimensionality()), this->df);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = this->df->getDistanceCount();

    }

    void runQuery(BasicArrayObject<T>* query, size_t k)
    {

        ans->clear();
        auto start = std::chrono::steady_clock::now();
        index->kNNInc(*query, k, index->getRoot(), ans, this->df);
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        this->time = elapsed.count();
        this->calcDist = this->df->getDistanceCount();
        this->diskAccess = index->getLeafNodeAccess();

    }

    std::string indexName()
    {

        return "VPTREE";

    }

};

#endif // VPEXPERIMENTS_H
