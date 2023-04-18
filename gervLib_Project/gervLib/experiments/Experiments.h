#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include <Hermes.h>
#include <Dataset.h>
#include <Pivots.h>
#include <fstream>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

template <class T>
class Experiments
{

protected:
    Dataset<T> *train, *test;
    DistanceFunction<BasicArrayObject<T>>* df;
    Pivot<T>* pvt;
    std::vector<size_t> kRange;
    std::string outputPath, dfName;
    size_t numPerLeaf, pivot_num, num_query, seed;
    size_t time, calcDist, diskAccess;
    double sample_size;
    const std::vector<std::string> names = {"RANDOM",
                                            "GNAT",
                                            "CONVEX",
                                            "KMEDOIDS",
                                            "MAXSEPARETED",
                                            "MAXVARIANCE",
                                            "SELECTION",
                                            "PCA",
                                            "SSS",
                                            "FFT",
                                            "HFI",
                                            "IS",
                                            "WDR",
                                            "BPP"};

public:
    Experiments()
    {

    }

    virtual ~Experiments()
    {

    }

    unsigned long long getTime()
    {

        return time;

    }

    size_t getCalcDist()
    {

        return calcDist;

    }

    size_t getDiskAccess()
    {

        return diskAccess;

    }

    virtual void buildIndex() = 0;
    virtual void runQuery(BasicArrayObject<T>* query, size_t k) = 0;
    virtual std::string indexName() = 0;
    virtual void runExperiment() = 0;
    virtual void runExperimentWithRepetitions(size_t rep) = 0;

    void setKRange(size_t min, size_t max, size_t step)
    {

        kRange.resize(3);
        kRange[0] = min;
        kRange[1] = max;
        kRange[2] = step;

    }

    void setOutputPath(std::string output_path)
    {

        outputPath = output_path;

    }

    void setTrainDataset(Dataset<T>* train_)
    {

        train = train_;

    }

    void setNumQuery(size_t num)
    {

        num_query = num;

    }

    void setTestDataset(Dataset<T>* test_)
    {

        test = test_;

    }

    void setDistanceFunctionName(std::string df_name)
    {

        dfName = df_name;

    }

    void setSeed(size_t seed_)
    {

        seed = seed_;

    }

    void setSampleSize(double smp_size)
    {

        sample_size = smp_size;

    }

    void setDistanceFunction(DistanceFunction<BasicArrayObject<T>>* df_)
    {

        df = df_;

    }

    void setPivotMethod(Pivot<T>* pvt_)
    {

        pvt = pvt_;

    }

    void modifySeed()
    {

        pvt->setSeed(seed);

    }

    void modifySampleSize()
    {

        pvt->setSampleSize(sample_size);

    }

    void setNumPerLeaf(size_t num)
    {

        numPerLeaf = num;

    }

};


#endif // EXPERIMENTS_H
