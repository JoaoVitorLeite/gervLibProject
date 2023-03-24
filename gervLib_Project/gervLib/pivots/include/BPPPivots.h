#ifndef BPPPIVOTS_H
#define BPPPIVOTS_H

#include <Pivot.h>

template <class DType>
class BPPPivots : public Pivot<DType>
{

private:
    size_t *flg, **rnk;
    double** rec;

    size_t tot(size_t pn);
    double cnt(size_t p, size_t pvtSize, size_t candSize);

public:
    BPPPivots();
    BPPPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots) : Pivot<DType>(){

        this->setPivotType(PIVOT_TYPE::BPP);
        generatePivots(sample, df, nPivots);

    }
    BPPPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t seed) : Pivot<DType>(){

        this->setPivotType(PIVOT_TYPE::BPP);
        this->setSeed(seed);
        generatePivots(sample, df, nPivots);

    }

    ~BPPPivots();

    void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
    void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());


};

template <class DType>
size_t BPPPivots<DType>::tot(size_t pn)
{

    size_t s = 0;
    for(size_t i = 0; i < pn; i++)
    {

        s += flg[i];

    }

    return s;

}

template <class DType>
double BPPPivots<DType>::cnt(size_t p, size_t pvtSize, size_t candSize)
{

    for(size_t i = 0; i < pvtSize; i++)
    {

        for(size_t k = 0; k <= p; k++)
        {

            rec[i][k] = 0;

        }

    }

    for(size_t k = 0; k < candSize; k++)
    {

        size_t j = 0;

        for(size_t i = 0; i < pvtSize; i++)
        {

            if(flg[rnk[k][i]] == 1)
            {

                j++;
                rec[rnk[k][i]][j]++;
                if(j == p) break;

            }

        }

    }

    for(size_t i = 0; i < pvtSize; i++)
    {

        if(flg[i])
        {

            for(size_t k = 1; k <= p; k++)
            {

                rec[i][0] += rec[i][k];

            }

            rec[i][0] /= 1.0*p;

        }

    }

    double s = 0;

    for(size_t i = 0; i < pvtSize; i++)
    {

        if(flg[i])
        {

            for(size_t k = 1; k <= p; k++)
            {

                s += (rec[i][k] - rec[i][0])*(rec[i][k] - rec[i][0])*1.0/p;

            }

        }

    }

    s = s*1.0/pvtSize;
    return s;

}

template <class DType>
BPPPivots<DType>::BPPPivots()
{

    this->setPivotType(PIVOT_TYPE::BPP);

}



template <class DType>
BPPPivots<DType>::~BPPPivots()
{

}

template <class DType>
void BPPPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t candSize = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2));
    size_t pvtSize = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2));
    size_t* pivot_index = uniqueRandomNumber(0, sample->getCardinality(), pvtSize, this->getSeed());
    size_t* cand_index = uniqueRandomNumber(0, sample->getCardinality(), candSize, this->getSeed());
    flg = new size_t[pvtSize];
    rnk = new size_t*[candSize];
    rec = new double*[pvtSize];
    double** dist = new double*[candSize];

    for(size_t i = 0; i < pvtSize; i++)
    {

        flg[i] = 1;
        rec[i] = new double[pvtSize];

    }

    for(size_t i = 0; i < candSize; i++)
    {

        rnk[i] = new size_t[pvtSize];
        dist[i] = new double[pvtSize];

        for(size_t j = 0; j < pvtSize; j++)
        {

            dist[i][j] = df->getDistance(sample->getFeatureVector(cand_index[i]), sample->getFeatureVector(pivot_index[j]));
            rnk[i][j] = j;

        }

    }

    for(size_t k = 0; k < candSize; k++)
    {

        for(size_t i = 0; i < pvtSize; i++)
        {

            for(size_t j = 0; j < pvtSize; j++)
            {

                if(dist[k][i] > dist[k][j])
                {

                    double q = dist[k][i];
                    dist[k][i] = dist[k][j];
                    dist[k][j] = q;
                    size_t p = rnk[k][i];
                    rnk[k][i] = rnk[k][j];
                    rnk[k][j] = p;

                }

            }

        }

    }

    size_t p = tot(pvtSize);

    while(tot(pvtSize) > nPivots)
    {

        if(p > 100) p = 100;
        double t1 = std::numeric_limits<double>::max();
        size_t min = -1;

        for(size_t i = 0; i < pvtSize; i++)
        {

            if(flg[i] == 1)
            {

                flg[i] = 0;
                double q = cnt(p, pvtSize, candSize);
                flg[i] = 1;

                if(q < t1)
                {

                    t1 = q;
                    min = i;

                }

            }

        }

        flg[min] = 0;

    }

    size_t pos = 0;
    for(size_t i = 0; i < pvtSize; i++)
    {

        if(flg[i])
        {

            this->setPivot(sample->getInstance(pivot_index[i]), pos++);

        }

        delete [] rec[i];

    }

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete [] pivot_index;
    delete [] cand_index;
    delete [] flg;

    for(size_t i = 0; i < candSize; i++)
    {

        delete [] rnk[i];
        delete [] dist[i];

    }

//    for(size_t i = 0; i < pvtSize; i++)
//    {

//        delete [] rec[i];

//    }

    delete [] rec;

    delete [] rnk;
    delete [] dist;

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());


}

template <class DType>
void BPPPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}


#endif // BPPPIVOTS_H
