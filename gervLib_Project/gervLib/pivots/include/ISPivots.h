#ifndef ISPIVOTS_H
#define ISPIVOTS_H

#include <Pivot.h>

template <class DType>
class ISPivots : public Pivot<DType>
{

    private:
        size_t** ppa;
        double* svg;
        size_t CandA = 300;
        double dfs(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t k);

    public:
        ISPivots();
        ISPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::IS);
            generatePivots(sample, function, nPivots);

        }

        ISPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::IS);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots);

        }

        ~ISPivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};


template <class DType>
double ISPivots<DType>::dfs(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t k)
{

    double sum = 0.0;
    for(size_t i = 0; i < CandA; i++)
    {

        double q = abs(df->getDistance(*sample->instance(k), *sample->instance(ppa[i][1])) - df->getDistance(*sample->instance(k), *sample->instance(ppa[i][0])));
        if(svg[i] > q) sum += svg[i];
        else sum += q;

    }

    return sum;

}

template <class DType>
ISPivots<DType>::ISPivots()
{

    this->setPivotType(PIVOT_TYPE::IS);

}

template <class DType>
ISPivots<DType>::~ISPivots()
{



}

template <class DType>
void ISPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

//    size_t CandA = 300;
    size_t NUUM = sample->getCardinality()/2;
    ppa = new size_t*[NUUM];

    for(size_t i = 0; i < NUUM; i++) ppa[i] = new size_t[2];

    size_t can[NUUM];
    size_t fgg[NUUM];
    size_t pvtIndex[nPivots];
    svg = new double[NUUM];

    for(size_t i = 0; i < CandA; i++)
    {

        size_t p = rand()%NUUM;
        ppa[i][0] = p;
        p = rand()%NUUM;
        ppa[i][1] = p;
        if(ppa[i][0] == ppa[i][1]) i--;
        svg[i] = 0;

    }

    double max = 0.0;
    size_t k = 0;

    for(size_t i = 0; i < CandA; i++) fgg[can[i]] = 0;
    for(size_t i = 0; i < CandA; i++)
    {

        size_t p = rand()%NUUM;
        if(fgg[p]) i--;
        else
        {

            can[i] = p;
            fgg[p] = 1;

        }

    }

    for(size_t i = 0; i < nPivots; i++)
    {

        max = 0.0;

        for(size_t j = 0; j < CandA; j++)
        {

            double tmp = dfs(sample, df, can[j]);
            if(tmp > max)
            {

                max = tmp;
                k = can[j];

            }

        }

        pvtIndex[i] = k;

        for(size_t i = 0; i < CandA; i++)
        {

            double q = abs(df->getDistance(*sample->instance(k), *sample->instance(ppa[i][1])) - df->getDistance(*sample->instance(k), *sample->instance(ppa[i][0])));
            if(svg[i] < q) svg[i] = q;

        }

    }

    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
        this->setPivot(sample->instance(pvtIndex[x]), x);

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

}

template <class DType>
void ISPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}












#endif // ISPIVOTS_H
