#ifndef ISPIVOTS_H
#define ISPIVOTS_H

#include <Pivot.h>

template <class DType>
class ISPivots : public Pivot<DType>
{

    private:
//        size_t** ppa;
//        double* svg;
//        size_t CandA = 300;
//        double dfs(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t k);

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


//template <class DType>
//double ISPivots<DType>::dfs(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t k)
//{

//    double sum = 0.0;
//    for(size_t i = 0; i < CandA; i++)
//    {

//        double q = abs(df->getDistance(*sample->instance(k), *sample->instance(ppa[i][1])) - df->getDistance(*sample->instance(k), *sample->instance(ppa[i][0])));
//        if(svg[i] > q) sum += svg[i];
//        else sum += q;

//    }

//    return sum;

//}

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

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

//    size_t pairsSize = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2)), candidatesSize = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2));
    //    size_t cand_size = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2)), pair_size = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2));
    size_t candidatesSize = std::min(2 * this->getNumberOfPivots(), sample->getCardinality());
    //size_t pair_size = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2));
    size_t pairsSize = 300;
    size_t* pvtIndex = new size_t[this->getNumberOfPivots()];
    bool* bitmap = new bool[candidatesSize];
    size_t** pairs = new size_t*[pairsSize];
    size_t* candidates = uniqueRandomNumber(0, sample->getCardinality(), candidatesSize, this->getSeed());
    double* svg = new double[pairsSize];
    size_t* aux;
    double max = std::numeric_limits<double>::min(), q = 0.0;
    size_t pos = 0, currentPivot = 0;

    for(size_t x = 0; x < candidatesSize; x++)
    {

        bitmap[x] = true;

    }

    for(size_t x = 0; x < pairsSize; x++)
    {

        svg[x] = 0.0;

        aux = uniqueRandomNumber(0, sample->getCardinality(), 2, this->getSeed()/2);
        pairs[x] = new size_t[2];
        pairs[x][0] = aux[0];
        pairs[x][1] = aux[1];

        delete [] aux;

    }

    for(size_t x = 0; x < this->getNumberOfPivots(); x++)
    {

        max = std::numeric_limits<double>::min();
        pos = 0;

        for(size_t y = 0; y < candidatesSize; y++)
        {

            if(bitmap[y])
            {

                double dfs = 0.0;
                double d = 0.0;
                for(size_t z = 0; z < pairsSize; z++)
                {

                    d = abs(df->getDistance(*sample->instance(candidates[y]), *sample->instance(pairs[z][1])) - df->getDistance(*sample->instance(candidates[y]), *sample->instance(pairs[z][0])));

                    if(svg[z] > d)
                    {

                        dfs += svg[z];

                    }
                    else
                    {

                        dfs += d;

                    }

                }

                if(dfs > max)
                {

                    max = dfs;
                    pos = y;

                }

            }

        }

        pvtIndex[currentPivot++] = candidates[pos];
        bitmap[pos] = false;

        for(size_t i = 0; i < pairsSize; i++)
        {

            q = abs(df->getDistance(*sample->instance(candidates[pos]), *sample->instance(pairs[i][1])) - df->getDistance(*sample->instance(candidates[pos]), *sample->instance(pairs[i][0])));

            if(svg[i] < q)
            {

                svg[i] = q;

            }

        }

    }

    for(size_t x = 0; x < this->getNumberOfPivots(); x++)
    {

        this->setPivot(sample->instance(pvtIndex[x]), x);

    }

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete [] pvtIndex;
    delete [] bitmap;
    delete [] candidates;
    delete [] svg;

    aux = nullptr;
    delete aux;

    for(size_t x = 0; x < pairsSize; x++)
    {

        delete [] pairs[x];

    }

    delete [] pairs;

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());


//**********************************************************************************************************************************
//    size_t CandA = 300;
//    size_t NUUM = sample->getCardinality()/2;
//    CandA = std::min(NUUM, CandA);
//    ppa = new size_t*[NUUM];

//    for(size_t i = 0; i < NUUM; i++) ppa[i] = new size_t[2];

//    size_t can[NUUM];
//    size_t fgg[NUUM];
//    size_t pvtIndex[nPivots];
//    svg = new double[NUUM];

//    for(size_t i = 0; i < CandA; i++)
//    {

//        size_t p = rand()%NUUM;
//        ppa[i][0] = p;
//        p = rand()%NUUM;
//        ppa[i][1] = p;
//        if(ppa[i][0] == ppa[i][1]) i--;
//        svg[i] = 0;

//        can[i] = i;

//    }

//    double max = 0.0;
//    size_t k = 0;

//    for(size_t i = 0; i < CandA; i++) fgg[can[i]] = 0;
//    for(size_t i = 0; i < CandA; i++)
//    {

//        size_t p = rand()%NUUM;
//        if(fgg[p]) i--;
//        else
//        {

//            can[i] = p;
//            fgg[p] = 1;

//        }

//    }

//    for(size_t i = 0; i < nPivots; i++)
//    {

//        max = 0.0;

//        for(size_t j = 0; j < CandA; j++)
//        {

//            double tmp = dfs(sample, df, can[j]);
//            if(tmp > max)
//            {

//                max = tmp;
//                k = can[j];

//            }

//        }

//        pvtIndex[i] = k;

//        for(size_t i = 0; i < CandA; i++)
//        {

//            double q = abs(df->getDistance(*sample->instance(k), *sample->instance(ppa[i][1])) - df->getDistance(*sample->instance(k), *sample->instance(ppa[i][0])));
//            if(svg[i] < q) svg[i] = q;

//        }

//    }

//    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
//        this->setPivot(sample->instance(pvtIndex[x]), x);

//    if(this->sample_size == -1.0)
//        sample = nullptr;

//    delete sample;

}

template <class DType>
void ISPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}












#endif // ISPIVOTS_H
