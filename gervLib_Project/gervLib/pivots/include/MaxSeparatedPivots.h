#ifndef MAXSEPARETEDPIVOTS_H
#define MAXSEPARETEDPIVOTS_H

#include <Pivot.h>

template <class DType>
class MaxSeparatedPivots : public Pivot<DType>
{

    public:

        MaxSeparatedPivots();
        MaxSeparatedPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::MAXSEPARATED);
            generatePivots(sample, function, nPivots);

        }
        MaxSeparatedPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::MAXSEPARATED);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots);

        }
        ~MaxSeparatedPivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};


template <class DType>
MaxSeparatedPivots<DType>::MaxSeparatedPivots()
{

    this->setPivotType(PIVOT_TYPE::MAXSEPARATED);

}

template <class DType>
MaxSeparatedPivots<DType>::~MaxSeparatedPivots()
{


}

template <class DType>
void MaxSeparatedPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t drop = 2, currentPivot = 0, pos = 0, p1;
    bool* bitmap = new bool[sample->getCardinality()];
    size_t* pvtIndex = new size_t[this->getNumberOfPivots()+drop];
    size_t* aux;
    double max = std::numeric_limits<double>::min(), dist, sum = 0;

    for(size_t x = 0; x < sample->getCardinality(); x++)
        bitmap[x] = false;

    for(size_t x = 0; x < (this->getNumberOfPivots()+drop); x++)
        pvtIndex[x] = 0;

    aux = uniqueRandomNumber(0, sample->getCardinality(), 1, this->getSeed());
    p1 = aux[0];
    bitmap[p1] = true;
    pvtIndex[currentPivot] = p1;
    currentPivot++;

    //std::cout << "AUX : " << p1 << std::endl;

    for(size_t x = 0; x < sample->getCardinality(); x++)
    {

        dist = df->getDistance(*sample->instance(p1), *sample->instance(x));

        if(dist > max)
        {

            max = dist;
            pos = x;

        }

    }

    bitmap[pos] = true;
    pvtIndex[currentPivot] = pos;
    currentPivot++;

    //std::cout << "AUX : " << pos << std::endl;

    while((currentPivot-drop) < this->getNumberOfPivots())
    {

        max = std::numeric_limits<double>::min();

        for(size_t x = 0; x < sample->getCardinality(); x++)
        {

            if(!bitmap[x])
            {

                sum = 0;

                for(size_t y = 0; y < currentPivot; y++)
                    sum += df->getDistance(*sample->instance(x), *sample->instance(pvtIndex[y]));

                if(sum > max)
                {

                    max = sum;
                    pos = x;

                }

            }

        }

        bitmap[pos] = true;
        pvtIndex[currentPivot] = pos;
        currentPivot++;

    }

    for(size_t x = drop; x < (this->getNumberOfPivots()+drop); x++)
        this->setPivot(sample->instance(pvtIndex[x]), x-drop);


    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete [] (bitmap);
    delete [] (pvtIndex);

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

}



template <class DType>
void MaxSeparatedPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}










#endif // MAXSEPARETEDPIVOTS_H
