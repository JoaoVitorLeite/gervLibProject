#ifndef SSSPIVOTS_H
#define SSSPIVOTS_H

#include <Pivot.h>

template <class DType>
class SSSPivots : public Pivot<DType>
{

    private:

        double alpha;

    public:

        SSSPivots();
        SSSPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, double alpha = 0.15) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::SSS);
            generatePivots(sample,function,nPivots,alpha);

        }
        SSSPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed, double alpha = 0.15) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::SSS);
            this->setSeed(seed);
            generatePivots(sample,function,nPivots,alpha);

        }
        ~SSSPivots();

        void setAlpha(double alpha);

        double getAlpha();

        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, double alpha, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, double alpha, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};


template <class DType>
SSSPivots<DType>::SSSPivots()
{

    setAlpha(0.35);
    this->setPivotType(PIVOT_TYPE::SSS);

}

template <class DType>
SSSPivots<DType>::~SSSPivots()
{


}

template <class DType>
void SSSPivots<DType>::setAlpha(double alpha)
{

    this->alpha = alpha;

}

template <class DType>
double SSSPivots<DType>::getAlpha()
{

    return alpha;

}

template <class DType>
void SSSPivots<DType>::generatePivots(Dataset<DType>* dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t drop = 2, currentPivot = 0, p1, index = 0, count = 0;
    size_t* pvtIndex = new size_t[this->getNumberOfPivots()+drop];
    size_t* aux;
    double max = std::numeric_limits<double>::min(), dist, threshold;

    for(size_t x = 0; x < (this->getNumberOfPivots()+drop); x++)
        pvtIndex[x] = 0;

    for(size_t x = 0; x < sample->getCardinality(); x++)
    {

        for(size_t y = 0; y < sample->getCardinality(); y++)
        {

            if(x < y)
            {

                dist = df->getDistance(*sample->instance(x), *sample->instance(y));

                if(dist > max)
                    max = dist;

            }

        }

    }

    threshold = max * getAlpha();

    aux = uniqueRandomNumber(0, sample->getCardinality(), 1, this->getSeed());
    p1 = aux[0];
    pvtIndex[currentPivot] = p1;
    currentPivot++;

    //std::cout << "AUX : " << p1 << std::endl;

    while(index < sample->getCardinality())
    {

        count = 0;

        for(size_t x = 0; x < currentPivot; x++)
        {

            dist = df->getDistance(*sample->instance(index), *sample->instance(pvtIndex[x]));

            if(dist >= threshold)
                count++;

        }

        if(count == currentPivot)
        {

            pvtIndex[currentPivot] = index;
            currentPivot++;

        }

        if((currentPivot-drop) == this->getNumberOfPivots())
            break;
        else
            index++;

    }

    //std::cout << "AUX : " << pvtIndex[1] << std::endl;

    for(size_t x = drop; x < (this->getNumberOfPivots()+drop); x++)
        this->setPivot(sample->instance(pvtIndex[x]), x-drop);

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete [] (aux);
    delete [] (pvtIndex);

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

}

template <class DType>
void SSSPivots<DType>::generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, double alpha, std::vector<std::string> args)
{

    setAlpha(alpha);
    generatePivots(sample, df, nPivots);

}

template <class DType>
void SSSPivots<DType>::generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, double alpha, size_t seed, std::vector<std::string> args)
{

    setAlpha(alpha);
    this->setSeed(seed);
    generatePivots(sample, df, nPivots);

}



#endif // SSSPIVOTS_H
