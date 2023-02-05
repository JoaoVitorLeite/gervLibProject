#ifndef FFTPIVOTS_H
#define FFTPIVOTS_H

#include <Pivot.h>

template <class DType>
class FFTPivots : public Pivot<DType>
{

    public:
        FFTPivots();
        FFTPivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::FFT);
            generatePivots(sample, function, nPivots);

        }

        FFTPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::FFT);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots);

        }

        ~FFTPivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};


template <class DType>
FFTPivots<DType>::FFTPivots()
{

    this->setPivotType(PIVOT_TYPE::FFT);

}

template <class DType>
FFTPivots<DType>::~FFTPivots()
{



}

template <class DType>
void FFTPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t p1, currentPivot = 0;
    size_t* pvtIndex = new size_t[this->getNumberOfPivots()];
    double *dist = new double[sample->getCardinality()];
    bool *bitmap = new bool[sample->getCardinality()];

    for(size_t x = 0; x < sample->getCardinality(); x++)
        bitmap[x] = false;

    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
        pvtIndex[x] = 0;

    size_t* aux = uniqueRandomNumber(0, sample->getCardinality(), 1, this->getSeed());
    p1 = aux[0];
    bitmap[p1] = true;
    pvtIndex[currentPivot] = p1;
    currentPivot++;

    for(size_t i = 0; i < sample->getCardinality(); i++)
    {

        dist[i] = df->getDistance(*sample->instance(p1), *sample->instance(i));

    }

    double max = 0.0;
    size_t pos = 0;

    for(size_t i = 0; i < this->getNumberOfPivots(); i++)
    {

        max = 0.0;
        pos = -1;

        for(size_t i = 0; i < sample->getCardinality(); i++)
        {

            if(!bitmap[i])
            {

                if(dist[i] > max)
                {

                    max = dist[i];
                    pos = i;

                }

            }

        }

        pvtIndex[currentPivot] = pos;
        currentPivot++;
        bitmap[pos] = true;

        for(size_t j = 0; j < sample->getCardinality(); j++)
            dist[j] = std::min(dist[j], max);

    }

    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
        this->setPivot(sample->instance(pvtIndex[x]), x);

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete [] (bitmap);
    delete [] (pvtIndex);
    delete [] (dist);

}

template <class DType>
void FFTPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}












#endif // FFTPIVOTS_H
