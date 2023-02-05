#ifndef RANDOMPIVOTS_H
#define RANDOMPIVOTS_H

#include <Pivot.h>

template <class DType>
class RandomPivots : public Pivot<DType>
{

    public:
        RandomPivots();
        RandomPivots(Dataset<DType> *sample, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::RANDOM);
            generatePivots(sample, NULL, nPivots);

        }
        RandomPivots(Dataset<DType> *sample, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::RANDOM);
            this->setSeed(seed);
            generatePivots(sample, NULL, nPivots);
        }

        ~RandomPivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());


};


template <class DType>
RandomPivots<DType>::~RandomPivots()
{



}

template <class DType>
RandomPivots<DType>::RandomPivots()
{

    this->setPivotType(PIVOT_TYPE::RANDOM);

}



template <class DType>
void RandomPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    //Gerar a lista de pivôs
    //Gerar uma lista de ids sem reposição de tamanha nPivots
    size_t *randomIds = uniqueRandomNumber(0, sample->getCardinality(), this->getNumberOfPivots(), this->getSeed());

    for(size_t x = 0; x < this->getNumberOfPivots(); x++)
    {

        this->setPivot(sample->instance(randomIds[x]), x);

    }


    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete[] randomIds;

}



template <class DType>
void RandomPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots);

}



#endif // RANDOMPIVOTS_H
