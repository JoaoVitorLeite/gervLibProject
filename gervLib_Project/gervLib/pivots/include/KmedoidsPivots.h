#ifndef KMEDOIDSPIVOTS_H
#define KMEDOIDSPIVOTS_H

#include <Pivot.h>
#include <Kmedoids.h>

template <class DType>
class KmedoidsPivots : public Pivot<DType>
{

    private:

        size_t nIterations;

    public:

        KmedoidsPivots();
        KmedoidsPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>()
        {

            setNumberOfIterations(500);
            this->setPivotType(PIVOT_TYPE::KMEDOIDS);
            generatePivots(sample, function, nPivots);

        }
        KmedoidsPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>()
        {

            this->setPivotType(PIVOT_TYPE::KMEDOIDS);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots);

        }
        KmedoidsPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed, size_t nIterations) : Pivot<DType>()
        {

            this->setPivotType(PIVOT_TYPE::KMEDOIDS);
            this->setSeed(seed);
            setNumberOfIterations(nIterations);
            generatePivots(sample, function, nPivots);

        }
        ~KmedoidsPivots();

        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t nItearations, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t nItearations, size_t seed, std::vector<std::string> args = std::vector<std::string>());

        void setNumberOfIterations(size_t value);

        size_t getNumberOfIterations();

};



template <class DType>
KmedoidsPivots<DType>::KmedoidsPivots()
{

    setNumberOfIterations(500);
    this->setPivotType(PIVOT_TYPE::KMEDOIDS);

}


template <class DType>
KmedoidsPivots<DType>::~KmedoidsPivots()
{

}

template <class DType>
void KmedoidsPivots<DType>::generatePivots(Dataset<DType>* dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    Kmedoids<DType>* aux = new Kmedoids<DType>();
    aux->setSeed(this->getSeed());
    aux->setNumberOfIterations(getNumberOfIterations());

    aux->run(sample, df, nPivots);

    for(size_t x = 0; x < this->getNumberOfPivots(); x++)
        this->setPivot(aux->getCluster(x)->getMedoid(), x);

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete (aux);

}

template <class DType>
void KmedoidsPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t nItearations, std::vector<std::string> args)
{

    setNumberOfIterations(nItearations);
    generatePivots(sample, df, nPivots);

}

template <class DType>
void KmedoidsPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t nItearations, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    setNumberOfIterations(nItearations);
    generatePivots(sample, df, nPivots);

}

template <class DType>
void KmedoidsPivots<DType>::setNumberOfIterations(size_t value)
{

    nIterations = value;

}

template <class DType>
size_t KmedoidsPivots<DType>::getNumberOfIterations()
{

    return nIterations;

}

















#endif // KMEDOIDSPIVOTS_H
