#ifndef SELECTIONPIVOTS_H
#define SELECTIONPIVOTS_H

#include <Pivot.h>

template <class DType>
class SelectionPivots : public Pivot<DType>
{

    private:

        size_t nGroups;

    public:

        SelectionPivots();
        SelectionPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t nGroups) : Pivot<DType>(){

            setNumberOfGroups(nGroups);
            this->setPivotType(PIVOT_TYPE::SELECTION);
            generatePivots(sample, function, nPivots, nGroups);

        }
        SelectionPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t nGroups, size_t seed) : Pivot<DType>(){

            setNumberOfGroups(nGroups);
            this->setPivotType(PIVOT_TYPE::SELECTION);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots, nGroups);

        }

        ~SelectionPivots();

        void setNumberOfGroups(const size_t &value);

        size_t getNumberOfGroups() const;

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t nGroups, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t nGroups, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};


template <class DType>
SelectionPivots<DType>::SelectionPivots()
{

    setNumberOfGroups(5);
    this->setPivotType(PIVOT_TYPE::SELECTION);


}

template <class DType>
SelectionPivots<DType>::~SelectionPivots()
{


}

template <class DType>
void SelectionPivots<DType>::setNumberOfGroups(const size_t &value)
{

    nGroups = value;

}

template <class DType>
size_t SelectionPivots<DType>::getNumberOfGroups() const
{

    return nGroups;

}

template <class DType>
void SelectionPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    auto start = std::chrono::steady_clock::now();

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    setNumberOfGroups(size_t(0.5 * sample->getCardinality()));
    this->setNumberOfPivots(nPivots);

    double max = std::numeric_limits<double>::lowest();
    double length = ((this->getNumberOfPivots()*(this->getNumberOfPivots()-1))*1.0)/2.0;
    double mean = 0, dist;
    std::vector<BasicArrayObject<DType>*> group = std::vector<BasicArrayObject<DType>*>();

    srand(this->getSeed());

    for(size_t g = 0; g < getNumberOfGroups(); g++){

        srand(this->getSeed());
        group = sample->sample(this->getNumberOfPivots(), false, rand());
        mean = 0;

        for(size_t x = 0; x < this->getNumberOfPivots(); x++)
        {

            for(size_t y = 0; y < this->getNumberOfPivots(); y++)
            {

                if(x < y)
                {

                    dist = df->getDistance(*group[x], *group[y]);
                    mean += dist/length;

                }

            }

        }

        if(mean > max)
        {

            max = mean;

            for(size_t k = 0; k < this->getNumberOfPivots(); k++)
            {

                this->setPivot(group[k], k);

            }

        }


    }

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

}

template <class DType>
void SelectionPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t nGroups, std::vector<std::string> args)
{

    setNumberOfGroups(nGroups);
    generatePivots(sample, df, nPivots);

}



template <class DType>
void SelectionPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t nGroups, size_t seed, std::vector<std::string> args)
{

    setNumberOfGroups(nGroups);
    this->setSeed(seed);
    generatePivots(sample, df, nPivots);

}







#endif // SELECTIONPIVOTS_H
