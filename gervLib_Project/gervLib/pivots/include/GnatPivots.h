#ifndef GNATPIVOTS_H
#define GNATPIVOTS_H

#include <Pivot.h>

template <class DType>
class GnatPivots : public Pivot<DType>
{

    public:

        GnatPivots();
        GnatPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::GNAT);
            generatePivots(sample, df, nPivots);

        }
        GnatPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::GNAT);
            this->setSeed(seed);
            generatePivots(sample, df, nPivots);

        }

        ~GnatPivots();

        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};



template <class DType>
GnatPivots<DType>::GnatPivots()
{

    this->setPivotType(PIVOT_TYPE::GNAT);

}



template <class DType>
GnatPivots<DType>::~GnatPivots()
{

}



template <class DType>
void GnatPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t drop = 2, currentPivot = 0, p1, pos = 0;
    bool* bitmap = new bool[sample->getCardinality()];
    size_t* pvtIndex = new size_t[this->getNumberOfPivots()+drop];
    size_t* aux;
    double max = std::numeric_limits<double>::min(), dist, currentDist;

    for(size_t x = 0; x < sample->getCardinality(); x++)
        bitmap[x] = false;

    for(size_t x = 0; x < (this->getNumberOfPivots()+drop); x++)
        pvtIndex[x] = 0;

    aux = uniqueRandomNumber(0, sample->getCardinality(), 1, this->getSeed());
    p1 = aux[0];
    bitmap[p1] = true;
    pvtIndex[currentPivot] = p1;
    currentPivot++;

//    std::cout << "AUX: " << p1 << std::endl;

    for(size_t x = 0; x < sample->getCardinality(); x++)
    {

        dist = df->getDistance(*sample->instance(pvtIndex[0]), *sample->instance(x));

        if(dist > max)
        {

            max = dist;
            pos = x;

        }

    }

    bitmap[pos] = true;
    pvtIndex[currentPivot] = pos;
    currentPivot++;

//    std::cout << "AUX: " << pos << std::endl;

    while((currentPivot - drop) < this->getNumberOfPivots())
    {

        pos = 0;
        max = std::numeric_limits<double>::min();

        for(size_t x = 0; x < sample->getCardinality(); x++)
        {

            if(!bitmap[x])
            {

                currentDist = std::numeric_limits<double>::max();

                for(size_t y = 0; y < currentPivot; y++)
                {

                    dist = df->getDistance(*sample->instance(x), *sample->instance(pvtIndex[y]));
                    currentDist = std::min(currentDist, dist);

                }

                if(currentDist > max)
                {

                    max = currentDist;
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
    delete [] (aux);
    delete [] (bitmap);
    delete [] (pvtIndex);

}



template <class DType>
void GnatPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}








#endif // GNATPIVOTS_H
