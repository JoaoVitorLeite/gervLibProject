#ifndef CONVEXPIVOTS_H
#define CONVEXPIVOTS_H

#include <Pivot.h>

template <class DType>
class ConvexPivots : public Pivot<DType>
{

    public:
        ConvexPivots();
        ConvexPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::CONVEX);
            generatePivots(sample, df, nPivots);

        }
        ConvexPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::CONVEX);
            this->setSeed(seed);
            generatePivots(sample, df, nPivots);

        }
        ~ConvexPivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};



template <class DType>
ConvexPivots<DType>::ConvexPivots()
{

    this->setPivotType(PIVOT_TYPE::CONVEX);

}



template <class DType>
ConvexPivots<DType>::~ConvexPivots()
{

}



template <class DType>
void ConvexPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t drop = 2;
    size_t currentPivot = 0, p1, pos = 0;
    size_t *aux;
    double max = std::numeric_limits<double>::min();

    bool *bitmap = new bool[sample->getCardinality()];
    size_t *pvtIndex = new size_t[this->getNumberOfPivots()+drop];

    aux = uniqueRandomNumber(0, sample->getCardinality(), 1, this->getSeed());
    p1 = aux[0];

    //std::cout << "AUX: " << p1 << std::endl;

    //Inicializa bitmap
    for (size_t x = 0; x < sample->getCardinality(); x++)
        bitmap[x] = false;

    //Incializa ponteiros de pivôs
    for (size_t x = 0; x < this->getNumberOfPivots()+drop; x++)
        pvtIndex[x] = 0;

    //Cria primeiro pivô aleatório
    bitmap[p1] = true;

    //Calcula a distância do primeiro pivô aleatório para toda a amostra
    for(size_t x = 0; x < sample->getCardinality(); x++)
    {

        double dist = df->getDistance(*sample->instance(x), *sample->instance(p1));

        if(dist > max)
        {

            max = dist;
            pos = x;

        }

    }

    //Insere elemento mais afastado na lista de pivôs
    bitmap[pos] = true;
    pvtIndex[currentPivot] = pos;
    currentPivot++;

    //std::cout << "AUX: " << pos << std::endl;

    max = std::numeric_limits<double>::min();
    pos = 0;

    //Calcula a distância do elemento mais afastado do elemento mais distante do primeiro pivô aleatório para toda a amostra
    for(size_t x = 0; x < sample->getCardinality(); x++)
    {

        double dist = df->getDistance(*sample->instance(x), *sample->instance(pvtIndex[0]));

        if((dist > max) && (!bitmap[x]))
        {

            max = dist;
            pos = x;

        }

    }

    //Insere elemento mais afastado mais distante do primeiro pivô na lista de pivôs
    bitmap[pos] = true;
    pvtIndex[currentPivot] = pos;
    currentPivot++;

    //std::cout << "AUX: " << pos << std::endl;

    //Definição da primeira aresta
    double edge = df->getDistance(*sample->instance(pvtIndex[0]), *sample->instance(pvtIndex[1]));

    while ((currentPivot-drop) < this->getNumberOfPivots())
    {

        double error = 0;
        double min = std::numeric_limits<double>::max();
        pos = 0;

        for(size_t x = 0; x < sample->getCardinality(); x++)
        {

            if(!bitmap[x])
            {

                error = 0;
                for(size_t y = 0; y < currentPivot; y++)
                    error += std::abs(edge - df->getDistance(*sample->instance(pvtIndex[y]), *sample->instance(x)));


                if(min > error)
                {

                    min = error;
                    pos = x;

                }

            }

        }

        bitmap[pos] = true;
        pvtIndex[currentPivot] = pos;
        currentPivot++;

    }

    //Transfere os pivôs que foram marcados no bitmap para a lista de pivôs herdada de Pivot
    for(size_t x = drop; x < (this->getNumberOfPivots()+drop); x++)
        this->setPivot(sample->instance(pvtIndex[x]), x-drop);


    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;
    delete[] bitmap;
    delete[] pvtIndex;
    delete[] aux;

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

}



template <class DType>
void ConvexPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}



#endif // CONVEXPIVOTS_H
