#ifndef KMEDOIDS_H
#define KMEDOIDS_H

#include <Cluster.h>
#include <Dataset.h>
#include <Hermes.h>

template <class DType>
class Kmedoids
{

    private:

        size_t nClusters, iterations, seed;
        Cluster<DType>** clusters;

    private:
        void initializeClusters(Dataset<DType>* sample);
        void reCenter(DistanceFunction<BasicArrayObject<DType>>* df);
        void assignment(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df);


    public:

        //Constructors and destructors
        Kmedoids();
        ~Kmedoids();

        //Setters
        void setNumberOfClusters(size_t nClusters);
        void setSeed(size_t seed);
        void setNumberOfIterations(size_t iterations);

        //Getters
        size_t getSeed();
        size_t getNumberOfClusters() const;
        size_t getNumberOfIterations();
        BasicArrayObject<DType>** getMedoids();
        Cluster<DType>* getCluster(size_t pos);

        //Public methods
        void run(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nClusters);

};

template <class DType>
void Kmedoids<DType>::initializeClusters(Dataset<DType>* sample)
{

    size_t* randomIndex = uniqueRandomNumber(0, sample->getCardinality(), getNumberOfClusters(), getSeed());

    for(size_t x = 0; x < getNumberOfClusters(); x++){
        clusters[x]->setMedoid(sample->instance(randomIndex[x]));
        //std::cout << "ini: " << clusters[x]->getMedoid()->GetOID() << std::endl;
    }

}

template <class DType>
void Kmedoids<DType>::reCenter(DistanceFunction<BasicArrayObject<DType>>* df)
{

    double min = std::numeric_limits<double>::max(), total = 0, dist = 0;
    size_t index = 0, size = 0;

    for(size_t z = 0; z < getNumberOfClusters(); z++)
    {

        size = getCluster(z)->size();
        min = std::numeric_limits<double>::max();

        for(size_t x = 0; x < size; x++)
        {

            total = 0;

            for(size_t y = 0; y < size; y++)
            {

                if(x != y)
                {

                    dist = df->getDistance(*getCluster(z)->getInstance(x), *getCluster(z)->getInstance(y));
                    total += dist/((size-1)*1.0);

                }

            }

            if(total < min)
            {

                min = total;
                index = x;

            }

        }

        getCluster(z)->setMedoid(getCluster(z)->getInstance(index));
        //std::cout << "C = " << getCluster(z)->getMedoid()->GetOID() << std::endl;

    }

}

template <class DType>
void Kmedoids<DType>::assignment(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df)
{

    double min, dist;
    size_t index;

    for(size_t x = 0; x < getNumberOfClusters(); x++)
        getCluster(x)->clearSet();

    for(size_t x = 0; x < sample->getCardinality(); x++)
    {

        min = std::numeric_limits<double>::max();
        index = 0;

        for(size_t y = 0; y < getNumberOfClusters(); y++){

            dist = df->getDistance(*sample->instance(x), *getCluster(y)->getMedoid());

            if(dist < min)
            {

                min = dist;
                index = y;

            }

        }

        getCluster(index)->addInstance(sample->instance(x));

    }

}

template <class DType>
Kmedoids<DType>::Kmedoids()
{

    nClusters = 0;
    iterations = 500;
    clusters = nullptr;
    seed = 0;

}

template <class DType>
Kmedoids<DType>::~Kmedoids()
{

    if(clusters != nullptr)
    {

        for(size_t x = 0; x < getNumberOfClusters(); x++)
            delete (getCluster(x));

        delete [] clusters;
    }

}

template <class DType>
void Kmedoids<DType>::setNumberOfClusters(size_t nClusters)
{

    this->nClusters = nClusters;
    clusters = new Cluster<DType>*[getNumberOfClusters()];

    for(size_t x = 0; x < getNumberOfClusters(); x++)
        clusters[x] = new Cluster<DType>();

}

template <class DType>
void Kmedoids<DType>::setSeed(size_t seed)
{

    this->seed = seed;

}

template <class DType>
void Kmedoids<DType>::setNumberOfIterations(size_t iterations)
{

    this->iterations = iterations;

}

template <class DType>
size_t Kmedoids<DType>::getSeed()
{

    return seed;

}

template <class DType>
size_t Kmedoids<DType>::getNumberOfClusters() const
{

    return nClusters;

}

template <class DType>
size_t Kmedoids<DType>::getNumberOfIterations()
{

    return iterations;

}

template <class DType>
BasicArrayObject<DType>** Kmedoids<DType>::getMedoids()
{

    BasicArrayObject<DType>** ans = nullptr;

    if(getNumberOfClusters() != 0)
    {

        ans = new BasicArrayObject<DType>*[getNumberOfClusters()];

        for(size_t x = 0; x < getNumberOfClusters(); x++)
            ans[x] = getCluster(x)->getMedoid();

    }
    else
        throw std::invalid_argument("Number of cluster equal to 0 !_!");

    return ans;

}

template <class DType>
Cluster<DType>* Kmedoids<DType>::getCluster(size_t pos)
{

    Cluster<DType>* ans = nullptr;

    if((getNumberOfClusters() != 0) && (pos < getNumberOfClusters()))
    {

        ans = clusters[pos];

    }
    else
        throw std::invalid_argument("Number of cluster equal to 0 !_!");

    return ans;

}


template <class DType>
void Kmedoids<DType>::run(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nClusters)
{

    setNumberOfClusters(nClusters);
    initializeClusters(sample);

    size_t ite = 0;
    size_t* aux1 = new size_t[getNumberOfClusters()];
    size_t* aux2 = new size_t[getNumberOfClusters()];

    for(size_t x = 0; x  < getNumberOfClusters(); x++)
    {

        aux1[x] = 0;
        aux2[x] = 0;

    }

    while(ite < getNumberOfIterations())
    {

        for(size_t x = 0; x < getNumberOfClusters(); x++)
            aux1[x] = getCluster(x)->getMedoid()->getOID();

        assignment(sample, df);
        reCenter(df);

        for(size_t y = 0; y < getNumberOfClusters(); y++)
            aux2[y] = getCluster(y)->getMedoid()->getOID();

        bool ans = true;

        for(size_t z = 0; z < getNumberOfClusters(); z++)
        {

            if(aux1[z] != aux2[z])
                ans = false;

        }

        if(ans)
            break;

        ite++;
    }

    delete [] aux1;
    delete [] aux2;

}










#endif // KMEDOIDS_H
