#ifndef CLUSTER_H
#define CLUSTER_H

#include <BasicArrayObject.h>

template <class DType>
class Cluster
{

    private:

        BasicArrayObject<DType>* medoid; //Medoid
        std::vector<BasicArrayObject<DType>*> set; //Set of elements inside the cluster(Medoid is also included)

    public:

        //Constructors and destructor
        Cluster();
        Cluster(BasicArrayObject<DType>* medoid);
        ~Cluster();

        //Setters
        void setMedoid(BasicArrayObject<DType>* medoid);
        void setSet(std::vector<BasicArrayObject<DType>*> vec);

        //Getters
        BasicArrayObject<DType>* getInstance(size_t pos);
        BasicArrayObject<DType>* getMedoid();
        size_t size();

        //Public method
        void addInstance(BasicArrayObject<DType>* inst);
        void clearSet();

        //Print
        void printCluster();

};



template <class DType>
Cluster<DType>::Cluster()
{

    medoid = nullptr;
    set = std::vector<BasicArrayObject<DType>*>();

}

template <class DType>
Cluster<DType>::Cluster(BasicArrayObject<DType>* medoid)
{

    this->medoid = medoid;
    set = std::vector<BasicArrayObject<DType>*>();

}

template <class DType>
Cluster<DType>::~Cluster()
{

    clearSet();

}

template <class DType>
void Cluster<DType>::setMedoid(BasicArrayObject<DType> *medoid)
{

    this->medoid = medoid;

}

template <class DType>
void Cluster<DType>::setSet(std::vector<BasicArrayObject<DType>*> vec)
{

    this->set = vec;

}

template <class DType>
BasicArrayObject<DType>* Cluster<DType>::getInstance(size_t pos)
{

    BasicArrayObject<DType>* aux = nullptr;

    if((pos >= set.size()) || (pos < 0))
    {

        throw std::invalid_argument("Index Out of Bound !_!");

    }
    else
    {

        aux = set[pos];

    }

    return aux;

}

template <class DType>
BasicArrayObject<DType>* Cluster<DType>::getMedoid()
{

    return medoid;

}

template <class DType>
size_t Cluster<DType>::size()
{

    return set.size();

}

template <class DType>
void Cluster<DType>::addInstance(BasicArrayObject<DType> *inst)
{

    set.push_back(inst);

}

template <class DType>
void Cluster<DType>::clearSet()
{

    set.clear();

}

template <class DType>
void Cluster<DType>::printCluster()
{

    std::cout << "Medoid = " << medoid->toStringWithOID() << std::endl;

    for(BasicArrayObject<DType>* aux : this->set)
    {

        std::cout << aux->toStringWithOID() << std::endl;

    }

}
















#endif // CLUSTER_H
