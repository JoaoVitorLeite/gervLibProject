#ifndef PCAPIVOTS_H
#define PCAPIVOTS_H

#include <Pivot.h>
#include <Eigenvalues>

template <class DType>
class PCAPivots : public Pivot<DType>
{

    public:

        PCAPivots();
        PCAPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::PCA);
            generatePivots(sample,function,nPivots);

        }
        PCAPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::PCA);
            this->setSeed(seed);
            generatePivots(sample,function,nPivots);

        }

        ~PCAPivots();

        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};

template <class DType>
PCAPivots<DType>::PCAPivots()
{

    this->setPivotType(PIVOT_TYPE::PCA);


}

template <class DType>
PCAPivots<DType>::~PCAPivots()
{


}

template <class DType>
void PCAPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    double max = std::numeric_limits<double>::min();
    double min = std::numeric_limits<double>::max();

    Eigen::MatrixXd A;
    A.resize(sample->getCardinality(), sample->getCardinality());

    for(size_t i = 0; i < sample->getCardinality(); i++)
    {

        for(size_t j = 0; j < sample->getCardinality(); j++)
        {

            if(i < j)
            {

                A(i,j) = df->getDistance(*sample->instance(i), *sample->instance(j));
                A(j,i) = A(i,j);

                if(A(i,j) > max) max = A(i,j);
                if(A(i,j) < min) min = A(i,j);

            }

            if(i == j) A(i,j) = 0.0;

        }

    }

    for(size_t i = 0; i < sample->getCardinality(); i++)
    {

        for(size_t j = 0; j < sample->getCardinality(); j++)
        {

            if(i < j)
            {

                A(i,j) = (A(i,j) - min)/(max - min);
                A(j,i) = A(i,j);

            }

        }

    }

    Eigen::EigenSolver<Eigen::MatrixXd> s(A);
    Eigen::VectorXd eigenValues = s.eigenvalues().real();

    std::vector<std::pair<double, long>> v;

    for(long i = 0; i < eigenValues.size(); i++)
    {

        std::pair<double, long> tupleAux(eigenValues[i], i);
        v.push_back(tupleAux);

    }

    std::sort(v.begin(), v.end(), [](const std::pair<double, long>& a, const std::pair<double, long>& b){ return std::get<0>(a) > std::get<0>(b); });

    for(size_t x = 0; x < this->getNumberOfPivots(); x++)
        this->setPivot(sample->instance(v[x].second), x);

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

}


template <class DType>
void PCAPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots);

}







#endif // PCAPIVOTS_H
