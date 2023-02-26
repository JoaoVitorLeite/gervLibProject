#ifndef MAXVARIANCEPIVOTS_H
#define MAXVARIANCEPIVOTS_H

#include <Pivot.h>

template <class DType>
class MaxVariancePivots : public Pivot<DType>
{

    public:

        MaxVariancePivots();
        MaxVariancePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::MAXVARIANCE);
            generatePivots(sample, function, nPivots);

        }
        MaxVariancePivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::MAXVARIANCE);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots);

        }
        ~MaxVariancePivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};



template <class DType>
MaxVariancePivots<DType>::MaxVariancePivots()
{

    this->setPivotType(PIVOT_TYPE::MAXVARIANCE);


}



template <class DType>
MaxVariancePivots<DType>::~MaxVariancePivots()
{


}



template <class DType>
void MaxVariancePivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    std::pair<std::vector<BasicArrayObject<DType>*>, std::vector<BasicArrayObject<DType>*>> split = sample->splitDataset(0.5,0.5, this->getSeed());

    double mean, sum;
    std::vector<std::pair<BasicArrayObject<DType>*, double>> v;

    for(size_t x = 0; x < split.first.size(); x++)
    {

        mean = 0, sum = 0;

        for(size_t y = 0; y < split.second.size(); y++)
            mean += df->getDistance(*std::get<0>(split)[x], *std::get<1>(split)[y]);

        mean /= split.second.size()*1.0;

        for(size_t z = 0; z < split.second.size(); z++)
        {

            double dist = df->getDistance(*std::get<0>(split)[x], *std::get<1>(split)[z]);
            sum += (dist-mean)*(dist-mean);

        }

        sum /= split.second.size()*1.0;

        v.push_back(std::make_pair(std::get<0>(split)[x], sum));

    }

    std::sort(v.begin(), v.end(), [](const std::pair<BasicArrayObject<DType>*, double> &p1, const std::pair<BasicArrayObject<DType>*, double> &p2){ return p1.second > p2.second; });

    for(size_t z = 0; z < this->getNumberOfPivots(); z++)
        this->setPivot(v[z].first, z);

    v.clear();
    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

}



template <class DType>
void MaxVariancePivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots);

}


#endif // MAXVARIANCEPIVOTS_H
