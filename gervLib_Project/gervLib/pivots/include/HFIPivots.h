#ifndef HFIPIVOTS_H
#define HFIPIVOTS_H

#include <Pivot.h>

template <class DType>
class HFIPivots : public Pivot<DType>
{


    private:
        double** maxPrunning(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t num, size_t num_cand, size_t* cand);

    public:
        HFIPivots();
        HFIPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::HFI);
            generatePivots(sample, function, nPivots);

        }
        HFIPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::HFI);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots);

        }

        ~HFIPivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

};



template <class DType>
double** HFIPivots<DType>::maxPrunning(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t num, size_t num_cand, size_t* cand)
{

    bool *indicator = new bool[num];

    for(size_t i = 0; i < num; i++)
    {

        indicator[i] = true;

    }

    size_t *idSet = new size_t[num_cand];

    double d = 0.0, t;
    size_t choose = 0;
    double **distMatrix = new double*[num];

    for(size_t i = 0; i < num; i++)
    {

        distMatrix[i] = new double[num_cand];

        for(size_t j = 0; j < num_cand; j++)
        {

            distMatrix[i][j] = 0.0;

        }

    }

    //Seleciona o primeiro pivo
    if(num_cand > 0)
    {

        for(size_t i = 1; i < num; i++)
        {

            t = df->getDistance(*sample->instance(i), *sample->instance(0));

            if(t > d)
            {

                d = t;
                choose = i;

            }

        }

        idSet[0] = choose;
        cand[0] = choose;
        indicator[choose] = false;

    }

    //Seleciona o segundo pivo
    if(num_cand > 1)
    {

        d = 0.0;

        for(size_t i = 0; i < num; i++)
        {

            if(indicator[i])
            {

                distMatrix[i][0] = df->getDistance(*sample->instance(cand[0]), *sample->instance(i));

                if(distMatrix[i][0] > d)
                {

                    d = distMatrix[i][0];
                    choose = i;

                }

            }

        }

        idSet[1] = choose;
        cand[1] = choose;
        indicator[choose] = false;

    }

    double edge = d;
    d = std::numeric_limits<double>::max();

    for(size_t i = 2; i < num_cand; i++)
    {

        d = std::numeric_limits<double>::max();

        for(size_t j = 0; j < num; j++)
        {

            if(indicator[j])
            {

                t = 0;
                for(size_t k = 0; k < i-1; k++)
                {

                    t += fabs(edge - distMatrix[j][k]);

                }

                distMatrix[j][i-1] = df->getDistance(*sample->instance(j), *sample->instance(cand[i-1]));
                t += fabs(edge - distMatrix[j][i-1]);

                if(t < d)
                {

                    d = t;
                    choose = j;

                }

            }

        }

        idSet[i] = choose;
        cand[i] = choose;
        indicator[choose] = false;

    }


    for(size_t i = 0; i < num; i++)
    {

        if(indicator[i])
        {

            distMatrix[i][num_cand - 1] = df->getDistance(*sample->instance(i), *sample->instance(cand[num_cand - 1]));

        }

    }

    for(size_t i = 0; i < num_cand; i++)
    {

        for(size_t j = i+1; j < num_cand; j++)
        {

            distMatrix[idSet[i]][j] = df->getDistance(*sample->instance(idSet[i]), *sample->instance(cand[j]));

        }

    }

    delete [] (indicator);
    delete [] (idSet);

    return distMatrix;

}


template <class DType>
HFIPivots<DType>::HFIPivots()
{

    this->setPivotType(PIVOT_TYPE::HFI);

}

template <class DType>
HFIPivots<DType>::~HFIPivots()
{



}

template <class DType>
void HFIPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t o_num = sample->getCardinality()/2, num_cand = 300, q_num = num_cand;
    size_t *cand = new size_t[num_cand];

    double** O_P_matrix = maxPrunning(sample, df, o_num, num_cand, cand);
    double** Q_O_matrix = new double*[q_num];
    double** Q_P_matrix = new double*[q_num];
    double** esti = new double*[q_num];

    for(size_t i = 0; i < q_num; i++)
    {

        Q_O_matrix[i] = new double[o_num];
        Q_P_matrix[i] = new double[num_cand];
        esti[i] = new double[o_num];

    }

    bool* indicator = new bool[num_cand];

    for(size_t i = 0; i < num_cand; i++)
    {

        indicator[i] = true;

    }

    for(size_t i = 0; i < q_num; i++)
    {

        for(size_t j = 0; j < o_num; j++)
        {

            Q_O_matrix[i][j] = df->getDistance(*sample->instance(i), *sample->instance(j));
            esti[i][j] = 0.0;

        }

        for(size_t j = 0; j < num_cand; j++)
        {

            Q_P_matrix[i][j] = df->getDistance(*sample->instance(i), *sample->instance(cand[j]));

        }

    }


    double d = 0.0, t = 0.0;
    size_t choose;
    size_t *pvtIndex = new size_t[nPivots];

    for(size_t i = 0; i < nPivots; i++)
    {

        choose = UINT_MAX;

        for(size_t j = 0; j < num_cand; j++)
        {

            if(indicator[j])
            {

                t = 0.0;

                for(size_t m = 0; m < q_num; m++)
                {

                    for(size_t n = 0; n < o_num; n++)
                    {

                        if(Q_O_matrix[m][n] != 0.0)
                        {

                            t += (std::max(fabs(Q_P_matrix[m][j] - O_P_matrix[n][j]), esti[m][n]))/Q_O_matrix[m][n];

                        }

                    }

                }

                t = t / (q_num*o_num);

                if(t > d)
                {

                    d = t;
                    choose = j;

                }

            }

        }

        if(choose == UINT_MAX)
            break;

        indicator[choose] = false;
        pvtIndex[i] = cand[choose];

        for(size_t m = 0; m < q_num; m++)
        {

            for(size_t n = 0; n < o_num; n++)
            {

                esti[m][n] = std::max(fabs(Q_P_matrix[m][choose] - O_P_matrix[n][choose]), esti[m][n]);

            }

        }

    }

    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
        this->setPivot(sample->instance(pvtIndex[x]), x);

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

}


template <class DType>
void HFIPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}








#endif // HFIPIVOTS_H
