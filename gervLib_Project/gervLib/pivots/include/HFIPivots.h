#ifndef HFIPIVOTS_H
#define HFIPIVOTS_H

#include <Pivot.h>
#include <ConvexPivots.h>

template <class DType>
class HFIPivots : public Pivot<DType>
{


//    private:
//        double** maxPrunning(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t num_cand, size_t* cand);

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



//template <class DType>
//double** HFIPivots<DType>::maxPrunning(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t num_cand, size_t* cand)
//{

//    bool* bitmap = new bool[sample->getCardinality()];
//    double** distMatrix = new double*[sample->getCardinality()];
//    size_t pos = 0;
//    double max = std::numeric_limits<double>::min(), min = std::numeric_limits<double>::max(), d = 0.0, edge = 0.0;

//    for(size_t x = 0; x < sample->getCardinality(); x++)
//    {

//        bitmap[x] = true;

//    }

//    for(size_t x = 0; x < sample->getCardinality(); x++)
//    {

//        distMatrix[x] = new double[num_cand];

//        for(size_t y = 0; y < num_cand; y++)
//        {

//            distMatrix[x][y] = 0.0;

//        }

//    }

//    //Primeiro pivo
//    for(size_t x = 1; x < sample->getCardinality(); x++)
//    {

//        d = df->getDistance(*sample->instance(x), *sample->instance(0));

//        if(d > max)
//        {

//            max = d;
//            pos = x;

//        }

//    }

//    cand[0] = pos;
//    bitmap[pos] = false;
//    //

//    //Segundo pivo
//    max = std::numeric_limits<double>::min();
//    pos = 0;

//    for(size_t x = 0; x < sample->getCardinality(); x++)
//    {

//        if(bitmap[x])
//        {

//            distMatrix[x][0] = df->getDistance(*sample->instance(x), *sample->instance(cand[0]));

//            if(distMatrix[x][0] > max)
//            {

//                max = distMatrix[x][0];
//                pos = x;

//            }

//        }

//    }

//    cand[1] = pos;
//    bitmap[pos] = false;
//    //

//    edge = d;

//    for(size_t x = 2; x < num_cand; x++)
//    {

//        min = std::numeric_limits<double>::max();
//        pos = 0;

//        for(size_t y = 0; y < sample->getCardinality(); y++)
//        {

//            if(bitmap[y])
//            {

//                d = 0.0;

//                for(size_t z = 0; z < x-1; z++)
//                {

//                    d += fabs(edge - distMatrix[y][z]);

//                }

//                distMatrix[y][x-1] = df->getDistance(*sample->instance(y), *sample->instance(cand[x-1]));
//                d += fabs(edge - distMatrix[y][x-1]);

//                if(d < min)
//                {

//                    min = d;
//                    pos = y;

//                }

//            }

//        }

//        cand[x] = pos;
//        bitmap[pos] = false;

//    }

//    for(size_t x = 0; x < sample->getCardinality(); x++)
//    {

//        if(bitmap[x])
//        {

//            distMatrix[x][num_cand-1] = df->getDistance(*sample->instance(x), *sample->instance(cand[num_cand-1]));

//        }

//    }

//    for(size_t x = 0; x < num_cand; x++)
//    {

//        for(size_t y = x+1; y < num_cand; y++)
//        {

//            distMatrix[cand[x]][y] = df->getDistance(*sample->instance(cand[x]), *sample->instance(cand[y]));

//        }

//    }

//    delete [] bitmap;

//    return distMatrix;

//    bool *indicator = new bool[num];

//    for(size_t i = 0; i < num; i++)
//    {

//        indicator[i] = true;

//    }

//    size_t *idSet = new size_t[num_cand];

//    double d = 0.0, t;
//    size_t choose = 0;
//    double **distMatrix = new double*[num];

//    for(size_t i = 0; i < num; i++)
//    {

//        distMatrix[i] = new double[num_cand];

//        for(size_t j = 0; j < num_cand; j++)
//        {

//            distMatrix[i][j] = 0.0;

//        }

//    }

//    //Seleciona o primeiro pivo
//    if(num_cand > 0)
//    {

//        for(size_t i = 1; i < num; i++)
//        {

//            t = df->getDistance(*sample->instance(i), *sample->instance(0));

//            if(t > d)
//            {

//                d = t;
//                choose = i;

//            }

//        }

//        idSet[0] = choose;
//        cand[0] = choose;
//        indicator[choose] = false;

//    }

//    //Seleciona o segundo pivo
//    if(num_cand > 1)
//    {

//        d = 0.0;

//        for(size_t i = 0; i < num; i++)
//        {

//            if(indicator[i])
//            {

//                distMatrix[i][0] = df->getDistance(*sample->instance(cand[0]), *sample->instance(i));

//                if(distMatrix[i][0] > d)
//                {

//                    d = distMatrix[i][0];
//                    choose = i;

//                }

//            }

//        }

//        idSet[1] = choose;
//        cand[1] = choose;
//        indicator[choose] = false;

//    }

//    double edge = d;
//    d = std::numeric_limits<double>::max();

//    for(size_t i = 2; i < num_cand; i++)
//    {

//        d = std::numeric_limits<double>::max();

//        for(size_t j = 0; j < num; j++)
//        {

//            if(indicator[j])
//            {

//                t = 0;
//                for(size_t k = 0; k < i-1; k++)
//                {

//                    t += fabs(edge - distMatrix[j][k]);

//                }

//                distMatrix[j][i-1] = df->getDistance(*sample->instance(j), *sample->instance(cand[i-1]));
//                t += fabs(edge - distMatrix[j][i-1]);

//                if(t < d)
//                {

//                    d = t;
//                    choose = j;

//                }

//            }

//        }

//        idSet[i] = choose;
//        cand[i] = choose;
//        indicator[choose] = false;

//    }


//    for(size_t i = 0; i < num; i++)
//    {

//        if(indicator[i])
//        {

//            distMatrix[i][num_cand - 1] = df->getDistance(*sample->instance(i), *sample->instance(cand[num_cand - 1]));

//        }

//    }

//    for(size_t i = 0; i < num_cand; i++)
//    {

//        for(size_t j = i+1; j < num_cand; j++)
//        {

//            distMatrix[idSet[i]][j] = df->getDistance(*sample->instance(cand[i]), *sample->instance(cand[j]));

//        }

//    }

//    for(size_t i = 0; i < num; i++)
//    {

//        for(size_t j = 0; j < num_cand; j++)
//        {

//            std::cout << distMatrix[i][j] << " ";

//        }
//        std::cout << "\n";

//    }

//    delete [] (indicator);
//    delete [] (idSet);

//    return distMatrix;

//}


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

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t pvtSize = std::min(2 * this->getNumberOfPivots(), sample->getCardinality());
    size_t cand = pvtSize;
    size_t pairSize = std::min((size_t)300, (size_t)sample->getCardinality()/2);
    size_t* pivot_index = new size_t[pvtSize];
    size_t** pairs_index = new size_t*[pairSize];
    size_t* aux;
    bool* bitmap = new bool[pvtSize];
    double* dist_pairs = new double[pairSize];
    double** dist_pivots_pairs = new double*[pvtSize];

    std::vector<BasicArrayObject<DType>> vec;
    size_t id = 0;

    for(size_t i = 0; i < sample->getCardinality(); i++)
    {

        vec.push_back(sample->getFeatureVector(i));
        vec[i].setOID(id++);

    }

    Dataset<DType> *data_aux = new Dataset<DType>(vec, vec.size(), vec[0].size());
    ConvexPivots<DType>* convex = new ConvexPivots<DType>();
    convex->setSeed(this->getSeed());
    convex->generatePivots(data_aux, df, pvtSize);

    for(size_t i = 0; i < pvtSize; i++)
    {

        pivot_index[i] = convex->get(i).getOID();
        //std::cout << "P ID = " << pivot_index[i] << std::endl;
        dist_pivots_pairs[i] = new double[pairSize];
        bitmap[i] = false;

    }

    delete convex;
    vec.clear();
    delete data_aux;

    for(size_t i = 0; i < pairSize; i++)
    {

        aux = uniqueRandomNumber(0, sample->getCardinality(), 2, (this->getSeed()*(i+1))%RAND_MAX);
        pairs_index[i] = new size_t[2];
        pairs_index[i][0] = aux[0];
        pairs_index[i][1] = aux[1];
        dist_pairs[i] = df->getDistance(*sample->instance(pairs_index[i][0]), *sample->instance(pairs_index[i][1]));

        for(size_t j = 0; j < pvtSize; j++)
        {

            dist_pivots_pairs[j][i] = fabs(df->getDistance(*sample->instance(pivot_index[j]), *sample->instance(pairs_index[i][0])) -
                                           df->getDistance(*sample->instance(pivot_index[j]), *sample->instance(pairs_index[i][1])));

        }

        delete [] aux;

    }

//    for(size_t i = 0; i < pairSize; i++)
//    {

//        std::cout << "CAND ID = " << pairs_index[i][0] << "/" << pairs_index[i][1] << std::endl;

//    }

//    for(size_t i = 0; i < pairSize; i++)
//    {

//        std::cout << dist_pairs[i] << "\t";

//    }

//    std::cout << std::endl << std::endl;

//    for(size_t i = 0; i < pvtSize; i++)
//    {

//        for(size_t j = 0; j < pairSize; j++)
//        {

//            std::cout << dist_pivots_pairs[i][j] << "\t";

//        }

//        std::cout << std::endl;

//    }

    size_t pvtIndex = 0;
    while(cand > this->getNumberOfPivots())
    {

        double max = std::numeric_limits<double>::min();
        size_t pos = 0;

        for(size_t x = 0; x < pvtSize; x++)
        {

            if(!bitmap[x])
            {

                bitmap[x] = true;

//                std::cout << x << std::endl;
//                for(size_t z = 0; z < pvtSize; z++) std::cout << bitmap[z] << "\t";
//                std::cout << "\n";

                double prec = 0.0;

                for(size_t i = 0; i < pairSize; i++)
                {

                    double max_p = std::numeric_limits<double>::min();
                    size_t max_pos = 0;

                    for(size_t j = 0; j < pvtSize; j++)
                    {

                        if((bitmap[j] == true) && (dist_pivots_pairs[j][i] > max_p))
                        {

                            max_p = dist_pivots_pairs[j][i];
                            max_pos = j;

                        }

                    }

                    prec += dist_pivots_pairs[max_pos][i]/dist_pairs[i];

                }

                prec /= pairSize;

//                std::cout << "PREC " << pivot_index[x] << " : " << prec << std::endl;

                if(prec > max)
                {

                    max = prec;
                    pos = x;

                }

                bitmap[x] = false;

            }

        }

        bitmap[pos] = true;
        this->setPivot(sample->instance(pivot_index[pos]), pvtIndex++);
        cand--;

    }

//    size_t j = 0;
//    for(size_t i = 0; i < pvtSize; i++)
//    {

//        if(bitmap[i])
//        {

//            this->setPivot(sample->instance(pivot_index[i]), j++);

//        }

//    }

    delete [] bitmap;
    delete [] pivot_index;
    delete [] dist_pairs;

    for(size_t i = 0; i < pairSize; i++)
    {

        delete [] pairs_index[i];

    }

    for(size_t i = 0; i < pvtSize; i++)
    {

        delete [] dist_pivots_pairs[i];

    }

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

    auto end = std::chrono::steady_clock::now();
    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());


//*******************************************************************************************************************************
//    auto start = std::chrono::steady_clock::now();

//    this->setNumberOfPivots(nPivots);

//    Dataset<DType>* sample = nullptr;
//    if(this->sample_size != -1.0)
//        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
//    else
//        sample = dataset;

////    size_t num_cand = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2)), pos = 0;
//    size_t num_cand = std::min(2 * this->getNumberOfPivots(), sample->getCardinality()), pos = 0;
//    size_t* cand = new size_t[num_cand];
//    double** O_P_matrix = maxPrunning(sample, df, num_cand, cand);
//    double** Q_O_matrix = new double*[num_cand];
//    double** Q_P_matrix = new double*[num_cand];
//    double** esti = new double*[num_cand];
//    bool* bitmap = new bool[num_cand];
//    double max = std::numeric_limits<double>::min(), d = 0.0;

//    for(size_t x = 0; x < num_cand; x++)
//    {

//        Q_O_matrix[x] = new double[sample->getCardinality()];
//        Q_P_matrix[x] = new double[num_cand];
//        esti[x] = new double[sample->getCardinality()];
//        bitmap[x] = true;

//    }

//    for(size_t x = 0; x < num_cand; x++)
//    {

//        for(size_t y = 0; y < sample->getCardinality(); y++)
//        {

//            Q_O_matrix[x][y] = df->getDistance(*sample->instance(x), *sample->instance(y));
//            esti[x][y] = 0.0;

//        }

//        for(size_t z = 0; z < num_cand; z++)
//        {

//            Q_P_matrix[x][z] = df->getDistance(*sample->instance(x), *sample->instance(cand[z]));

//        }

//    }

//    for(size_t i = 0; i < nPivots; i++)
//    {

//        pos = UINT_MAX;
//        for(size_t j = 0; j < num_cand; j++)
//        {

//            if(bitmap[j])
//            {

//                d = 0.0;
//                for(size_t m = 0; m < num_cand; m++)
//                {

//                    for(size_t n = 0; n < sample->getCardinality(); n++)
//                    {

//                        if(Q_O_matrix[m][n] != 0.0)
//                        {

//                            d += (std::max(fabs(Q_P_matrix[m][j] - O_P_matrix[n][j]), esti[m][n]))/Q_O_matrix[m][n];

//                        }

//                    }

//                }

//                d = d/(num_cand*sample->getCardinality());

//                if(d > max)
//                {

//                    max = d;
//                    pos = j;


//                }

//            }

//        }

//        if(pos == UINT_MAX)
//            break;

//        bitmap[pos] = false;
//        this->setPivot(sample->instance(cand[pos]), i);

//        for(size_t m = 0; m < num_cand; m++)
//        {

//            for(size_t n = 0; n < sample->getCardinality(); n++)
//            {

//                esti[m][n] = std::max(fabs(Q_P_matrix[m][pos] - O_P_matrix[n][pos]), esti[m][n]);

//            }

//        }

//    }

//    for(size_t x = 0; x < sample->getCardinality(); x++)
//    {

//        delete [] O_P_matrix[x];

//    }

//    delete [] O_P_matrix;

//    if(this->sample_size == -1.0)
//        sample = nullptr;

//    delete sample;
//    delete [] cand;
//    delete [] bitmap;

//    for(size_t x = 0; x < num_cand; x++)
//    {

//        delete [] Q_O_matrix[x];
//        delete [] Q_P_matrix[x];
//        delete [] esti[x];

//    }

//    delete [] Q_O_matrix;
//    delete [] Q_P_matrix;
//    delete [] esti;

//    auto end = std::chrono::steady_clock::now();
//    this->setElapsedTime(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

//*******************************************************************************************************************************
//    size_t o_num = sample->getCardinality()/2, num_cand = std::min(o_num, (size_t)300), q_num = num_cand;
//    size_t *cand = new size_t[num_cand];

//    double** O_P_matrix = maxPrunning(sample, df, o_num, num_cand, cand);
//    double** Q_O_matrix = new double*[q_num];
//    double** Q_P_matrix = new double*[q_num];
//    double** esti = new double*[q_num];

//    for(size_t i = 0; i < q_num; i++)
//    {

//        Q_O_matrix[i] = new double[o_num];
//        Q_P_matrix[i] = new double[num_cand];
//        esti[i] = new double[o_num];

//    }

//    bool* indicator = new bool[num_cand];

//    for(size_t i = 0; i < num_cand; i++)
//    {

//        indicator[i] = true;

//    }

//    for(size_t i = 0; i < q_num; i++)
//    {

//        for(size_t j = 0; j < o_num; j++)
//        {

//            Q_O_matrix[i][j] = df->getDistance(*sample->instance(i), *sample->instance(j));
//            esti[i][j] = 0.0;

//        }

//        for(size_t j = 0; j < num_cand; j++)
//        {

//            Q_P_matrix[i][j] = df->getDistance(*sample->instance(i), *sample->instance(cand[j]));

//        }

//    }


//    double d = 0.0, t = 0.0;
//    size_t choose;
//    size_t *pvtIndex = new size_t[nPivots];

//    for(size_t i = 0; i < nPivots; i++)
//    {

//        choose = UINT_MAX;

//        for(size_t j = 0; j < num_cand; j++)
//        {

//            if(indicator[j])
//            {

//                t = 0.0;

//                for(size_t m = 0; m < q_num; m++)
//                {

//                    for(size_t n = 0; n < o_num; n++)
//                    {

//                        if(Q_O_matrix[m][n] != 0.0)
//                        {

//                            t += (std::max(fabs(Q_P_matrix[m][j] - O_P_matrix[n][j]), esti[m][n]))/Q_O_matrix[m][n];

//                        }

//                    }

//                }

//                t = t / (q_num*o_num);

//                if(t > d)
//                {

//                    d = t;
//                    choose = j;

//                }

//            }

//        }

//        if(choose == UINT_MAX)
//            break;

//        indicator[choose] = false;
//        pvtIndex[i] = cand[choose];

//        for(size_t m = 0; m < q_num; m++)
//        {

//            for(size_t n = 0; n < o_num; n++)
//            {

//                esti[m][n] = std::max(fabs(Q_P_matrix[m][choose] - O_P_matrix[n][choose]), esti[m][n]);

//            }

//        }

//    }

//    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
//        this->setPivot(sample->instance(pvtIndex[x]), x);

//    if(this->sample_size == -1.0)
//        sample = nullptr;

//    delete sample;

}


template <class DType>
void HFIPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}








#endif // HFIPIVOTS_H
