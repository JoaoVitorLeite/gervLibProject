#ifndef WDRPIVOTS_H
#define WDRPIVOTS_H

#include <Pivot.h>
#include <map>


template <class DType>
class WDRPivots : public Pivot<DType>
{

    private:
//        size_t CandA = 300, ctop = CandA;
//        size_t stop = 0, top = 0;
//        double distsim1[200000];
//        double distsim[400][200000];
//        size_t s_p[3000];
//        double alpha = 2;
//*************************************************************
//        double* distsim1; //distancia entre os pares
//        double** distsim; //distancia dos pares para os pivos
//        size_t* s_p;
//        double alpha = 2;
//*************************************************************
        double* dist_pair;
        double** dist_pair_cand;
        size_t* pivot_index;
        bool* bitmap;
        double alpha = 2;

    private:
//        double objr(size_t ctop, size_t stop, size_t num_cand);
        double func(size_t cand_size, size_t pair_size, size_t p_num);

    public:
        WDRPivots();

        WDRPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::WDR);
            generatePivots(sample, function, nPivots);

        }

        WDRPivots(Dataset<DType>* sample, DistanceFunction<BasicArrayObject<DType>>* function, size_t nPivots, size_t seed) : Pivot<DType>(){

            this->setPivotType(PIVOT_TYPE::WDR);
            this->setSeed(seed);
            generatePivots(sample, function, nPivots);

        }

        ~WDRPivots();

        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>());
        void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args = std::vector<std::string>());

        //void setAlpha(size_t alpha_);

};


template <class DType>
double WDRPivots<DType>::func(size_t cand_size, size_t pair_size, size_t p_num)
{

    double sum = 0, _dist_pair, _dist_pair_pivot, max = std::numeric_limits<double>::min();

    for(size_t i = 0; i < pair_size; i++)
    {

        _dist_pair = dist_pair[i];

        for(size_t j = 0; j < p_num; j++)
        {

            _dist_pair_pivot = dist_pair_cand[i][pivot_index[j]];

            if(_dist_pair_pivot > max)
            {

                max = _dist_pair_pivot;

            }

        }

        sum += _dist_pair * pow(1-max/_dist_pair, alpha)/cand_size;

    }

    return sum;

}


//template <class DType>
//double WDRPivots<DType>::objr(size_t ctop, size_t stop, size_t num_cand)
//{


////**********************************************************************************************
////    double sum = 0.0;

////    for(size_t i = 0; i < ctop; i++)
////    {

////        double max = 0.0, q, p;
////        p=distsim1[i];

////        for(size_t k = 0; k < stop; k++)
////        {

////            q = distsim[s_p[k]][i];

////            if(q > max)
////            {

////                max = q;

////            }

////        }

//////        std::cout << "MAX = " << max << std::endl;

////        sum += p * pow(1-max/p, alpha)/num_cand;

////        //std::cout << "SUM = " << sum << std::endl;

////    }

////    //std::cout << "SUM = " << sum << std::endl;

////    return sum;
////*************************************************************************************************

////    double sum=0;
////    for(int i=0;i<ctop;i++){
////        double max=0,q,p;
////        p=distsim1[i];
////        for(int k=0;k<stop;k++){
////            q=distsim[s_p[k]][i];
////            if(q>max) max=q;
////        }                    bitmap[i] = true;

////        //	cout<<"this is "<<p<<" "<<max<<endl;
////        sum+=p*pow(1-max/p,alpha)/CandA;
////    }
////    return sum;

//}



template <class DType>
WDRPivots<DType>::WDRPivots()
{

    this->setPivotType(PIVOT_TYPE::WDR);

}

template <class DType>
WDRPivots<DType>::~WDRPivots()
{



}




template <class DType>
void WDRPivots<DType>::generatePivots(Dataset<DType> *dataset, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args)
{

    auto start = std::chrono::steady_clock::now();

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;


//    size_t cand_size = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2)), pair_size = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2));
    size_t cand_size = std::min(2 * this->getNumberOfPivots(), sample->getCardinality());;
    //size_t pair_size = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2));
    size_t pair_size = std::min((size_t)300, (size_t)sample->getCardinality()/2);
    size_t *cand = uniqueRandomNumber(0, sample->getCardinality(), cand_size, this->getSeed());
    size_t **pair = new size_t*[pair_size];
    size_t *aux;
    bitmap = new bool[cand_size];
    pivot_index = new size_t[this->getNumberOfPivots()];

    for(size_t i = 0; i < cand_size; i++)
    {

        bitmap[i] = true;

    }

    for(size_t i = 0; i < pair_size; i++)
    {

        aux = uniqueRandomNumber(0, sample->getCardinality(), 2, (this->getSeed()*(i+1))%RAND_MAX);
        pair[i] = new size_t[2];
        pair[i][0] = aux[0];
        pair[i][1] = aux[1];
        delete [] aux;

    }

    dist_pair = new double[pair_size];
    dist_pair_cand = new double*[pair_size];

    for(size_t i = 0; i < pair_size; i++)
    {

        dist_pair_cand[i] = new double[cand_size];
        dist_pair[i] = df->getDistance(*sample->getInstance(pair[i][0]), *sample->getInstance(pair[i][1]));

        for(size_t j = 0; j < cand_size; j++)
        {

            dist_pair_cand[i][j] = fabs(df->getDistance(*sample->getInstance(cand[j]), *sample->getInstance(pair[i][0])) - df->getDistance(*sample->getInstance(cand[j]), *sample->getInstance(pair[i][1])));

        }

    }

    size_t p_num = 0;
    double min = std::numeric_limits<double>::max(), f;
    size_t pos = 0;
    size_t iter = std::max((size_t)500, this->getNumberOfPivots());
    bool flag = true;

    while(iter--)
    {

        //std::cout << "ITER = " << iter << std::endl;

        flag = true;

        if(p_num <= this->getNumberOfPivots())
        {

            min = std::numeric_limits<double>::max();

            for(size_t i = 0; i < cand_size; i++)
            {

                if(bitmap[i])
                {

                    pivot_index[p_num] = i;
                    f = func(cand_size, pair_size, p_num+1);

                    if(f < min)
                    {

                        min = f;
                        pos = i;

                    }

                }

            }

            pivot_index[p_num] = pos;
            bitmap[pos] = false;
            p_num++;

        }

        double temp = func(cand_size, pair_size, p_num);

        for(size_t i = 0; i < p_num; i++)
        {

            pos = pivot_index[i];

            for(size_t j = 0; j < cand_size; j++)
            {

                if(bitmap[j])
                {

                    pivot_index[i] = j;
                    f = func(cand_size, pair_size, p_num);

                    if(temp > f)
                    {

                        bitmap[i] = false;
                        temp = f;
                        pos = j;
                        flag = false;

                    }

                }

            }

            pivot_index[i] = pos;

        }

        if(flag && p_num == this->getNumberOfPivots())
        {

            break;

        }

    }

    for(size_t i = 0; i < p_num; i++)
    {

//        std::cout << cand[pivot_index[i]] << std::endl;
        this->setPivot(sample->getInstance(cand[pivot_index[i]]), i);

    }

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

    delete [] cand;

    for(size_t i = 0; i < pair_size; i++)
    {

        delete [] pair[i];
        delete [] dist_pair_cand[i];

    }

    delete [] pair;
    delete [] bitmap;
    delete [] pivot_index;
    delete [] dist_pair;
    delete [] dist_pair_cand;

    aux = nullptr;
    delete aux;

    auto end = std::chrono::steady_clock::now();
    this->elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();


//***************************************************************************************************************************************************************************************************
//    auto start = std::chrono::steady_clock::now();

//    this->setNumberOfPivots(nPivots);

//    Dataset<DType>* sample = nullptr;
//    if(this->sample_size != -1.0)
//        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
//    else
//        sample = dataset;

//    size_t top = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2)), ctop = std::max(nPivots, (size_t)std::ceil(sample->getCardinality()/2)), stop = 1;
//    double min = std::numeric_limits<double>::max();
//    size_t** cstack = new size_t*[ctop];
//    size_t* aux;
//    size_t* stack = uniqueRandomNumber(0, sample->getCardinality(), top, this->getSeed());
//    double* recc = new double[nPivots+1];

//    for(size_t x = 0; x < ctop; x++)
//    {

//        cstack[x] = new size_t[2];
//        aux = uniqueRandomNumber(0, sample->getCardinality(), 2, this->getSeed());
//        cstack[x][0] = aux[0];
//        cstack[x][1] = aux[1];

//        delete [] aux;


//    }

//    distsim1 = new double[ctop];
//    distsim = new double*[top];

//    for(size_t x = 0; x < top; x++)
//    {

//        distsim[x] = new double[ctop];

//    }1

//    s_p = new size_t[nPivots+1];

//    for(size_t i = 0; i < ctop; i++)
//    {

//        distsim1[i] = df->getDistance(*sample->instance(cstack[i][0]), *sample->instance(cstack[i][1]));

//        for(size_t j = 0; j < top; j++)
//        {

//            distsim[j][i] = abs(df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][0])) - df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][1])));

//        }

//    }

//    stop = 1;
//    min = std::numeric_limits<double>::max();
//    size_t pos = 0;

//    for(size_t j = 0; j < top; j++)
//    {

//        s_p[stop-1] = j;
//        double tt = objr(ctop, stop, top);

//        //std::cout << "MIN PVT 1 = " << tt << " / " << j << std::endl;

//        if(tt < min)
//        {

//            pos = j;
//            min = tt;

//        }

//    }

//    s_p[0] = pos;
//    //std::cout << "p1 = " << pos << "\n";
//    stop = 1;
//    size_t iter = 50;
//    bool* bitmap = new bool[top];
//    for(size_t i = 0; i < top; i++)
//        bitmap[i] = true;

//    while(iter--)
//    {

//        if(stop < nPivots)
//        {

//            stop++;
//            min = std::numeric_limits<double>::max();
//            pos = UINT_MAX;

//            for(size_t j = 0; j < top; j++)
//            {

//                s_p[stop-1] = j;
//                double tt = objr(ctop, stop, top);
//                //std::cout << "MIN PVT=  " << tt << " / " << j << std::endl;

//                if(tt < min && tt != min && bitmap[j])
//                {

//                    pos = j;
//                    min = tt;

//                }

//            }

//            s_p[stop-1] = pos;
//            bitmap[pos] = false;
//            //std::cout << "pi = " << pos << std::endl;

//            if(pos == UINT_MAX)
//                break;

//            continue;

//        }

//        //min = objr(ctop, stop, top);
//        bool flag = false;

//        for(size_t i = 0; i < top; i++)
//        {

//            for(size_t j = 0; j < stop; j++)
//            {

//                size_t q = s_p[j];
//                s_p[j] = i;
//                double tt = objr(ctop, stop, top);

//                if(tt < min && tt != min)
//                {

//                    //std::cout << "OI\n";
//                    flag = true;
//                    q = i;

//                }

//                s_p[j] = q;

//            }

//        }
    //        double* distsim1; //distancia entre os pares
    //        double** distsim; //distancia dos pares para os pivos
    //        size_t* s_p;
    //        double alpha = 2;
//        if(!flag)
//            break;

//    }

//    stop = nPivots;
//    stop = 1;

//    for(size_t j = 0; j < nPivots; j++)
//    {

//        s_p[nPivots] = s_p[0];
//        s_p[0] = s_p[j];
//        recc[j] = objr(ctop, stop, top);
//        s_p[0] = s_p[nPivots];

//    }

//    for(size_t i = 0; i < nPivots; i++)
//    {

//        for(size_t j = i; j < nPivots; j++)
//        {

//            if(recc[i] > recc[j])
//            {

//                recc[nPivots] = recc[i];
//                recc[i] = recc[j];
//                recc[j] = recc[nPivots];
//                s_p[nPivots] = s_p[i];
//                s_p[i] = s_p[j];
//                s_p[j] = s_p[nPivots];

//            }

//        }

//    }

//    stop = nPivots;

//    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
//        this->setPivot(sample->instance(stack[s_p[x]]), x);

//    delete [] distsim1;

//    for(size_t x = 0; x < top; x++)
//    {

//        delete [] distsim[x];

//    }

//    delete [] distsim;
//    delete [] s_p;
//    delete [] bitmap;

//    for(size_t x = 0; x < ctop; x++)
//    {

//        delete [] cstack[x];

//    }

//    delete [] cstack;
//    delete [] stack;
//    delete [] recc;

//    aux = nullptr;
//    delete aux;

//    if(this->sample_size == -1.0)
//        sample = nullptr;

//    delete sample;

//    auto end = std::chrono::steady_clock::now();
//    this->elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

//**********************************************************************************************************************************
//    size_t NUUM = sample->getCardinality()/2;
//    size_t cstack[NUUM][2];
//    size_t stack[NUUM], fa[NUUM];
//    size_t num = NUUM;
//    size_t pnum = nPivots;
//    CandA = std::min(NUUM, CandA);
//    ctop = CandA;

//    top=0;
//    for(int i=0;i<ctop;i++){
//        cstack[i][0]=rand()%num;
//        cstack[i][1]=rand()%num;
//        if(cstack[i][0]==cstack[i][1]) i--;
//    }
//    for(int i=0;i<CandA;i++) fa[i]=-1;
//    while (CandA>top){
//        int maxi=rand()%num;
//        stack[top++]=maxi;
//        if(fa[maxi]==1) top--;
//        fa[maxi]=1;
//    }

//    for(int i=0;i<ctop;i++)		{
//        //	cout<<"firtst "<<dis(41,184)<<endl;
//        distsim1[i]=df->getDistance(*sample->instance(cstack[i][0]), *sample->instance(cstack[i][1]));
//        for(int j=0;j<top;j++)	distsim[j][i]=abs(df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][0])) - df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][1])));
//        //		cout<<num<<" "<<cstack[i][0]<<" "<<cstack[i][1]<<"test " <<dis(cstack[i][0],cstack[i][1])<<endl;
//    }

//    stop=1;
//    double min=10000000;
//    int pos=-1;
//    for(int j=0;j<top;j++)
//    {
//        s_p[stop-1]=j;
//        double tt=objr(alpha);
//        if(tt<min-0.00000001){
//            pos=j;
//            min=tt;
//        }
//    }
//    s_p[0]=pos;
//    stop=1;
//    int i=50,fgg;
//    while(i--){
//        //cout<<"doing : "<<50-i<<" Cand: ";
//        //                for(int j=0;j<stop;j++){
//        //                    cout<<stack[s_p[j]]<<" ";
//        //                }
//        //cout<<endl;
//        if(stop<pnum){
//            //cout<<" ans= "<<objr(alpha)<<endl;
//            stop++;
//            min=10000000;
//            pos=-1;
//            for(int j=0;j<top;j++)
//            {
//                s_p[stop-1]=j;
//                double tt=objr(alpha);
//                if(tt<min-0.00000001){
//                    pos=j;
//                    min=tt;
//                }
//            }
//            s_p[stop-1]=pos;
//            if(pos==-1) break;
//            continue;
//        }
//        min=objr(alpha);
//        int f=0;
//        //cout<<top<<" is top and stop "<< stop<<" min "<<min<<endl;
//        for(int i=0;i<top;i++)
//            for(int j=0;j<stop;j++){
//                int q=s_p[j];
//                s_p[j]=i;
//                double tt=objr(alpha);
//                if(tt<min-0.00000001){
//                    f=1;
//                    q=i;
//                }
//                s_p[j]=q;
//            }
//        if(f==0) break;
//    }
//    stop=pnum;

//    stop=1;
//    double recc[100];
//    for(int j=0;j<pnum;j++){
//        s_p[pnum]=s_p[0];
//        s_p[0]=s_p[j];
//        recc[j]=objr(alpha);
//        s_p[0]=s_p[pnum];
//    }
//    for(int i=0;i<pnum;i++)
//        for(int j=i;j<pnum;j++)
//            if(recc[i]>recc[j]){
//                recc[pnum]=recc[i];
//                recc[i]=recc[j];
//                recc[j]=recc[pnum];
//                s_p[pnum]=s_p[i];
//                s_p[i]=s_p[j];
//                s_p[j]=s_p[pnum];
//            }
//    stop=pnum;

//    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
//        this->setPivot(sample->instance(stack[s_p[x]]), x);

//    if(this->sample_size == -1.0)
//        sample = nullptr;

//    delete sample;

}

template <class DType>
void WDRPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}




#endif // WDRPIVOTS_H
