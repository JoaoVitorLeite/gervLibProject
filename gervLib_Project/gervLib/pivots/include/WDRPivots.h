#ifndef WDRPIVOTS_H
#define WDRPIVOTS_H

#include <Pivot.h>

template <class DType>
class WDRPivots : public Pivot<DType>
{

    private:
        size_t CandA = 300, ctop = CandA;
        size_t stop = 0, top = 0;
        double distsim1[200000];
        double distsim[400][200000];
        size_t s_p[3000];
        double alpha = 2;


    private:
        double objr(double alpha);

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

};



template <class DType>
double WDRPivots<DType>::objr(double alpha)
{

    double sum=0;
    for(int i=0;i<ctop;i++){
        double max=0,q,p;
        p=distsim1[i];
        for(int k=0;k<stop;k++){
            q=distsim[s_p[k]][i];
            if(q>max) max=q;
        }
        //	cout<<"this is "<<p<<" "<<max<<endl;
        sum+=p*pow(1-max/p,alpha)/CandA;
    }
    return sum;

}



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

    this->setNumberOfPivots(nPivots);

    Dataset<DType>* sample = nullptr;
    if(this->sample_size != -1.0)
        sample = dataset->sampleDataset(std::ceil(dataset->getCardinality()*this->sample_size), false, this->getSeed());
    else
        sample = dataset;

    size_t NUUM = sample->getCardinality()/2;
    size_t cstack[NUUM][2];
    size_t stack[NUUM], fa[NUUM];
    size_t num = NUUM;
    size_t pnum = nPivots;
    CandA = std::min(NUUM, CandA);
    ctop = CandA;

    top=0;
    for(int i=0;i<ctop;i++){
        cstack[i][0]=rand()%num;
        cstack[i][1]=rand()%num;
        if(cstack[i][0]==cstack[i][1]) i--;
    }
    for(int i=0;i<CandA;i++) fa[i]=-1;
    while (CandA>top){
        int maxi=rand()%num;
        stack[top++]=maxi;
        if(fa[maxi]==1) top--;
        fa[maxi]=1;
    }

    for(int i=0;i<ctop;i++)		{
        //	cout<<"firtst "<<dis(41,184)<<endl;
        distsim1[i]=df->getDistance(*sample->instance(cstack[i][0]), *sample->instance(cstack[i][1]));
        for(int j=0;j<top;j++)	distsim[j][i]=abs(df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][0])) - df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][1])));
        //		cout<<num<<" "<<cstack[i][0]<<" "<<cstack[i][1]<<"test " <<dis(cstack[i][0],cstack[i][1])<<endl;
    }

    stop=1;
    double min=10000000;
    int pos=-1;
    for(int j=0;j<top;j++)
    {
        s_p[stop-1]=j;
        double tt=objr(alpha);
        if(tt<min-0.00000001){
            pos=j;
            min=tt;
        }
    }
    s_p[0]=pos;
    stop=1;
    int i=50,fgg;
    while(i--){
        //cout<<"doing : "<<50-i<<" Cand: ";
        //                for(int j=0;j<stop;j++){
        //                    cout<<stack[s_p[j]]<<" ";
        //                }
        //cout<<endl;
        if(stop<pnum){
            //cout<<" ans= "<<objr(alpha)<<endl;
            stop++;
            min=10000000;
            pos=-1;
            for(int j=0;j<top;j++)
            {
                s_p[stop-1]=j;
                double tt=objr(alpha);
                if(tt<min-0.00000001){
                    pos=j;
                    min=tt;
                }
            }
            s_p[stop-1]=pos;
            if(pos==-1) break;
            continue;
        }
        min=objr(alpha);
        int f=0;
        //cout<<top<<" is top and stop "<< stop<<" min "<<min<<endl;
        for(int i=0;i<top;i++)
            for(int j=0;j<stop;j++){
                int q=s_p[j];
                s_p[j]=i;
                double tt=objr(alpha);
                if(tt<min-0.00000001){
                    f=1;
                    q=i;
                }
                s_p[j]=q;
            }
        if(f==0) break;
    }
    stop=pnum;

    stop=1;
    double recc[100];
    for(int j=0;j<pnum;j++){
        s_p[pnum]=s_p[0];
        s_p[0]=s_p[j];
        recc[j]=objr(alpha);
        s_p[0]=s_p[pnum];
    }
    for(int i=0;i<pnum;i++)
        for(int j=i;j<pnum;j++)
            if(recc[i]>recc[j]){
                recc[pnum]=recc[i];
                recc[i]=recc[j];
                recc[j]=recc[pnum];
                s_p[pnum]=s_p[i];
                s_p[i]=s_p[j];
                s_p[j]=s_p[pnum];
            }
    stop=pnum;
    //                for(int i=0;i<pnum;i++)
    //                    cout<<"rec pivot "<<i<<" rec "<<recc[i]<<endl;

//****************************************************************************

//    for(size_t i = 0; i < ctop; i++)
//    {

//        cstack[i][0] = rand()%NUUM;
//        cstack[i][1] = rand()%NUUM;
//        if(cstack[i][0] == cstack[i][1]) i--;

//    }

//    while(CandA > top)
//    {

//        size_t maxi = rand()%NUUM;
//        stack[top++] = maxi;
//        if(fa[maxi] == 1) top--;
//        fa[maxi] = 1;

//    }

//    for(size_t i = 0; i < ctop; i++)
//    {

//        distsim1[i] = df->getDistance(*sample->instance(cstack[i][0]), *sample->instance(cstack[i][1]));
//        for(size_t j = 0; j < top; j++) distsim[i][j] = abs(df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][0])) - df->getDistance(*sample->instance(stack[j]), *sample->instance(cstack[i][1])));

//    }

//    stop=1;
//    double min=10000000;
//    int pos=-1;

//    for(size_t j = 0; j < top; j++)
//    {

//        s_p[stop-1]=j;
//        double tt=objr(alpha);
//        //std::cout << "pos = " << j << " / objr = " << tt << "\n";
//        if(tt<min-0.00000001){
//            pos=j;
//            min=tt;
//        }

//    }

//    s_p[0]=pos;
//    stop=1;
//    int i=50,fgg;

//    while(i--){
//        std::cout<<"doing : "<<50-i<<" Cand: ";
//        for(int j=0;j<stop;j++){
//            std::cout<<stack[s_p[j]]<<" ";
//        }
//        std::cout<<std::endl;
//        if(stop<nPivots){
//            //cout<<" ans= "<<objr(alpha)<<endl;
//            objr(alpha);
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
////            bitmap[pos] = false;
//            if(pos==-1) break;
//            continue;
//        }
//        min=objr(alpha);
//        int f=0;
////        std::cout<<top<<" is top and stop "<< stop<<" min "<<min<<std::endl;
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
//    stop=nPivots;

//    stop=1;
//    double recc[100];
//    for(int j=0;j<nPivots;j++){
//        s_p[nPivots]=s_p[0];
//        s_p[0]=s_p[j];
//        recc[j]=objr(alpha);
//        s_p[0]=s_p[nPivots];
//    }
//    for(int i=0;i<nPivots;i++)
//        for(int j=i;j<nPivots;j++)
//            if(recc[i]>recc[j]){
//                recc[nPivots]=recc[i];
//                recc[i]=recc[j];
//                recc[j]=recc[nPivots];
//                s_p[nPivots]=s_p[i];
//                s_p[i]=s_p[j];
//                s_p[j]=s_p[nPivots];
//            }
//    stop=nPivots;
////    for(int i=0;i<nPivots;i++)
////        cout<<"rec pivot "<<i<<" rec "<<recc[i]<<endl;

    for(size_t x = 0; x < (this->getNumberOfPivots()); x++)
        this->setPivot(sample->instance(stack[s_p[x]]), x);

    if(this->sample_size == -1.0)
        sample = nullptr;

    delete sample;

}

template <class DType>
void WDRPivots<DType>::generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, size_t seed, std::vector<std::string> args)
{

    this->setSeed(seed);
    generatePivots(sample, df, nPivots, args);

}




#endif // WDRPIVOTS_H
