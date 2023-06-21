#ifndef MVPTREE_H
#define MVPTREE_H

#include <Dataset.h>
#include <Hermes.h>
#include <Pivots.h>
#include <MVPData.h>
#include <MVPNode.h>
#include <queue>
#include <config_spb.h>
#include <MemoryManagerUtils.h>

template <class T>
struct MVPPartition
{

    MVPNode<T>* node;
    double min, max;

    MVPPartition()
    {



    }

    MVPPartition(MVPNode<T>* _node, double _min, double _max)
    {

        node = _node;
        min = _min;
        max = _max;

    }

};

template <class T>
struct CompareMVPPartition
{

    bool operator()(MVPPartition<T> const& a, MVPPartition<T> const& b)
    {

        bool ans;

        if(a.min != b.min)
        {

            ans = a.min > b.min;

        }
        else
        {

            ans = a.max > b.max;

        }

        return ans;

    }

};

template <class T>
struct KnnEntryMVP
{

    BasicArrayObject<T> element;
    double distance;

    KnnEntryMVP(BasicArrayObject<T> &element_, double distance_)
    {

        element = element_;
        distance = distance_;

    }

    KnnEntryMVP()
    {



    }

    bool operator<(const KnnEntryMVP<T>& item) const
    {

        return distance < item.distance;

    }

    bool operator>(const KnnEntryMVP<T>& item) const
    {

        return distance > item.distance;

    }

};

template<class Type, class Comp>
std::vector<Type> dequeueInOrderMVP_Results(std::priority_queue<Type, std::vector<Type>, Comp> pq)
{

    std::vector<Type> ans;
    std::priority_queue<Type, std::vector<Type>, Comp> pqClone = pq;

    while(!pqClone.empty())
    {

        ans.push_back(pqClone.top());
        pqClone.pop();

    }

    return ans;

}


template <class T>
class MVPTree
{

private:
    MVPNode<T>* root;
    size_t BF, PL, LC, LPN, FO, NS;
    DistanceFunction<BasicArrayObject<T>>* df;
    Pivot<T>* pvt;
    enum cases{R0, R1, R2, R3, R4};
    size_t leafNodeAccess;

private:
    void build(std::vector<datapoint_t<T>> &vec);
    std::vector<vp_t<T>> selectVantagePoints(std::vector<datapoint_t<T>> &vec);
    void markDistances(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec, size_t level);
    std::vector<double> calculatePivotDistances(vp_t<T> &vp, std::vector<datapoint_t<T>> &points);
    void calculateSplits(MVPNode<T>* node, std::vector<double> &dists, size_t n, size_t split_index);
    bool compareDistance(double a, double b, bool less);
    std::vector<datapoint_t<T>>* cullPoints(std::vector<datapoint_t<T>> &list, std::vector<double> &dists, double split, bool less);
    std::vector<std::vector<datapoint_t<T>>*> splitNode(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec);
    double minDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3);
    double minDist2(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3);
    double maxDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3);
    void initDisk();
    void saveLeafNode(MVPNode<T>* curr);
    Dataset<T>* readLeafNode(MVPNode<T>* curr);


public:
    MVPTree(Dataset<T>* dataset, DistanceFunction<BasicArrayObject<T>>* _df, Pivot<T>* _pvt, size_t bf, size_t pl, size_t lc, size_t lpn, size_t fo, size_t ns);

    MVPNode<T>* getRoot();
    void knn(BasicArrayObject<T> query, size_t k, std::vector<KnnEntryMVP<T>> &ans);
    size_t getLeafNodeAccess();

};


template <class T>
void MVPTree<T>::build(std::vector<datapoint_t<T>> &vec)
{

    std::queue<std::tuple<MVPNode<T>*, std::vector<datapoint_t<T>>, size_t>> b_queue;
    root = new MVPNode<T>(LPN, FO, NS);
    b_queue.push({root, vec, 0});
    MVPNode<T>* curr = nullptr;
    std::tuple<MVPNode<T>*, std::vector<datapoint_t<T>>, size_t> tuple;

    while(!b_queue.empty())
    {

        tuple = b_queue.front();
        b_queue.pop();
        curr = std::get<0>(tuple);
        std::vector<datapoint_t<T>> vec = std::get<1>(tuple);
        size_t level = std::get<2>(tuple);
        std::vector<vp_t<T>> vps = selectVantagePoints(vec);

        for(size_t i = 0; i < vps.size(); i++)
        {

            curr->setPivot(i, vps[i]);

        }

        vps.clear();

        if(vec.size() <= (LC + LPN))
        {

            curr->pushData(vec);
            curr->initializeLeafDists();

            for(size_t i = 0; i < curr->getDataSize(); i++)
            {

                for(size_t j = 0; j < LPN; j++)
                {

                    curr->setLeafDists(i, j, df->getDistance(*curr->getData(i).key, *curr->getVP(j).key));

                }

            }

        }
        else
        {

            markDistances(curr, vec, level);

            std::vector<std::vector<datapoint_t<T>>*> childrensVec = splitNode(curr, vec);

            for(size_t i = 0; i < FO; i++)
            {

                if(childrensVec[i] != nullptr)
                {

                    curr->setChildren(i, new MVPNode<T>(LPN, FO, NS));
                    b_queue.push({curr->getChild(i), *childrensVec[i], level + LPN});

                }

            }

        }

        vec.clear();

    }

    vec.clear();

}

template <class T>
std::vector<vp_t<T>> MVPTree<T>::selectVantagePoints(std::vector<datapoint_t<T>> &vec)
{

    std::vector<BasicArrayObject<T>> data;
    size_t dim = 0;

    for(size_t i = 0; i < vec.size(); i++)
    {

        dim = std::max(dim, vec[i].key->size());
        data.push_back(*vec[i].key);

    }

    size_t nPivots = std::min(LPN, data.size());
    Dataset<T>* dataset = new Dataset<T>(data, data.size(), dim);
    pvt->generatePivots(dataset, df, nPivots);
    std::vector<vp_t<T>> ans;

    for(size_t i = 0; i < nPivots; i++)
    {

        vp_t<T> vp = vp_t<T>(pvt->getPivot(i));
        ans.push_back(vp);

    }

    delete dataset;
    data.clear();

    return ans;

}

template <class T>
void MVPTree<T>::markDistances(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec, size_t level)
{

    for(size_t p = 0; p < LPN; p++)
    {

        if((level + p) < PL)
        {

            for(size_t i = 0; i < vec.size(); i++)
            {

                vec[i].dists[level + p] = df->getDistance(*node->getVP(p).key, *vec[i].key);

            }

        }

    }

}

template <class T>
std::vector<double> MVPTree<T>::calculatePivotDistances(vp_t<T> &vp, std::vector<datapoint_t<T>> &points)
{

    std::vector<double> ans;

    for(size_t i = 0; i < points.size(); i++)
    {

        ans.push_back(df->getDistance(*vp.key, *points[i].key));

    }

    return ans;

}

template <class T>
void MVPTree<T>::calculateSplits(MVPNode<T>* node, std::vector<double> &dists, size_t n, size_t split_index)
{

    size_t len = BF - 1;

    if(dists.size() > 0)
    {

        if(node->getSplit(n, split_index * len) == -1.0)
        {

            std::vector<double> tmp = dists;
            std::sort(tmp.begin(), tmp.end());
            double factor = (tmp.size() * 1.0)/BF;

            for(size_t i = 0; i < len; i++)
            {

                size_t pos = (i + 1) * factor;
                size_t lo = floor(pos);
                size_t hi = (pos <= tmp.size()) ? ceil(pos) : 0;
                node->setSplit(n, split_index * len + i, (tmp[lo] + tmp[hi])/2.0);

            }

        }

    }

}

template <class T>
bool MVPTree<T>::compareDistance(double a, double b, bool less)
{

    return (less) ? (a <= b) : (a > b);

}

template <class T>
std::vector<datapoint_t<T>>* MVPTree<T>::cullPoints(std::vector<datapoint_t<T>> &list, std::vector<double> &dists, double split, bool less)
{

    std::vector<datapoint_t<T>>* ans = new std::vector<datapoint_t<T>>();
    auto list_iter = list.begin();
    auto dist_iter = dists.begin();

    while((list_iter != list.end()) && (dist_iter != dists.end()))
    {

        if(compareDistance(*dist_iter, split, less))
        {

            ans->push_back(*list_iter);
            list_iter = list.erase(list_iter);
            dist_iter = dists.erase(dist_iter);

        }
        else
        {

            list_iter++;
            dist_iter++;

        }

    }

    if(ans->size() > 0)
    {

        return ans;

    }

    delete ans;
    return nullptr;

}

template <class T>
std::vector<std::vector<datapoint_t<T>>*> MVPTree<T>::splitNode(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec)
{

    std::map<long long, std::vector<datapoint_t<T>>*> pnts, pnts2;
    pnts[0] = &vec;
    std::vector<std::vector<datapoint_t<T>>*> ans = std::vector<std::vector<datapoint_t<T>>*>(FO, nullptr);

    size_t len = BF - 1, n = 0;
    double m;
    std::vector<datapoint_t<T>>* culledPts = nullptr;

    do{

        for(auto it = pnts.begin(); it != pnts.end(); it++)
        {

            long long node_index = it->first;
            std::vector<datapoint_t<T>>* list = it->second;
            vp_t<T> vp = node->getVP(n);
            std::vector<double> dists = calculatePivotDistances(vp, *list);

            if(dists.size() > 0)
            {

                calculateSplits(node, dists, n, node_index);

                for(size_t j = 0; j < len; j++)
                {

                    m = node->getSplit(n, node_index * len + j);
                    culledPts = cullPoints(*list, dists, m, true);

                    if(culledPts != nullptr)
                    {

                        pnts2[node_index * BF + j] = culledPts;

                    }

                }

                m = node->getSplit(n, node_index * len + len - 1);
                culledPts = cullPoints(*list, dists, m, false);

                if(culledPts != nullptr)
                {

                    pnts2[node_index * BF + BF - 1] = culledPts;

                }

            }

            if(list->size() > 0)
            {

                throw std::runtime_error("not fully collated !_!");

            }

            if(n > 0)
            {

                delete list;

            }

        }

        pnts = std::move(pnts2);
        n++;

    } while (n < LPN);

    for(auto it = pnts.begin(); it != pnts.end(); it++)
    {

        long long index = it->first;
        std::vector<datapoint_t<T>>* list = it->second;

        if(list != nullptr)
        {

            ans[index] = list;

        }

    }

    return ans;

}

template <class T>
double MVPTree<T>::minDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3)
{

    double ans;

    //std::cout << "MIN : " << sqCase << "\t" << nodeCase << "\n";

//    enum formatPartition{form_mu2_1, form_mu2_2, form_mu3_1, form_mu3_2, form_mu3_3};
//    enum formatMu2_1{mu3_sob, mu3_ext_mu1};
//    enum formatMu2_2{mu3_ext_mu1, mu3_menor_mu2, mu3_maior_mu2, mu3_contain_mu1};
//    enum formatMu3_1{mu2_maior_mu3};
//    enum formatMu3_2{mu3_menor_mu2, mu3_maior_mu2, mu3_contain_mu2};
//    enum formatMu3_3{mu3_contain_mu1_contain_mu2, mu3_maior_mu2};

//    formatPartition partitionInterno, partitionExterno;

//    if(d_sq_p2 <= mu2 && (d_sq_p2 + mu2) <= fabs(d_sq_p1 - mu1))
//    {

//        partitionInterno = formatPartition::form_mu2_1;

//    }
//    else
//    {

//        partitionInterno = formatPartition::form_mu2_2;

//    }

//    if(d_sq_p2 > mu3 && (fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) < (mu1 + mu3))
//    {

//        partitionExterno = formatPartition::form_mu3_1;

//    }
//    else
//    {

//        if(d_sq_p1 + mu1 <= fabs(d_sq_p2 - mu3))
//        {

//            partitionExterno = formatPartition::form_mu3_3;

//        }
//        else
//        {

//            partitionExterno = formatPartition::form_mu3_1;

//        }

//    }

    if (sqCase == R0)
    {

        if(nodeCase == R0)
            ans = 0.0;
        else if(nodeCase == R1)
            ans = fabs(d_sq_p2 - mu2);
        else if(nodeCase == R2)
        {

//            if(d_sq_p2 > mu3 && ((fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) > (mu1 + mu3))) //Externo
//            {

//                ans = fabs(d_sq_p2 - mu3);

//            }
//            else
//            {
//                ans = fabs(d_sq_p1 - mu1);

//            }

            if(d_sq_p2 > mu3 && ((fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) > (mu1 + mu3))) //bolas disjuntas em mu1 e mu3
            {
                ans = fabs(d_sq_p2 - mu3);
            }
            else
            {
                ans = fabs(d_sq_p1 - mu1);

            }

//            if(partitionInterno == formatPartition::form_mu2_1)
//            {

//                if((d_sq_p2 + mu2) <= (fabs(d_sq_p2 - mu3))) //Externo
//                {

//                    ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//                }
//                else //Sob
//                {

//                    ans = fabs(d_sq_p1 - mu1);

//                }

//            }
//            else
//            {

//                //mu3_ext_mu1, mu3_menor_mu2, mu3_maior_mu2, mu3_contain_mu1

//                //mu3_ext_mu1
//                if(d_sq_p2 > mu3 && (fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) < (mu1 + mu3)) //Externo
//                {

//                    ans = fabs(d_sq_p2 - mu3);

//                }
//                else
//                {

//                    //mu3_menor_mu2
//                    if(fabs(d_sq_p2 - mu3) > fabs(d_sq_p2 - mu2))
//                    {

//                        ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//                    }
//                    else
//                    {

//                        ans = fabs(d_sq_p1 - mu1);

//                    }

//                }

//            }


        }
        else if(nodeCase == R3)
        {

//            if((d_sq_p1 + mu1) <= (mu3 - fabs(d_sq_p2 - mu3))) //Externo
//            {
//                ans = fabs(d_sq_p2 - mu3);
//            }
//            else
//            {
//                ans = fabs(d_sq_p1 - mu1);

//            }

            if ((d_sq_p2 + d_sq_p1 + mu1) <= mu3) //mu3 cobre m1
            {
                ans = fabs(d_sq_p2 - mu3);
            }
            else
            {
                if ((d_sq_p2 + d_sq_p1 + mu2) <= mu3) //mu3 cobre mu2
                {
                    ans = fabs(d_sq_p2 - mu3);
                } else {
                    ans = fabs(d_sq_p1 - mu1);
                }
            }

//            if(partitionInterno == formatPartition::form_mu2_1)
//            {

//                ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//            }
//            else
//            {

//                //mu3_ext_mu1, mu3_menor_mu2, mu3_maior_mu2, mu3_contain_mu1

//                //mu3_ext_mu1
//                if(d_sq_p2 > mu3 && (fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) < (mu1 + mu3)) //Externo
//                {

//                    ans = fabs(d_sq_p1 - mu1);

//                }
//                else
//                {

//                    //mu3_maior_mu2
//                    if(fabs(d_sq_p2 - mu3) < fabs(d_sq_p2 - mu2))
//                    {

//                        ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//                    }
//                    else
//                    {

//                        ans = fabs(d_sq_p2 - mu3);

//                    }

//                }

//            }

        }
        else
            throw std::runtime_error("non-existent case");

    }
    else if(sqCase == R1)
    {

        if(nodeCase == R0)
            ans = fabs(d_sq_p2 - mu2);
        else if(nodeCase == R1)
            ans = 0.0;
        else if(nodeCase == R2)
        {
            //ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));
//            if(d_sq_p2 > mu3 && ((fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) > (mu1 + mu3))) //mu3 disjunto mu1
//            {

//                ans = fabs(d_sq_p2 - mu3);

//            }
//            else
//            {
//                if(((d_sq_p1 + mu1) > (mu3 - fabs(d_sq_p2 - mu3))) && (mu3 <= mu2)) //mu3 não cobre mu1
//                {
//                    ans = fabs(d_sq_p2 - mu3);
//                }
//                else
//                {
//                    ans = fabs(d_sq_p1 - mu1);

//                }

//            }

            if(d_sq_p2 > mu3 && ((fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) > (mu1 + mu3))) //bolas disjuntas em mu1 e mu3
            {
                ans = fabs(d_sq_p2 - mu3);
            }
            else
            {
                ans = fabs(d_sq_p1 - mu1);
            }

        }
        else if(nodeCase == R3)
        {
            //ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));

//            if((d_sq_p1 + mu1) <= (mu3 - fabs(d_sq_p2 - mu3))) //Externo
//            {
//                ans = fabs(d_sq_p2 - mu3);
//            }
//            else
//            {
//                ans = fabs(d_sq_p1 - mu1);

//            }

            if ((d_sq_p2 + d_sq_p1 + mu1) <= mu3) //mu3 cobre mu1
            {
                ans = fabs(d_sq_p2 - mu3);
            }
            else
            {
                ans = fabs(d_sq_p1 - mu1);

            }

        }
        else
            throw std::runtime_error("non-existent case");

    }
    else if(sqCase == R2)
    {

        if(nodeCase == R0)
        {
            //ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));

//            if((d_sq_p2 + mu2) <= (mu1 - fabs(d_sq_p1 - mu1))) //mu1 cobre mu2
//            {
//                ans = fabs(d_sq_p2 - mu2);
//            }
//            else
//            {
//                ans = fabs(d_sq_p1 - mu1);

//            }

            //ans = std::max(fabs(d_sq_p2 - mu2), fabs(d_sq_p1 - mu1));

//            if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){
//                ans = fabs(d_sq_p2 - mu2);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }

            if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){ //mu2 está contido em mu1
                ans = fabs(d_sq_p2 - mu2);
            } else {
                ans = fabs(d_sq_p1 - mu1);
            }

        }
        else if(nodeCase == R1)
        {
            //ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));

            ans = fabs(d_sq_p1 - mu1);

        }
        else if(nodeCase == R2)
            ans = 0.0;
        else if(nodeCase == R3)
        {
            //ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3));

            ans = fabs(d_sq_p2 - mu3);

        }
        else
            throw std::runtime_error("non-existent case");

    }
    else if(sqCase == R3)
    {

        if(nodeCase == R0)
        {
            //ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
//            if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){
//                ans = fabs(d_sq_p2 - mu2);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }

            if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){ //mu2 está contido em mu1
                ans = fabs(d_sq_p2 - mu2);
            } else {
                ans = fabs(d_sq_p1 - mu1);
            }

        }
        else if(nodeCase == R1){
            //ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
            //ans = fabs(d_sq_p1 - mu1);
            ans = fabs(d_sq_p1 - mu1);
        }
        else if(nodeCase == R2){
            //ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3));
            //ans = fabs(d_sq_p2 - mu3);
            ans = fabs(d_sq_p2 - mu3);
        }
        else if(nodeCase == R3)
            ans = 0.0;
        else
            throw std::runtime_error("non-existent case");

    }
    else if(sqCase == R4)
    {

        if(nodeCase == R0)
        {
            //ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2));
//            if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){
//                ans = fabs(d_sq_p2 - mu2);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }
            if ((d_sq_p1 > d_sq_p2) && ((d_sq_p2 + d_sq_p1 + mu2) <= mu1)){ //mu2 está contido em mu1
                ans = fabs(d_sq_p2 - mu2);
            } else {
                ans = fabs(d_sq_p1 - mu1);
            }
        }
        else if(nodeCase == R1){
            //ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2));
            //ans = fabs(d_sq_p1 - mu1);
            ans = fabs(d_sq_p1 - mu1);
        }
        else if(nodeCase == R2){
            //ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            //ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
        }
        else if(nodeCase == R3){
            //ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            ans = fabs(d_sq_p1 - M);

        }
        else
            throw std::runtime_error("non-existent case");

    }

//    switch(sqCase)
//    {

//    case R0:{

//        switch (nodeCase){

//        case R0:{
//            ans = 0.0;
//            break;
//        }

//        case R1:{
//            std::cout << "OLA\n";
//            ans = fabs(d_sq_p2 - mu2);
//            break;
//        }

//        case R2:{

//            if(partitionType == formatPartition::form_mu2_1)
//            {

//                if((d_sq_p2 + mu2) <= (fabs(d_sq_p2 - mu3))) //Externo
//                {

//                    ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//                }
//                else //Sob
//                {

//                    ans = fabs(d_sq_p1 - mu1);

//                }

//            }
//            else
//            {

//                //mu3_ext_mu1, mu3_menor_mu2, mu3_maior_mu2, mu3_contain_mu1

//                //mu3_ext_mu1
//                if(d_sq_p2 > mu3 && (fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) < (mu1 + mu3)) //Externo
//                {

//                    ans = fabs(d_sq_p2 - mu3);

//                }
//                else
//                {

//                    //mu3_menor_mu2
//                    if(fabs(d_sq_p2 - mu3) > fabs(d_sq_p2 - mu2))
//                    {

//                        ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//                    }
//                    else
//                    {

//                        ans = fabs(d_sq_p1 - mu1);

//                    }

//                }

//            }

//            break;
//        }

//        case R3:{

//            if(partitionType == formatPartition::form_mu2_1)
//            {

//                ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//            }
//            else
//            {

//                //mu3_ext_mu1, mu3_menor_mu2, mu3_maior_mu2, mu3_contain_mu1

//                //mu3_ext_mu1
//                if(d_sq_p2 > mu3 && (fabs(d_sq_p1 - mu1) + fabs(d_sq_p2 - mu3)) < (mu1 + mu3)) //Externo
//                {

//                    ans = fabs(d_sq_p1 - mu1);

//                }
//                else
//                {

//                    //mu3_maior_mu2
//                    if(fabs(d_sq_p2 - mu3) < fabs(d_sq_p2 - mu2))
//                    {

//                        ans = std::max(fabs(d_sq_p2 - mu3), fabs(d_sq_p1 - mu1));

//                    }
//                    else
//                    {

//                        ans = fabs(d_sq_p2 - mu3);

//                    }

//                }

//            }

//            break;
//        }


//        break;

//    }

//    }

//    case R1:{

//    switch(nodeCase)
//    {

//    case R0:{
//            ans = fabs(d_sq_p2 - mu2);
//            break;
//    }

//    case R1:{
//            ans = 0.0;
//            break;
//    }

//    case R2:{
//            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));
//            break;
//    }

//    case R3:{
//            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));
//            break;
//    }

//    case R4:{
//            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
//            break;
//    }

//    }
//    break;

//    }

//    case R2:{

//    switch(nodeCase)
//    {

//    case R0:{
//            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));//
//            break;
//    }

//    case R1:{
//            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
//            break;
//    }

//    case R2:{
//            ans = 0.0;
//            break;
//    }

//    case R3:{
//            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3));
//            break;
//    }

//    case R4:{
//            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
//            break;
//    }

//    }
//    break;

//    }

//    case R3:{

//    switch(nodeCase)
//    {

//    case R0:{
//            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
//            break;
//    }

//    case R1:{
//            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
//            break;
//    }

//    case R2:{
//            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3));
//            break;
//    }

//    case R3:{
//            ans = 0.0;
//            break;
//    }

//    case R4:{
//            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
//            break;
//    }

//    }
//    break;

//    }

//    case R4:{

//    switch(nodeCase)
//    {

//    case R0:{
//            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2));
//            break;
//    }

//    case R1:{
//            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2));
//            break;
//    }

//    case R2:{
//            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
//            break;
//    }

//    case R3:{
//            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
//            break;
//    }

//    case R4:{
//            ans = 0.0;
//            break;
//    }

//    }
//    break;

//    }

//    }


//    switch(sqCase)
//    {

//    case R0:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = 0.0;
//            break;
//        }

//        case R1:{
//            ans = fabs(d_sq_p2 - mu2);
//            break;
//        }

//        case R2:{

//            if ((d_sq_p2 >  mu3) && (fabs(d_sq_p2 - d_sq_p1) > fabs(mu1 - mu3))){
//                ans = fabs(d_sq_p2 - mu3);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }
//            break;
//        }

//        case R3:{
//            if (fabs(d_sq_p2 - mu3) > (d_sq_p1 + mu1)){
//                ans = fabs(d_sq_p2 - mu3);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }
//            break;
//        }

//        }
//        break;
//    }

//    case R1:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = fabs(d_sq_p2 - mu2);
//            break;
//        }

//        case R1:{
//            ans = 0.0;
//            break;
//        }

//        case R2:{

//            if ((d_sq_p2 >  mu3) && (fabs(d_sq_p2 - d_sq_p1) > fabs(mu1 - mu3))){
//                ans = fabs(d_sq_p2 - mu3);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }
//            break;
//        }

//        case R3:{
//            if (fabs(d_sq_p2 - mu3) > (d_sq_p1 + mu1)){
//                ans = fabs(d_sq_p2 - mu3);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }
//            break;
//        }

//        }
//        break;
//    }

//    case R2:{

//        switch(nodeCase)
//        {
//        case R0:{
//            if ((d_sq_p2 + mu2) < fabs(d_sq_p1 - mu1)){
//                ans = fabs(d_sq_p2 - mu2);
//            } else {
//                ans = fabs(d_sq_p1 - mu1);
//            }
//            break;
//        }

//        case R1:{
//            if ((d_sq_p2 + mu2) < fabs(d_sq_p1 - mu1)){
//                ans = fabs(d_sq_p1 - mu1);
//            } else {
//                ans = fabs(d_sq_p2 - mu2);
//            }
//            break;
//        }

//        case R2:{
//            ans = 0.0;
//            break;
//        }

//        case R3:{
//            ans = fabs(d_sq_p2 - mu3);
//            break;
//        }

//        }
//        break;
//    }

//    case R3:{

//        switch(nodeCase)
//        {

//        case R0:{
//            if ((d_sq_p2 + mu2) < fabs(d_sq_p1 - mu1)){
//                ans = fabs(d_sq_p2 - mu2);
//            } else {
//                if (fabs(d_sq_p2 - mu2) > fabs(d_sq_p2 - mu3)){
//                    ans = fabs(d_sq_p1 - mu1);
//                } else {
//                    ans = fabs(d_sq_p2 - mu2);
//                }
//            }
//            break;
//        }

//        case R1:{
//            if ((d_sq_p2 + mu2) < fabs(d_sq_p1 - mu1)){
//                ans = fabs(d_sq_p1 - mu1);
//            } else {
//                if (fabs(d_sq_p2 - mu2) > fabs(d_sq_p2 - mu3)){
//                    ans = fabs(d_sq_p2 - mu2);
//                } else {
//                    ans = fabs(d_sq_p1 - mu1);
//                }
//            }
//            break;
//        }

//        case R2:{

//            ans = fabs(d_sq_p2 - mu3);
//            break;
//        }

//        case R3:{
//            ans = 0.0;
//            break;
//        }

//        }
//        break;
//    }
//    case R4:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2));
//            break;
//        }

//        case R1:{
//            ans = fabs(d_sq_p1 - mu1);
//            break;
//        }

//        case R2:{
//            ans = std::max(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
//            break;
//        }

//        case R3:{
//            ans = fabs(d_sq_p1 - M);
//            break;
//        }

//        }
//        break;

//    }

//    }

    //std::cout << "DIST = " << ans << "\n";
    return ans;

}


template <class T>
double MVPTree<T>::minDist2(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3)
{

    double ans;

    switch(sqCase)
    {

    case R0:{

        switch(nodeCase)
        {

        case R0:{
            ans = 0.0;
            break;
        }

        case R1:{
            ans = fabs(d_sq_p2 - mu2);
            break;
        }

        case R2:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));
            break;
        }

        case R3:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));
            break;
        }

        case R4:{
            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            break;
        }

        }
        break;

    }

    case R1:{

        switch(nodeCase)
        {

        case R0:{
            ans = fabs(d_sq_p2 - mu2);
            break;
        }

        case R1:{
            ans = 0.0;
            break;
        }

        case R2:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));
            break;
        }

        case R3:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3)));
            break;
        }

        case R4:{
            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            break;
        }

        }
        break;

    }

    case R2:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));//
            break;
        }

        case R1:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
            break;
        }

        case R2:{
            ans = 0.0;
            break;
        }

        case R3:{
            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3));
            break;
        }

        case R4:{
            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            break;
        }

        }
        break;

    }

    case R3:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
            break;
        }

        case R1:{
            ans = std::max(fabs(d_sq_p1 - mu1), std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2)));
            break;
        }

        case R2:{
            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu3));
            break;
        }

        case R3:{
            ans = 0.0;
            break;
        }

        case R4:{
            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            break;
        }

        }
        break;

    }

    case R4:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2));
            break;
        }

        case R1:{
            ans = std::min(fabs(d_sq_p1 - mu1), fabs(d_sq_p2 - mu2));
            break;
        }

        case R2:{
            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            break;
        }

        case R3:{
            ans = std::min(fabs(d_sq_p1 - M), fabs(d_sq_p2 - mu3));
            break;
        }

        case R4:{
            ans = 0.0;
            break;
        }

        }
        break;

    }

    }

    return ans;

}

template <class T>
double MVPTree<T>::maxDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double mu1, double M, double mu2, double mu3)
{

    double ans;

    //std::cout << "MAX : " << sqCase << "\t" << nodeCase << "\n";

    switch(sqCase)
    {

    case R0:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
            break;
        }

        case R1:{
            ans = d_sq_p1 + mu1;
            break;
        }

        case R2:{
            ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
            break;
        }

        case R3:{
            ans = d_sq_p1 + M;
            break;
        }

        }
        break;

    }

    case R1:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
            break;
        }

        case R1:{
            ans = d_sq_p1 + mu1;
            break;
        }

        case R2:{
            ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
            break;
        }

        case R3:{
            ans = d_sq_p1 + M;
            break;
        }

        }
        break;

    }

    case R2:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
            break;
        }

        case R1:{
            ans = d_sq_p1 + mu1;
            break;
        }

        case R2:{
            ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
            break;
        }

        case R3:{
            ans = d_sq_p1 + M;
            break;
        }

        }
        break;

    }

    case R3:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
            break;
        }

        case R1:{
            ans = d_sq_p1 + mu1;
            break;
        }

        case R2:{
            ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
            break;
        }

        case R3:{
            ans = d_sq_p1 + M;
            break;
        }


        }
        break;

    }

    case R4:{

        switch(nodeCase)
        {

        case R0:{
            ans = std::min(d_sq_p1 + mu1, d_sq_p2 + mu2);
            break;
        }

        case R1:{
            ans = d_sq_p1 + mu1;
            break;
        }

        case R2:{
            ans = std::min(d_sq_p1 + M, d_sq_p2 + mu3);
            break;
        }

        case R3:{
            ans = d_sq_p1 + M;
            break;
        }

        }
        break;

    }

    }

    //std::cout << "DIST = " << ans << "\n";
    return ans;

}

template <class T>
void MVPTree<T>::initDisk()
{

    baseFilePath = "../mvptree/mvp_files";
    setBaseFilePath("MVPfiles");

    std::queue<MVPNode<T>*> queue;
    queue.push(getRoot());
    MVPNode<T>* curr = nullptr;
    size_t id = 0;

    while(!queue.empty())
    {

        curr = queue.front();
        queue.pop();

        if(curr->isLeaf())
        {

            curr->setPageID(id++);
            saveLeafNode(curr);

        }
        else
        {

            for(size_t i = 0; i < FO; i++)
            {

                if(curr->getChild(i) != nullptr) queue.push(curr->getChild(i));

            }

        }

    }

}

template <class T>
void MVPTree<T>::saveLeafNode(MVPNode<T>* curr)
{

    std::vector<BasicArrayObject<T>> data = curr->purgeData();
    Dataset<T>* datasetLeaf = new Dataset<T>(data, data.size(), data[0].size());
    write_dataset_to_disk(datasetLeaf, curr->getPageID());
    delete datasetLeaf;
    curr->clear();
    data.clear();

}

template <class T>
Dataset<T>* MVPTree<T>::readLeafNode(MVPNode<T>* curr)
{

    Dataset<T>* datasetLeaf = new Dataset<T>();
    read_dataset_from_disk(datasetLeaf, curr->getPageID());
    return datasetLeaf;

}

template <class T>
MVPTree<T>::MVPTree(Dataset<T>* dataset, DistanceFunction<BasicArrayObject<T>>* _df, Pivot<T>* _pvt, size_t bf, size_t pl, size_t lc, size_t lpn, size_t fo, size_t ns)
{

    std::vector<datapoint_t<T>> vec;
    df = _df;
    pvt = _pvt;
    BF = bf;
    PL = pl;
    LC = lc;
    LPN = lpn;
    FO = fo;
    NS = ns;
    leafNodeAccess = 0;

    for(size_t i = 0; i < dataset->getCardinality(); i++)
    {

        datapoint_t<T> dp = datapoint_t<T>(dataset->getInstance(i), pl);
        vec.push_back(dp);

    }

    build(vec);
    initDisk();

}

template <class T>
MVPNode<T>* MVPTree<T>::getRoot()
{

    return root;

}

template <class T>
void MVPTree<T>::knn(BasicArrayObject<T> query, size_t k, std::vector<KnnEntryMVP<T>> &ans)
{

    df->resetStatistics();
    leafNodeAccess = 0;
    std::priority_queue<MVPPartition<T>, std::vector<MVPPartition<T>>, CompareMVPPartition<T>> nodeQueue;
    std::priority_queue<KnnEntryMVP<T>, std::vector<KnnEntryMVP<T>>, std::greater<KnnEntryMVP<T>>> candidatesQueue;
    std::priority_queue<KnnEntryMVP<T>, std::vector<KnnEntryMVP<T>>, std::less<KnnEntryMVP<T>>> resultQueue;
    nodeQueue.push(MVPPartition<T>(getRoot(), 0.0, std::numeric_limits<double>::max()));
    MVPNode<T>* node = nullptr;
    MVPPartition<T> partition;
    std::vector<cases> casesNodeVec = {R0, R1, R2, R3};
    Dataset<T>* datasetLeaf;

    while(!nodeQueue.empty() || candidatesQueue.size() > 0)
    {

        if(candidatesQueue.size() == 0)
        {

            partition = nodeQueue.top();
            node = partition.node;
            nodeQueue.pop();

            if(node->isLeaf())
            {

                leafNodeAccess++;

                datasetLeaf = readLeafNode(node);

                for(size_t i = 0; i < datasetLeaf->getCardinality(); i++)
                {

                    candidatesQueue.push(KnnEntryMVP<T>(datasetLeaf->getFeatureVector(i), df->getDistance(datasetLeaf->getFeatureVector(i), query)));

                }

                delete datasetLeaf;

//                std::vector<BasicArrayObject<T>> data = node->purgeData();

//                for(size_t i = 0; i < data.size(); i++)
//                {

//                    candidatesQueue.push(KnnEntryMVP<T>(data[i], df->getDistance(query, data[i])));

//                }

            }
            else
            {

                cases sqCase;
                double d_sq_p1 = df->getDistance(query, *node->getVP(0).key);
                double d_sq_p2 = df->getDistance(query, *node->getVP(1).key);

                if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 0))
                    sqCase = R0;
                else if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 0))
                    sqCase = R1;
                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 1))
                    sqCase = R2;
                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 1))
                    sqCase = R3;
                else
                    sqCase = R4;

//                std::cout << "\n\n\n\n";
//                std::cout << "sqCase: " << sqCase << std::endl;
//                std::cout << "sq: " << query.toString(",") << std::endl;
//                std::cout << "p1: " << node->getVP(0).key->toString(",") << std::endl;
//                std::cout << "p2: " << node->getVP(1).key->toString(",") << std::endl;
//                std::cout << "d_sq_p1 = " << d_sq_p1 << std::endl;
//                std::cout << "d_sq_p2 = " << d_sq_p2 << std::endl;
//                std::cout << "u1 = " << node->getSplit(0, 0) << std::endl;
//                std::cout << "u2 = " << node->getSplit(1, 0) << std::endl;
//                std::cout << "u3 = " << node->getSplit(1, 1) << std::endl;
//                std::cout << "M = " << partition.max << std::endl;

                for(size_t i = 0; i < FO; i++)
                {

                    if(node->getChild(i) != nullptr)
                    {

                        //std::cout << "filho " << i+1 << ": \n";
                        nodeQueue.push(MVPPartition<T>(node->getChild(i),
                                                       minDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1)),
                                                       maxDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1))));

                    }

                }


            }

        }
        else if((nodeQueue.size() > 0) && nodeQueue.top().min < candidatesQueue.top().distance)
        {

            partition = nodeQueue.top();
            node = partition.node;
            nodeQueue.pop();

            if(node->isLeaf())
            {

                leafNodeAccess++;

                datasetLeaf = readLeafNode(node);

                for(size_t i = 0; i < datasetLeaf->getCardinality(); i++)
                {

                    candidatesQueue.push(KnnEntryMVP<T>(datasetLeaf->getFeatureVector(i), df->getDistance(datasetLeaf->getFeatureVector(i), query)));

                }

                delete datasetLeaf;

//                std::vector<BasicArrayObject<T>> data = node->purgeData();

//                for(size_t i = 0; i < data.size(); i++)
//                {

//                    candidatesQueue.push(KnnEntryMVP<T>(data[i], df->getDistance(query, data[i])));

//                }

            }
            else
            {

                cases sqCase;
                double d_sq_p1 = df->getDistance(query, *node->getVP(0).key);
                double d_sq_p2 = df->getDistance(query, *node->getVP(1).key);

                if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 0))
                    sqCase = R0;
                else if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 0))
                    sqCase = R1;
                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 1))
                    sqCase = R2;
                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 1))
                    sqCase = R3;
                else
                    sqCase = R4;

//                std::cout << "\n\n\n\n";
//                std::cout << "sqCase: " << sqCase << std::endl;
//                std::cout << "sq: " << query.toString(",") << std::endl;
//                std::cout << "p1: " << node->getVP(0).key->toString(",") << std::endl;
//                std::cout << "p2: " << node->getVP(1).key->toString(",") << std::endl;
//                std::cout << "d_sq_p1 = " << d_sq_p1 << std::endl;
//                std::cout << "d_sq_p2 = " << d_sq_p2 << std::endl;
//                std::cout << "u1 = " << node->getSplit(0, 0) << std::endl;
//                std::cout << "u2 = " << node->getSplit(1, 0) << std::endl;
//                std::cout << "u3 = " << node->getSplit(1, 1) << std::endl;
//                std::cout << "M = " << partition.max << std::endl;

                for(size_t i = 0; i < FO; i++)
                {

                    if(node->getChild(i) != nullptr)
                    {

                        //std::cout << "filho " << i+1 << ": \n";
                        nodeQueue.push(MVPPartition<T>(node->getChild(i),
                                                       minDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1)),
                                                       maxDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1))));

                    }

                }

            }

        }
        else
        {

            if(!resultQueue.empty() && !candidatesQueue.empty() && resultQueue.size() >= k && candidatesQueue.top().distance > resultQueue.top().distance)
            {

                break;

            }

            resultQueue.push(candidatesQueue.top());
            candidatesQueue.pop();

            while(resultQueue.size() > k)
            {

                resultQueue.pop();

            }

        }

    }

    ans = dequeueInOrderMVP_Results(resultQueue);
    std::reverse(ans.begin(), ans.end());

    while(!candidatesQueue.empty())
    {

        candidatesQueue.pop();

    }

    while(!resultQueue.empty())
    {

        resultQueue.pop();

    }

    while(!nodeQueue.empty())
    {

        nodeQueue.pop();

    }

}

template <class T>
size_t MVPTree<T>::getLeafNodeAccess()
{

    return leafNodeAccess;

}


#endif // MVPTREE_H

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////ARQUIVO ANTIGO/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#ifndef MVPTREE_H
//#define MVPTREE_H

//#include <Dataset.h>
//#include <Hermes.h>
//#include <Pivots.h>
//#include <MVPData.h>
//#include <MVPNode.h>
//#include <queue>
//#include <config_spb.h>
//#include <MemoryManagerUtils.h>

//template <class T>
//struct MVPPartition
//{

//    MVPNode<T>* node;
//    double min, max;

//    MVPPartition()
//    {



//    }

//    MVPPartition(MVPNode<T>* _node, double _min, double _max)
//    {

//        node = _node;
//        min = _min;
//        max = _max;

//    }

//};

//template <class T>
//struct CompareMVPPartition
//{

//    bool operator()(MVPPartition<T> const& a, MVPPartition<T> const& b)
//    {

//        bool ans;

//        if(a.min != b.min)
//        {

//            ans = a.min > b.min;

//        }
//        else
//        {

//            ans = a.max > b.max;

//        }

//        return ans;

//    }

//};

//template <class T>
//struct KnnEntryMVP
//{

//    BasicArrayObject<T> element;
//    double distance;

//    KnnEntryMVP(BasicArrayObject<T> &element_, double distance_)
//    {

//        element = element_;
//        distance = distance_;

//    }

//    KnnEntryMVP()
//    {



//    }

//    bool operator<(const KnnEntryMVP<T>& item) const
//    {

//        return distance < item.distance;

//    }

//    bool operator>(const KnnEntryMVP<T>& item) const
//    {

//        return distance > item.distance;

//    }

//};

//template<class Type, class Comp>
//std::vector<Type> dequeueInOrderMVP_Results(std::priority_queue<Type, std::vector<Type>, Comp> pq)
//{

//    std::vector<Type> ans;
//    std::priority_queue<Type, std::vector<Type>, Comp> pqClone = pq;

//    while(!pqClone.empty())
//    {

//        ans.push_back(pqClone.top());
//        pqClone.pop();

//    }

//    return ans;

//}


//template <class T>
//class MVPTree
//{

//private:
//    MVPNode<T>* root;
//    size_t BF, PL, LC, LPN, FO, NS;
//    DistanceFunction<BasicArrayObject<T>>* df;
//    Pivot<T>* pvt;
//    enum cases{R0, R1, R2, R3, R4};
//    size_t leafNodeAccess;

//private:
//    void build(std::vector<datapoint_t<T>> &vec);
//    std::vector<vp_t<T>> selectVantagePoints(std::vector<datapoint_t<T>> &vec);
//    void markDistances(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec, size_t level);
//    std::vector<double> calculatePivotDistances(vp_t<T> &vp, std::vector<datapoint_t<T>> &points);
//    void calculateSplits(MVPNode<T>* node, std::vector<double> &dists, size_t n, size_t split_index);
//    bool compareDistance(double a, double b, bool less);
//    std::vector<datapoint_t<T>>* cullPoints(std::vector<datapoint_t<T>> &list, std::vector<double> &dists, double split, bool less);
//    std::vector<std::vector<datapoint_t<T>>*> splitNode(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec);
//    double minDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2);
//    double maxDist(cases nodeCase, cases sqCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2);
//    void initDisk();
//    void saveLeafNode(MVPNode<T>* curr);
//    Dataset<T>* readLeafNode(MVPNode<T>* curr);


//public:
//    MVPTree(Dataset<T>* dataset, DistanceFunction<BasicArrayObject<T>>* _df, Pivot<T>* _pvt, size_t bf, size_t pl, size_t lc, size_t lpn, size_t fo, size_t ns);

//    MVPNode<T>* getRoot();
//    void knn(BasicArrayObject<T> query, size_t k, std::vector<KnnEntryMVP<T>> &ans);
//    size_t getLeafNodeAccess();

//};


//template <class T>
//void MVPTree<T>::build(std::vector<datapoint_t<T>> &vec)
//{

//    std::queue<std::tuple<MVPNode<T>*, std::vector<datapoint_t<T>>, size_t>> b_queue;
//    root = new MVPNode<T>(LPN, FO, NS);
//    b_queue.push({root, vec, 0});
//    MVPNode<T>* curr = nullptr;
//    std::tuple<MVPNode<T>*, std::vector<datapoint_t<T>>, size_t> tuple;

//    while(!b_queue.empty())
//    {

//        tuple = b_queue.front();
//        b_queue.pop();
//        curr = std::get<0>(tuple);
//        std::vector<datapoint_t<T>> vec = std::get<1>(tuple);
//        size_t level = std::get<2>(tuple);
//        std::vector<vp_t<T>> vps = selectVantagePoints(vec);

//        for(size_t i = 0; i < vps.size(); i++)
//        {

//            curr->setPivot(i, vps[i]);

//        }

//        vps.clear();

//        if(vec.size() <= (LC + LPN))
//        {

//            curr->pushData(vec);
//            curr->initializeLeafDists();

//            for(size_t i = 0; i < curr->getDataSize(); i++)
//            {

//                for(size_t j = 0; j < LPN; j++)
//                {

//                    curr->setLeafDists(i, j, df->getDistance(*curr->getData(i).key, *curr->getVP(j).key));

//                }

//            }

//        }
//        else
//        {

//            markDistances(curr, vec, level);

//            std::vector<std::vector<datapoint_t<T>>*> childrensVec = splitNode(curr, vec);

//            for(size_t i = 0; i < FO; i++)
//            {

//                if(childrensVec[i] != nullptr)
//                {

//                    curr->setChildren(i, new MVPNode<T>(LPN, FO, NS));
//                    b_queue.push({curr->getChild(i), *childrensVec[i], level + LPN});

//                }

//            }

//        }

//        vec.clear();

//    }

//    vec.clear();

//}

//template <class T>
//std::vector<vp_t<T>> MVPTree<T>::selectVantagePoints(std::vector<datapoint_t<T>> &vec)
//{

//    std::vector<BasicArrayObject<T>> data;
//    size_t dim = 0;

//    for(size_t i = 0; i < vec.size(); i++)
//    {

//        dim = std::max(dim, vec[i].key->size());
//        data.push_back(*vec[i].key);

//    }

//    size_t nPivots = std::min(LPN, data.size());
//    Dataset<T>* dataset = new Dataset<T>(data, data.size(), dim);
//    pvt->generatePivots(dataset, df, nPivots);
//    std::vector<vp_t<T>> ans;

//    for(size_t i = 0; i < nPivots; i++)
//    {

//        vp_t<T> vp = vp_t<T>(pvt->getPivot(i));
//        ans.push_back(vp);

//    }

//    delete dataset;
//    data.clear();

//    return ans;

//}

//template <class T>
//void MVPTree<T>::markDistances(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec, size_t level)
//{

//    for(size_t p = 0; p < LPN; p++)
//    {

//        if((level + p) < PL)
//        {

//            for(size_t i = 0; i < vec.size(); i++)
//            {

//                vec[i].dists[level + p] = df->getDistance(*node->getVP(p).key, *vec[i].key);

//            }

//        }

//    }

//}

//template <class T>
//std::vector<double> MVPTree<T>::calculatePivotDistances(vp_t<T> &vp, std::vector<datapoint_t<T>> &points)
//{

//    std::vector<double> ans;

//    for(size_t i = 0; i < points.size(); i++)
//    {

//        ans.push_back(df->getDistance(*vp.key, *points[i].key));

//    }

//    return ans;

//}

//template <class T>
//void MVPTree<T>::calculateSplits(MVPNode<T>* node, std::vector<double> &dists, size_t n, size_t split_index)
//{

//    size_t len = BF - 1;

//    if(dists.size() > 0)
//    {

//        if(node->getSplit(n, split_index * len) == -1.0)
//        {

//            std::vector<double> tmp = dists;
//            std::sort(tmp.begin(), tmp.end());
//            double factor = (tmp.size() * 1.0)/BF;

//            for(size_t i = 0; i < len; i++)
//            {

//                size_t pos = (i + 1) * factor;
//                size_t lo = floor(pos);
//                size_t hi = (pos <= tmp.size()) ? ceil(pos) : 0;
//                node->setSplit(n, split_index * len + i, (tmp[lo] + tmp[hi])/2.0);

//            }

//        }

//    }

//}

//template <class T>
//bool MVPTree<T>::compareDistance(double a, double b, bool less)
//{

//    return (less) ? (a <= b) : (a > b);

//}

//template <class T>
//std::vector<datapoint_t<T>>* MVPTree<T>::cullPoints(std::vector<datapoint_t<T>> &list, std::vector<double> &dists, double split, bool less)
//{

//    std::vector<datapoint_t<T>>* ans = new std::vector<datapoint_t<T>>();
//    auto list_iter = list.begin();
//    auto dist_iter = dists.begin();

//    while((list_iter != list.end()) && (dist_iter != dists.end()))
//    {

//        if(compareDistance(*dist_iter, split, less))
//        {

//            ans->push_back(*list_iter);
//            list_iter = list.erase(list_iter);
//            dist_iter = dists.erase(dist_iter);

//        }
//        else
//        {

//            list_iter++;
//            dist_iter++;

//        }

//    }

//    if(ans->size() > 0)
//    {

//        return ans;

//    }

//    delete ans;
//    return nullptr;

//}

//template <class T>
//std::vector<std::vector<datapoint_t<T>>*> MVPTree<T>::splitNode(MVPNode<T>* node, std::vector<datapoint_t<T>> &vec)
//{

//    std::map<long long, std::vector<datapoint_t<T>>*> pnts, pnts2;
//    pnts[0] = &vec;
//    std::vector<std::vector<datapoint_t<T>>*> ans = std::vector<std::vector<datapoint_t<T>>*>(FO, nullptr);

//    size_t len = BF - 1, n = 0;
//    double m;
//    std::vector<datapoint_t<T>>* culledPts = nullptr;

//    do{

//        for(auto it = pnts.begin(); it != pnts.end(); it++)
//        {

//            long long node_index = it->first;
//            std::vector<datapoint_t<T>>* list = it->second;
//            vp_t<T> vp = node->getVP(n);
//            std::vector<double> dists = calculatePivotDistances(vp, *list);

//            if(dists.size() > 0)
//            {

//                calculateSplits(node, dists, n, node_index);

//                for(size_t j = 0; j < len; j++)
//                {

//                    m = node->getSplit(n, node_index * len + j);
//                    culledPts = cullPoints(*list, dists, m, true);

//                    if(culledPts != nullptr)
//                    {

//                        pnts2[node_index * BF + j] = culledPts;

//                    }

//                }

//                m = node->getSplit(n, node_index * len + len - 1);
//                culledPts = cullPoints(*list, dists, m, false);

//                if(culledPts != nullptr)
//                {

//                    pnts2[node_index * BF + BF - 1] = culledPts;

//                }

//            }

//            if(list->size() > 0)
//            {

//                throw std::runtime_error("not fully collated !_!");

//            }

//            if(n > 0)
//            {

//                delete list;

//            }

//        }

//        pnts = std::move(pnts2);
//        n++;

//    } while (n < LPN);

//    for(auto it = pnts.begin(); it != pnts.end(); it++)
//    {

//        long long index = it->first;
//        std::vector<datapoint_t<T>>* list = it->second;

//        if(list != nullptr)
//        {

//            ans[index] = list;

//        }

//    }

//    return ans;

//}

//template <class T>
//double MVPTree<T>::minDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2)
//{

//    double ans;

//    switch(sqCase)
//    {

//    case R0:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = 0.0;
//            break;
//        }

//        case R1:{
//            ans = fabs(d_sq_p2 - p2_ct1);
//            break;
//        }

//        case R2:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
//            break;
//        }

//        case R3:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
//            break;
//        }

//        case R4:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        }
//        break;

//    }

//    case R1:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = fabs(d_sq_p2 - p2_ct1);
//            break;
//        }

//        case R1:{
//            ans = 0.0;
//            break;
//        }

//        case R2:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
//            break;
//        }

//        case R3:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
//            break;
//        }

//        case R4:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        }
//        break;

//    }

//    case R2:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));//
//            break;
//        }

//        case R1:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));
//            break;
//        }

//        case R2:{
//            ans = 0.0;
//            break;
//        }

//        case R3:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        case R4:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        }
//        break;

//    }

//    case R3:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));
//            break;
//        }

//        case R1:{
//            ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));
//            break;
//        }

//        case R2:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        case R3:{
//            ans = 0.0;
//            break;
//        }

//        case R4:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        }
//        break;

//    }

//    case R4:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1));
//            break;
//        }

//        case R1:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1));
//            break;
//        }

//        case R2:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        case R3:{
//            ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
//            break;
//        }

//        case R4:{
//            ans = 0.0;
//            break;
//        }

//        }
//        break;

//    }

//    }

//    return ans;

//}

//template <class T>
//double MVPTree<T>::maxDist(cases nodeCase, cases sqCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2)
//{

//    double ans;

//    switch(sqCase)
//    {

//    case R0:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
//            break;
//        }

//        case R1:{
//            ans = d_sq_p1 + p1_ct1;
//            break;
//        }

//        case R2:{
//            ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
//            break;
//        }

//        case R3:{
//            ans = d_sq_p1 + p1_ct2;
//            break;
//        }

//        case R4:{
//            ans = -1.0; //!!!!!
//            break;
//        }

//        }
//        break;

//    }

//    case R1:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
//            break;
//        }

//        case R1:{
//            ans = d_sq_p1 + p1_ct1;
//            break;
//        }

//        case R2:{
//            ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
//            break;
//        }

//        case R3:{
//            ans = d_sq_p1 + p1_ct2;
//            break;
//        }

//        case R4:{
//            ans = -1.0; //!!!!!
//            break;
//        }

//        }
//        break;

//    }

//    case R2:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
//            break;
//        }

//        case R1:{
//            ans = d_sq_p1 + p1_ct1;
//            break;
//        }

//        case R2:{
//            ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
//            break;
//        }

//        case R3:{
//            ans = d_sq_p1 + p1_ct2;
//            break;
//        }

//        case R4:{
//            ans = -1.0; //!!!!!
//            break;
//        }

//        }
//        break;

//    }

//    case R3:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
//            break;
//        }

//        case R1:{
//            ans = d_sq_p1 + p1_ct1;
//            break;
//        }

//        case R2:{
//            ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
//            //                    ans = d_sq_p1 + p1_ct2;
//            break;
//        }

//        case R3:{
//            ans = d_sq_p1 + p1_ct2;
//            break;
//        }

//        case R4:{
//            ans = -1.0; //!!!!!
//            break;
//        }

//        }
//        break;

//    }

//    case R4:{

//        switch(nodeCase)
//        {

//        case R0:{
//            ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
//            break;
//        }

//        case R1:{
//            ans = d_sq_p1 + p1_ct1;
//            break;
//        }

//        case R2:{
//            ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
//            break;
//        }

//        case R3:{
//            ans = d_sq_p1 + p1_ct2;
//            break;
//        }

//        case R4:{
//            ans = -1.0; //!!!!!
//            break;
//        }

//        }
//        break;

//    }

//    }

//    return ans;

//}

//template <class T>
//void MVPTree<T>::initDisk()
//{

//    baseFilePath = "../mvptree/mvp_files";
//    setBaseFilePath("MVPfiles");

//    std::queue<MVPNode<T>*> queue;
//    queue.push(getRoot());
//    MVPNode<T>* curr = nullptr;
//    size_t id = 0;

//    while(!queue.empty())
//    {

//        curr = queue.front();
//        queue.pop();

//        if(curr->isLeaf())
//        {

//            curr->setPageID(id++);
//            saveLeafNode(curr);

//        }
//        else
//        {

//            for(size_t i = 0; i < FO; i++)
//            {

//                if(curr->getChild(i) != nullptr) queue.push(curr->getChild(i));

//            }

//        }

//    }

//}

//template <class T>
//void MVPTree<T>::saveLeafNode(MVPNode<T>* curr)
//{

//    std::vector<BasicArrayObject<T>> data = curr->purgeData();
//    Dataset<T>* datasetLeaf = new Dataset<T>(data, data.size(), data[0].size());
//    write_dataset_to_disk(datasetLeaf, curr->getPageID());
//    delete datasetLeaf;
//    curr->clear();
//    data.clear();

//}

//template <class T>
//Dataset<T>* MVPTree<T>::readLeafNode(MVPNode<T>* curr)
//{

//    Dataset<T>* datasetLeaf = new Dataset<T>();
//    read_dataset_from_disk(datasetLeaf, curr->getPageID());
//    return datasetLeaf;

//}

//template <class T>
//MVPTree<T>::MVPTree(Dataset<T>* dataset, DistanceFunction<BasicArrayObject<T>>* _df, Pivot<T>* _pvt, size_t bf, size_t pl, size_t lc, size_t lpn, size_t fo, size_t ns)
//{

//    std::vector<datapoint_t<T>> vec;
//    df = _df;
//    pvt = _pvt;
//    BF = bf;
//    PL = pl;
//    LC = lc;
//    LPN = lpn;
//    FO = fo;
//    NS = ns;
//    leafNodeAccess = 0;

//    for(size_t i = 0; i < dataset->getCardinality(); i++)
//    {

//        datapoint_t<T> dp = datapoint_t<T>(dataset->getInstance(i), pl);
//        vec.push_back(dp);

//    }

//    build(vec);
//    initDisk();

//}

//template <class T>
//MVPNode<T>* MVPTree<T>::getRoot()
//{

//    return root;

//}

//template <class T>
//void MVPTree<T>::knn(BasicArrayObject<T> query, size_t k, std::vector<KnnEntryMVP<T>> &ans)
//{

//    df->resetStatistics();
//    leafNodeAccess = 0;
//    std::priority_queue<MVPPartition<T>, std::vector<MVPPartition<T>>, CompareMVPPartition<T>> nodeQueue;
//    std::priority_queue<KnnEntryMVP<T>, std::vector<KnnEntryMVP<T>>, std::greater<KnnEntryMVP<T>>> candidatesQueue;
//    std::priority_queue<KnnEntryMVP<T>, std::vector<KnnEntryMVP<T>>, std::less<KnnEntryMVP<T>>> resultQueue;
//    nodeQueue.push(MVPPartition<T>(getRoot(), 0.0, std::numeric_limits<double>::max()));
//    MVPNode<T>* node = nullptr;
//    MVPPartition<T> partition;
//    std::vector<cases> casesNodeVec = {R0, R1, R2, R3};
//    Dataset<T>* datasetLeaf;

//    while(!nodeQueue.empty() || candidatesQueue.size() > 0)
//    {

//        if(candidatesQueue.size() == 0)
//        {

//            partition = nodeQueue.top();
//            node = partition.node;
//            nodeQueue.pop();

//            if(node->isLeaf())
//            {

//                leafNodeAccess++;

//                datasetLeaf = readLeafNode(node);

//                for(size_t i = 0; i < datasetLeaf->getCardinality(); i++)
//                {

//                    candidatesQueue.push(KnnEntryMVP<T>(datasetLeaf->getFeatureVector(i), df->getDistance(datasetLeaf->getFeatureVector(i), query)));

//                }

//                delete datasetLeaf;

//                //                std::vector<BasicArrayObject<T>> data = node->purgeData();

//                //                for(size_t i = 0; i < data.size(); i++)
//                //                {

//                //                    candidatesQueue.push(KnnEntryMVP<T>(data[i], df->getDistance(query, data[i])));

//                //                }

//            }
//            else
//            {

//                cases sqCase;
//                double d_sq_p1 = df->getDistance(query, *node->getVP(0).key);
//                double d_sq_p2 = df->getDistance(query, *node->getVP(1).key);

//                if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 0))
//                    sqCase = R0;
//                else if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 0))
//                    sqCase = R1;
//                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 1))
//                    sqCase = R2;
//                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 1))
//                    sqCase = R3;
//                else
//                    sqCase = R4;


//                for(size_t i = 0; i < FO; i++)
//                {

//                    if(node->getChild(i) != nullptr)
//                    {

//                        nodeQueue.push(MVPPartition<T>(node->getChild(i),
//                                                       minDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1)),
//                                                       maxDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1))));

//                    }

//                }


//            }

//        }
//        else if((nodeQueue.size() > 0) && nodeQueue.top().min < candidatesQueue.top().distance)
//        {

//            partition = nodeQueue.top();
//            node = partition.node;
//            nodeQueue.pop();

//            if(node->isLeaf())
//            {

//                leafNodeAccess++;

//                datasetLeaf = readLeafNode(node);

//                for(size_t i = 0; i < datasetLeaf->getCardinality(); i++)
//                {

//                    candidatesQueue.push(KnnEntryMVP<T>(datasetLeaf->getFeatureVector(i), df->getDistance(datasetLeaf->getFeatureVector(i), query)));

//                }

//                delete datasetLeaf;

//                //                std::vector<BasicArrayObject<T>> data = node->purgeData();

//                //                for(size_t i = 0; i < data.size(); i++)
//                //                {

//                //                    candidatesQueue.push(KnnEntryMVP<T>(data[i], df->getDistance(query, data[i])));

//                //                }

//            }
//            else
//            {

//                cases sqCase;
//                double d_sq_p1 = df->getDistance(query, *node->getVP(0).key);
//                double d_sq_p2 = df->getDistance(query, *node->getVP(1).key);

//                if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 0))
//                    sqCase = R0;
//                else if(d_sq_p1 <= node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 0))
//                    sqCase = R1;
//                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 <= node->getSplit(1, 1))
//                    sqCase = R2;
//                else if(d_sq_p1 > node->getSplit(0, 0) && d_sq_p2 > node->getSplit(1, 1))
//                    sqCase = R3;
//                else
//                    sqCase = R4;

//                for(size_t i = 0; i < FO; i++)
//                {

//                    if(node->getChild(i) != nullptr)
//                    {

//                        nodeQueue.push(MVPPartition<T>(node->getChild(i),
//                                                       minDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1)),
//                                                       maxDist(sqCase, casesNodeVec[i], d_sq_p1, d_sq_p2, node->getSplit(0, 0), partition.max, node->getSplit(1, 0), node->getSplit(1, 1))));

//                    }

//                }

//            }

//        }
//        else
//        {

//            if(!resultQueue.empty() && !candidatesQueue.empty() && resultQueue.size() >= k && candidatesQueue.top().distance > resultQueue.top().distance)
//            {

//                break;

//            }

//            resultQueue.push(candidatesQueue.top());
//            candidatesQueue.pop();

//            while(resultQueue.size() > k)
//            {

//                resultQueue.pop();

//            }

//        }

//    }

//    ans = dequeueInOrderMVP_Results(resultQueue);
//    std::reverse(ans.begin(), ans.end());

//    while(!candidatesQueue.empty())
//    {

//        candidatesQueue.pop();

//    }

//    while(!resultQueue.empty())
//    {

//        resultQueue.pop();

//    }

//    while(!nodeQueue.empty())
//    {

//        nodeQueue.pop();

//    }

//}

//template <class T>
//size_t MVPTree<T>::getLeafNodeAccess()
//{

//    return leafNodeAccess;

//}


//#endif // MVPTREE_H
