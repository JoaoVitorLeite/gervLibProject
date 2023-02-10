#ifndef MVPTREE_H
#define MVPTREE_H

//#include <mvpnode.h>

//namespace mvp
//{

//    template <typename T, class F, class PVT, size_t BF = 2, size_t PL = 8, size_t LC = 30, size_t LPN = 2, size_t FO = 4, size_t NS = 2>
//    class MVPTree
//    {

//        private:
//            std::vector<datapoint_t<T,PL>> m_arrivals;
//            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* m_top;
//            size_t n_sync;

//            void LinkNodes(std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &nodes,
//                           std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes)const
//            {

//                for(auto iter = nodes.begin(); iter != nodes.end(); iter++)
//                {

//                    size_t i = iter->first;
//                    MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvpnode = iter->second;

//                    if(mvpnode != NULL)
//                    {

//                        for(size_t j = 0; j < FO; j++)
//                        {

//                            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* child = childnodes[i*FO+j];

//                            if(child != NULL)
//                            {

//                                mvpnode->SetChildNode(j, child);

//                            }

//                        }

//                    }

//                }

//            }


//            void ExpandNode(MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node,
//                            std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes, const size_t index)const
//            {

//                if(node != NULL)
//                {

//                    for(size_t i = 0; i < FO; i++)
//                    {

//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* child = node->GetChildNode(i);

//                        if(child != NULL)
//                        {

//                            childnodes[index*FO+i] = child;

//                        }

//                    }

//                }

//            }


//            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* ProcessNode(const size_t level,
//                                                             const size_t index,
//                                                             MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node,
//                                                             std::vector<datapoint_t<T,PL>> &points,
//                                                             std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
//                                                             std::map<size_t, std::vector<datapoint_t<T,PL>>*> &childpoints)
//            {

//                MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* retnode = node;

//                if(node == NULL)
//                {

//                    retnode = MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CreateNode(points, childpoints, level, index);

//                }
//                else
//                {

//                    retnode = node->AddDataPoints(points, childpoints, level, index);

//                }

//                if(retnode == NULL)
//                {

//                    throw std::runtime_error("failure to assign node");

//                }

//                return retnode;

//            }


//        public:
//            MVPTree()
//            {

//                m_top = NULL;
//                n_sync = 100;

//            }

//            ~MVPTree()
//            {



//            }

//            void setDistanceFunction(F* df_)
//            {

//                m_top->setDistanceFunction(df_);

//            }

//            void setPivotMethod(PVT* pvt_)
//            {

//                m_top->setPivotMethod(pvt_);

//            }

//            void Add(datapoint_t<T,PL> &item)
//            {

//                m_arrivals.push_back({item.id, item.key});

//                if(m_arrivals.size() >= n_sync)
//                {

//                    Add(m_arrivals);

//                }

//            }

//            void Add(std::vector<datapoint_t<T,PL>> &points)
//            {

//                if(points.empty())
//                {

//                    return;

//                }

//                std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> prevnodes, currnodes, childnodes;

//                if(m_top != NULL)
//                {

//                    currnodes[0] = m_top;

//                }

//                std::map<size_t, std::vector<datapoint_t<T,PL>>*> pnts, pnts2;
//                pnts[0] = &points;

//                size_t n = 0;

//                do{

//                    for(auto iter = pnts.begin(); iter != pnts.end(); iter++)
//                    {

//                        size_t index = iter->first;
//                        std::vector<datapoint_t<T,PL>>* list = iter->second;
//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvpnode = currnodes[index];
//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* newnode = ProcessNode(n, index, mvpnode, *list, childnodes, pnts2);

//                        if(newnode != mvpnode)
//                        {

//                            if(mvpnode != NULL)
//                            {

//                                delete mvpnode;

//                            }

//                            if(n == 0)
//                            {

//                                m_top = newnode;

//                            }

//                        }

//                        currnodes[index] = newnode;
//                        ExpandNode(newnode, childnodes, index);

//                        if(n > 0)
//                        {

//                            delete list;

//                        }

//                    }

//                    if(!prevnodes.empty())
//                    {

//                        LinkNodes(prevnodes, currnodes);

//                    }

//                    prevnodes = std::move(currnodes);
//                    currnodes = std::move(childnodes);
//                    pnts = std::move(pnts2);
//                    childnodes.clear();
//                    n += LPN;

//                } while(!pnts.empty());


//            }

//            void Sync()
//            {

//                if(m_arrivals.size() > 0)
//                {

//                    Add(m_arrivals);

//                }

//            }

//            const size_t Size()const
//            {

//                size_t count = 0;

//                std::queue<MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> nodes;

//                if(m_top != NULL)
//                {

//                    nodes.push(m_top);

//                }

//                while(!nodes.empty())
//                {

//                    MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* current = nodes.front();

//                    count += current->Size();

//                    for(size_t i = 0; i < FO; i++)
//                    {

//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* child = current->GetChildNode(i);

//                        if(child != NULL)
//                        {

//                            nodes.push(child);

//                        }

//                    }

//                    nodes.pop();

//                }

//                return count;

//            }

//            void Clear()
//            {

//                std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;

//                if(m_top != NULL)
//                {

//                    currnodes[0] = m_top;

//                }

//                size_t n = 0;

//                do{

//                    size_t n_nodes = pow(BF, n);
//                    size_t n_childnodes = pow(BF, n+LPN);

//                    for(auto iter = currnodes.begin(); iter != currnodes.end(); iter++)
//                    {

//                        size_t index = iter->first;
//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvpnode = iter->second;

//                        if(mvpnode != NULL)
//                        {

//                            ExpandNode(mvpnode, childnodes, index);
//                            delete mvpnode;

//                        }

//                    }

//                    currnodes = std::move(childnodes);
//                    n += LPN;

//                } while(!currnodes.empty());

//                m_top = NULL;

//            }

//            const std::vector<item_t<T>> Query(const BasicArrayObject<T> &target, const double radius)const
//            {

//                std::vector<item_t<T>> results;
//                std::map<size_t, std::vector<double>*> tdistances, tdistances2;
//                std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;

//                if(m_top != NULL)
//                {

//                    currnodes[0] = m_top;

//                }

//                size_t n = 0;

//                do{

//                    size_t n_nodes = pow(BF, n);
//                    size_t n_childnodes = pow(BF, n+LPN);

//                    for(auto const &iter : currnodes)
//                    {

//                        size_t node_index = iter.first;
//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvpnode = iter.second;
//                        std::vector<double>* tpath = tdistances[node_index];

//                        if(mvpnode)
//                        {

//                            mvpnode->TraverseNode(target, radius, childnodes, tpath, tdistances2, node_index, n, false, results);

//                        }

//                    }

//                    currnodes = std::move(childnodes);
//                    tdistances = std::move(tdistances2);
//                    n += LPN;

//                } while(!currnodes.empty());

//                return results;

//            }

//            const size_t DeletePoint(const BasicArrayObject<T> &target)
//            {

//                std::vector<item_t<T>> results;
//                std::map<size_t, std::vector<double>*> tdistances, tdistances2;
//                std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;

//                if(m_top != NULL)
//                {

//                    currnodes[0] = m_top;

//                }

//                size_t n = 0;

//                while(!currnodes.empty())
//                {

//                    size_t n_nodes = pow(BF, n);
//                    size_t n_childnodes = pow(BF, n+LPN);

//                    for(auto const &iter : currnodes)
//                    {

//                        size_t node_index = iter.first;
//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvpnode = iter.second;
//                        std::vector<double>* tpath = tdistances[node_index];

//                        if(mvpnode)
//                        {

//                            mvpnode->TraverseNode(target, 0.0, childnodes, tpath, tdistances2, node_index, n, true, results);

//                        }

//                    }

//                    currnodes = std::move(childnodes);
//                    tdistances = std::move(tdistances2);
//                    n += LPN;

//                }

//                return results.size();

//            }

//            void Print(std::ostream &ostrm)const //METODO NN ESTA FINALIZADO
//            {

//                std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;

//                if(m_top != NULL)
//                {

//                    currnodes[0] = m_top;

//                }

//                size_t n = 0;

//                do{

//                    ostrm << std::dec << "level = " << n << " (" << currnodes.size() << " nodes" << std::endl;;
//                    size_t n_nodes = pow(BF, n);
//                    size_t n_childnodes = pow(BF, n+LPN);
//                    for (auto const &iter : currnodes)
//                    {

//                        size_t node_index = iter.first;
//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;

//                        if (mvpnode)
//                        {

//                            ExpandNode(mvpnode, childnodes, node_index);

//                        }

//                    }

//                    ostrm << std::endl;
//                    currnodes = std::move(childnodes);
//                    n+= LPN;

//                } while(!currnodes.empty());

//                return;

//            }

//            size_t MemoryUsage()const
//            {

//                std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;

//                if(m_top != NULL)
//                {

//                    currnodes[0] = m_top;

//                }

//                size_t n_internal = 0, n_leaf = 0;

//                do{

//                    for(auto iter = currnodes.begin(); iter != currnodes.end(); iter++)
//                    {

//                        size_t node_index = iter->first;
//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node = iter->second;

//                        if (typeid(*node).hash_code() == typeid(MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>).hash_code())
//                        {

//                            n_internal++;

//                        }
//                        else
//                        {

//                            n_leaf++;

//                        }

//                        ExpandNode(node, childnodes, index);

//                    }

//                    currnodes = std::move(childnodes);

//                } while(!currnodes.empty());

//                return n_internal*sizeof(MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>)
//                        + n_leaf*sizeof(MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>)
//                        + sizeof(mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>);

//            }

//            void SetSync(size_t n)
//            {

//                n_sync = n;

//            }

//            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* getRoot()
//            {

//                return m_top;

//            }

//    };


//}

//**********************************************************************************************************************************

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <queue>
#include <ostream>
#include "mvpnode.h"
#include "datapoint.h"


/**
 * template arguments
 *  T  - key type defined in  key.hpp
 *  BF - branchfactor, no. of units to split space at each level of tree
 *  PL - pathlength, no. vantage point distances to maintain for each datapoint
 *  LC - leafcapacity, no. points in each leaf node
 *  LPN - levelspernode, no. levels per internal node in which the space is successively split into BF
 *  FO - fanout, no childnodes per internal node BF^LPN
 *  NS - numsplits, max no. split points for last level of internal node BF^(LPN-1)
 **/

namespace mvp {

    template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
    struct MVPPartition
    {


        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node;
        double min, max;

        MVPPartition()
        {



        }

        MVPPartition(MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node_, double min_, double max_)
        {

            node = node_;
            min = min_;
            max = max_;

        }

    };

    template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
    struct CompareMVPPartition
    {

        bool operator()(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> const& a, MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> const& b)
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

        T element;
        double distance;

        KnnEntryMVP(T& element_, double distance_)
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

    enum cases{R0, R1, R2, R3, R4};

    template<typename T,class F,class PVT, class DT,int BF=2, int PL=8, int LC=30, int LPN=2, int FO=4, int NS=2>
    class MVPTree {
    private:
        std::vector<datapoint_t<T,PL>> m_arrivals;

        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *m_top;

        F* df;

        size_t leafNodeAccess;

        int n_sync;

        void LinkNodes(std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &nodes,
                       std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes)const;
        void ExpandNode(MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node,
                        std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes, const int index)const;
        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* ProcessNode(const int level,
                                                   const int index,
                                                   MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node, std::vector<datapoint_t<T,PL>> &points,
                                                   std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                                                   std::map<int, std::vector<datapoint_t<T,PL>>*> &childpoints);
    public:
        MVPTree(F* distanceFunction)/*:m_top(NULL),n_sync(100)*/{
            n_sync = 100;
            m_top = new MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();
            df = distanceFunction;
            //m_top->setDistanceFunction(distanceFunction);
            //m_top->setPivotMethod(pivotMethod);
        };

        void Add(datapoint_t<T,PL> &item);

        void Add(std::vector<datapoint_t<T,PL>> &items);

        void Sync();

        const size_t Size()const;

        void Clear();

        const std::vector<item_t<T>> Query(const T &target, const double radius) const;

        const int DeletePoint(const T &target);

        void Print(std::ostream &ostrm)const;

        size_t MemoryUsage()const;

        void SetSync(int n);

        //JOAO
        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* getRoot(){ return m_top; }

        //JOAO
        F* getDistanceFunction()
        {

            return m_top->getDistanceFunction();

        }

        size_t getLeafNodeAccess(){ return leafNodeAccess; }

        double minDist(cases nodeCase, cases sqCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2);

        double maxDist(cases nodeCase, cases sqCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2);

        void knn(T query, size_t k, std::vector<KnnEntryMVP<T>> &ans);

    };

//    struct DistBrowsingEntry
//    {

//        double ct1, ct2;
//        int pivot_index;

//        DistBrowsingEntry()
//        {



//        }

//        DistBrowsingEntry(int pivot_index_, double ct1_, double ct2_)
//        {

//            pivot_index = pivot_index_;
//            ct1 = ct1_;
//            ct2 = ct2_;

//        }

//        bool operator<(const DistBrowsingEntry& item) const
//        {

//            bool ans;

//            if(ct1 != item.ct1)
//            {

//                ans = ct1 < item.ct1;

//            }
//            else
//            {

//                ans = ct2 < item.ct2;

//            }

//            return ans;

//        }

//        bool operator>(const DistBrowsingEntry& item) const
//        {

//            bool ans;

//            if(ct1 != item.ct1)
//            {

//                ans = ct1 > item.ct1;

//            }
//            else
//            {

//                ans = ct2 > item.ct2;

//            }

//            return ans;

//        }

//    };



}

/** -----------------------------------
 *
 *     MVPTree implementation
 *
 **/

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::LinkNodes(std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &nodes,
                                            std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes)const{
    for (auto iter=nodes.begin();iter!=nodes.end();iter++){
        int i = iter->first;
        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *mvpnode = iter->second;
        if (mvpnode != NULL){
            for (int j=0;j < FO;j++){
                MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *child = childnodes[i*FO+j];
                if (child != NULL) mvpnode->SetChildNode(j, child);
            }
        }
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::ExpandNode(MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node,
                                             std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                                             const int index)const{
    if (node != NULL){
        for (int i=0;i < FO;i++){
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *child = node->GetChildNode(i);
            if (child != NULL) childnodes[index*FO+i] = child;
        }
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::ProcessNode(const int level,
                                                                          const int index,
                                                                          MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node,
                                                                          std::vector<datapoint_t<T,PL>> &points,
                                                                          std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                                                                          std::map<int, std::vector<datapoint_t<T,PL>>*> &childpoints){
    MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *retnode = node;
    if (node == NULL){ // create new node
        retnode = MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CreateNode(points, childpoints, level, index);
    } else {           // node exists
        retnode = node->AddDataPoints(points, childpoints, level, index);
    }

    if (retnode == NULL)
        throw std::runtime_error("failure to assign node");

    return retnode;
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Add(datapoint_t<T,PL> &dp){
    m_arrivals.push_back({ dp.id, dp.key });
    if ((int)m_arrivals.size() >= n_sync){
        Add(m_arrivals);
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Add(std::vector<datapoint_t<T,PL>> &points){
    if (points.empty()) return;

    std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> prevnodes, currnodes, childnodes;
    if (m_top != NULL)
        currnodes[0] = m_top;

    std::map<int, std::vector<datapoint_t<T,PL>>*> pnts, pnts2;
    pnts[0] = &points;

    int n = 0;
    do {
        for (auto iter=pnts.begin();iter!=pnts.end();iter++){
            int index = iter->first;
            std::vector<datapoint_t<T,PL>> *list = iter->second;
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *mvpnode = currnodes[index];
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *newnode = ProcessNode(n, index, mvpnode, *list, childnodes, pnts2);
            if (newnode != mvpnode){
                if (mvpnode != NULL){
                    delete mvpnode;
                }
                if (n == 0) m_top = newnode;
            }
            currnodes[index] = newnode;
            ExpandNode(newnode, childnodes, index);

            if (n > 0) delete list;
        }

        if (!prevnodes.empty()) {
            LinkNodes(prevnodes, currnodes);
        }
        prevnodes = move(currnodes);
        currnodes = move(childnodes);
        pnts = move(pnts2);
        childnodes.clear();
        n += LPN;
    } while (!pnts.empty());

}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Sync(){
    if (m_arrivals.size() > 0) {
        Add(m_arrivals);
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
const size_t mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Size()const{
    size_t count = 0;

    std::queue<MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> nodes;
    if (m_top != NULL) nodes.push(m_top);

    while (!nodes.empty()){
        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *current = nodes.front();

        count += current->Size();

        for (int i=0;i < FO;i++){
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *child = current->GetChildNode(i);
            if (child != NULL) nodes.push(child);
        }

        nodes.pop();
    }
    return count;
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Clear(){
    std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
    if (m_top != NULL) currnodes[0] = m_top;

    int n = 0;
    do {
        int n_nodes = pow(BF, n);
        int n_childnodes = pow(BF, n+LPN);
        for (auto iter=currnodes.begin();iter!=currnodes.end();iter++){
            int index = iter->first;
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *mvpnode = iter->second;

            if (mvpnode != NULL){
                ExpandNode(mvpnode, childnodes, index);
                delete mvpnode;
            }
        }
        currnodes = move(childnodes);
        n += LPN;

    } while (!currnodes.empty());
    m_top = NULL;
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
const std::vector<mvp::item_t<T>> mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Query(const T &target, const double radius) const{
    std::vector<item_t<T>> results;
    std::map<int, std::vector<double>*> tdistances, tdistances2;
    std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
    if (m_top != NULL) currnodes[0] = m_top;

    int n = 0;
    do {
        int n_nodes = pow(BF, n);
        int n_childnodes = pow(BF, n+LPN);

        for (auto const &iter : currnodes){
            int node_index = iter.first;
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;
            std::vector<double> *tpath = tdistances[node_index];
            if (mvpnode)
                mvpnode->TraverseNode(target, radius, childnodes,  tpath, tdistances2, node_index, n, false, results);
        }
        currnodes = move(childnodes);
        tdistances = move(tdistances2);
        n += LPN;
    } while (!currnodes.empty());

    return results;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const int mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::DeletePoint(const T &target){
    std::vector<item_t<T>> results;

    std::map<int, std::vector<double>*> tdistances, tdistances2;
    std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
    if (m_top != NULL) currnodes[0] = m_top;

    int n = 0;
    while (!currnodes.empty()) {
        int n_nodes = pow(BF, n);
        int n_childnodes = pow(BF, n+LPN);
        for (auto const &iter : currnodes){
            int node_index = iter.first;
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;
            std::vector<double> *tpath = tdistances[node_index];
            if (mvpnode){
                mvpnode->TraverseNode(target, 0, childnodes, tpath, tdistances2, node_index, n, true, results);
            }
        }
        currnodes = move(childnodes);
        tdistances = move(tdistances2);
        n += LPN;
    };

    return (int)results.size();
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Print(std::ostream &ostrm)const {
    std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
    if (m_top != NULL) currnodes[0] = m_top;

    int n=0;
    do {

        ostrm << std::dec << "level = " << n << " (" << currnodes.size() << " nodes" << std::endl;;
        int n_nodes = pow(BF, n);
        int n_childnodes = pow(BF, n+LPN);
        for (auto const &iter : currnodes){
            int node_index = iter.first;
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *mvpnode = iter.second;
            if (mvpnode){
                ExpandNode(mvpnode, childnodes, node_index);
            }
        }
        ostrm << std::endl;
        currnodes = move(childnodes);
        n+= LPN;
    } while (!currnodes.empty());

    return;
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
size_t mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::MemoryUsage()const{
    std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> currnodes, childnodes;
    if (m_top != NULL) currnodes[0] = m_top;

    int n_internal=0, n_leaf=0;
    do {
        for (auto iter=currnodes.begin();iter!=currnodes.end();iter++){
            int index = iter->first;
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node = iter->second;

            if (typeid(*node).hash_code() == typeid(MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>).hash_code())
                n_internal++;
            else
                n_leaf++;
            ExpandNode(node, childnodes, index);
        }
        currnodes = move(childnodes);
    } while (!currnodes.empty());

    return  n_internal*sizeof(MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>)
        + n_leaf*sizeof(MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>)
        + sizeof(mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>);
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::SetSync(int n){
    n_sync = n;
}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
double mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::minDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2)
{

//    return 0.0;
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
                    ans = fabs(d_sq_p2 - p2_ct1);
                    break;
                }

                case R2:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
                    break;
                }

                case R3:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
                    break;
                }

                case R4:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
                    break;
                }

            }
            break;

        }

        case R1:{

            switch(nodeCase)
            {

                case R0:{
                    ans = fabs(d_sq_p2 - p2_ct1);
                    break;
                }

                case R1:{
                    ans = 0.0;
                    break;
                }

                case R2:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
                    break;
                }

                case R3:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2)));
                    break;
                }

                case R4:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
                    break;
                }

            }
            break;

        }

        case R2:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));//
                    break;
                }

                case R1:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));
                    break;
                }

                case R2:{
                    ans = 0.0;
                    break;
                }

                case R3:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2));
                    break;
                }

                case R4:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
                    break;
                }

            }
            break;

        }

        case R3:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));
                    break;
                }

                case R1:{
                    ans = std::max(fabs(d_sq_p1 - p1_ct1), std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1)));
                    break;
                }

                case R2:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct2));
                    break;
                }

                case R3:{
                    ans = 0.0;
                    break;
                }

                case R4:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
                    break;
                }

            }
            break;

        }

        case R4:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1));
                    break;
                }

                case R1:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct1), fabs(d_sq_p2 - p2_ct1));
                    break;
                }

                case R2:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
                    break;
                }

                case R3:{
                    ans = std::min(fabs(d_sq_p1 - p1_ct2), fabs(d_sq_p2 - p2_ct2));
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

    //std::cout << "MIN DIST = " << ans << " / " << sqCase << " / " << nodeCase << " / " << d_sq_p1 << "/" << d_sq_p2 << "/" << p1_ct1 << "/" << p1_ct2 << "/" << p2_ct1 << "/" << p2_ct2 << "\n";
    return ans;

}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
double mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::maxDist(cases sqCase, cases nodeCase, double d_sq_p1, double d_sq_p2, double p1_ct1, double p1_ct2, double p2_ct1, double p2_ct2)
{

//    return 0.0;
    double ans;

    switch(sqCase)
    {

        case R0:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
                    break;
                }

                case R1:{
                    ans = d_sq_p1 + p1_ct1;
                    break;
                }

                case R2:{
                    ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
                    break;
                }

                case R3:{
                    ans = d_sq_p1 + p1_ct2;
                    break;
                }

                case R4:{
                    ans = -1.0; //!!!!!
                    break;
                }

            }
            break;

        }

        case R1:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
                    break;
                }

                case R1:{
                    ans = d_sq_p1 + p1_ct1;
                    break;
                }

                case R2:{
                    ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
                    break;
                }

                case R3:{
                    ans = d_sq_p1 + p1_ct2;
                    break;
                }

                case R4:{
                    ans = -1.0; //!!!!!
                    break;
                }

            }
            break;

        }

        case R2:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
                    break;
                }

                case R1:{
                    ans = d_sq_p1 + p1_ct1;
                    break;
                }

                case R2:{
                    ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
                    break;
                }

                case R3:{
                    ans = d_sq_p1 + p1_ct2;
                    break;
                }

                case R4:{
                    ans = -1.0; //!!!!!
                    break;
                }

            }
            break;

        }

        case R3:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
                    break;
                }

                case R1:{
                    ans = d_sq_p1 + p1_ct1;
                    break;
                }

                case R2:{
                    ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
//                    ans = d_sq_p1 + p1_ct2;
                    break;
                }

                case R3:{
                    ans = d_sq_p1 + p1_ct2;
                    break;
                }

                case R4:{
                    ans = -1.0; //!!!!!
                    break;
                }

            }
            break;

        }

        case R4:{

            switch(nodeCase)
            {

                case R0:{
                    ans = std::min(d_sq_p1 + p1_ct1, d_sq_p2 + p2_ct1);
                    break;
                }

                case R1:{
                    ans = ans = d_sq_p1 + p1_ct1;
                    break;
                }

                case R2:{
                    ans = std::min(d_sq_p1 + p1_ct2, d_sq_p2 + p2_ct2);
                    break;
                }

                case R3:{
                    ans = d_sq_p1 + p1_ct2;
                    break;
                }

                case R4:{
                    ans = -1.0; //!!!!!
                    break;
                }

            }
            break;

        }

    }

    //std::cout << "MAX DIST = " << ans << " / " << sqCase << " / " << nodeCase << " / " << d_sq_p1 << "/" << d_sq_p2 << "/" << p1_ct1 << "/" << p1_ct2 << "/" << p2_ct1 << "/" << p2_ct2 << "\n";
    return ans;

}

template<typename T,class F,class PVT,class DT,int BF,int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPTree<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::knn(T query, size_t k, std::vector<KnnEntryMVP<T>> &ans)
{

    df->resetStatistics();
    leafNodeAccess = 0;
//    size_t cnt = 0;

    std::priority_queue<MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>, std::vector<MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>>, CompareMVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>> nodeQueue;
    std::priority_queue<KnnEntryMVP<T>, std::vector<KnnEntryMVP<T>>, std::greater<KnnEntryMVP<T>>> candidatesQueue;
    std::priority_queue<KnnEntryMVP<T>, std::vector<KnnEntryMVP<T>>, std::less<KnnEntryMVP<T>>> resultQueue;
    nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(getRoot(), 0.0, std::numeric_limits<double>::max()));
    MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node = nullptr;
    MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> partition;
    MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* leaf;
    MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* internal;

    while(!nodeQueue.empty() || candidatesQueue.size() > 0)
    {

        //std::cout << "CAND QUEUE SIZE = " << candidatesQueue.size() << std::endl;
//        if(!nodeQueue.empty() && candidatesQueue.size() > 0 && !resultQueue.empty())
//            std::cout << "AUX = " << nodeQueue.top().min << " / " << candidatesQueue.top().distance << " / " << resultQueue.size() << "\n";

        if(candidatesQueue.size() == 0)
        {

            partition = nodeQueue.top();
            node = partition.node;
            nodeQueue.pop();

            if(node->IsLeaf())
            {

                //std::cout << "LEAF\n";
                leafNodeAccess++;
                leaf = static_cast<MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*>(node);
                std::vector<datapoint_t<T,PL>> datapointsLeaf = leaf->GetDataPoints();
                std::vector<vp_t<T>> datavpsLeaf = leaf->GetVantagePoints();

//                cnt += datapointsLeaf.size();
//                cnt += datavpsLeaf.size();

//                std::cout << datapointsLeaf.size() + datavpsLeaf.size() << "\n";

                for(size_t x = 0; x < datapointsLeaf.size(); x++)
                {

                    candidatesQueue.push(KnnEntryMVP<T>(datapointsLeaf[x].key, df->getDistance(query, datapointsLeaf[x].key)));

                }

                for(size_t x = 0; x < datavpsLeaf.size(); x++)
                {

                    candidatesQueue.push(KnnEntryMVP<T>(datavpsLeaf[x].key, df->getDistance(query, datavpsLeaf[x].key)));

                }

            }
            else
            {

                //std::cout << "INTERNAL\n";
                internal = static_cast<MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*>(node);
                std::array<std::array<double,NS>,LPN> splits = internal->getSplits();
                std::vector<vp_t<T>> datavpsInternal = internal->GetVantagePoints();

//                for(size_t x = 0; x < datavpsInternal.size(); x++)
//                {

//                    candidatesQueue.push(KnnEntryMVP<T>(datavpsInternal[x].key, df->getDistance(query, datavpsInternal[x].key)));

//                }

                cases sqCase;
                double d_sq_p1 = df->getDistance(query, datavpsInternal[0].key);
                double d_sq_p2 = df->getDistance(query, datavpsInternal[1].key);

//                cnt += datavpsInternal.size();

                candidatesQueue.push(KnnEntryMVP<T>(datavpsInternal[0].key, d_sq_p1));
                candidatesQueue.push(KnnEntryMVP<T>(datavpsInternal[1].key, d_sq_p2));

                if(d_sq_p1 <= splits[0][0] && d_sq_p2 <= splits[1][0])
                    sqCase = R0;
                else if(d_sq_p1 <= splits[0][0] && d_sq_p2 > splits[1][0])
                    sqCase = R1;
                else if(d_sq_p1 > splits[0][0] && d_sq_p2 <= splits[1][1])
                    sqCase = R2;
                else if(d_sq_p1 > splits[0][0] && d_sq_p2 > splits[1][1])
                    sqCase = R3;
                else
                    sqCase = R4;

                if(node->GetChildNode(0) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(0), minDist(sqCase, R0, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R0, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));
                if(node->GetChildNode(1) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(1), minDist(sqCase, R1, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R1, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));
                if(node->GetChildNode(2) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(2), minDist(sqCase, R2, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R2, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));
                if(node->GetChildNode(3) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(3), minDist(sqCase, R3, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R3, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));

            }

        }
        else if((nodeQueue.size() > 0) && nodeQueue.top().min < candidatesQueue.top().distance)
        {

            partition = nodeQueue.top();
            node = partition.node;
            nodeQueue.pop();

            if(node->IsLeaf())
            {

                leafNodeAccess++;
                leaf = static_cast<MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*>(node);
                std::vector<datapoint_t<T,PL>> datapointsLeaf = leaf->GetDataPoints();
                std::vector<vp_t<T>> datavpsLeaf = leaf->GetVantagePoints();

//                cnt += datapointsLeaf.size();
//                cnt += datavpsLeaf.size();

                //std::cout << datapointsLeaf.size() + datavpsLeaf.size() << "\n";

                for(size_t x = 0; x < datapointsLeaf.size(); x++)
                {

                    candidatesQueue.push(KnnEntryMVP<T>(datapointsLeaf[x].key, df->getDistance(query, datapointsLeaf[x].key)));

                }

                for(size_t x = 0; x < datavpsLeaf.size(); x++)
                {

                    candidatesQueue.push(KnnEntryMVP<T>(datavpsLeaf[x].key, df->getDistance(query, datavpsLeaf[x].key)));

                }

            }
            else
            {

                internal = static_cast<MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*>(node);
                std::array<std::array<double,NS>,LPN> splits = internal->getSplits();
                std::vector<vp_t<T>> datavpsInternal = internal->GetVantagePoints();

//                for(size_t x = 0; x < datavpsInternal.size(); x++)
//                {

//                    candidatesQueue.push(KnnEntryMVP<T>(datavpsInternal[x].key, df->getDistance(query, datavpsInternal[x].key)));

//                }
//                cnt += datavpsInternal.size();


                cases sqCase;
                double d_sq_p1 = df->getDistance(query, datavpsInternal[0].key);
                double d_sq_p2 = df->getDistance(query, datavpsInternal[1].key);

                candidatesQueue.push(KnnEntryMVP<T>(datavpsInternal[0].key, d_sq_p1));
                candidatesQueue.push(KnnEntryMVP<T>(datavpsInternal[1].key, d_sq_p2));

                if(d_sq_p1 <= splits[0][0] && d_sq_p2 <= splits[1][0])
                    sqCase = R0;
                else if(d_sq_p1 <= splits[0][0] && d_sq_p2 > splits[1][0])
                    sqCase = R1;
                else if(d_sq_p1 > splits[0][0] && d_sq_p2 <= splits[1][1])
                    sqCase = R2;
                else if(d_sq_p1 > splits[0][0] && d_sq_p2 > splits[1][1])
                    sqCase = R3;
                else
                    sqCase = R4;

                if(node->GetChildNode(0) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(0), minDist(sqCase, R0, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R0, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));
                if(node->GetChildNode(1) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(1), minDist(sqCase, R1, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R1, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));
                if(node->GetChildNode(2) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(2), minDist(sqCase, R2, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R2, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));
                if(node->GetChildNode(3) != NULL) nodeQueue.push(MVPPartition<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>(node->GetChildNode(3), minDist(sqCase, R3, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1]), maxDist(sqCase, R3, d_sq_p1, d_sq_p2, splits[0][0], partition.max, splits[1][0], splits[1][1])));

            }

        }
        else
        {

            //std::cout << "RES\n";

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

//    std::cout << "DIST IN = " << cnt << "\n";

}


#endif // MVPTREE_H
