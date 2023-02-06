#ifndef MVPNODE_H
#define MVPNODE_H

//#include <datapoint.h>
//#include <Hermes.h>
//#include <Pivots.h>
//#include <Dataset.h>
//#include <array>

///**
// * template arguments
// *  T   - BasicArrayObject type (double, std::vector<char>, float, ...)
// *  F   - Distance Function
// *  PVT - Pivot method
// *  BF  - branchfactor, no. of units to split space at each level of tree
// *  PL  - pathlength, no. vantage point distances to maintain for each datapoint
// *  LC  - leafcapacity, no. points in each leaf node
// *  LPN - levelspernode, no. levels per internal node in which the space is successively split into BF
// *  FO  - fanout, no childnodes per internal node BF^LPN
// *  NS  - numsplits, max no. split points for last level of internal node BF^(LPN-1)
// **/

//namespace mvp
//{

//    template <typename T, class F, class PVT, size_t BF = 2, size_t PL = 8, size_t LC = 30, size_t LPN = 2, size_t FO = 4, size_t NS = 2>
//    class MVPNode
//    {

//        protected:
//            size_t m_nvps;
//            F* df = new F();
//            PVT* pvt = new PVT();
//            std::array<vp_t<T>, LPN> m_vps;

//            void SelectVantagePoints(std::vector<datapoint_t<T,PL>> &points)
//            {


////                if(this->m_nvps < LPN && points.size() > 0){

////                    std::vector<BasicArrayObject<T>> vec;
////                    size_t id = 0;
//////                    std::transform(points.begin(), points.end(), vec.begin(), [&id](datapoint_t<T,PL> dp){
//////                        BasicArrayObject<T> b = dp.key;
//////                        b.setOID(id++);
//////                        return b;
//////                    });
////                    for(datapoint_t<T,PL> &dp : points)
////                    {
////                        BasicArrayObject<T> b = dp.key;
////                        b.setOID(id++);
////                        vec.push_back(b);

////                    }

////                    Dataset<T>* dataset = new Dataset<T>(vec, vec.size(), vec[0].size());
////                    pvt->generatePivots(dataset, this->df, LPN);

////                    this->m_nvps = 0;
////                    while(this->m_nvps < LPN)
////                    {

////                        size_t pos = pvt->getPivot(this->m_nvps)->getOID();
////                        this->m_vps[this->m_nvps++] = {points[pos].id, points[pos].key};

////                        mvp::datapoint_t<T,PL> temp = points[points.size() - 1];
////                        points[pos] = points[points.size() - 1];
////                        points[points.size() - 1] = temp;
////                        points.pop_back();

////                    }

////                }

//                /////

//                if (this->m_nvps == 0 && points.size() > 0){
//                    std::cout << "VANTAGE = " << points.back().id << "\n";
//                    this->m_vps[this->m_nvps++] = { points.back().id, points.back().key };
//                    points.pop_back();
//                }

//                int slimit = 10;
//                while (this->m_nvps < LPN && points.size() > 0){
//                    int max_pos = 0;
//                    double maxd = 0;
//                    slimit = ((int)points.size() <= slimit) ? points.size() : slimit;
//                    for (int i=0;i < slimit;i++){
//                        double d = df->getDistance(points[i].key, this->m_vps[this->m_nvps-1].key);
//                        //double d = points[i].distance(this->m_vps[this->m_nvps-1]);
//                        if (d > maxd){
//                            max_pos = i;
//                            maxd = d;
//                        }
//                    }
//                    this->m_vps[this->m_nvps++] = { points[max_pos].id, points[max_pos].key };
//                    std::cout << "VANTAGE = " << points[max_pos].id << "\n";
//                    //JOAO
//                    mvp::datapoint_t<T,PL> temp = points[points.size()-1];
//                    points[max_pos] = points[points.size()-1];
//                    points[points.size()-1] = temp;
//                    //
//                    points.pop_back();
//                }

//            }

//            void MarkPointDistances(std::vector<datapoint_t<T,PL>> &points, const size_t level)
//            {

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    if(level + i < PL)
//                    {

//                        for(mvp::datapoint_t<T,PL> &pnt : points)
//                        {

//                            pnt.dists[level + i] = df->getDistance(pnt.key, this->m_vps[i].key);

//                        }

//                    }

//                }

//            }

//        public:
//            MVPNode()
//            {

//                m_nvps = 0;

//            }

//            virtual ~MVPNode(){};

//            void setDistanceFunction(F* df_)
//            {

//                df = df_;

//            }

//            void setPivotMethod(PVT* pvt_)
//            {

//                pvt = pvt_;

//            }

//            static MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* CreateNode(std::vector<datapoint_t<T,PL>> &points,
//                                                                   std::map<size_t, std::vector<datapoint_t<T,PL>>*> &childpoints,
//                                                                   size_t level,
//                                                                   size_t index);
////            {

////                MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node = NULL;

////                if(points.size() <= (LC + LPN))
////                {

////                    node = new MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();

////                }
////                else
////                {

////                    node = new MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();

////                }

////                node = node->AddDataPoints(points, childpoints, level, index);

////                if(node == NULL)
////                {

////                    throw std::runtime_error("unable to create node");

////                }

////                return node;

////            }

//            virtual MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
//                                                                       std::map<size_t, std::vector<datapoint_t<T,PL>>*> &childpoints,
//                                                                       const size_t level,
//                                                                       const size_t index) = 0;

//            virtual const size_t nbytes()const = 0;


//            virtual const size_t Size()const = 0;


//            virtual const bool IsLeaf()const = 0;

//            virtual void SetChildNode(const size_t n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node) = 0;

//            virtual MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* GetChildNode(const size_t n)const = 0;

//            virtual const std::vector<vp_t<T>> GetVantagePoints()const = 0;

//            virtual const std::vector<datapoint_t<T,PL>> GetDataPoints()const = 0;

//            virtual void FilterDataPoints(const BasicArrayObject<T> &target,
//                                          std::vector<double> &tpath,
//                                          const double radius,
//                                          const size_t level,
//                                          const bool delete_points,
//                                          std::vector<item_t<T>> &results) = 0;

//            virtual void TraverseNode(const BasicArrayObject<T> &target,
//                                      const double radius,
//                                      std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
//                                      std::vector<double> *tpath,
//                                      std::map<size_t, std::vector<double>*> &child_tdistances,
//                                      const size_t index,
//                                      const size_t level,
//                                      const bool delete_points,
//                                      std::vector<item_t<T>> &results) = 0;


//            virtual const std::vector<datapoint_t<T,PL>> PurgeDataPoints() = 0;

//    };

//    bool CompareDistance(const double a, const double b, const bool less)
//    {

//        return (less) ? (a <= b) : (a > b);

//    }

//    template <typename T, class F, class PVT, size_t BF = 2, size_t PL = 8, size_t LC = 30, size_t LPN = 2, size_t FO = 4, size_t NS = 2>
//    class MVPInternal : public MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>
//    {

//        private:
//            std::array<MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*, FO> m_childnodes;

//            //LPN X NS
//            std::array<std::array<double, NS>, LPN> m_splits;

//            void CalcSplitPoints(const std::vector<double> &dists, size_t n, size_t split_index)
//            {

//                size_t lengthM = BF - 1;

//                if(dists.size() > 0)
//                {

//                    if(m_splits[n][split_index*lengthM] == -1)
//                    {

//                        std::vector<double> tmpdists = dists;
//                        std::sort(tmpdists.begin(), tmpdists.end()); //TROCAR POR nth_element
//                        double factor = (double)tmpdists.size()/(double)BF;

//                        for(size_t i = 0; i < lengthM; i++)
//                        {

//                            double pos = (i+1)*factor;
//                            int lo = floor(pos);
//                            int hi = (pos <= tmpdists.size()-1) ? ceil(pos) : 0;
//                            m_splits[n][split_index*lengthM+i] = (tmpdists[lo] + tmpdists[hi])/2.0;

//                        }

//                    }

//                }

//            }

//            std::vector<double> CalcPointDistances(vp_t<T> &vp, std::vector<datapoint_t<T,PL>> &points)
//            {

//                std::vector<double> results;

//                for(datapoint_t<T,PL> &dp : points)
//                {

//                    results.push_back(this->df->getDistance(vp.key, dp.key));

//                }

//                return results;

//            }

//            std::vector<datapoint_t<T,PL>>* CullPoints(std::vector<datapoint_t<T,PL>> &list, std::vector<double> &dists, double split, bool less)
//            {

//                std::vector<datapoint_t<T,PL>>* results = new std::vector<datapoint_t<T,PL>>();
//                auto list_iter = list.begin();
//                auto dist_iter = dists.begin();

//                while(list_iter != list.end() && dist_iter != dists.end())
//                {

//                    if(CompareDistance(*dist_iter, split, less))
//                    {

//                        results->push_back(*list_iter);
//                        list_iter = list.erase(list_iter);
//                        dist_iter = dists.erase(dist_iter);

//                    }
//                    else
//                    {

//                        list_iter++;
//                        dist_iter++;

//                    }

//                }

//                if(results->size() > 0)
//                {

//                    return results;

//                }

//                delete results;
//                return NULL;

//            }

//            void CollatePoints(std::vector<datapoint_t<T,PL>> &points, std::map<size_t, std::vector<datapoint_t<T,PL>>*> &childpoints, const size_t level, const size_t index)
//            {

//                std::map<size_t, std::vector<datapoint_t<T,PL>>*> pnts, pnts2;
//                pnts[0] = &points;

//                size_t lengthM = BF - 1, n = 0;

//                do{

//                    for(auto iter = pnts.begin(); iter != pnts.end(); iter++)
//                    {

//                        size_t node_index = iter->first;
//                        std::vector<datapoint_t<T,PL>>* list = iter->second;

//                        std::vector<double> dists = CalcPointDistances((this->m_vps[n]), *list);

//                        if(dists.size() > 0)
//                        {

//                            CalcSplitPoints(dists, n, node_index);

//                            double m;
//                            std::vector<datapoint_t<T,PL>>* culledpts = NULL;

//                            for(size_t j = 0; j < lengthM; j++)
//                            {

//                                m = m_splits[n][node_index*lengthM+j];

//                                culledpts = CullPoints(*list, dists, m, true);

//                                if(culledpts != NULL)
//                                {

//                                    pnts2[node_index*BF+j] = culledpts;

//                                }

//                            }

//                            m = m_splits[n][node_index*lengthM+lengthM-1];
//                            culledpts = CullPoints(*list, dists, m, false);

//                            if(culledpts != NULL)
//                            {

//                                pnts2[node_index*BF+BF-1] = culledpts;

//                            }

//                        }

//                        if(list->size() > 0)
//                        {

//                            throw std::length_error("not fully collated");

//                        }

//                        if(n > 0)
//                        {

//                            delete list;

//                        }

//                    }

//                    pnts = std::move(pnts2);
//                    n++;

//                } while(n < LPN);

//                for(auto iter = pnts.begin(); iter != pnts.end(); iter++)
//                {

//                    size_t i = iter->first;
//                    std::vector<datapoint_t<T,PL>>* list = iter->second;

//                    if(list != NULL)
//                    {

//                        childpoints[index*FO+i] = list;

//                    }

//                }

//            }

//        public:
//            MVPInternal()
//            {

//                for(size_t i = 0; i < FO; i++)
//                {

//                    m_childnodes[i] = NULL;

//                }

//                for(size_t i = 0; i < LPN; i++)
//                {

//                    for(size_t j = 0; j < NS; j++)
//                    {

//                        m_splits[i][j] = -1.0;

//                    }

//                }

//            }

//            ~MVPInternal()
//            {



//            }

//            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
//                                                               std::map<size_t, std::vector<datapoint_t<T,PL>>*> &childpoints,
//                                                               const size_t level,
//                                                               const size_t index)
//            {

//                this->SelectVantagePoints(points);
//                this->MarkPointDistances(points, level);

//                if(this->m_nvps < LPN)
//                {

//                    throw std::invalid_argument("too few points for internal node");

//                }

//                CollatePoints(points, childpoints, level, index);
//                points.clear();
//                return this;

//            }


//            const size_t nbytes()const
//            {

//                return sizeof(MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>);

//            }

//            const size_t Size()const
//            {

//                size_t total = 0;

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    if(this->m_vps[i].active)
//                    {

//                        total++;

//                    }

//                }

//                return total;

//            }

//            const bool IsLeaf()const
//            {

//                return false;

//            }

//            void SetChildNode(const size_t n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node)
//            {

//                if(n < 0 || n >= FO)
//                {

//                    throw std::invalid_argument("index out of range");

//                }

//                m_childnodes[n] = node;

//            }

//            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* GetChildNode(const size_t n)const
//            {

//                if(n < 0 || n >= FO)
//                {

//                    throw std::invalid_argument("index out of range");

//                }

//                return m_childnodes[n];

//            }

//            const std::vector<vp_t<T>> GetVantagePoints()const
//            {

//                std::vector<vp_t<T>> results;

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    results.push_back(this->m_vps[i]);

//                }

//                return results;

//            }

//            const std::vector<datapoint_t<T,PL>> GetDataPoints()const
//            {

//                std::vector<datapoint_t<T,PL>> results;
//                return results;

//            }

//            void FilterDataPoints(const BasicArrayObject<T> &target,
//                                  std::vector<double> &tpath,
//                                  const double radius,
//                                  const size_t level,
//                                  const bool delete_points,
//                                  std::vector<item_t<T>> &results)
//            {

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    double d = this->df->getDistance(this->m_vps[i].key, target);
//                    tpath.push_back(d);

//                    if(this->m_vps[i].active && d <= radius)
//                    {

//                        vp_t<T> vp = this->m_vps[i];

//                        if(delete_points)
//                        {

//                            this->m_vps[i].active = false;

//                        }

//                        results.push_back({vp.id, vp.key});

//                    }

//                }

//            }

//            void TraverseNode(const BasicArrayObject<T> &target,
//                              const double radius,
//                              std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
//                              std::vector<double>* tpath,
//                              std::map<size_t, std::vector<double>*> &child_tdistances,
//                              const size_t index,
//                              const size_t level,
//                              const bool delete_points,
//                              std::vector<item_t<T>> &results)
//            {

//                size_t lengthM = BF - 1, n = 0;
//                bool *currnodes = new bool[1];
//                currnodes[0] = true;

//                std::vector<double> *path = (tpath != NULL) ? tpath : new std::vector<double>();

//                FilterDataPoints(target, *path, radius, level, delete_points, results);

//                do{

//                    size_t n_nodes = pow(BF, n), n_childnodes = pow(BF, n+1);
//                    bool *nextnodes = new bool[n_childnodes];

//                    for(size_t i = 0; i < n_childnodes; i++)
//                    {

//                        nextnodes[i] = false;

//                    }

//                    double d = path->at(level+n);

//                    for(size_t node_index = 0; node_index < n_nodes; node_index++)
//                    {

//                        if(currnodes[node_index])
//                        {

//                            if(m_splits[n][node_index*lengthM] >= 0)
//                            {

//                                double m = m_splits[n][node_index*lengthM];

//                                for(size_t j = 0; j < lengthM; j++)
//                                {

//                                    m = m_splits[n][node_index*lengthM+j];

//                                    if(d <= m + radius)
//                                    {

//                                        nextnodes[node_index*BF+j] = true;

//                                    }

//                                }

//                                if(d > m - radius)
//                                {

//                                    nextnodes[node_index*BF+BF-1] = true;

//                                }

//                            }

//                        }

//                    }

//                    delete [] (currnodes);
//                    currnodes = nextnodes;
//                    n++;
//                } while(n < LPN);

//                for(size_t i = 0; i < FO; i++)
//                {

//                    if(currnodes[i])
//                    {

//                        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* child = GetChildNode(i);

//                        if(child != NULL)
//                        {

//                            childnodes[index*FO+i] = child;

//                            std::vector<double>* tmppath = new std::vector<double>();

//                            for(size_t j = 0; j < path->size(); j++)
//                            {

//                                tmppath->push_back(path->at(j));

//                            }

//                            child_tdistances[index*FO+i] = tmppath;

//                        }

//                    }

//                }

//                delete [] (currnodes);
//                delete path;

//            }

//            const std::vector<datapoint_t<T,PL>> PurgeDataPoints()
//            {

//                std::vector<datapoint_t<T,PL>> results;

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    if(this->m_vps[i].active)
//                    {

//                        vp_t<T> vp = this->m_vps[i];
//                        results.push_back({vp.id, vp.key});

//                    }

//                }

//                return results;

//            }

//    };


//    template <typename T, class F, class PVT, size_t BF = 2, size_t PL = 8, size_t LC = 30, size_t LPN = 2, size_t FO = 4, size_t NS = 2>
//    class MVPLeaf : public MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>
//    {

//        private:
//            std::array<std::array<double, LC>, LPN> m_pdists;
//            size_t m_npoints;
//            std::array<datapoint_t<T, PL>, LC> m_points;

//            void MarkLeafDistances(std::vector<datapoint_t<T,PL>> &points)
//            {

//                if(m_npoints + points.size() > LC)
//                {

//                    throw std::invalid_argument("no. points exceed leaf capacity");

//                }

//                for(size_t m = 0; m < this->m_nvps; m++)
//                {

//                    size_t index = m_npoints;
//                    vp_t<T> vp = this->m_vps[m];

//                    for(datapoint_t<T,PL> &dp : points)
//                    {

//                        m_pdists[m][index++] = this->df->getDistance(vp.key, dp.key);

//                    }

//                }

//            }

//        public:
//            MVPLeaf()
//            {

//                m_npoints = 0;

//                for(size_t i = 0; i < LPN; i++)
//                {

//                    for(size_t j = 0; j < LC; j++)
//                    {

//                        m_pdists[i][j] = -1.0;

//                    }

//                }

//            }

//            ~MVPLeaf()
//            {



//            }

//            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
//                                                               std::map<size_t, std::vector<datapoint_t<T,PL>>*> &childpoints,
//                                                               const size_t level,
//                                                               const size_t index)
//            {

//                this->SelectVantagePoints(points);
//                MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* retnode = this;

//                if(m_npoints + points.size() <= LC)
//                {

//                    MarkLeafDistances(points);
//                    this->MarkPointDistances(points, level);

//                    for(datapoint_t<T,PL> &dp : points)
//                    {

//                        m_points[m_npoints++] = dp;

//                    }

//                    points.clear();

//                }
//                else
//                {

//                    std::vector<datapoint_t<T,PL>> pts = PurgeDataPoints();

//                    for(datapoint_t<T,PL> &dp : pts)
//                    {

//                        points.push_back(dp);

//                    }

//                    if(points.size() <= LPN + LC)
//                    {

//                        this->m_npoints = 0;
//                        this->m_nvps = 0;
//                        this->SelectVantagePoints(points);
//                        MarkLeafDistances(points);
//                        this->MarkPointDistances(points, level);

//                        for(datapoint_t<T,PL> &dp : points)
//                        {

//                            m_points[m_npoints++] = dp;

//                        }

//                        points.clear();

//                    }
//                    else
//                    {

//                        retnode = new MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();
//                        retnode = retnode->AddDataPoints(points, childpoints, level, index);

//                    }

//                }

//                if(retnode == NULL)
//                {

//                    throw std::runtime_error("unable to create node");

//                }

//            }


//            const size_t nbytes()const
//            {

//                return sizeof(MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>);

//            }

//            const size_t Size()const
//            {

//                size_t total = 0;

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    if(this->m_vps[i].active)
//                    {

//                        total++;

//                    }

//                }

//                for(size_t i = 0; i < this->m_npoints; i++)
//                {

//                    if(m_points[i].active)
//                    {

//                        total++;

//                    }

//                }

//                return total;

//            }

//            const bool IsLeaf()const
//            {

//                return true;

//            }

//            void SetChildNode(const size_t n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node)
//            {



//            }


//            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* GetChildNode(const size_t n)const
//            {

//                return NULL;

//            }

//            const std::vector<vp_t<T>> GetVantagePoints()const
//            {

//                std::vector<vp_t<T>> results;

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    results.push_back(this->m_vps[i]);

//                }

//                return results;

//            }

//            const std::vector<datapoint_t<T,PL>> GetDataPoints()const
//            {

//                std::vector<datapoint_t<T,PL>> results;

//                for(size_t i = 0; i < m_npoints; i++)
//                {

//                    results.push_back(m_points[i]);

//                }

//                return results;

//            }

//            void FilterDataPoints(const BasicArrayObject<T> &target,
//                                  std::vector<double> &tpath,
//                                  const double radius,
//                                  const size_t level,
//                                  const bool delete_points,
//                                  std::vector<item_t<T>> &results)
//            {

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    double d = this->df->getDistance(this->m_vps[i].key, target);
//                    tpath.push_back(d);

//                    if(this->m_vps[i].active && d <= radius)
//                    {

//                        vp_t<T> vp = this->m_vps[i];

//                        if(delete_points)
//                        {

//                            vp.active = false;

//                        }

//                        results.push_back({vp.id, vp.key});

//                    }

//                }

//                size_t pathlimit = (tpath.size() <= PL) ? tpath.size() : PL;

//                for(size_t index = 0; index < m_npoints; index++)
//                {

//                    bool skip = false;
//                    datapoint_t<T,PL> dp = m_points[index];

//                    if(dp.active)
//                    {

//                        for(size_t i = 0; i < this->m_nvps; i++)
//                        {

//                            if(!(m_pdists[i][index] >= tpath[level+i] - radius) && (m_pdists[i][index] <= tpath[level+i] + radius))
//                            {

//                                skip = true;
//                                break;

//                            }

//                        }

//                        if(!skip)
//                        {

//                            for(size_t i = 0; i < pathlimit; i++)
//                            {

//                                if(!(dp.dists[i] >= tpath[i] - radius && dp.dists[i] <= tpath[i] + radius))
//                                {

//                                    skip = true;
//                                    break;

//                                }

//                            }

//                        }

//                        if(!skip)
//                        {

//                            if(this->df->getDistance(dp.key, target) <= radius)
//                            {

//                                if(delete_points)
//                                {

//                                    m_points[index].active = false;

//                                }

//                                results.push_back({dp.id, dp.key});

//                            }

//                        }

//                    }

//                }

//            }

//            void TraverseNode(const BasicArrayObject<T> &target,
//                              const double radius,
//                              std::map<size_t, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
//                              std::vector<double> *tpath,
//                              std::map<size_t, std::vector<double>*> &child_tdistances,
//                              const size_t index,
//                              const size_t level,
//                              const bool delete_points,
//                              std::vector<item_t<T>> &results)
//            {

//                std::vector<double>* path = (tpath) ? tpath : new std::vector<double>();
//                FilterDataPoints(target, *path, radius, level, delete_points, results);
//                delete path;

//            }


//            const std::vector<datapoint_t<T,PL>> PurgeDataPoints()
//            {

//                std::vector<datapoint_t<T,PL>> results;

//                for(size_t i = 0; i < this->m_nvps; i++)
//                {

//                    vp_t<T> vp = this->m_vps[i];

//                    if(vp.active)
//                    {

//                        results.push_back({vp.id, vp.key});

//                    }

//                }

//                for(size_t i = 0; i < m_npoints; i++)
//                {

//                    if(m_points[i].active)
//                    {

//                        results.push_back(m_points[i]);

//                    }

//                }

//                return results;

//            }

//    };

//}

//template <typename T, class F, class PVT, size_t BF, size_t PL, size_t LC, size_t LPN, size_t FO, size_t NS>
//mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CreateNode(std::vector<datapoint_t<T,PL>> &points,
//                                                                                               std::map<size_t, std::vector<datapoint_t<T,PL>>*> &childpoints,
//                                                                                               size_t level,
//                                                                                               size_t index)
//{

//    MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node = NULL;

//    if(points.size() <= (LC + LPN))
//    {

//        node = new MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();

//    }
//    else
//    {

//        node = new MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();

//    }

//    node = node->AddDataPoints(points, childpoints, level, index);

//    if(node == NULL)
//    {

//        throw std::runtime_error("unable to create node");

//    }

//    return node;

//}

//*******************************************************************************************************************************

#include <map>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "datapoint.h"
#include <iostream>
#include <BasicArrayObject.h>
#include <Dataset.h>

namespace mvp {

    template <typename Type, typename Container>
    void erase_indices(
            const Container& indices_to_erase,
            std::vector<Type>& vec) {
        typedef typename Container::value_type IndexType;
        static_assert(std::is_same<IndexType, std::size_t>::value,
            "Indices to be erased have to be of type std::size_t");
        std::vector<bool> erase_index(vec.size(), false);
        for (const IndexType idx_erase: indices_to_erase) {
            erase_index[idx_erase] = true;
        }
        std::vector<bool>::const_iterator it_to_erase = erase_index.cbegin();
        typename std::vector<Type>::iterator it_erase_from = std::remove_if(
            vec.begin(), vec.end(),
            [&it_to_erase](const Type&) -> bool {
              return *it_to_erase++ == true;
            }
        );
        vec.erase(it_erase_from, vec.end());
    }

    template<typename T, class F, class PVT, class DT, int BF=2, int PL=8, int LC=30, int LPN=2, int FO=4, int NS=2>
    class MVPNode {
    protected:

        F* df = new F();
        //PVT* pvt = new PVT();
        int m_nvps;
        std::array<vp_t<T>,LPN>  m_vps;

        void SelectVantagePoints(std::vector<datapoint_t<T,PL>> &points);

        void MarkPointDistances(std::vector<datapoint_t<T,PL>> &points, const int level);

    public:
        MVPNode():m_nvps(0){}

        virtual ~MVPNode(){};

        void setDistanceFunction(F* df_)
        {

            df = df_;

        }

//        void setPivotMethod(PVT* pvt_)
//        {

//            pvt = pvt_;

//        }

        F* getDistanceFunction()
        {

            return df;

        }

        static MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* CreateNode(std::vector<datapoint_t<T,PL>> &points,
                                                         std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
                                                         int level,
                                                         int index);

        virtual MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
                                                             std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
                                                             const int level,
                                                             const int index) = 0;

        virtual const size_t nbytes()const = 0;

        virtual const int Size()const = 0;

        virtual const bool IsLeaf()const = 0;

        virtual void SetChildNode(const int n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node) = 0;

        virtual MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* GetChildNode(int n)const = 0;

        virtual const std::vector<vp_t<T>> GetVantagePoints()const = 0;

        virtual const std::vector<datapoint_t<T,PL>> GetDataPoints()const = 0;

        virtual void FilterDataPoints(const T &target,
                                      std::vector<double> &tpath,
                                      const double radius,
                                      const int level,
                                      const bool delete_points,
                                      std::vector<item_t<T>> &results) = 0;

        virtual void TraverseNode(const T &target,
                                  const double radius,
                                  std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                                  std::vector<double> *tpath,
                                  std::map<int, std::vector<double>*> &child_tdistances,
                                  const int index,
                                  const int level,
                                  const bool delete_points,
                                  std::vector<item_t<T>> &results) = 0;

        virtual const std::vector<datapoint_t<T,PL>> PurgeDataPoints()=0;

    };

    template<typename T, class F, class PVT, class DT, int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
    class MVPInternal : public MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> {
    private:

        std::array<MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*,FO> m_childnodes;

        // LPN x NS;
        std::array<std::array<double,NS>,LPN> m_splits;

        void CalcSplitPoints(const std::vector<double> &dists, int n, int split_index);

        std::vector<double> CalcPointDistances(vp_t<T> &vp, std::vector<datapoint_t<T,PL>> &points);

        std::vector<datapoint_t<T,PL>>* CullPoints(std::vector<datapoint_t<T,PL>> &list, std::vector<double> &dists,
                                                   double split, bool less);

        void CollatePoints(std::vector<datapoint_t<T,PL>> &points,
                           std::map<int, std::vector<datapoint_t<T,PL>>*> &childpoints,
                           const int level, const int index);

    public:
        MVPInternal();
        ~MVPInternal(){};

        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
                                                     std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
                                                     const int level,
                                                     const int index);

        const size_t nbytes()const;

        const int Size()const;

        const bool IsLeaf() const;

        void SetChildNode(const int n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node);

        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

        const std::vector<vp_t<T>> GetVantagePoints()const;

        const std::vector<datapoint_t<T,PL>> GetDataPoints()const;

        void FilterDataPoints(const T &target,
                              std::vector<double> &tpath,
                              const double radius,
                              const int level,
                              const bool delete_points,
                              std::vector<item_t<T>> &results);

        void TraverseNode(const T &target,const double radius,
                          std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                          std::vector<double> *tpath,
                          std::map<int, std::vector<double>*> &child_tdistances,
                          const int index,
                          const int level,
                          const bool delete_points,
                          std::vector<item_t<T>> &results);

        const std::vector<datapoint_t<T,PL>> PurgeDataPoints();

        //JOAO
        std::array<std::array<double,NS>,LPN> getSplits(){ return m_splits; }

    };

    template<typename T, class F, class PVT, class DT, int BF=2,int PL=8,int LC=30,int LPN=2, int FO=4, int NS=2>
    class MVPLeaf : public MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> {
    private:

        // PN x LC;
        std::array<std::array<double,LC>,LPN> m_pdists;
        int m_npoints;
        std::array<datapoint_t<T,PL>,LC>  m_points;

        void MarkLeafDistances(std::vector<datapoint_t<T,PL>> &points);

    public:
        MVPLeaf();
        ~MVPLeaf(){};
        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* AddDataPoints(std::vector<datapoint_t<T,PL>> &points,
                                                     std::map<int,std::vector<datapoint_t<T,PL>>*> &childpoints,
                                                     const int level, const int index);

        const size_t nbytes()const;

        const int Size()const;

        const bool IsLeaf()const;

        void SetChildNode(const int n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node);

        MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* GetChildNode(const int n)const;

        const std::vector<vp_t<T>> GetVantagePoints()const;

        const std::vector<datapoint_t<T,PL>> GetDataPoints()const;

        void FilterDataPoints(const T &target,
                              std::vector<double> &tpath,
                              const double radius,
                              const int level,
                              const bool delete_points,
                              std::vector<item_t<T>> &results);

        void TraverseNode(const T &target,const double radius,
                          std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                          std::vector<double> *tpath,
                          std::map<int, std::vector<double>*> &child_tdistances,
                          const int index,
                          const int level,
                          const bool delete_points,
                          std::vector<item_t<T>> &results);

        const std::vector<datapoint_t<T,PL>> PurgeDataPoints();
    };


    /**
     *  MVPNode implementation
     *
     **/

    bool CompareDistance(const double a, const double b, const bool less){
        return (less) ? (a <= b) : (a > b);
    }

}

/********** MVPNode methods *********************/

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CreateNode(std::vector<mvp::datapoint_t<T,PL>> &points,
                                                                                   std::map<int, std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
                                                                                   int level,
                                                                                   int index){
    MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node = NULL;
    if (points.size() <= LC + LPN){
        node = new MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();
    } else {
        node = new MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();
    }

    node = node->AddDataPoints(points, childpoints, level, index);
    if (node == NULL) throw std::runtime_error("unable to create node");
    return node;
}


template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::SelectVantagePoints(std::vector<mvp::datapoint_t<T,PL>> &points){

    if(this->m_nvps < LPN && points.size() > 0)
    {

        int nPivots = std::min((size_t)LPN, points.size());

        std::vector<T> vec;
        size_t id = 0;

        for(datapoint_t<T,PL> &dp : points)
        {

            T b = *dp.key.clone();
            b.setOID(id++);
            vec.push_back(b);

            //std::cout << "CANDIDATES = " << dp.id << " / " << dp.key.toStringWithOID() << "\n";

        }


        DT* dataset = new DT(vec, vec.size(), vec[0].size());
        PVT* pvt = new PVT();
        pvt->generatePivots(dataset, this->df, nPivots);

        std::vector<size_t> remove_points_index;

        this->m_nvps = nPivots;
        for(int i = 0; i < nPivots; i++)
        {

            int pos = pvt->getPivot(i)->getOID();
            remove_points_index.push_back(pos);
            this->m_vps[i] = {points[pos].id, points[pos].key};

            //std::cout << "PIVO = " << this->m_vps[i].id << " / " << this->m_vps[i].key.toStringWithOID() << "\n";

//            mvp::datapoint_t<T,PL> temp = points[points.size()-1];
//            points[pos] = points[points.size()-1];
//            points[points.size()-1] = temp;
//            points.pop_back();
//            points.erase(points.begin() + pos);

        }

        erase_indices(remove_points_index, points);

//        for(int pos : remove_points_index)
//        {
////            mvp::datapoint_t<T,PL> temp = points[points.size()-1];
////            points[pos] = points[points.size()-1];
////            points[points.size()-1] = temp;
////            points.pop_back();

//            points.erase(points.begin() + pos);
//        }

        delete dataset;
        delete pvt;



//        DT dataset = DT(vec, vec.size(), vec[0].size());
//        pvt->generatePivots(&dataset, this->df, LPN);
//        std::vector<vp_t<T>> vp_aux;

//        this->m_nvps = LPN;
//        for(int i = 0; i < LPN; i++)
//        {

//            size_t pos = pvt->getPivot(i)->getOID();
//            vp_aux.push_back({points[pos].id, points[pos].key});
//            std::cout << "PIVO = " << vp_aux[i].id << " / " << vp_aux[i].key.toStringWithOID() << " \/ POS = " << pos << "\n";
//            this->m_vps[i] = { points[pos].id, points[pos].key };
//            mvp::datapoint_t<T,PL> temp = points[points.size()-1];
//            points[pos] = points[points.size()-1];
//            points[points.size()-1] = temp;
//            points.pop_back();

//        }

    }

    //PRECISA DAS POSIÇÕES DOS PIVOS

//    if (this->m_nvps == 0 && points.size() > 0){
//        //std::cout << "VANTAGE = " << points.back().id << "\n";
//        this->m_vps[this->m_nvps++] = { points.back().id, points.back().key };
//        points.pop_back();
//    }

//    int slimit = 10;
//    while (this->m_nvps < LPN && points.size() > 0){
//        int max_pos = 0;
//        double maxd = 0;
//        slimit = ((int)points.size() <= slimit) ? points.size() : slimit;
//        for (int i=0;i < slimit;i++){
//            double d = this->df->getDistance(points[i].key, this->m_vps[this->m_nvps-1].key);
////            double d = points[i].distance(this->m_vps[this->m_nvps-1]);
//            if (d > maxd){
//                max_pos = i;
//                maxd = d;
//            }
//        }
//        this->m_vps[this->m_nvps++] = { points[max_pos].id, points[max_pos].key };
//        //std::cout << "VANTAGE = " << points[max_pos].id << "\n";
//        //JOAO
//        mvp::datapoint_t<T,PL> temp = points[points.size()-1];
//        points[max_pos] = points[points.size()-1];
//        points[points.size()-1] = temp;
//        //
//        points.pop_back();
//    }

//    std::cout << "\n\n\n\n";

}

template<typename T,class F,class PVT,class DT,int BF, int PL, int LC, int LPN, int FO, int NS>
void mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::MarkPointDistances(std::vector<mvp::datapoint_t<T,PL>> &points, const int level){
    for (int i=0;i < this->m_nvps;i++){
        if (level + i < PL){
            for (mvp::datapoint_t<T,PL> &pnt : points){
                pnt.dists[level+i] = this->df->getDistance(pnt.key, m_vps[i].key);
//                pnt.dists[level+i] = pnt.distance(m_vps[i]);
            }
        }
    }
}

/********** MVPInternal methods *******************/
template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::MVPInternal(){
    for (int i=0;i < FO;i++) m_childnodes[i] = NULL;
    for (int i=0;i < LPN;i++){
        for (int j=0;j < NS;j++){
            m_splits[i][j] = -1.0;
        }
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CalcSplitPoints(const std::vector<double> &dists, int n, int split_index){
    int lengthM = BF - 1;
    if (dists.size() > 0){
        if (m_splits[n][split_index*lengthM] == -1){
            std::vector<double> tmpdists = dists;
            sort(tmpdists.begin(), tmpdists.end());
            double factor = (double)tmpdists.size()/(double)BF;
            for (int i=0;i<lengthM;i++){
                double pos = (i+1)*factor;
                int lo = floor(pos);
                int hi = (pos <= tmpdists.size()-1) ? ceil(pos) : 0;
                m_splits[n][split_index*lengthM+i] = (tmpdists[lo] + tmpdists[hi])/2.0;
            }
        }
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
std::vector<double> mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CalcPointDistances(vp_t<T> &vp, std::vector<mvp::datapoint_t<T,PL>> &points){
    std::vector<double> results;
    for (mvp::datapoint_t<T,PL> &dp : points){
        //results.push_back(dp.distance(vp));
        results.push_back(this->df->getDistance(dp.key, vp.key));
    }
    return results;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
std::vector<mvp::datapoint_t<T,PL>>* mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CullPoints(std::vector<mvp::datapoint_t<T,PL>> &list,
                                                                         std::vector<double> &dists,
                                                                         double split,
                                                                         bool less){
    std::vector<mvp::datapoint_t<T,PL>> *results = new std::vector<mvp::datapoint_t<T,PL>>();

    auto list_iter = list.begin();
    auto dist_iter = dists.begin();
    while (list_iter != list.end() && dist_iter != dists.end()){
        if (CompareDistance(*dist_iter,split,less)){
            results->push_back(*list_iter);
            list_iter = list.erase(list_iter);
            dist_iter = dists.erase(dist_iter);
        } else {
            list_iter++;
            dist_iter++;
        }
    }
    if (results->size() > 0){
        return results;
    }
    delete results;
    return NULL;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::CollatePoints(std::vector<mvp::datapoint_t<T,PL>> &points,
                                std::map<int, std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
                                const int level, const int index){
    std::map<int, std::vector<mvp::datapoint_t<T,PL>>*> pnts, pnts2;
    pnts[0] = &points;

    int lengthM = BF - 1;
    int n = 0;
    do {
        for (auto iter=pnts.begin();iter!=pnts.end();iter++){
            int node_index = iter->first;
            std::vector<mvp::datapoint_t<T,PL>> *list = iter->second;

            //int list_size = list->size();

            std::vector<double> dists = CalcPointDistances((this->m_vps[n]), *list);

//            for(double d : dists)
//                std::cout << d << " ";
//            std::cout << "\n";

            if (dists.size() > 0){
                CalcSplitPoints(dists, n, node_index);

                double m;
                std::vector<mvp::datapoint_t<T,PL>> *culledpts = NULL;
                for (int j=0;j<lengthM;j++){
                    m = m_splits[n][node_index*lengthM+j];

                    culledpts = CullPoints(*list, dists, m, true);
                    if (culledpts != NULL){
                        pnts2[node_index*BF+j] = culledpts;
                    }
                }
                m = m_splits[n][node_index*lengthM+lengthM-1];
                culledpts = CullPoints(*list, dists, m, false);
                if (culledpts != NULL){
                    pnts2[node_index*BF+BF-1] = culledpts;
                }
            }
            if (list->size() > 0)
                throw std::length_error("not fully collated");

            if (n > 0) delete list;
        }
        pnts = move(pnts2);
        n++;
    } while (n < LPN);

    for (auto iter=pnts.begin();iter!=pnts.end();iter++){
        int i=iter->first;
        std::vector<mvp::datapoint_t<T,PL>> *list = iter->second;
        if (list != NULL)
            childpoints[index*FO+i] = list;
    }

}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::AddDataPoints(std::vector<mvp::datapoint_t<T,PL>> &points,
                                                                                std::map<int,std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
                                                                                const int level,
                                                                                const int index){
    this->SelectVantagePoints(points);
    this->MarkPointDistances(points, level);
    if (this->m_nvps < LPN) throw std::invalid_argument("too few points for internal node");
    CollatePoints(points, childpoints, level, index);
    points.clear();
    return this;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const size_t mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::nbytes()const{
    return sizeof(MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>);
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const int mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Size()const{
    int total = 0;
    for (int i=0;i < this->m_nvps;i++){
        if (this->m_vps[i].active) total++;
    }
    return total;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const bool mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::IsLeaf()const{
    return false;
}


template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::SetChildNode(const int n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *node){
    if (n < 0 || n >= FO) throw std::invalid_argument("index out of range");
    m_childnodes[n] = node;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::GetChildNode(const int n)const{
    if (n < 0 || n >= FO) throw std::invalid_argument("index out of range");
    return m_childnodes[n];
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::vp_t<T>> mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
    std::vector<mvp::vp_t<T>> results;
    for (int i=0;i < this->m_nvps;i++) results.push_back(this->m_vps[i]);
    return results;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
    std::vector<mvp::datapoint_t<T,PL>> results;
    return results;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const T &target,
                                                         std::vector<double> &tpath,
                                                         const double radius,
                                                         const int level,
                                                         const bool delete_points,
                                                         std::vector<item_t<T>> &results){
    for (int i=0;i < this->m_nvps;i++){
//        double d = this->m_vps[i].distance(target);
        double d = this->df->getDistance(this->m_vps[i].key, target);
        tpath.push_back(d);
        if (this->m_vps[i].active && d <= radius){
            mvp::vp_t<T> vp = this->m_vps[i];
            if (delete_points) this->m_vps[i].active = false;
            results.push_back({vp.id, vp.key});
        }
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::TraverseNode(const T &target,
                                                     const double radius,
                                                     std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                                                     std::vector<double> *tpath,
                                                     std::map<int, std::vector<double>*> &child_tdistances,
                                                     const int index,
                                                     const int level,
                                                     const bool delete_points,
                                                     std::vector<item_t<T>> &results){
    int lengthM = BF - 1;
    int n = 0;
    bool *currnodes  = new bool[1];
    currnodes[0] = true;


    std::vector<double> *path = (tpath != NULL) ? tpath : new std::vector<double>();

    FilterDataPoints(target, *path, radius, level, delete_points, results);

    do {
        int n_nodes = pow(BF, n);
        int n_childnodes = pow(BF, n+1);
        bool *nextnodes = new bool[n_childnodes];

        for (int i=0;i<n_childnodes;i++)
            nextnodes[i] = false;

        double d = path->at(level+n);

        //int lengthMn = lengthM*n_nodes;
        for (int node_index=0;node_index < n_nodes;node_index++){
            if (currnodes[node_index]){
                if (m_splits[n][node_index*lengthM] >= 0){
                    double m = m_splits[n][node_index*lengthM];
                    for (int j=0;j < lengthM;j++){
                        m = m_splits[n][node_index*lengthM+j];
                        if (d <= m + radius) nextnodes[node_index*BF+j] = true;
                    }
                    if (d > m - radius) nextnodes[node_index*BF+BF-1] = true;
                }
            }
        }

        delete[] currnodes;
        currnodes = nextnodes;
        n++;
    } while (n < LPN);

    for (int i=0;i < FO;i++){
        if (currnodes[i]){
            MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *child = GetChildNode(i);
            if (child != NULL) {
                childnodes[index*FO+i] = child;

                std::vector<double> *tmppath = new std::vector<double>();
                for (int j=0;j < (int)path->size();j++){
                    tmppath->push_back(path->at(j));
                }
                child_tdistances[index*FO+i] = tmppath;
            }
        }
    }

    delete[] currnodes;
    delete path;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
    std::vector<mvp::datapoint_t<T,PL>> results;
    for (int i=0;i < this->m_nvps;i++){
        if (this->m_vps[i].active){
            vp_t<T> vp = this->m_vps[i];
            results.push_back({ vp.id, vp.key} );
        }
    }

    return results;
}

/********* MVPLeaf methods **********************/

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::MVPLeaf(){
    m_npoints = 0;
    for (int i=0;i < LPN;i++)
        for (int j=0;j<LC;j++)
            m_pdists[i][j] = -1.0;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::MarkLeafDistances(std::vector<mvp::datapoint_t<T,PL>> &points){
    if (m_npoints + points.size() > LC)
        throw std::invalid_argument("no. points exceed leaf capacity");

    for (int m = 0; m < this->m_nvps;m++){
        int index = m_npoints;
        vp_t<T> vp = this->m_vps[m];
        for (mvp::datapoint_t<T,PL> &dp : points){
//            m_pdists[m][index++] = dp.distance(vp);
            m_pdists[m][index++] = this->df->getDistance(dp.key, vp.key);
        }
    }
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::AddDataPoints(std::vector<mvp::datapoint_t<T,PL>> &points,
                                std::map<int,std::vector<mvp::datapoint_t<T,PL>>*> &childpoints,
                                const int level, const int index){
    this->SelectVantagePoints(points);
    MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS> *retnode = this;
    if (m_npoints + points.size() <= LC){
        // add points to existing leaf
        MarkLeafDistances(points);
        this->MarkPointDistances(points, level);
        for (mvp::datapoint_t<T,PL> &dp : points){
            m_points[m_npoints++] = dp;
        }
        points.clear();
    } else {  // create new internal node

        // get existing points, purge inactive poins
        std::vector<mvp::datapoint_t<T,PL>> pts = PurgeDataPoints();

        // merge points
        for (mvp::datapoint_t<T,PL> &dp : pts)
            points.push_back(dp);

        if (points.size() <= LPN + LC){
            // clear out points
            this->m_npoints = 0;
            this->m_nvps = 0;
            this->SelectVantagePoints(points);
            MarkLeafDistances(points);
            this->MarkPointDistances(points, level);
            for (mvp::datapoint_t<T,PL> &dp : points){
                m_points[m_npoints++] = dp;
            }
            points.clear();
        } else {
            retnode = new mvp::MVPInternal<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>();
            retnode = retnode->AddDataPoints(points, childpoints, level, index);
        }
    }
    if (retnode == NULL) throw std::runtime_error("unable to create node");

    return retnode;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const size_t mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::nbytes()const{
    return sizeof(mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>);
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const int mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::Size()const{
    int total = 0;
    for (int i=0;i < this->m_nvps;i++){
        if (this->m_vps[i].active) total++;
    }
    for (int i=0;i < this->m_npoints;i++){
        if (m_points[i].active) total++;
    }

    return total;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const bool mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::IsLeaf()const{
    return true;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::SetChildNode(const int n, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* node){}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
mvp::MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>* mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::GetChildNode(const int n)const{return NULL;}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::vp_t<T>> mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::GetVantagePoints()const{
    std::vector<vp_t<T>> results;
    for (int i=0;i < this->m_nvps;i++)
        results.push_back(this->m_vps[i]);

    return results;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::GetDataPoints()const{
    std::vector<mvp::datapoint_t<T,PL>> results;
    for (int i=0;i<m_npoints;i++)
        results.push_back(m_points[i]);
    return results;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::FilterDataPoints(const T &target,
                                                     std::vector<double> &tpath,
                                                     const double radius,
                                                     const int level,
                                                     const bool delete_points,
                                                     std::vector<item_t<T>> &results){
    for (int i=0;i < this->m_nvps;i++){
//        double d = this->m_vps[i].distance(target);
        double d = this->df->getDistance(this->m_vps[i].key, target);
        tpath.push_back(d);
        if (this->m_vps[i].active && d <= radius){
            vp_t<T> vp = this->m_vps[i];
            if (delete_points){
                vp.active = false;
            }
            results.push_back({vp.id, vp.key });
        }
    }

    int pathlimit = (tpath.size() <= PL) ? tpath.size() : PL;
    for (int index=0;index < (int)m_npoints;index++){
        bool skip = false;
        mvp::datapoint_t<T,PL> dp = m_points[index];

        if (dp.active){
            // filter using precomputed distances in node
            for (int i=0;i < this->m_nvps;i++){
                if (!(m_pdists[i][index] >=  tpath[level+i] - radius) && (m_pdists[i][index] <= tpath[level+i] + radius)){
                    skip = true;
                    break;
                }
            }

            if (!skip){
                // filter using precomputed path distances between target and vantage points
                // from top down to current place in tree
                for (int i=0;i < pathlimit;i++){
                    if (!(dp.dists[i] >= tpath[i] - radius && dp.dists[i] <= tpath[i] + radius)){
                        skip = true;
                        break;
                    }
                }
            }

            if (!skip){
                // still not ruled out
                if (/*dp.distance(target)*/this->df->getDistance(dp.key, target) <= radius){
                    if (delete_points) {
                        m_points[index].active = false;
                    }
                    results.push_back({ dp.id, dp.key });
                }
            }
        }
    }
}


template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
void mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::TraverseNode(const T &target,
                                                 const double radius,
                                                 std::map<int, MVPNode<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>*> &childnodes,
                                                 std::vector<double> *tpath,
                                                 std::map<int, std::vector<double>*> &child_tdistances,
                                                 const int index,
                                                 const int level,
                                                 const bool delete_points,
                                                 std::vector<item_t<T>> &results){

    std::vector<double> *path = (tpath) ? tpath : new std::vector<double>();
    FilterDataPoints(target, *path, radius, level, delete_points, results);
    delete path;
}

template<typename T,class F,class PVT,class DT,int BF,int PL,int LC,int LPN,int FO,int NS>
const std::vector<mvp::datapoint_t<T,PL>> mvp::MVPLeaf<T,F,PVT,DT,BF,PL,LC,LPN,FO,NS>::PurgeDataPoints(){
    std::vector<mvp::datapoint_t<T,PL>> results;
    for (int i=0;i < this->m_nvps;i++){
        vp_t<T> vp = this->m_vps[i];
        if (vp.active){
            results.push_back({ vp.id, vp.key });
        }
    }
    for (int i=0;i < m_npoints;i++){
        if (m_points[i].active)
            results.push_back(m_points[i]);
    }
    return results;
}


#endif // MVPNODE_H
