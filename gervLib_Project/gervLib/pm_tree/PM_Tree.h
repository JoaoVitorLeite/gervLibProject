#ifndef PM_TREE_H
#define PM_TREE_H

#include <Dataset.h>
#include <Hermes.h>
#include <Pivot.h>
#include <vector>
#include <set>
#include <Pivots.h>
#include <algorithm>
#include <KdTree.h>

//enum PromoteFunc_E { RANDOM_e,m_RAD_2 };

template <class DType>
struct PM_Node
{

    PM_Node* parent_node;
    size_t node_category; //0: represent routing node, 1: leaf node, 2: data entry
    BasicArrayObject<DType> feature_val;
    double dist_to_parent;
    size_t level;

    size_t id;
    std::vector<double> pivot_distance;

    double range;
    std::vector<PM_Node<DType>*> ptr_sub_tree;
    std::vector<std::pair<double, double>> hyper_rings;

    PM_Node()
    {



    }

    PM_Node(PM_Node* parent_node_, size_t node_category_, double dist_to_parent_, size_t id_)
    {

        parent_node = parent_node_;
        node_category = node_category_;
        dist_to_parent = dist_to_parent_;
        id = id_;


    }

};

template <class DType>
struct PM_Partition
{

    PM_Node<DType>* node;
    double min, max;

    PM_Partition(PM_Node<DType>* node_, double min_, double max_)
    {

        node = node_;
        min = min_;
        max = max_;

    }

    PM_Partition()
    {



    }

};


template <class DType>
struct Compare_PM_Nodes
{

    bool operator()(PM_Partition<DType> const& a, PM_Partition<DType> const& b)
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


template <class DType>
struct KnnEntry
{

    BasicArrayObject<DType> element;
    double distance;

    KnnEntry(BasicArrayObject<DType>& element_, double distance_)
    {

        element = element_;
        distance = distance_;

    }

    KnnEntry()
    {



    }

    bool operator<(const KnnEntry<DType>& item) const
    {

        return distance < item.distance;

    }

    bool operator>(const KnnEntry<DType>& item) const
    {

        return distance > item.distance;

    }

};

template<class Type, class Comp>
std::vector<Type> dequeueInOrderPM_Results(std::priority_queue<Type, std::vector<Type>, Comp> pq)
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


template <class DType>
class PM_Tree
{

    private:
        size_t M, HR, PD, Pivot_Num, Max_Level;
        Dataset<DType> pivot_vec;
        PM_Node<DType>* root;
        size_t Pivot_Filter_Num = 0, Cal_Distance_Num = 0, Via_subNode_Num, Via_Node_Num;
        size_t leafNodeAccess = 0;
        DistanceFunction<BasicArrayObject<DType>>* df;
        Pivot<DType>* pvt;

    private:
        void update_pivot_distance(std::vector<double>& pivot_distance_, BasicArrayObject<DType>& data_feature_val_);
        double cal_cover_radius(const PM_Node<DType>* node_);
        double cal_cover_radius(const PM_Node<DType>* node_, const PM_Node<DType>* cmp_node_);
        double cal_dist_to_parent(const PM_Node<DType>* node_);
        bool is_leaf_node(const PM_Node<DType>* node_);
        void update_hyper_rings(std::vector<std::pair<double,double>>& hyper_rings_, BasicArrayObject<DType>& data_feature_val_);
        void insert(PM_Node<DType> ** cur_node_address_, PM_Node<DType> ** insert_node_address_);
        bool is_full(PM_Node<DType>* node_);
        PM_Node<DType>** get_next_node_and_update_range(PM_Node<DType> ** cur_node_address_, PM_Node<DType> ** insert_node_address_);
        void split(PM_Node<DType> ** cur_node_ptr_address_, PM_Node<DType> ** insert_node_address_);
        bool is_root_node(const PM_Node<DType>* node_);
        void merge_subNode_HR(PM_Node<DType>* cur_node_);
        void assign_node_all_value(PM_Node<DType>* cur_node_, PM_Node<DType>* new_node_);
        void partition(std::vector<PM_Node<DType>*>& entries_, PM_Node<DType>* node_1_, PM_Node<DType>* node_2_);
        void promote(std::vector<PM_Node<DType>*>& entries_, PM_Node<DType>* node_1_, PM_Node<DType>* node_2_);
        void subRange_search(PM_Node<DType>** cur_node_address_, BasicArrayObject<DType>& q_feature_val_, double search_range_, std::vector<std::pair<double, size_t>>& res_vec_, double& dist_parent_q_, std::vector<double>& dist_q_pivot_);
        bool is_pivot_filter(std::vector<double>& dist_q_pivot_, double range_, PM_Node<DType>* cur_node_);
        bool is_data_node(const PM_Node<DType>* node_);
        bool is_rounting_node(const PM_Node<DType>* node_);
        void sub_traverse_get_volumem(PM_Node<DType>* cur_node_, double & volume_);
        double cal_hypersphere_volume(double radius_, size_t dim_);
        double minDistNode(PM_Node<DType>* cur_node_, BasicArrayObject<DType>& element, std::vector<double> dist_to_query);
        double maxDistNode(PM_Node<DType>* cur_node_, BasicArrayObject<DType>& element, std::vector<double> dist_to_query);
//        double minLimInf(PM_Node<DType>* cur_node_);
//        double maxLimSup(PM_Node<DType>* cur_node_);

    public:
        PM_Tree(Dataset<DType>* dataset, DistanceFunction<BasicArrayObject<DType>>* df_, Pivot<DType>* pvt, size_t m_, size_t pivot_num_);
        ~PM_Tree();

        void set_pivot(size_t pivot_num_, const Dataset<DType>& pivot_vec_);
        void insert(BasicArrayObject<DType>& feature_val, int id_);
        size_t get_pivot_filter_num();
        size_t get_cal_distance_num();
        size_t get_via_subNode_num();
        size_t get_via_node_num();
        void range_search(BasicArrayObject<DType>& q_feature_val_, double search_range_, std::vector<std::pair<double, size_t>>& res_vec_);
        void update_level(PM_Node<DType>* cur_node_, size_t level_);
        double traverse_get_volume();
        PM_Node<DType>* get_root();
        void output_node_info(std::ofstream& out_, const PM_Node<DType>* cur_node_);
        void traverse_bread_tree(std::string out_file_path);
        size_t getMaxLevel();
        //void test();
        void kNN(BasicArrayObject<DType> query, size_t k, std::vector<KnnEntry<DType>>& ansVec);
        size_t getLeafNodeAccess();
        std::vector<PM_Node<DType>*> leafsNodes();

};



template <class DType>
void PM_Tree<DType>::update_pivot_distance(std::vector<double>& pivot_distance_, BasicArrayObject<DType>& data_feature_val_)
{

    if(pivot_distance_.size() == 0)
    {

        pivot_distance_.resize(PD, 0);

    }

    for(size_t i = 0; i < PD; ++i)
    {

        double dist = df->getDistance(data_feature_val_, *pivot_vec.instance(i));
        pivot_distance_[i] = dist;

    }

}



template <class DType>
double PM_Tree<DType>::cal_cover_radius(const PM_Node<DType>* node_)
{

    if(node_->ptr_sub_tree.size() == 0)
    {

        return 0;

    }

//MUDANCAS JOAO
//    double max_range = std::numeric_limits<double>::min();
    double max_range = 0.0;


    if(is_leaf_node(node_))
    {

        for(size_t i = 0; i < node_->ptr_sub_tree.size(); ++i)
        {

            double dist = df->getDistance(node_->ptr_sub_tree[i]->feature_val, node_->feature_val);

            if(dist > max_range)
            {

                max_range = dist;

            }

        }

    }
    else
    {

        for(size_t i = 0; i < node_->ptr_sub_tree.size(); ++i)
        {

            double dist = df->getDistance(node_->ptr_sub_tree[i]->feature_val, node_->feature_val);
            dist += node_->ptr_sub_tree[i]->range;

            if(dist > max_range)
            {

                max_range = dist;

            }

        }

    }

    return max_range;

}


template <class DType>
double PM_Tree<DType>::cal_cover_radius(const PM_Node<DType>* node_, const PM_Node<DType>* cmp_node_)
{

    double dist;

    if(is_leaf_node(node_))
    {

        dist = df->getDistance(cmp_node_->feature_val, node_->feature_val);

    }
    else
    {

        dist = df->getDistance(cmp_node_->feature_val, node_->feature_val);
        dist += cmp_node_->range;

    }

    if(node_->range < dist)
    {

        return dist;

    }
    else
    {

        return node_->range;

    }

}


template <class DType>
double PM_Tree<DType>::cal_dist_to_parent(const PM_Node<DType>* node_)
{

    double dist;

    if(!(node_ == nullptr || node_->parent_node == nullptr))
    {

        dist = df->getDistance(node_->feature_val, node_->parent_node->feature_val);

    }
    else
    {

        throw std::invalid_argument("Error in cal_dist_to_parent \n");
        return -1;

    }

    return dist;

}


template <class DType>
bool PM_Tree<DType>::is_leaf_node(const PM_Node<DType>* node_)
{

    return node_->node_category == 1;

}


template <class DType>
void PM_Tree<DType>::update_hyper_rings(std::vector<std::pair<double,double>>& hyper_rings_, BasicArrayObject<DType>& data_feature_val_)
{

    if(hyper_rings_.size() == 0)
    {

        hyper_rings_.resize(HR, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::min()));

    }

    for(size_t i = 0; i < HR; ++i)
    {

        double dist = df->getDistance(data_feature_val_, *pivot_vec.instance(i));

        if(dist < hyper_rings_[i].first)
        {

            hyper_rings_[i].first = dist;

        }

        if(dist > hyper_rings_[i].second)
        {

            hyper_rings_[i].second = dist;

        }
    }

}



template <class DType>
void PM_Tree<DType>::insert(PM_Node<DType> ** cur_node_address_, PM_Node<DType> ** insert_node_address_)
{

    PM_Node<DType>* cur_node_ = *cur_node_address_;
    PM_Node<DType> * insert_node_ = *insert_node_address_;

    if(is_leaf_node(cur_node_))
    {

        if(is_full(cur_node_))
        {

            split(&cur_node_, &insert_node_);

        }
        else
        {

            cur_node_->range = cal_cover_radius(cur_node_, insert_node_);
            cur_node_->ptr_sub_tree.push_back(insert_node_);

            insert_node_->parent_node = cur_node_;
            insert_node_->dist_to_parent = cal_dist_to_parent(insert_node_);
            update_hyper_rings(cur_node_->hyper_rings, insert_node_->feature_val);

            //MUDANCAS JOAO
            PM_Node<DType>* update_ptr = cur_node_->parent_node;

            while(update_ptr != nullptr)
            {

                update_ptr->range = cal_cover_radius(update_ptr);
                merge_subNode_HR(update_ptr);
                update_ptr = update_ptr->parent_node;

            }
            //MUDANCAS JOAO

        }

    }
    else
    {

        PM_Node<DType>** next_node = get_next_node_and_update_range(&cur_node_, &insert_node_);
        update_hyper_rings((*next_node)->hyper_rings, insert_node_->feature_val);
        insert(next_node, &insert_node_);

    }

}


template <class DType>
bool PM_Tree<DType>::is_full(PM_Node<DType>* node_)
{

    return node_->ptr_sub_tree.size() >= M;

}

template <class DType>
PM_Node<DType>** PM_Tree<DType>::get_next_node_and_update_range(PM_Node<DType> ** cur_node_address_, PM_Node<DType> ** insert_node_address_)
{

    PM_Node<DType>* cur_node_ = *cur_node_address_;
    PM_Node<DType>* insert_node_ = *insert_node_address_;

    std::vector<std::pair<double, PM_Node<DType>**>> data_vec;
    std::vector<std::pair<double, PM_Node<DType>**>> data2_vec;

    for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
    {

        double dist = df->getDistance(cur_node_->ptr_sub_tree[i]->feature_val, insert_node_->feature_val);
        //std::cout << "range = " << cur_node_->ptr_sub_tree[i]->range << "\n";

        if(dist <= cur_node_->ptr_sub_tree[i]->range)
        {

            //std::cout << "DIS in range= "  << cur_node_->ptr_sub_tree[i]->range - dist << "\n";
            data_vec.push_back(std::make_pair(dist, &(cur_node_->ptr_sub_tree[i])));

        }

        if(data_vec.size() == 0)
        {

            dist -= cur_node_->ptr_sub_tree[i]->range;
            //std::cout << "DIS out range= "  << dist << "\n";
            data2_vec.push_back(std::make_pair(dist, &(cur_node_->ptr_sub_tree[i])));

        }

    }

    if(data_vec.size() != 0)
    {

        std::sort(data_vec.begin(), data_vec.end());
        return data_vec[0].second;

    }

    std::sort(data2_vec.begin(), data2_vec.end());
    (*(data2_vec[0].second))->range = data2_vec[0].first + (*(data2_vec[0].second))->range;
    //std::cout << (*data2_vec[0].second)->feature_val.getOID() << "\n";
    return data2_vec[0].second;

}



template <class DType>
void PM_Tree<DType>::split(PM_Node<DType> ** cur_node_ptr_address_, PM_Node<DType> ** insert_node_address_)
{

    //std::cout << "SPLIT\n";
    PM_Node<DType> * cur_node_ = *cur_node_ptr_address_;
    PM_Node<DType> * insert_node_ = *insert_node_address_;

    PM_Node<DType>* new_Mnode_1 = new PM_Node<DType>(nullptr, 0, -1, -1);
    PM_Node<DType>* new_Mnode_2 = new PM_Node<DType>(nullptr, 0, -1, -1);

    std::vector<PM_Node<DType>*> entries = cur_node_->ptr_sub_tree;
    entries.push_back(insert_node_);

    promote(entries, new_Mnode_1, new_Mnode_2);

    partition(entries, new_Mnode_1, new_Mnode_2);

    if(is_root_node(cur_node_))
    {

        PM_Node<DType>* new_root = new PM_Node<DType>(nullptr, 0, -1, -1);
        new_root->feature_val = new_Mnode_1->feature_val;

        new_Mnode_1->parent_node = new_root;
        new_Mnode_1->dist_to_parent = cal_dist_to_parent(new_Mnode_1);
        new_Mnode_2->parent_node = new_root;
        new_Mnode_2->dist_to_parent = cal_dist_to_parent(new_Mnode_2);
        merge_subNode_HR(new_Mnode_1);
        merge_subNode_HR(new_Mnode_2);

        new_root->ptr_sub_tree.push_back(new_Mnode_1);
        new_root->ptr_sub_tree.push_back(new_Mnode_2);
        new_root->range = cal_cover_radius(new_root);
        merge_subNode_HR(new_root);
        root = new_root;

    }
    else
    {

        new_Mnode_1->parent_node = cur_node_->parent_node;
        new_Mnode_1->dist_to_parent = cal_dist_to_parent(new_Mnode_1);
        new_Mnode_2->parent_node = cur_node_->parent_node;
        new_Mnode_2->dist_to_parent = cal_dist_to_parent(new_Mnode_2);
        merge_subNode_HR(new_Mnode_1); //JOAO
        assign_node_all_value(cur_node_, new_Mnode_1);
        //merge_subNode_HR(cur_node_); //JOAO: PROVAVEL COMANDO SEM UTILIDADE
        merge_subNode_HR(new_Mnode_2);

        delete (new_Mnode_1);

        if(is_full(new_Mnode_2->parent_node))
        {

            split(&(new_Mnode_2->parent_node), &new_Mnode_2);

        }
        else
        {

            new_Mnode_2->parent_node->ptr_sub_tree.push_back(new_Mnode_2);

        }

        //MUDANCAS JOAO
        PM_Node<DType>* update_ptr = cur_node_->parent_node;

        while(update_ptr != nullptr)
        {

            update_ptr->range = cal_cover_radius(update_ptr);
            merge_subNode_HR(update_ptr);
            update_ptr = update_ptr->parent_node;

        }
        //MUDANCAS JOAO


//        std::cout << "DPS 1 = " << root->range << std::endl;
//        root->range = cal_cover_radius(root);
//        std::cout << "DPS 2 = " << root->range << "\n\n";

    }

}



template <class DType>
bool PM_Tree<DType>::is_root_node(const PM_Node<DType>* node_)
{

    return node_->parent_node == nullptr;

}


template <class DType>
void PM_Tree<DType>::merge_subNode_HR(PM_Node<DType>* cur_node_)
{

    if(is_leaf_node(cur_node_))
    {

        for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
        {

            update_hyper_rings(cur_node_->hyper_rings, cur_node_->ptr_sub_tree[i]->feature_val);

        }

    }
    else
    {

        if(cur_node_->hyper_rings.size() == 0)
        {

            cur_node_->hyper_rings.resize(HR, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::min()));

        }

        for(size_t i = 0; i < HR; i++)
        {

            for(size_t j = 0; j < cur_node_->ptr_sub_tree.size(); ++j)
            {

                double sub_min;
                double sub_max;

                if(cur_node_->ptr_sub_tree[j]->ptr_sub_tree.size() != 0)
                {

                    sub_min = cur_node_->ptr_sub_tree[j]->hyper_rings[i].first;
                    sub_max = cur_node_->ptr_sub_tree[j]->hyper_rings[i].second;

                }
                else
                    continue;
                double cur_min = cur_node_->hyper_rings[i].first;
                double cur_max = cur_node_->hyper_rings[i].second;

                if(cur_min > sub_min)
                {

                    cur_node_->hyper_rings[i].first = cur_node_->ptr_sub_tree[j]->hyper_rings[i].first;

                }

                if(cur_max < sub_max)
                {

                    cur_node_->hyper_rings[i].second = cur_node_->ptr_sub_tree[j]->hyper_rings[i].second;

                }

            }

        }

    }

}



template <class DType>
void PM_Tree<DType>::assign_node_all_value(PM_Node<DType>* cur_node_, PM_Node<DType>* new_node_)
{

    cur_node_->dist_to_parent = new_node_->dist_to_parent;
    cur_node_->feature_val = new_node_->feature_val;
    cur_node_->id = new_node_->id;
    cur_node_->node_category = new_node_->node_category;
    cur_node_->parent_node = new_node_->parent_node;
    cur_node_->ptr_sub_tree = new_node_->ptr_sub_tree;
    cur_node_->range = new_node_->range;
    cur_node_->pivot_distance = new_node_->pivot_distance;
    //cur_node_->hyper_rings = cur_node_->hyper_rings; //DUVIDA
    cur_node_->hyper_rings = new_node_->hyper_rings; //DUVIDA

    for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
    {

        cur_node_->ptr_sub_tree[i]->parent_node = cur_node_;

    }

}

template <class DType>
void PM_Tree<DType>::partition(std::vector<PM_Node<DType>*>& entries_, PM_Node<DType>* node_1_, PM_Node<DType>* node_2_)
{

    if(entries_[0]->node_category == 1)
    {

        node_1_->node_category = node_2_->node_category = 0;

    }
    else if(entries_[0]->node_category == 2)
    {

        node_1_->node_category = node_2_->node_category = 1;

    }
    else if(entries_[0]->node_category == 0)
    {

        node_1_->node_category = node_2_->node_category = 0;

    }
    else
    {

        throw std::invalid_argument("Error in partition");

    }

    for(size_t i = 0; i < entries_.size(); i++)
    {

        if(df->getDistance(entries_[i]->feature_val, node_1_->feature_val) <= df->getDistance(entries_[i]->feature_val, node_2_->feature_val))
        {

            node_1_->ptr_sub_tree.push_back(entries_[i]);
            entries_[i]->parent_node = node_1_;
            //std::cout << "NODE 1 - PARTITION\n";

        }
        else
        {

            node_2_->ptr_sub_tree.push_back(entries_[i]);
            entries_[i]->parent_node = node_2_;
            //std::cout << "NODE 2 - PARTITION\n";

        }

        entries_[i]->dist_to_parent = cal_dist_to_parent(entries_[i]);

    }

    node_1_->range = cal_cover_radius(node_1_);
    node_2_->range = cal_cover_radius(node_2_);

}

template <class DType>
void PM_Tree<DType>::promote(std::vector<PM_Node<DType>*>& entries_, PM_Node<DType>* node_1_, PM_Node<DType>* node_2_)
{


    std::vector<BasicArrayObject<DType>> v;

    for(size_t i = 0; i < entries_.size(); i++)
    {

        v.push_back(entries_[i]->feature_val);

    }

    Dataset<DType>* data = new Dataset<DType>(v, entries_.size(), entries_[0]->feature_val.size());

    pvt->generatePivots(data, df, 2);

    node_1_->feature_val = *pvt->getPivot(0);
    node_2_->feature_val = *pvt->getPivot(1);

    //std::cout << "PROMOTE : " << node_1_->feature_val.toStringWithOID() << "\n";
    //std::cout << "PROMOTE : " << node_2_->feature_val.toStringWithOID() << "\n\n";

    v.clear();
    delete (data);

}

template <class DType>
void PM_Tree<DType>::subRange_search(PM_Node<DType>** cur_node_address_, BasicArrayObject<DType>& q_feature_val_, double search_range_, std::vector<std::pair<double, size_t>>& res_vec_, double& dist_parent_q_, std::vector<double>& dist_q_pivot_)
{

    double dist = 0;
    PM_Node<DType>* cur_node_ = *cur_node_address_;
    Via_Node_Num++;

    if(is_leaf_node(cur_node_))
    {

        for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
        {

            if(fabs(dist_parent_q_ - cur_node_->ptr_sub_tree[i]->dist_to_parent) <= search_range_ && is_pivot_filter(dist_q_pivot_, search_range_, cur_node_->ptr_sub_tree[i]))
            {

                Via_subNode_Num++;
                dist = df->getDistance(cur_node_->ptr_sub_tree[i]->feature_val, q_feature_val_);
                Cal_Distance_Num++;

                if(dist <= search_range_)
                {

                    res_vec_.push_back(std::make_pair(dist, cur_node_->ptr_sub_tree[i]->id));

                }


            }

        }

    }
    else
    {

        for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
        {

            if(fabs(dist_parent_q_ - cur_node_->ptr_sub_tree[i]->dist_to_parent) <= search_range_ + cur_node_->ptr_sub_tree[i]->range)
            {

                Via_subNode_Num++;
                dist = df->getDistance(cur_node_->ptr_sub_tree[i]->feature_val, q_feature_val_);

                if(dist <= search_range_ + cur_node_->ptr_sub_tree[i]->range && is_pivot_filter(dist_q_pivot_, search_range_, cur_node_->ptr_sub_tree[i]))
                {

                    subRange_search(&(cur_node_->ptr_sub_tree[i]), q_feature_val_, search_range_, res_vec_, dist, dist_q_pivot_);

                }

            }

        }

    }

}


template <class DType>
bool PM_Tree<DType>::is_pivot_filter(std::vector<double>& dist_q_pivot_, double range_, PM_Node<DType>* cur_node_)
{

    if(is_data_node(cur_node_))
    {

        double PD_i = 0;

        for(size_t i = 0; i < PD; ++i)
        {

            PD_i = cur_node_->pivot_distance[i];

            if(fabs(dist_q_pivot_[i] - PD_i) > range_)
            {

                Pivot_Filter_Num++;
                return false;

            }

        }

    }
    else
    {

        for(size_t i = 0; i < HR; ++i)
        {

            if((dist_q_pivot_[i] - range_ > cur_node_->hyper_rings[i].second) || (dist_q_pivot_[i] + range_) < cur_node_->hyper_rings[i].first)
            {

                Pivot_Filter_Num++;
                return false;

            }

        }

    }

    return true;

}

template <class DType>
bool PM_Tree<DType>::is_data_node(const PM_Node<DType>* node_)
{

    return node_->node_category == 2;

}

template <class DType>
bool PM_Tree<DType>::is_rounting_node(const PM_Node<DType>* node_)
{

    return node_->node_category == 0;

}

template <class DType>
void PM_Tree<DType>::sub_traverse_get_volumem(PM_Node<DType>* cur_node_, double & volume_)
{

    if(is_leaf_node(cur_node_))
    {

        return;

    }

    for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
    {

        double v = cal_hypersphere_volume(cur_node_->ptr_sub_tree[i]->range, cur_node_->ptr_sub_tree[i]->feature_val.size());
        volume_ += v;
        sub_traverse_get_volumem(cur_node_->ptr_sub_tree[i], volume_);

    }

}


template <class DType>
double PM_Tree<DType>::cal_hypersphere_volume(double radius_, size_t dim_)
{

    double a = 2 * pow(M_PI, double(dim_) / 2.0);
    double b = double(dim_) * tgamma(double(dim_) / 2.0);
    double c = pow(radius_, double(dim_));

    return (a*c)/b;

}

template <class DType>
double PM_Tree<DType>::minDistNode(PM_Node<DType>* cur_node_, BasicArrayObject<DType>& element, std::vector<double> dist_to_query)
{

//    double dist = df->getDistance(cur_node_->feature_val, element), ans;
//    std::cout << cur_node_->feature_val.toStringWithOID() << "\n";
//    std::cout << element.toStringWithOID() << "\n";
//    std::cout << "MINDIST = " << dist << " / RANGE = " << cur_node_->range << "\n\n\n";

//    if(dist <= cur_node_->range) //Dentro da bola
//    {

//        ans = 0.0;

//    }
//    else //Fora da bola
//    {

//        ans = dist - cur_node_->range;

//    }

//    return ans;

    double d_hr_max = 0.0, d_hr_min = 0.0, dist_p_i_q;

    for(size_t i = 0; i < pivot_vec.getCardinality(); i++)
    {

//        dist_p_i_q = df->getDistance(pivot_vec.getFeatureVector(i), element);
        dist_p_i_q = dist_to_query[i];
        d_hr_max = std::max(d_hr_max, dist_p_i_q - cur_node_->hyper_rings[i].second);
        d_hr_min = std::max(d_hr_min, cur_node_->hyper_rings[i].first - dist_p_i_q);

    }

    return std::max({0.0, df->getDistance(cur_node_->feature_val, element) - cur_node_->range, d_hr_max, d_hr_min});

}

template <class DType>
double PM_Tree<DType>::maxDistNode(PM_Node<DType>* cur_node_, BasicArrayObject<DType>& element, std::vector<double> dist_to_query)
{

    //return df->getDistance(cur_node_->feature_val, element) + cur_node_->range;

    double d_hr = std::numeric_limits<double>::max();

    for(size_t i = 0; i < pivot_vec.getCardinality(); i++)
    {

//        d_hr = std::min(d_hr, df->getDistance(pivot_vec.getFeatureVector(i), element) + cur_node_->hyper_rings[i].second);
        d_hr = std::min(d_hr, dist_to_query[i] + cur_node_->hyper_rings[i].second);

    }

    return std::min(df->getDistance(cur_node_->feature_val, element) + cur_node_->range, d_hr);

}

//template <class DType>
//double PM_Tree<DType>::minLimInf(PM_Node<DType>* cur_node_)
//{

//    double ans = std::numeric_limits<double>::max();

//    for(size_t i = 0; i < pivot_vec.getCardinality(); i++)
//    {

//        if(cur_node_->hyper_rings[i].first < ans)
//        {

//            ans = cur_node_->hyper_rings[i].first;

//        }

//    }

//    return ans;

//}

//template <class DType>
//double PM_Tree<DType>::maxLimSup(PM_Node<DType>* cur_node_)
//{

//    double ans = std::numeric_limits<double>::min();

//    for(size_t i = 0; i < pivot_vec.getCardinality(); i++)
//    {

//        if(cur_node_->hyper_rings[i].second > ans)
//        {

//            ans = cur_node_->hyper_rings[i].second;

//        }

//    }

//    return ans;

//}

template <class DType>
PM_Tree<DType>::PM_Tree(Dataset<DType>* dataset, DistanceFunction<BasicArrayObject<DType>>* df_, Pivot<DType>* pvt_, size_t m_, size_t pivot_num_)
{

    M = m_;
//    promotefunc_e = RANDOM_e;
    root = nullptr;
    Max_Level = 0;

    df = df_;

    Pivot_Num = pivot_num_;
    HR = pivot_num_;
    PD = pivot_num_;

    pvt = pvt_;
    pvt->generatePivots(dataset, df, pivot_num_);

    pivot_vec = Dataset<DType>();
    pivot_vec.setCardinality(pivot_num_);
    pivot_vec.setDimensionality(dataset->getDimensionality());

    for(size_t x = 0; x < pivot_num_; x++)
    {

        pivot_vec.push_back(*pvt->getPivot(x));

    }

    for(size_t x = 0; x < dataset->getCardinality(); x++)
    {

        insert(dataset->getFeatureVector(x), x);

    }

}

template <class DType>
PM_Tree<DType>::~PM_Tree()
{


}

template <class DType>
void PM_Tree<DType>::set_pivot(size_t pivot_num_, const Dataset<DType>& pivot_vec_)
{

    Pivot_Num = pivot_num_;
    HR = pivot_num_;
    PD = pivot_num_;
    pivot_vec = pivot_vec_;

}


template <class DType>
void PM_Tree<DType>::insert(BasicArrayObject<DType>& feature_val_, int id_)
{

    //std::cout << "INSERT: " << feature_val_.toStringWithOID() << "\n\n\n";
    PM_Node<DType>* new_node = new PM_Node<DType>(nullptr, 2, -1, id_);
    new_node->feature_val = feature_val_;
    new_node->range = -1;
    new_node->id = id_; //DUVIDA
    update_pivot_distance(new_node->pivot_distance, feature_val_);

    if(root == nullptr)
    {

        PM_Node<DType>* new_root = new PM_Node<DType>(nullptr, 1, -1, -1);
        new_root->feature_val = feature_val_;
        new_root->ptr_sub_tree.push_back(new_node);
        new_root->range = cal_cover_radius(new_root);
        update_hyper_rings(new_root->hyper_rings, feature_val_);
        root = new_root;

        new_node->parent_node = new_root;
        new_node->dist_to_parent = cal_dist_to_parent(new_node);

    }
    else
    {

        insert(&root, &new_node);

    }

}


template <class DType>
size_t PM_Tree<DType>::get_pivot_filter_num()
{

    return Pivot_Filter_Num;

}

template <class DType>
size_t PM_Tree<DType>::get_cal_distance_num()
{

    return df->getDistanceCount(); //return Cal_Distance_Num;

}

template <class DType>
size_t PM_Tree<DType>::get_via_subNode_num()
{

    return Via_subNode_Num;

}

template <class DType>
size_t PM_Tree<DType>::get_via_node_num()
{

    return Via_Node_Num;

}

template <class DType>
void PM_Tree<DType>::range_search(BasicArrayObject<DType>& q_feature_val_, double search_range_, std::vector<std::pair<double, size_t>>& res_vec_)
{

    Via_Node_Num = 0;
    Via_subNode_Num = 0;
    Pivot_Filter_Num = 0;
    Cal_Distance_Num = 0;

    std::vector<double> dist_q_pivot((std::max(PD,HR)));

    for(size_t i = 0; i < pivot_vec.getSize(); ++i)
    {

        dist_q_pivot[i] = df->getDistance(q_feature_val_, *pivot_vec.instance(i));

    }

    double dist = df->getDistance(root->feature_val, q_feature_val_);
    Via_subNode_Num++;
    Via_Node_Num++;
    subRange_search(&(root), q_feature_val_, search_range_, res_vec_, dist, dist_q_pivot);
    std::sort(res_vec_.begin(), res_vec_.end());

}


template <class DType>
void PM_Tree<DType>::update_level(PM_Node<DType>* cur_node_, size_t level_)
{

    if(level_ > Max_Level)
    {

        Max_Level = level_;

    }

    for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); ++i)
    {

        cur_node_->ptr_sub_tree[i]->level = level_;

        if(!(is_data_node(cur_node_->ptr_sub_tree[i])))
        {

            update_level(cur_node_->ptr_sub_tree[i], level_ + 1);

        }

    }

}



template <class DType>
double PM_Tree<DType>::traverse_get_volume()
{

    double volume = cal_hypersphere_volume(root->range, root->feature_val.size());
    sub_traverse_get_volumem(root, volume);
    return volume;

}


template <class DType>
PM_Node<DType>* PM_Tree<DType>::get_root()
{

    return root;

}

template <class DType>
void PM_Tree<DType>::output_node_info(std::ofstream& out_, const PM_Node<DType>* cur_node_)
{

    out_ << cur_node_->level << " ";

    if(Pivot_Num == 0)
    {

        if(is_data_node(cur_node_))
        {

            out_ << 0 << " " << 1 << "     ";

        }
        else
        {

            out_ << cur_node_->range << " " << cur_node_->ptr_sub_tree.size() << "     ";

        }

    }
    else
    {

        if(is_data_node(cur_node_))
        {

            out_ << 0 << " " << 1 << "     ";

            for(size_t i = 0; i < cur_node_->pivot_distance.size(); ++i)
            {

                out_ << cur_node_->pivot_distance[i] << "  " << cur_node_->pivot_distance[i] << "  ";

            }

        }
        else
        {

            out_ << cur_node_->range << " " << cur_node_->ptr_sub_tree.size() << "     ";

            for(size_t i = 0; i < cur_node_->hyper_rings.size(); ++i)
            {

                out_ << cur_node_->hyper_rings[i].first << "  " << cur_node_->hyper_rings[i].second << "  ";

            }

        }

    }

    out_ << std::endl;

}


template <class DType>
void PM_Tree<DType>::traverse_bread_tree(std::string out_file_path)
{

    std::ofstream outfile(out_file_path);
    std::queue<PM_Node<DType>*> my_queue;
    my_queue.push(root);
    size_t cur_level = 0;

    if(outfile.is_open())
    {

        while(!my_queue.empty())
        {

            PM_Node<DType>* cur_node_ = my_queue.front();
            my_queue.pop();

            if(cur_node_->level - cur_level >= 2)
            {

                throw std::invalid_argument("Error in traverse_bread_tree\n");

            }

            if(cur_node_->level != cur_level)
            {

                cur_level = cur_node_->level;

            }

            output_node_info(outfile, cur_node_);

            for(size_t i = 0; i < cur_node_->ptr_sub_tree.size(); i++)
            {

                my_queue.push(cur_node_->ptr_sub_tree[i]);

            }

        }

    }
    else
    {

        std::invalid_argument("Error in traverse_bread_tree\n");

    }

}



template <class DType>
size_t PM_Tree<DType>::getMaxLevel()
{

    return Max_Level;

}


//template <class DType>
//void PM_Tree<DType>::test()
//{

//    BasicArrayObject<DType> element = BasicArrayObject<DType>(15, {19.0,4.0});
//    std::priority_queue<PM_Partition<DType>, std::vector<PM_Partition<DType>>, Compare_PM_Nodes<DType>> queue;
//    queue.push(PM_Partition(get_root(), minDistNode(get_root(), element), maxDistNode(get_root(), element)));
//    PM_Partition<DType> ft;
//    PM_Node<DType>* curr = nullptr;

//    while(!queue.empty())
//    {

//        ft = queue.top();
//        curr = ft.node;
//        queue.pop();

//        for(size_t i = 0; i < curr->ptr_sub_tree.size(); i++)
//        {

//            if(!is_leaf_node(curr->ptr_sub_tree[i]))
//            {

//                queue.push(PM_Partition(curr->ptr_sub_tree[i], minDistNode(curr->ptr_sub_tree[i], element), maxDistNode(curr->ptr_sub_tree[i], element)));

//            }

//        }

//        std::cout << "FEATURE = " << curr->feature_val.toStringWithOID() << "\n";
//        std::cout << "MIN DIST = " << ft.min << "\n";
//        std::cout << "MAX DIST = " << ft.max << "\n\n\n";

//    }

////    std::cout << "MIN LIM INF = " << minLimInf(get_root()->ptr_sub_tree[2]) << "\n";
////    std::cout << "MAX LIM INF = " << maxLimSup(get_root()->ptr_sub_tree[2]) << "\n";

//}


template <class DType>
void PM_Tree<DType>::kNN(BasicArrayObject<DType> query, size_t k, std::vector<KnnEntry<DType>>& ansVec)
{

    Via_Node_Num = 0;
    Via_subNode_Num = 0;
    Pivot_Filter_Num = 0;
    Cal_Distance_Num = 0;
    df->resetStatistics();
    leafNodeAccess = 0;

    std::priority_queue<PM_Partition<DType>, std::vector<PM_Partition<DType>>, Compare_PM_Nodes<DType>> nodeQueue;
    std::priority_queue<KnnEntry<DType>, std::vector<KnnEntry<DType>>, std::greater<KnnEntry<DType>>> candidatesQueue;
    std::priority_queue<KnnEntry<DType>, std::vector<KnnEntry<DType>>, std::less<KnnEntry<DType>>> resultQueue;
    PM_Node<DType>* node = nullptr;
    PM_Partition<DType> partition;
    std::vector<double> query_to_pivot;
    update_pivot_distance(query_to_pivot, query);
    nodeQueue.push(PM_Partition<DType>(get_root(), minDistNode(get_root(), query, query_to_pivot), maxDistNode(get_root(), query, query_to_pivot)));

    while(!nodeQueue.empty() || candidatesQueue.size() > 0)
    {

        if(candidatesQueue.size() == 0)
        {

            partition = nodeQueue.top();
            node = partition.node;
            nodeQueue.pop();

            if(is_leaf_node(node))
            {

                leafNodeAccess++;

                for(size_t i = 0; i < node->ptr_sub_tree.size(); i++)
                {

                    candidatesQueue.push(KnnEntry<DType>(node->ptr_sub_tree[i]->feature_val, df->getDistance(query, node->ptr_sub_tree[i]->feature_val)));

                }

            }
            else
            {

                for(size_t i = 0; i < node->ptr_sub_tree.size(); i++)
                {

                    nodeQueue.push(PM_Partition<DType>(node->ptr_sub_tree[i], minDistNode(node->ptr_sub_tree[i], query, query_to_pivot), maxDistNode(node->ptr_sub_tree[i], query, query_to_pivot)));

                }

            }

        }
        else if((nodeQueue.size() > 0) && nodeQueue.top().min < candidatesQueue.top().distance)
        {

            partition = nodeQueue.top();
            node = partition.node;
            nodeQueue.pop();

            if(is_leaf_node(node))
            {

                leafNodeAccess++;

                for(size_t i = 0; i < node->ptr_sub_tree.size(); i++)
                {

                    candidatesQueue.push(KnnEntry<DType>(node->ptr_sub_tree[i]->feature_val, df->getDistance(query, node->ptr_sub_tree[i]->feature_val)));

                }

            }
            else
            {

                for(size_t i = 0; i < node->ptr_sub_tree.size(); i++)
                {

                    nodeQueue.push(PM_Partition<DType>(node->ptr_sub_tree[i], minDistNode(node->ptr_sub_tree[i], query, query_to_pivot), maxDistNode(node->ptr_sub_tree[i], query, query_to_pivot)));

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

    ansVec = dequeueInOrderPM_Results(resultQueue);
    std::reverse(ansVec.begin(), ansVec.end());

    while(!resultQueue.empty())
    {

        resultQueue.pop();

    }

    while(!nodeQueue.empty())
    {

        nodeQueue.pop();

    }

    while(!candidatesQueue.empty())
    {

        candidatesQueue.pop();

    }

}

template <class DType>
size_t PM_Tree<DType>::getLeafNodeAccess()
{

    return leafNodeAccess;

}

template <class DType>
std::vector<PM_Node<DType>*> PM_Tree<DType>::leafsNodes()
{

    std::queue<PM_Node<DType>*> queue;
    queue.push(get_root());
    PM_Node<DType>* curr = nullptr;
    std::vector<PM_Node<DType>*> ans;

    while(!queue.empty())
    {

        curr = queue.front();
        queue.pop();

        if(is_leaf_node(curr))
        {

            ans.push_back(curr);

        }
        else
        {

            for(size_t i = 0; i < curr->ptr_sub_tree.size(); i++)
            {

                queue.push(curr->ptr_sub_tree[i]);

            }

        }


    }

    return ans;

}


#endif // PM_TREE_H
