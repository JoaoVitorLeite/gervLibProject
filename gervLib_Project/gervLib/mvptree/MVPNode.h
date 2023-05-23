#ifndef MVPNODE_H
#define MVPNODE_H

#include <MVPData.h>

template <class T>
class MVPNode
{

private:
    std::vector<vp_t<T>> vps;
    std::vector<MVPNode<T>*> childrens;
    std::vector<std::vector<double>> splits;
    std::vector<datapoint_t<T>> data;
    std::vector<std::vector<double>> leaf_dists;
    size_t pageID;
    bool leaf;

public:
    MVPNode();
    MVPNode(size_t fo);
    MVPNode(size_t lpn, size_t fo);
    MVPNode(size_t lpn, size_t fo, size_t ns);

    void setChildren(size_t pos, MVPNode<T>* child);
    void setPivot(size_t pos, vp_t<T> vp);
    void setSplit(size_t r, size_t c, double d);
    void setLeafDists(size_t r, size_t c, double d);
    void setPageID(size_t page);

    MVPNode<T>* getChild(size_t pos);
    vp_t<T> getVP(size_t pos);
    double getSplit(size_t r, size_t c);
    double getLeafDists(size_t r, size_t c);
    datapoint_t<T> getData(size_t pos);
    size_t getPageID();

    bool isLeaf();
    void pushData(std::vector<datapoint_t<T>> dp);
    void pushData(datapoint_t<T> dp);
    void initializeLeafDists();
    size_t getDataSize();
    std::vector<BasicArrayObject<T>> purgeData();
    void clear();

};



template <class T>
MVPNode<T>::MVPNode()
{

    pageID = std::numeric_limits<size_t>::max();
    leaf = false;

}

template <class T>
MVPNode<T>::MVPNode(size_t fo)
{

    childrens = std::vector<MVPNode<T>*>(fo, nullptr);
    pageID = std::numeric_limits<size_t>::max();
    leaf = false;

}

template <class T>
MVPNode<T>::MVPNode(size_t lpn, size_t fo)
{

    childrens = std::vector<MVPNode<T>*>(fo, nullptr);
    vps = std::vector<vp_t<T>>(lpn, {nullptr});
    pageID = std::numeric_limits<size_t>::max();
    leaf = false;

}

template <class T>
MVPNode<T>::MVPNode(size_t lpn, size_t fo, size_t ns)
{

    childrens = std::vector<MVPNode<T>*>(fo, nullptr);
    vps = std::vector<vp_t<T>>(lpn, {nullptr});
    splits = std::vector<std::vector<double>>(lpn, std::vector<double>(ns, -1.0));
    pageID = std::numeric_limits<size_t>::max();
    leaf = false;

}


template <class T>
void MVPNode<T>::setChildren(size_t pos, MVPNode<T>* child)
{

    if((pos < 0) || (pos >= childrens.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    childrens[pos] = child;

}

template <class T>
void MVPNode<T>::setPivot(size_t pos, vp_t<T> vp)
{

    if((pos < 0) || (pos >= vps.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    vps[pos] = vp;

}

template <class T>
void MVPNode<T>::setSplit(size_t r, size_t c, double d)
{

    if((r < 0) || (r > splits.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    if((c < 0) || (c > splits[0].size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    splits[r][c] = d;

}

template <class T>
void MVPNode<T>::setLeafDists(size_t r, size_t c, double d)
{

    if((r < 0) || (r >= leaf_dists.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    if((c < 0) || (c >= vps.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    leaf_dists[r][c] = d;

}

template <class T>
void MVPNode<T>::setPageID(size_t page)
{

    pageID = page;

}

template <class T>
MVPNode<T>* MVPNode<T>::getChild(size_t pos)
{

    if((pos < 0) || (pos >= childrens.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    return childrens[pos];

}

template <class T>
vp_t<T> MVPNode<T>::getVP(size_t pos)
{

    if((pos < 0) || (pos >= vps.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    return vps[pos];

}

template <class T>
double MVPNode<T>::getSplit(size_t r, size_t c)
{

    if((r < 0) || (r > splits.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    if((c < 0) || (c > splits[0].size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    return splits[r][c];

}

template <class T>
double MVPNode<T>::getLeafDists(size_t r, size_t c)
{

    if((r < 0) || (r >= leaf_dists.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    if((c < 0) || (c >= vps.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    return leaf_dists[r][c];

}

template <class T>
datapoint_t<T> MVPNode<T>::getData(size_t pos)
{

    if((pos < 0) || (pos >= data.size()))
    {

        throw std::invalid_argument("Index out of bounds !_!");

    }

    return data[pos];

}

template <class T>
size_t MVPNode<T>::getPageID()
{

    return pageID;

}

template <class T>
bool MVPNode<T>::isLeaf()
{

    return (!data.empty() || leaf);

}

template <class T>
void MVPNode<T>::pushData(std::vector<datapoint_t<T>> dp)
{

    data.insert(data.end(), dp.begin(), dp.end());

}

template <class T>
void MVPNode<T>::pushData(datapoint_t<T> dp)
{

    data.push_back(dp);

}

template <class T>
void MVPNode<T>::initializeLeafDists()
{

    leaf_dists = std::vector<std::vector<double>>(data.size(), std::vector<double>(vps.size(), 0.0));

}

template <class T>
size_t MVPNode<T>::getDataSize()
{

    return data.size();

}

template <class T>
std::vector<BasicArrayObject<T>> MVPNode<T>::purgeData()
{

    std::vector<BasicArrayObject<T>> ans;

    for(size_t i = 0; i < data.size(); i++)
    {

        ans.push_back(*data[i].key);

    }

    return ans;

}

template <class T>
void MVPNode<T>::clear()
{

    data.clear();

    for(auto vec : leaf_dists)
        vec.clear();

    leaf_dists.clear();

    leaf = true;

}

#endif // MVPNODE_H
