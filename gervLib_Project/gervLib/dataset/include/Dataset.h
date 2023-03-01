#ifndef DATASET_H
#define DATASET_H

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <BasicArrayObject.h>
#include <regex>
#include <vector>
#include <DistanceFunction.h>
#include <Util.h>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>

template <class DType>
class Dataset{

    private:
        std::vector<BasicArrayObject<DType>> data;
        size_t cardinality;
        size_t dimensionality;
        size_t seed;

    private:
        std::vector<double> split(std::string str, std::string delimiter="[//s+,;]");
        std::vector<std::string> splitString(std::string str, std::string delimiter="[//s+,;]");
        size_t serializedSize();

    public:
        static void loadNumericDataset(Dataset<DType> *dataset, std::string filePath, size_t cardinality, size_t dimensionality);
        static void loadNumericDataset(Dataset<DType> *dataset, std::string filePath, std::string separator);
        static void loadTextDataset(Dataset<DType> *dataset, std::string filePath, size_t cardinality, size_t dimensionality);
        static void loadTextDataset(Dataset<DType> *dataset, std::string filePath, std::string separator="[//s+,;]");

    public:
        //Constructors and Destructors
        Dataset();
        Dataset(std::vector<BasicArrayObject<DType>> elements, size_t cardinality, size_t dimensionality);
        Dataset(std::vector<BasicArrayObject<DType>*> elements, size_t cardinality, size_t dimensionality);
        Dataset(std::vector<BasicArrayObject<DType>> elements, size_t cardinality, size_t dimensionality, size_t seed);
        ~Dataset();

        //Operator Overloading
        BasicArrayObject<DType>& operator[](const size_t pos);
        const BasicArrayObject<DType>& operator[](const size_t pos) const;
        bool operator==(Dataset<DType>* other);

        //Setters
        void setInstance(BasicArrayObject<DType> *instance, size_t pos);
        //void setInstanceSize();
        void setSeed(const size_t &value);
        void setCardinality(const size_t &value);
        void setDimensionality(const size_t &value);
        void setData(std::vector<BasicArrayObject<DType>> data);

        //Getters
        BasicArrayObject<DType>* instance(size_t pos);
        BasicArrayObject<DType>* getInstance(size_t pos);
        BasicArrayObject<DType> &getFeatureVector(const size_t pos);
        const BasicArrayObject<DType> &getFeatureVector(const size_t pos) const;
        size_t getSerializedSize();
        std::vector<BasicArrayObject<DType>> getElements();
        size_t getDimensionality() const;
        size_t getCardinality() const;
        size_t getSeed() const;
        size_t getSize() const;
        void clear();

        //Public Methods
        void reserve(size_t n);
        void push_back(BasicArrayObject<DType> instance);
        std::vector<BasicArrayObject<DType>*> sample(size_t size, bool reposition);
        Dataset<DType>* sampleDataset(size_t size, bool reposition);
        std::vector<BasicArrayObject<DType>*> sample(size_t size, bool reposition, size_t seed);
        Dataset<DType>* sampleDataset(size_t size, bool reposition, size_t seed);
        void swap(size_t pos1, size_t pos2);
        void shuffle();
        void shuffle(size_t seed);
        bool isEqual(Dataset* other);
        void saveDataset(std::string path, std::string separator = ",");
        void printDataset(std::string separator=",");
        std::pair<std::vector<BasicArrayObject<DType>*>, std::vector<BasicArrayObject<DType>*>> splitDataset(double p1, double p2);
        std::pair<std::vector<BasicArrayObject<DType>*>, std::vector<BasicArrayObject<DType>*>> splitDataset(double p1, double p2, size_t seed);
        unsigned char* serialize();
        void unserialize(unsigned char* dataIn);
        void saveToFile(std::string fileName);
        void loadFromFile(std::string fileName);

        //VPTREE
        template <class DistanceFunction>
        Dataset<DType> rangeQuery(const BasicArrayObject<DType> &qElement, const double radius, DistanceFunction df) const;

        template <class DistanceFunction>
        Dataset<DType> kNN(const BasicArrayObject<DType> &qElement, const size_t k, DistanceFunction df) const;

        template<class DistanceFunction>
        void getDistanceVector(DistanceFunction df, std::vector<double> *distances) const;

        template<class DistanceFunction>
        void getDistanceVector(const BasicArrayObject<DType> &element, DistanceFunction df, std::vector<double> *distances) const;

        template <class DistanceFunction>
        double medianDistance(const BasicArrayObject<DType> &pivot, DistanceFunction df);

        template <class DistanceFunction>
        double meanDistance(const BasicArrayObject<DType> &pivot, DistanceFunction df);

        template <class DistanceFunction>
        static double getMax(const Dataset<DType> &dataset, DistanceFunction df);

        template <class DistanceFunction>
        static double getMax(const Dataset<DType> &dataset, const BasicArrayObject<DType> &element, DistanceFunction df);

        template <class DistanceFunction>
        static double getMax(const Dataset<DType> &dataset, const BasicArrayObject<DType> &element, double mu, DistanceFunction df);

        void insert(uint_fast32_t pos, const BasicArrayObject<DType> &fv);

        bool contains(const BasicArrayObject<DType> &fv);

        void erase(uint_fast32_t position);

        void erase(const BasicArrayObject<DType> &fv);

        void resize(size_t n);


};

//**************************************************************************************************************************
//**************************************************************************************************************************
//**************************************************************************************************************************
//**************************************************************************************************************************



template <class DType>
std::vector<double> Dataset<DType>::split(std::string str, std::string delimiter)
{

    std::regex rgx(delimiter);
    std::sregex_token_iterator iter(str.begin(),str.end(),rgx,-1);
    std::sregex_token_iterator end;
    std::vector<double> v;
    double d;

    for(; iter != end; ++iter)
    {

        d = stod(*iter);
        v.push_back(d);

    }

    return v;

}



template <class DType>
std::vector<std::string> Dataset<DType>::splitString(std::string str, std::string delimiter)
{

    std::regex rgx(delimiter);
    std::sregex_token_iterator iter(str.begin(),str.end(),rgx,-1);
    std::sregex_token_iterator end;
    std::vector<std::string> v;

    for(; iter != end; ++iter)
    {

        v.push_back(*iter);

    }

    return v;

}



template <class DType>
size_t Dataset<DType>::serializedSize()
{

    size_t total = 3*sizeof(size_t);

    for(size_t x = 0; x < getCardinality(); x++)
    {

        total += sizeof(size_t) + data[x].getSerializedSize();

    }

    return total;

}



template <class DType>
void Dataset<DType>::loadNumericDataset(Dataset<DType> *dataset, std::string filePath, size_t cardinality, size_t dimensionality)
{

    dataset->setSeed(0);
    dataset->setData(std::vector<BasicArrayObject<DType>>());
    dataset->setCardinality(cardinality);
    dataset->setDimensionality(dimensionality);
    dataset->reserve(cardinality);

    std::ifstream file(filePath);
    double readin;

    if (file.is_open())
    {

        for (size_t x = 0; x < dataset->getCardinality(); x++)
        {

            BasicArrayObject<double> aux;
            aux.setOID(x);
            for (size_t y = 0; y < dataset->getDimensionality(); y++)
            {

                file >> readin;
                aux.set(readin);

            }
            dataset->push_back(aux);

        }

    }
    else
    {

        throw std::runtime_error("File not found : " + filePath);

    }

    file.close();

}



template <class DType>
void Dataset<DType>::loadTextDataset(Dataset<DType> *dataset, std::string filePath, size_t cardinality, size_t dimensionality)
{

    dataset->setSeed(0);
    dataset->setCardinality(cardinality);
    dataset->setDimensionality(dimensionality);
    dataset->setData(std::vector<BasicArrayObject<DType>>());
    dataset->reserve(cardinality);

    std::ifstream file;

    file.open(filePath);

    if (file.is_open())
    {

        for(size_t x = 0; x < dataset->getCardinality(); x++)
        {

            BasicArrayObject<std::vector<char>> aux;
            aux.setOID(x);

            for(size_t y = 0; y < dataset->getDimensionality(); y++)
            {

                std::string s;
                file >> s;

                std::vector<char> charVec(s.begin(), s.end());
                aux.set(charVec);

            }

            dataset->push_back(aux);

        }

    }
    else
    {

        throw std::runtime_error("File not found : " + filePath);

    }

    file.close();

}


template <class DType>
void Dataset<DType>::loadTextDataset(Dataset<DType> *dataset, std::string filePath, std::string separator)
{

    dataset->setSeed(0);
    dataset->setCardinality(0);
    dataset->setDimensionality(0);
    dataset->setSeed(0);
    dataset->setData(std::vector<BasicArrayObject<DType>>());

    std::ifstream file;
    file.open(filePath);

    if (file.is_open())
    {

        std::string first, line;
        size_t oid = 0;
        std::getline(file, first);
        BasicArrayObject<std::vector<char>> aux;

        aux.setOID(oid++);

        std::vector<std::string> vecST = dataset->splitString(first, separator);

        dataset->setDimensionality(vecST.size());

        for(size_t x = 0; x < dataset->getDimensionality(); x++)
        {

            std::string s = vecST[x];
            std::vector<char> charVec(s.begin(), s.end());
            aux.set(charVec);

        }

        dataset->push_back(aux);

        while(std::getline(file, line))
        {

            BasicArrayObject<std::vector<char>> aux;
            aux.setOID(oid++);

            std::regex rgx(separator);
            std::sregex_token_iterator iter(line.begin(),line.end(),rgx,-1);
            std::sregex_token_iterator end;

            for(; iter != end; ++iter)
            {

                std::string st = *iter;
                std::vector<char> charVec(st.begin(), st.end());
                aux.set(charVec);

            }

            dataset->push_back(aux);

        }

    }
    else
    {

        throw std::runtime_error("File not found : " + filePath);

    }

    file.close();

}



template <class DType>
void Dataset<DType>::loadNumericDataset(Dataset<DType> *dataset, std::string filePath, std::string separator)
{

    dataset->setSeed(0);
    dataset->setCardinality(0);
    dataset->setDimensionality(0);
    dataset->setData(std::vector<BasicArrayObject<DType>>());
    std::ifstream file;
    std::string line;
    std::vector<double> aux;

    file.open(filePath);

    if (file.is_open()){

            size_t oid = 1;

            while(std::getline(file, line) && oid)
            {

                aux = dataset->split(line, separator);

                BasicArrayObject<double> inst = BasicArrayObject<double>(oid++-1,aux);
                dataset->push_back(inst);

            }

            dataset->setCardinality(dataset->getElements().size());
            dataset->setDimensionality(aux.size());


        }
        else
        {

            throw std::runtime_error("File not found : " + filePath);

        }

        file.close();

}



template <class DType>
Dataset<DType>::Dataset()
{

    setSeed(0);
    setCardinality(0);
    setDimensionality(0);
    setData(std::vector<BasicArrayObject<DType>>());


}



template <class DType>
Dataset<DType>::Dataset(std::vector<BasicArrayObject<DType>> elements, size_t cardinality, size_t dimensionality)
{

    setSeed(0);
    setCardinality(cardinality);
    setDimensionality(dimensionality);
    setData(elements);

}

template <class DType>
Dataset<DType>::Dataset(std::vector<BasicArrayObject<DType>*> elements, size_t cardinality, size_t dimensionality)
{

    setSeed(0);
    setCardinality(cardinality);
    setDimensionality(dimensionality);
    data = std::vector<BasicArrayObject<DType>>(cardinality);

    BasicArrayObject<DType>* b;

    for(size_t x = 0; x < getCardinality(); x++){
        b = elements[x];
        data[x] = *b;
    }

}



template <class DType>
Dataset<DType>::Dataset(std::vector<BasicArrayObject<DType>> elements, size_t cardinality, size_t dimensionality, size_t seed)
{

    setSeed(seed);
    setCardinality(cardinality);
    setDimensionality(dimensionality);
    setData(elements);

}



template <class DType>
Dataset<DType>::~Dataset()
{

    clear();
    //data.erase(data.begin(), data.end());

}



template <class DType>
BasicArrayObject<DType>& Dataset<DType>::operator[](const size_t pos)
{

    BasicArrayObject<DType> ans = BasicArrayObject<DType>();

    if(((pos < getCardinality()) && (pos >= 0)))
    {

        ans = data[pos];

    }
    else
    {

        throw std::invalid_argument("Out of Bound Exception !_!");

    }

    return ans;

}



template <class DType>
const BasicArrayObject<DType>& Dataset<DType>::operator[](const size_t pos) const
{

    BasicArrayObject<DType> ans = BasicArrayObject<DType>();

    if(((pos < getCardinality()) && (pos >= 0)))
    {

        ans = data[pos];

    }
    else
    {

        throw std::invalid_argument("Out of Bound Exception !_!");

    }

    return ans;

}



template <class DType>
bool Dataset<DType>::operator==(Dataset<DType>* other)
{

    return isEqual(other);

}



template <class DType>
void Dataset<DType>::setInstance(BasicArrayObject<DType> *instance, size_t pos)
{

    if(((pos < getCardinality()) & (pos >= 0)))
    {

        data[pos] = *instance;

    }
    else
    {

        throw std::invalid_argument("Out of Bound Exception !_!");

    }

}



//template <class DType>
//void Dataset<DType>::setInstanceSize()
//{

//    instanceSize = instance(0)->getSerializedSize();

//}



template <class DType>
void Dataset<DType>::setSeed(const size_t &value)
{

    this->seed = value;

}



template <class DType>
void Dataset<DType>::setCardinality(const size_t &value)
{

    this->cardinality = value;

}



template <class DType>
void Dataset<DType>::setDimensionality(const size_t &value)
{

    this->dimensionality = value;

}



template <class DType>
void Dataset<DType>::setData(std::vector<BasicArrayObject<DType>> data)
{

    this->data = data;

}



template <class DType>
BasicArrayObject<DType>* Dataset<DType>::instance(size_t pos)
{

    BasicArrayObject<DType>* ans = nullptr;

    if(((pos < getCardinality()) & (pos >= 0)))
    {

        ans = &data[pos];

    }
    else
    {

        throw std::invalid_argument("Out of Bound Exception !_!");

    }

    return ans;

}



template <class DType>
BasicArrayObject<DType> *Dataset<DType>::getInstance(size_t pos)
{

    return instance(pos);

}



template <class DType>
BasicArrayObject<DType> &Dataset<DType>::getFeatureVector(const size_t pos){ return data[pos]; }

template <class DType>
const BasicArrayObject<DType> &Dataset<DType>::getFeatureVector(const size_t pos) const{ return data[pos]; }

template <class DType>
size_t Dataset<DType>::getSerializedSize()
{

    return serializedSize();

}



template <class DType>
std::vector<BasicArrayObject<DType>> Dataset<DType>::getElements()
{

    return data;

}



template <class DType>
size_t Dataset<DType>::getDimensionality() const
{

    return dimensionality;

}



template <class DType>
size_t Dataset<DType>::getCardinality() const
{

    return cardinality;

}



template <class DType>
size_t Dataset<DType>::getSeed() const
{

    return seed;

}



template <class DType>
size_t Dataset<DType>::getSize() const
{

    return cardinality;

}



template <class DType>
void Dataset<DType>::clear()
{

    data.clear();

}



template <class DType>
void Dataset<DType>::reserve(size_t n)
{

    data.reserve(n);

}



template <class DType>
void Dataset<DType>::push_back(BasicArrayObject<DType> instance)
{

    data.push_back(instance);

    if (data.size() > cardinality)
    {

            ++cardinality;

    }

}



template <class DType>
std::vector<BasicArrayObject<DType>*> Dataset<DType>::sample(size_t size, bool reposition)
{

    std::vector<BasicArrayObject<DType>*> ans;
    ans.reserve(size);
    size_t* aux;

    if(reposition) aux = randomNumber(0, getCardinality(), size, getSeed());
    else aux = uniqueRandomNumber(0, getCardinality(), size, getSeed());

    for(size_t x = 0; x < size; x++)
        ans[x] = instance(aux[x]);

    delete[] aux;

    return ans;

}

template <class DType>
Dataset<DType>* Dataset<DType>::sampleDataset(size_t size, bool reposition)
{

    std::vector<BasicArrayObject<DType>> ans;
    size_t* aux;

    if(reposition) aux = randomNumber(0, getCardinality(), size, getSeed());
    else aux = uniqueRandomNumber(0, getCardinality(), size, getSeed());

    for(size_t x = 0; x < size; x++)
        ans.push_back(data[aux[x]]);

    delete[] aux;

    return new Dataset<DType>(ans, ans.size(), ans[0].size());

}

template <class DType>
std::vector<BasicArrayObject<DType>*> Dataset<DType>::sample(size_t size, bool reposition, size_t seed)
{

    std::vector<BasicArrayObject<DType>*> ans;
    ans.reserve(size);
    size_t* aux;

    if(reposition) aux = randomNumber(0, getCardinality(), size, seed);
    else aux = uniqueRandomNumber(0, getCardinality(), size, seed);

    for(size_t x = 0; x < size; x++)
        ans[x] = instance(aux[x]);

    delete[] aux;

    return ans;

}

template <class DType>
Dataset<DType>* Dataset<DType>::sampleDataset(size_t size, bool reposition, size_t seed)
{

    std::vector<BasicArrayObject<DType>> ans;
    size_t* aux;

    if(reposition) aux = randomNumber(0, getCardinality(), size, seed);
    else aux = uniqueRandomNumber(0, getCardinality(), size, seed);

    for(size_t x = 0; x < size; x++)
        ans.push_back(data[aux[x]]);


    delete[] aux;

    return new Dataset<DType>(ans, ans.size(), ans[0].size());

}

template <class DType>
void Dataset<DType>::swap(size_t pos1, size_t pos2)
{

    BasicArrayObject<DType>* temp = instance(pos1);
    setInstance(instance(pos2), pos1);
    setInstance(temp, pos2);

}



template <class DType>
void Dataset<DType>::shuffle()
{

    srand(getSeed());

    for(size_t z = getCardinality() - 1; z > 0; z--)
    {

        size_t aux = rand() % (z+1);
        swap(z, aux);

    }

}



template <class DType>
void Dataset<DType>::shuffle(size_t seed)
{

    srand(seed);

    for(size_t z = getCardinality() - 1; z > 0; z--)
    {

        size_t aux = rand() % (z+1);
        swap(z, aux);

    }

}



template <class DType>
bool Dataset<DType>::isEqual(Dataset<DType>* other)
{

    bool equal = true;

    if((getCardinality() == other->getCardinality()) && (getDimensionality() == other->getDimensionality()))
    {

        for(size_t x = 0; x < getCardinality(); x++)
        {

            if(!(instance(x)->isEqual(other->instance(x))))
            {

                equal = false;
                break;

            }

        }

    }

    return equal;

}



template <class DType>
void Dataset<DType>::saveDataset(std::string path, std::string separator)
{

    std::ofstream mFile(path);

    for(size_t x = 0; x < getCardinality(); x++)
    {

        mFile << instance(x)->toString(separator) << "\n";

    }

    mFile.close();

}



template <class DType>
void Dataset<DType>::printDataset(std::string separator)
{

    for (size_t x = 0; x < getCardinality(); x++)
    {

        std::cout << data[x].toStringWithOID(separator) << std::endl;

    }

}



template <class DType>
std::pair<std::vector<BasicArrayObject<DType>*>, std::vector<BasicArrayObject<DType>*>> Dataset<DType>::splitDataset(double p1, double p2)
{

    size_t pos, sizeTrain = getCardinality() * p1, sizeTest = std::max(getCardinality() - sizeTrain*1.0, getCardinality() * p2);
    std::vector<BasicArrayObject<DType>*> ansTrain = std::vector<BasicArrayObject<DType>*>(sizeTrain);
    std::vector<BasicArrayObject<DType>*> ansTest = std::vector<BasicArrayObject<DType>*>(sizeTest);

    size_t* aux;

    aux = shuffleIndex(0, getCardinality(), getSeed());

    for(size_t x = 0; x < sizeTrain; x++)
        ansTrain[x] = instance(aux[x]);

    pos = 0;

    for(size_t y = sizeTrain; y < getCardinality(); y++)
        ansTest[pos++] = instance(aux[y]); //Pode ser trocado

    return std::make_pair(ansTrain, ansTest);

}



template <class DType>
std::pair<std::vector<BasicArrayObject<DType>*>, std::vector<BasicArrayObject<DType>*>> Dataset<DType>::splitDataset(double p1, double p2, size_t seed)
{

    size_t pos, sizeTrain = getCardinality() * p1, sizeTest = std::max(getCardinality() - sizeTrain*1.0, getCardinality() * p2);
    std::vector<BasicArrayObject<DType>*> ansTrain = std::vector<BasicArrayObject<DType>*>(sizeTrain);
    std::vector<BasicArrayObject<DType>*> ansTest = std::vector<BasicArrayObject<DType>*>(sizeTest);

    size_t* aux;

    aux = shuffleIndex(0, getCardinality(), seed);

    for(size_t x = 0; x < sizeTrain; x++)
        ansTrain[x] = instance(aux[x]);

    pos = 0;

    for(size_t y = sizeTrain; y < getCardinality(); y++)
        ansTest[pos++] = instance(aux[y]); //Pode ser trocado

    return std::make_pair(ansTrain, ansTest);

}



template <class DType>
unsigned char* Dataset<DType>::serialize()
{

    //CARDINALITY
    //DIMENSIONALITY
    //SEED
    //INSTANCE SIZE + BASIC

    unsigned char* dataOut = new unsigned char[serializedSize()];

    size_t total = 0;

    memcpy(dataOut + total, &cardinality, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(dataOut + total, &dimensionality, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(dataOut + total, &seed, sizeof(size_t));
    total += sizeof(size_t);


    BasicArrayObject<DType> b = data[0];
    size_t inst = b.getSerializedSize();

    for(size_t x = 0; x < getCardinality(); x++)
    {

        b = data[x];
        inst = b.getSerializedSize();

        memcpy(dataOut + total, &inst, sizeof(size_t));
        total += sizeof(size_t);

        memcpy(dataOut + total, b.serialize(), inst);
        total += inst;

    }

    return dataOut;

}



template <class DType>
void Dataset<DType>::unserialize(unsigned char* dataIn)
{

    size_t total = 0;

    memcpy(&cardinality, dataIn + total, sizeof (size_t));
    total += sizeof(size_t);
    memcpy(&dimensionality, dataIn + total, sizeof (size_t));
    total += sizeof(size_t);
    memcpy(&seed, dataIn + total, sizeof (size_t));
    total += sizeof(size_t);

    BasicArrayObject<DType> b;
    for(size_t x = 0; x < getCardinality(); x++)
    {

        /*BasicArrayObject<DType>*/ b = BasicArrayObject<DType>();
        size_t inst;

        memcpy(&inst, dataIn + total, sizeof(size_t));
        total += sizeof(size_t);

        unsigned char* aux = new unsigned char[inst];

        memcpy(aux, dataIn + total, inst);
        total += inst;

        b.unserialize(aux);

        data.push_back(b);

        delete [] aux;

    }

}



template <class DType>
void Dataset<DType>::saveToFile(std::string fileName)
{

    std::ofstream file(fileName, std::ios::out | std::ios::binary);


    size_t sizeData = serializedSize();
    unsigned char* ans = new unsigned char[sizeof(size_t)];

    memcpy(ans, &sizeData, sizeof(size_t));

    file.write((char *)ans, sizeof(size_t));

    ans = serialize();

    file.write((char* )ans, serializedSize());
    file.close();

}



template <class DType>
void Dataset<DType>::loadFromFile(std::string fileName)
{

//    std::ifstream file(fileName, std::ios::in | std::ios::binary);
//    unsigned char* dataIn = new unsigned char[4*sizeof (size_t)];
//    file.read((char *)dataIn, 4*sizeof (size_t));

//    memcpy(&cardinality, dataIn, sizeof (size_t));
//    memcpy(&dimensionality, dataIn + sizeof (size_t), sizeof (size_t));
//    memcpy(&seed, dataIn + sizeof (size_t) + sizeof (size_t), sizeof (size_t));
//    memcpy(&instanceSize, dataIn + sizeof(size_t)*3, sizeof(size_t));

//    unsigned char* ans = new unsigned char[serializedSize()];
//    file.seekg(0);
//    file.read(ans, serializedSize());
//    unserialize(ans);

//    file.close();

//    delete [] (dataIn);
//    delete [] (ans);


    std::ifstream file(fileName, std::ios::in | std::ios::binary);
    size_t size = 0;
    unsigned char* dataIn = new unsigned char[sizeof(size_t)];

    file.read((char *)dataIn, sizeof(size_t));

    memcpy(&size, dataIn, sizeof(size_t));

    dataIn = new unsigned char[size];

    file.read((char *)dataIn, size);

    unserialize(dataIn);


}



template <class DType>
template <class DistanceFunction>
Dataset<DType> Dataset<DType>::rangeQuery(const BasicArrayObject<DType> &qElement, const double radius, DistanceFunction df) const
{

    Dataset<DType> answer;
    answer.reserve(getCardinality());
    answer.setCardinality(getCardinality());
    answer.setDimensionality(getDimensionality());
    size_t ansSize = 0;

    for (size_t i = 0; i < getCardinality(); i++){

        BasicArrayObject<DType> fv = operator [](i);

        double dist = df->getDistance(qElement, fv);

        if (dist <= radius)
        {

            answer.push_back(fv);
            ansSize++;

        }

    }

    answer.setCardinality(ansSize);

    if (ansSize == 0)
    {

        answer.setDimensionality(0);

    }

    return answer;

}



template <class DType>
template <class DistanceFunction>
Dataset<DType> Dataset<DType>::kNN(const BasicArrayObject<DType> &qElement, const size_t k, DistanceFunction df) const
{

    Dataset<DType> answer;

    answer.reserve(k);
    answer.setCardinality(k);
    answer.setDimensionality(getDimensionality());

    std::vector<double> distances;
    getDistanceVector(qElement, df, &distances);

    for (size_t i = 0; i < k; i++)
    {

        size_t min = i;

        for (size_t j = i+1; j < getCardinality(); j++)
        {

            if (distances[j] < distances[min])
            {

                min = j;

            }

        }

        answer.push_back(data[min]);
        std::swap(distances[i], distances[min]);

    }

    return answer;

}



template <class DType>
template <class DistanceFunction>
void Dataset<DType>::getDistanceVector(DistanceFunction df, std::vector<double> *distances) const
{

    distances->reserve(getCardinality()/2);

    for (size_t i = 0; i < getCardinality(); i++)
    {

        BasicArrayObject<DType> fv = operator [](i);

        for (size_t j = i+1; j < getCardinality(); j++)
        {

            distances->push_back(df->getDistance(fv, operator [](j)));

        }

    }

}



template <class DType>
template <class DistanceFunction>
void Dataset<DType>::getDistanceVector(const BasicArrayObject<DType> &element, DistanceFunction df, std::vector<double> *distances) const
{

    distances->reserve(getCardinality());

    for (size_t i = 0; i < getCardinality(); i++)
    {

        distances->push_back(df->getDistance(element, data[i]));

    }

}



template <class DType>
template <class DistanceFunction>
double Dataset<DType>::medianDistance(const BasicArrayObject<DType> &pivot, DistanceFunction df)
{

    std::vector<double> distances;
    getDistanceVector(pivot, df, &distances);

    const size_t size = distances.size();

    if ((size % 2) == 0)
    {

        std::nth_element(distances.begin(), distances.begin()+(size/2)-1, distances.end());
        double firstValue = distances[(size/2)-1];

        std::nth_element(distances.begin(), distances.begin()+(size/2), distances.end());
        double secondValue = distances[(size/2)];

        return (firstValue + secondValue)/2;

    }
    else
    {

        std::nth_element(distances.begin(), distances.begin()+(size/2), distances.end());
        return distances[size/2];

    }

}



template <class DType>
template <class DistanceFunction>
double Dataset<DType>::meanDistance(const BasicArrayObject<DType> &pivot, DistanceFunction df)
{

    std::vector<double> distances;
    getDistanceVector(pivot, df, &distances);

    return std::accumulate(distances.begin(), distances.end(), 0.0)/distances.size();

}



template <class DType>
template <class DistanceFunction>
double Dataset<DType>::getMax(const Dataset<DType> &dataset, DistanceFunction df)
{

    std::vector<double> distances;
    dataset.getDistanceVector(df, &distances);

    return *std::max_element(distances.begin(), distances.end());

}



template <class DType>
template <class DistanceFunction>
double Dataset<DType>::getMax(const Dataset<DType> &dataset, const BasicArrayObject<DType> &element, DistanceFunction df)
{

    std::vector<double> distances;
    dataset.getDistanceVector(element, df, &distances);
    return *std::max_element(distances.begin(), distances.end());

}



template <class DType>
template <class DistanceFunction>
double Dataset<DType>::getMax(const Dataset<DType> &dataset, const BasicArrayObject<DType> &element, double mu, DistanceFunction df)
{

    double max = 0.0;
    for (size_t i = 0; i < dataset.getCardinality(); i++){

        double dist = df->getDistance(dataset[i], element);

        if ((dist > max) && (dist <= mu))
            max = dist;
    }

    return max;

}



template <class DType>
void Dataset<DType>::insert(uint_fast32_t pos, const BasicArrayObject<DType> &fv)
{

    data.insert(data.begin() + pos, fv);
    setCardinality(data.size());

}



template <class DType>
bool Dataset<DType>::contains(const BasicArrayObject<DType> &fv)
{

    return (std::find_if(data.begin(), data.end(), [&fv](BasicArrayObject<DType> &datasetFv){
                return datasetFv.getOID() == fv.getOID();
            }) != data.end());

}



template <class DType>
void Dataset<DType>::erase(uint_fast32_t position)
{

    data.erase(data.begin() + position);
    data.shrink_to_fit();
    setCardinality(data.size());

}



template <class DType>
void Dataset<DType>::erase(const BasicArrayObject<DType> &fv)
{

    for (size_t i = 0; i < getCardinality(); i++)
    {

        if (fv.getOID() == data[i].getOID())
        {

            erase(data.begin() + i);
            break;

        }

    }

}



template <class DType>
void Dataset<DType>::resize(size_t n)
{

    data.resize(n);
    setCardinality(n);

}



#endif // DATASET_H
