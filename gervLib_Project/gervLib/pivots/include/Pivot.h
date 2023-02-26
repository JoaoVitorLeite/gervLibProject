#ifndef PIVOT_H
#define PIVOT_H

#include <Dataset.h>
#include <limits.h>
#include <Util.h>
#include <Hermes.h>
#include <chrono>

enum PIVOT_TYPE{RANDOM, GNAT, CONVEX, KMEDOIDS, MAXSEPARETED, MAXVARIANCE, SELECTION, PCA, SSS, FFT, HFI, IS, WDR};

template <class DType>
class Pivot{

    protected:
        size_t currentPivot;
        size_t seed;
        size_t nPivots;
        PIVOT_TYPE pivotType;
        Dataset<DType>* pivots;
        long long elapsedTime;
        std::string path;

    private:
        void clearPivots();

        //Public functions
    public:
        static double sample_size;
        //Constructors and destructors
        Pivot();
        virtual ~Pivot();

        //Public metods
        BasicArrayObject<DType>* getNextPivot();

        //Setters
        void setNumberOfPivots(const size_t &value);
        void setSeed(const size_t &seed);
        void setPivot(BasicArrayObject<DType> *value, const size_t &pos);
        void setPivots(std::vector<BasicArrayObject<DType>> value);
        void setPivotType(PIVOT_TYPE pvtType);
        void setSampleSize(double sample_size_);
        void setElapsedTime(long long _elapsedTime);
        void setPath(std::string _path);

        //Getters
        BasicArrayObject<DType> get(const size_t pos){
            return pivots->getFeatureVector(pos);
        }

        size_t getNumberOfPivots() const;
        size_t getSeed() const;
        BasicArrayObject<DType> *getPivot(const size_t pos);
//        BasicArrayObject<DType> getPivot2(const size_t pos);
        //            const BasicArrayObject<DType> &getPivot(const size_t pos) const;
        //            std::vector<BasicArrayObject<DType>> &getPivots() const;
        PIVOT_TYPE getPivotType() const;
        size_t getSerializedSize();
        double getSampleSize(){ return sample_size; }
        long long getElapsedTime(){ return elapsedTime; }
        std::string getPath(){ return path; }

        //Métodos virtuais públicos
        virtual void generatePivots(Dataset<DType> *sample, DistanceFunction<BasicArrayObject<DType>> *df, size_t nPivots, std::vector<std::string> args = std::vector<std::string>()) = 0;

        bool isEqual(Pivot* pivot);

        unsigned char* serialize();
        void unserialize(unsigned char* dataIn);

        void saveToFile(std::string fileName);
        void loadFromFile(std::string fileName);
        void writePivotsToFile();
        void writePivotsIdsToFile();

        void printPivots();

};

//------------------------------------------------------------------------------------------------------------

template <class DType>
double Pivot<DType>::sample_size = -1.0;

template <class DType>
void Pivot<DType>::clearPivots()
{

    pivots->clear();

}



template <class DType>
Pivot<DType>::Pivot()
{

    this->currentPivot = 0;
    this->nPivots = 0;
    this->pivotType = PIVOT_TYPE::RANDOM;
    this->pivots = new Dataset<DType>();
    this->seed = 0;

}



template <class DType>
Pivot<DType>::~Pivot()
{

    clearPivots();

}



template <class DType>
void Pivot<DType>::setNumberOfPivots(const size_t &value)
{

    this->nPivots = value;
    this->pivots = new Dataset<DType>();
    this->pivots->setData(std::vector<BasicArrayObject<DType>>(nPivots));
    pivots->setCardinality(nPivots);

}



template <class DType>
void Pivot<DType>::setSeed(const size_t &seed)
{

    this->seed = seed;

}



template <class DType>
void Pivot<DType>::setPivot(BasicArrayObject<DType> *value, const size_t &pos)
{

    if((pos < getNumberOfPivots()) & (pos >= 0))
    {

        pivots->setInstance(value, pos);

    }
    else
    {

        throw std::invalid_argument("Out of bounds pivot position !_!");

    }

    if(pivots->getDimensionality() != value->getSize())
    {

        pivots->setDimensionality(value->getSize());

    }

}



template <class DType>
void Pivot<DType>::setPivots(std::vector<BasicArrayObject<DType>> value)
{

    pivots->setCardinality(value.size());
    pivots->setDimensionality(value[0].size());
    pivots = new Dataset<DType>(value, pivots->getCardinality(), pivots->getDimensionality());

}



template <class DType>
void Pivot<DType>::setPivotType(PIVOT_TYPE pvtType)
{

    this->pivotType = pvtType;

}

template <class DType>
void Pivot<DType>::setSampleSize(double sample_size_)
{

    sample_size = sample_size_;

}

template <class DType>
void Pivot<DType>::setElapsedTime(long long _elapsedTime)
{

    elapsedTime = _elapsedTime;

}

template <class DType>
void Pivot<DType>::setPath(std::string _path)
{

    path = _path;

}

template <class DType>
size_t Pivot<DType>::getNumberOfPivots() const
{

    return nPivots;

}



template <class DType>
size_t Pivot<DType>::getSeed() const
{

    return seed;

}



template <class DType>
BasicArrayObject<DType>* Pivot<DType>::getPivot(const size_t pos)
{

    BasicArrayObject<DType>* ans = nullptr;

    if((pos < getNumberOfPivots()) & (pos >= 0))
    {

        ans = pivots->instance(pos);

    }
    else
    {

        throw std::invalid_argument("Out of bounds pivot position !_!");

    }

    return ans;

}

//template <class DType>
//BasicArrayObject<DType> Pivot<DType>::getPivot2(const size_t pos)
//{

//    BasicArrayObject<DType> ans = BasicArrayObject<DType>();

//    if((pos < getNumberOfPivots()) & (pos >= 0))
//    {

//        ans = pivots[pos];

//    }
//    else
//    {

//        throw std::invalid_argument("Out of bounds pivot position !_!");

//    }

//    return ans;

//}



//template <class DType>
//const BasicArrayObject<DType>& Pivot<DType>::getPivot(const size_t pos) const
//{

//    if((pos < getNumberOfPivots()) & (pos >= 0))
//    {

//        return pivots[pos];

//    }
//    else
//    {

//        throw std::invalid_argument("Out of bounds pivot position !_!");

//    }

//}



//template <class DType>
//std::vector<BasicArrayObject<DType>>& Pivot<DType>::getPivots() const
//{

//    return &this->pivots.getElements();

//}



template <class DType>
PIVOT_TYPE Pivot<DType>::getPivotType() const
{

    return this->pivotType;

}



template <class DType>
size_t Pivot<DType>::getSerializedSize()
{

    return (4*sizeof(size_t) + sizeof(PIVOT_TYPE) + pivots->getSerializedSize());

}



template <class DType>
bool Pivot<DType>::isEqual(Pivot<DType>* pivot)
{

    bool answer = true;

    if((getNumberOfPivots() != pivot->getNumberOfPivots()) || (getPivotType() != pivot->getPivotType()))
        answer = false;
    else
    {

        for(size_t x = 0; x < getNumberOfPivots(); x++)
        {

            if(!getPivot(x)->isEqual(pivot->getPivot(x)))
            {

                answer = false;
                break;

            }


        }

    }

    return answer;

}



template <class DType>
unsigned char* Pivot<DType>::serialize()
{

    /*
    nPivots
    seed
    Dimensionality
    InstanceSize
    PivotType
    Instance**
    */

//    char* serialized = new char[getSerializedSize()];
//    size_t aux = getPivot(0)->getSize();
//    size_t size = getNumberOfPivots();

//    memcpy(serialized, &nPivots, sizeof(size_t));
//    memcpy(serialized + sizeof(size_t), &seed, sizeof(size_t));
//    memcpy(serialized + sizeof(size_t)*2, &aux, sizeof(size_t));
//    memcpy(serialized + sizeof(size_t)*3, &instanceSize, sizeof(size_t));
//    memcpy(serialized + sizeof(size_t)*4, &pivotType, sizeof(PIVOT_TYPE));

//    for(size_t x = 0; x < size; x++)
//    {

//        memcpy(serialized + sizeof(size_t)*4 + sizeof(PIVOT_TYPE) + instanceSize*x, getPivot(x)->serialize(), instanceSize);

//    }

//    return serialized;

    size_t datasetSize = pivots->getSerializedSize(), total = 0;
    unsigned char* serialized = new unsigned char[getSerializedSize()];


    memcpy(serialized + total, &currentPivot, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(serialized + total, &seed, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(serialized + total, &nPivots, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(serialized + total, &pivotType, sizeof(PIVOT_TYPE));
    total += sizeof(PIVOT_TYPE);
    memcpy(serialized + total, &datasetSize, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(serialized + total, pivots->serialize(), datasetSize);
    total += datasetSize;

    return serialized;

}



template <class DType>
void Pivot<DType>::unserialize(unsigned char* dataIn)
{

    size_t datasetSize, total = 0;

    memcpy(&currentPivot, dataIn + total, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(&seed, dataIn + total, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(&nPivots, dataIn + total, sizeof(size_t));
    total += sizeof(size_t);
    memcpy(&pivotType, dataIn + total, sizeof(PIVOT_TYPE));
    total += sizeof(PIVOT_TYPE);
    memcpy(&datasetSize, dataIn + total, sizeof(size_t));
    total += sizeof(size_t);

    unsigned char* aux = new unsigned char[datasetSize];

    memcpy(aux, dataIn + total, datasetSize);
    total += datasetSize;

    pivots = new Dataset<DType>();
    pivots->unserialize(aux);

}



template <class DType>
void Pivot<DType>::saveToFile(std::string fileName)
{

    std::ofstream file(fileName, std::ios::out | std::ios::binary);

    size_t sizeData = getSerializedSize();
    unsigned char* dataSize = new unsigned char[sizeof(size_t)];

    memcpy(dataSize, &sizeData, sizeof(size_t));

    file.write((char *)dataSize, sizeof(size_t));

    unsigned char* dataOut = new unsigned char[sizeData];

    memcpy(dataOut, serialize(), sizeData);

    file.write((char* )dataOut, sizeData);

    file.close();

}



template <class DType>
void Pivot<DType>::loadFromFile(std::string fileName)
{

    std::ifstream file(fileName, std::ios::in | std::ios::binary);

    size_t size = 0;
    unsigned char* dataSize = new unsigned char[sizeof(size_t)];

    file.read((char *)dataSize, sizeof(size_t));

    memcpy(&size, dataSize, sizeof(size_t));

    unsigned char* dataIn = new unsigned char[size];

    file.read((char *)dataIn, size);

    unserialize(dataIn);

    delete [] dataSize;
    delete [] dataIn;

}

template <class DType>
void Pivot<DType>::writePivotsToFile()
{

    std::fstream file;
    file.open(path, std::fstream::in | std::fstream::out | std::fstream::app);
    std::string aux;

    for(size_t x = 0; x < pivots->getCardinality(); x++)
    {

        aux = pivots->getFeatureVector(x).toString();
        file << aux.substr(1, aux.size() - 2) << "\n";

    }

//    if(file)
//    {

//        for(size_t x = 0; x < pivots->getCardinality(); x++)
//        {

//            aux = pivots->getFeatureVector(x).toString();
//            file << aux.substr(1, aux.size() - 2) << "\n";

//        }

//    }
//    else
//        throw std::runtime_error("Could not open the file !_!");

    file.close();

}

template <class DType>
void Pivot<DType>::writePivotsIdsToFile()
{

    std::fstream file;
    file.open(path, std::fstream::in | std::fstream::out | std::fstream::app);

    for(size_t x = 0; x < pivots->getCardinality(); x++)
    {

        file << pivots->getFeatureVector(x).getOID() << "\n";

    }

//    if(file)
//    {

//        for(size_t x = 0; x < pivots->getCardinality(); x++)
//        {

//            file << pivots->getFeatureVector(x).getOID() << "\n";

//        }

//    }
//    else
//        throw std::runtime_error("Could not open the file !_!");

    file.close();

}

template <class DType>
void Pivot<DType>::printPivots()
{

    for(size_t x = 0; x < this->getNumberOfPivots(); x++)
    {

        std::cout << this->getPivot(x)->toStringWithOID(",") << std::endl;

    }

}



#endif // PIVOT_H
