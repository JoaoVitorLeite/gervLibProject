#ifndef MEMORYMANAGERUTILS_H
#define MEMORYMANAGERUTILS_H

#include "config_spb.h"
#include <Dataset.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <experimental/filesystem>

using fs = std::experimental::filesystem;

void setBaseFilePath(std::string auxName = "SPBfiles")
{

    if(!fs::exists(baseFilePath))
    {

        throw std::runtime_error("Base file path not exists !_!");

    }

    std::vector<int> sub_folders;
    std::string auxPath, subName;
    //size_t num;
    for(const auto & entry : fs::directory_iterator(baseFilePath))
    {

        auxPath = entry.path().string();

        //std::cout << auxPath << "\n";
        if(entry.is_directory() && (entry.path().string().find(auxName) != std::string::npos))
        {

//            std::cout << "SUB = " << auxPath.substr(auxPath.find_last_of(std::filesystem::path::preferred_separator) + 1) << "\n";
            subName = auxPath.substr(auxPath.find_last_of(fs::path::preferred_separator) + 1);
            sub_folders.push_back(std::stoi(subName.substr(subName.find("_") + 1)));

        }

    }

    if(!sub_folders.empty())
    {

        sort(sub_folders.begin(), sub_folders.end());

    }

    std::stringstream ss;

    if(sub_folders.empty())
    {

        ss << baseFilePath << fs::path::preferred_separator << auxName << "_0";

    }
    else
    {

//        size_t pos = sub_folders.back().find("_");
//        std::cout << sub_folders.back().substr(pos+1) << std::endl;
//        std::stringstream ssAux(sub_folders.back().substr(pos+1));
//        ssAux >> num;
//        num++;
//        ss.str("");
        ss << baseFilePath << fs::path::preferred_separator << auxName << "_" << (++sub_folders.back());
//        ss << "SPBfiles_" << num;

    }

    std::string sub_folderPath = ss.str();

    //std::cout << sub_folderPath << "\n";

    if (!fs::is_directory(sub_folderPath) || !fs::exists(sub_folderPath))
    {

        fs::create_directory(sub_folderPath);

    }

    baseFilePath = sub_folderPath;
    sub_folders.clear();

}

std::string genDiskFileName(size_t pageID)
{

    std::stringstream ss;

    ss << baseFilePath << fs::path::preferred_separator << "tmp_" << pageID << ".dat";

    return ss.str();

}

template <class type>
void write_to_disk(Dataset<type>* dataset, std::vector<size_t> ids, size_t pageCost, size_t pageSize, size_t pageID)
{

    IOwrite++;

    std::ofstream file(genDiskFileName(pageID), std::ios::out | std::ios::binary);

    unsigned char* data = new unsigned char[pageCost + sizeof(size_t)];
    size_t total = 0, sizeElem;

    memcpy(data + total, &pageSize, sizeof(size_t));
    total += sizeof(size_t);


    for(size_t i = 0; i < ids.size(); i++)
    {

        sizeElem = dataset->getFeatureVector(ids[i]).getSerializedSize();

        memcpy(data + total, &sizeElem, sizeof(size_t));
        total += sizeof(size_t);

        memcpy(data + total, dataset->getFeatureVector(ids[i]).serialize(), sizeElem);
        total += sizeElem;

    }

    file.write((char*)data, pageCost + sizeof(size_t));
    file.close();
    delete [] data;

}

template <class type>
void write_dataset_to_disk(Dataset<type>* dataset, size_t pageID)
{

    IOwrite++;

    std::ofstream file(genDiskFileName(pageID), std::ios::out | std::ios::binary);
    unsigned char* data = new unsigned char[dataset->getSerializedSize() + sizeof(size_t)];
    size_t sz = dataset->getSerializedSize();
    memcpy(data, &sz, sizeof(size_t));
    memcpy(data + sizeof(size_t), dataset->serialize(), sz);
    file.write((char*)data, dataset->getSerializedSize() + sizeof(size_t));
    file.close();
    delete [] data;

}

template <class type>
void read_dataset_from_disk(Dataset<type>* dataset, size_t pageID)
{

    IOread++;

    std::ifstream file(genDiskFileName(pageID), std::ios::in | std::ios::binary);
    size_t size;
    unsigned char* sizeChar = new unsigned char[sizeof(size_t)];
    file.read((char*)sizeChar, sizeof(size_t));
    memcpy(&size, sizeChar, sizeof(size_t));
    unsigned char* data = new unsigned char[size];
    file.read((char*)data, size);
    dataset->unserialize(data);
    file.close();
    delete [] sizeChar;
    delete [] data;

}

template <class type>
void read_from_disk(std::vector<BasicArrayObject<type>*>& read, size_t pageID)
{

    IOread++;

    std::ifstream file(genDiskFileName(pageID), std::ios::in | std::ios::binary);

    size_t size;
//    BasicArrayObject<type> b;
    unsigned char* sizeChar = new unsigned char[sizeof(size_t)];

    file.read((char*)sizeChar, sizeof(size_t));
    //total += sizeof(size_t);
    memcpy(&size, sizeChar, sizeof(size_t));

    size_t pageSize = size;

    for(size_t i = 0; i < pageSize; i++)
    {

        file.read((char*)sizeChar, sizeof(size_t));
        //total += sizeof(size_t);
        memcpy(&size, sizeChar, sizeof(size_t));

        unsigned char* basicChar = new unsigned char[size];
        file.read((char*)basicChar, size);
        //total += size;

        BasicArrayObject<type>* b = new BasicArrayObject<type>();
        b->unserialize(basicChar);
        read.push_back(b);

        delete [] basicChar;

    }

    delete [] sizeChar;

}

void checkRAFDisk()
{

    std::string path = baseFilePath;
    std::vector<BasicArrayObject<double>*> ans;

    for (const auto & entry : fs::directory_iterator(path))
    {

        if (entry.path().extension() == ".dat")
        {

            std::string p = entry.path().stem().string();
            std::stringstream ss(p.substr(p.find("_") + 1));
            size_t pageID;
            ss >> pageID;
            read_from_disk(ans, pageID);

        }
    }

    for(auto b : ans)
        std::cout << b->toStringWithOID() << std::endl;

    std::cout << ans.size() << std::endl;

    for(auto i : ans)
        delete i;
    ans.clear();

}

#endif // MEMORYMANAGERUTILS_H
