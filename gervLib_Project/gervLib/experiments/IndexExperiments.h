#ifndef INDEXEXPERIMENTS_H
#define INDEXEXPERIMENTS_H

#include <Experiments.h>

template <class T>
class IndexExperiments : public Experiments<T>
{

public:

    IndexExperiments()
    {

    }

    ~IndexExperiments()
    {

    }

    void runExperiment() override
    {

        std::string filePath = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator)) +
                std::filesystem::path::preferred_separator + this->indexName() + "_" + this->names[this->pvt->getPivotType()] + ".csv";
        std::ofstream file(filePath, std::fstream::app);

        if(std::filesystem::file_size(filePath) == 0)
        {

            file << this->indexName()
                 << ","
                 << this->dfName
                 << ","
                 << this->names[this->pvt->getPivotType()]
                 << ","
                 << this->pvt->getNumberOfPivots()
                 << ","
                 << this->num_query
                 << ","
                 << this->seed
                 << "\n";

            file << "k,time,count,disk\n";

        }

        for(size_t i = 0; i < this->num_query; i++)
        {

            for(size_t k = this->kRange[0]; k < this->kRange[1]; k+= this->kRange[2])
            {

                runQuery(this->test->getInstance(i), k);

                file << k
                     << ","
                     << this->getTime()
                     << ","
                     << this->getCalcDist()
                     << ","
                     << this->getDiskAccess()
                     << "\n";

            }

        }

        file.close();

    }

    void runExperimentWithRepetitions(size_t rep) override
    {

        for(size_t r = 0; r < rep; r++)
        {

            std::string filePath = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator)) +
                    std::filesystem::path::preferred_separator + this->indexName() + "_" + this->names[this->pvt->getPivotType()]
                    + "_rep_" + std::to_string(r) + ".csv";
            std::ofstream file(filePath, std::fstream::app);

            if(std::filesystem::file_size(filePath) == 0)
            {

                file << this->indexName()
                     << ","
                     << this->dfName
                     << ","
                     << this->names[this->pvt->getPivotType()]
                     << ","
                     << this->pvt->getNumberOfPivots()
                     << ","
                     << this->num_query
                     << ","
                     << this->seed
                     << "\n";

                file << "k,time,count,disk\n";

            }

            for(size_t i = 0; i < this->num_query; i++)
            {

                for(size_t k = this->kRange[0]; k < this->kRange[1]; k+= this->kRange[2])
                {

                    runQuery(this->test->getInstance(i), k);

                    file << k
                         << ","
                         << this->getTime()
                         << ","
                         << this->getCalcDist()
                         << ","
                         << this->getDiskAccess()
                         << "\n";

                }

            }

            file.close();

        }

    }

};

#endif // INDEXEXPERIMENTS_H
