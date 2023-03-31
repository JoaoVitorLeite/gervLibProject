#ifndef PIVOTEXPERIMENTS_H
#define PIVOTEXPERIMENTS_H

#include <Experiments.h>

template <class T>
class PivotExperiments : public Experiments<T>
{

private:
    std::vector<size_t> pvtNum;
    bool savePivots;

public:
    PivotExperiments()
    {

    }

    ~PivotExperiments()
    {


    }

    void setPivotNum(std::vector<size_t>& pvt)
    {

        pvtNum = pvt;

    }

    void setPivotNum(std::vector<size_t> pvt)
    {

        pvtNum = pvt;

    }

    void setPivotNum(size_t start, size_t end, size_t step)
    {

        pvtNum.clear();

        for(size_t i = start; i < end; i+= step)
        {

            pvtNum.push_back(i);

        }

    }

    void buildIndex() override
    {

    }

    void setSavePivot(bool b)
    {

        savePivots = b;

    }

    void runQuery(BasicArrayObject<T>* query, size_t k) override
    {

    }

    std::string indexName() override
    {

    }

    void runExperiment() override
    {

        std::string filePath = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator)) +
                std::filesystem::path::preferred_separator + this->names[this->pvt->getPivotType()] + ".csv";
        std::ofstream file(filePath, std::fstream::app);

        if(std::filesystem::file_size(filePath) == 0)
        {

            file << this->dfName
                 << ","
                 << this->names[this->pvt->getPivotType()]
                 << ","
                 << this->sample_size
                 << ","
                 << this->seed
                 << "\n";

            file << "p,time,count\n";

        }

        for(size_t &p : pvtNum)
        {

            this->df->resetStatistics();
            auto start = std::chrono::steady_clock::now();
            this->pvt->generatePivots(this->train, this->df, p);
            auto end = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
            this->time = elapsed.count();
            this->calcDist = this->df->getDistanceCount();

            file << p
                 << ","
                 << this->getTime()
                 << ","
                 << this->getCalcDist()
                 << std::endl;

            if(savePivots)
            {

                std::string filePath2 = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator))
                        + std::filesystem::path::preferred_separator + this->names[this->pvt->getPivotType()] + "_PIVOTS_ID_" + std::to_string(p) + ".csv";
                std::string filePath3 = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator))
                        + std::filesystem::path::preferred_separator + this->names[this->pvt->getPivotType()] + "_PIVOTS_VALUE_" + std::to_string(p) + ".csv";
                std::ofstream file2(filePath2, std::fstream::app);
                std::ofstream file3(filePath3, std::fstream::app);

                for(size_t i = 0; i < p; i++)
                {

                    file2 << this->pvt->get(i).getOID() << std::endl;
                    std::string aux = this->pvt->get(i).toString(" ");
                    file3 << aux.substr(1, aux.size()-2)  << std::endl;

                }

            }

        }

    }

    void runExperimentWithRepetitions(size_t rep) override
    {

        bool repSavePivots = true;

        for(size_t r = 0; r < rep; r++)
        {

            std::string filePath = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator)) +
                    std::filesystem::path::preferred_separator + this->names[this->pvt->getPivotType()] + "_rep_"
                    + std::to_string(r) + ".csv";
            std::ofstream file(filePath, std::fstream::app);

            if(std::filesystem::file_size(filePath) == 0)
            {

                file << this->dfName
                     << ","
                     << this->names[this->pvt->getPivotType()]
                     << ","
                     << this->sample_size
                     << ","
                     << this->seed
                     << "\n";

                file << "p,time,count\n";

            }

            for(size_t &p : pvtNum)
            {

                auto start = std::chrono::steady_clock::now();
                this->pvt->generatePivots(this->train, this->df, p);
                auto end = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
                this->time = elapsed.count();
                this->calcDist = this->df->getDistanceCount();

                file << p
                     << ","
                     << this->getTime()
                     << ","
                     << this->getCalcDist()
                     << "\n";

                if(savePivots && repSavePivots)
                {

                    std::string filePath2 = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator))
                            + std::filesystem::path::preferred_separator + this->names[this->pvt->getPivotType()] + "_PIVOTS_REP_ID_" + std::to_string(p) + ".csv";
                    std::string filePath3 = this->outputPath.substr(0, this->outputPath.find_last_of(std::filesystem::path::preferred_separator))
                            + std::filesystem::path::preferred_separator + this->names[this->pvt->getPivotType()] + "_PIVOTS_REP_VALUE_" + std::to_string(p) + ".csv";
                    std::ofstream file2(filePath2, std::fstream::app);
                    std::ofstream file3(filePath3, std::fstream::app);

                    for(size_t i = 0; i < p; i++)
                    {

                        file2 << this->pvt->get(i).getOID() << std::endl;
                        std::string aux = this->pvt->get(i).toString(" ");
                        file3 << aux.substr(1, aux.size()-2)  << std::endl;

                    }

                }

            }

            repSavePivots = false;

        }

    }

};

#endif // PIVOTEXPERIMENTS_H
