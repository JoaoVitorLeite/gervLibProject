#ifndef EQUIDEPTH_H
#define EQUIDEPTH_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <Dataset.h>
#include <fstream>
#include <config_spb.h>

using std::numeric_limits, std::vector, std::sort, std::cout, std::endl, std::ofstream, std::ifstream, std::pair, std::make_pair;

class Bins
{
public:
    long id;
    double min, max;

    Bins()
    {

        id = 0;
        min = numeric_limits<double>::max();
        max = numeric_limits<double>::lowest();

    }

    Bins(long id_)
    {

        id = id_;
        min = numeric_limits<double>::max();
        max = numeric_limits<double>::lowest();

    }

    Bins(long id_, double min_, double max_)
    {

        id = id_;
        min = min_;
        max = max_;

    }

    ~Bins(){

    }


    bool isInterval(double test)
    {

        return ((test >= min) && (test < max));

    }

};

template <class T>
class EquiDepth
{
private:
    long num_bins;
    size_t pivot_num;
    bool toINF;
    vector<vector<Bins>> bins;

public:
    EquiDepth()
    {

        num_bins = -1;
        pivot_num = -1;
        toINF = true;
        bins = vector<vector<Bins>>();

    }

    EquiDepth(long num_bins_, size_t pivot_num_, bool toINF_ = true)
    {

        num_bins = num_bins_;
        pivot_num = pivot_num_;
        toINF = toINF_;
        bins = std::vector<std::vector<Bins>>(pivot_num, std::vector<Bins>(num_bins));

    }

    void clear()
    {

        for(auto v : bins)
        {

            v.clear();

        }
        bins.clear();

    }

    ~EquiDepth()
    {

        clear();

    }

    void build(vector<vector<double>>& aux)
    {

        std::vector<double> v(aux.size());

        for(size_t i = 0; i < pivot_num; i++)
        {

            for(size_t k = 0; k < aux.size(); k++)
                v[k] = aux[k][i];

            std::sort(v.begin(), v.end());
            size_t index = 0;
            size_t size = v.size()/num_bins;
            size_t sz = 0;
            size_t j = 0;

            //cout << size << endl;
            double last_value = -1.0;

            for(j = 0; j < v.size(); j++)
            {

                if(sz >= size && v[j] != last_value)
                {

                    if(index != (num_bins - 1))
                    {

                        bins[i][index+1].min = bins[i][index].max;
                        bins[i][index].max = std::nextafter(bins[i][index].max, std::numeric_limits<double>::lowest());
                        index++;

                    }

                    sz = 0;

                }

                bins[i][index].min = std::min(bins[i][index].min, v[j]);
                bins[i][index].max = std::max(bins[i][index].max, v[j]);
                sz++;
                last_value = v[j];

            }

            if(toINF)
            {

                bins[i][num_bins-1].max = std::numeric_limits<double>::max();

            }


        }

    }

    void print()
    {

        cout << "PIVOT NUM: " << bins.size() << endl;
        cout << "NUM BINS: " << bins[0].size() << endl << endl;

        for(size_t i = 0; i < bins.size(); i++)
        {

            for(size_t j = 0; j < bins[0].size(); j++)
            {

                cout << "[" << bins[i][j].min << ", " << bins[i][j].max << ")" << endl;

            }

            cout << "\n\n";

        }

    }

    void saveToFile()
    {

        ofstream file(baseFilePath + std::filesystem::path::preferred_separator + "equi_depth.txt");
        file << bins.size() << " "  << bins[0].size() << endl;

        for(size_t i = 0; i < bins.size(); i++)
        {

            for(size_t j = 0; j < bins[0].size(); j++)
            {

                file << bins[i][j].min << " " << bins[i][j].max << endl;

            }

        }

        file.close();

    }

    void readFromFile()
    {

        ifstream file(baseFilePath + std::filesystem::path::preferred_separator + "equi_depth.txt");
        size_t pivot_num;
        double minV, maxV;

        file >> pivot_num;
        file >> num_bins;

        bins.resize(pivot_num, vector<Bins>(num_bins));

        for(size_t i = 0; i < pivot_num; i++)
        {

            for(long j = 0; j < num_bins; j++)
            {

                file >> minV;
                file >> maxV;

                bins[i][j] = Bins(j, minV, maxV);

            }

        }

        file.close();

    }

    long getNumberOfBins()
    {

        return num_bins;

    }

    long getBin(size_t pivot, double value)
    {

        for(long i = 0; i < num_bins; i++)
        {

            if(bins[pivot][i].isInterval(value))
            {

                return i;

            }

        }

        throw std::runtime_error("Bins definition error !_!");

    }

    pair<double, double> getInterval(size_t pivot, long id)
    {

        return make_pair(bins[pivot][id].min, bins[pivot][id].max);

    }

};


#endif // EQUIDEPTH_H
