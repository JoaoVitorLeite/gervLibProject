#include <PivotExperiments.h>

using namespace std;

typedef std::vector<char> str;

/*Keywords

    -DATASET_TRAIN => Caminho do dataset de treino
    ** -DATASET_TRAIN_CARDINALITY => Cardinalidade do dataset de treino
    ** -DATASET_TRAIN_DIMENSIONALITY => Dimensionalidade do dataset de treino
    ** -DATASET_TRAIN_SEPARATOR => Caracter separador do dataset de treino !DEFAULT = COMMA
    -DISTANCE_FUNCTION => Função de distância !DEFAULT = EUCLIDEAN DISTANCE
    -PIVOT_TYPE => Pivot
    ** -SAMPLE_SIZE_PIVOT => Tamanho da amostra para gerar os pivôs
    ** -SEED => Seed
    -REP => Numero de repeticoes !DEFAULT = 1
    -PATH_SAVE_RESULTS => Caminho para salvar os arquivos gerados !DEFAULT = ../results/

*/

int main(int argc, char *argv[])
{

    if((argc-1) % 2 != 0)
    {

        throw std::invalid_argument("Invalid number of arguments !_!");

    }
    else
    {

        Dataset<double>* train = new Dataset<double>();
        DistanceFunction<BasicArrayObject<double>>* df = nullptr;
        Pivot<double>* pvt = nullptr;
        std::string names[] = {"RANDOM", "GNAT", "CONVEX", "KMEDOIDS", "MAXSEPARATED", "MAXVARIANCE", "SELECTION", "PCA", "SSS", "FFT", "HFI", "IS", "WDR", "BPP"};

        std::string *dataset_train = nullptr,
            *dataset_train_separator = nullptr,
                *distanceFunction = nullptr,
                    *pivot_type = nullptr,
                        *pivot_optional = nullptr,
                            *path_save_results = nullptr;

        size_t *dataset_train_cardinality = nullptr,
            *dataset_train_dimensionality = nullptr,
                *seed = nullptr,
                    *num = nullptr,
                        *rep = nullptr;

        double *sample_size_pivot = nullptr;

        for(int x = 1; x < argc; x += 2)
        {

            std::string key = argv[x];
            for(size_t x = 0; x < key.size(); x++)
                key[x] = std::toupper(key[x]);

            std::string value = argv[x+1];

            if(key == "-DATASET_TRAIN")
            {

                dataset_train = new std::string(value);

            }
            else if(key == "-DATASET_TRAIN_CARDINALITY")
            {

                dataset_train_cardinality = new size_t(std::stoi(value));

            }
            else if(key == "-DATASET_TRAIN_DIMENSIONALITY")
            {

                dataset_train_dimensionality = new size_t(std::stoi(value));

            }
            else if(key == "-DATASET_TRAIN_SEPARATOR")
            {

                dataset_train_separator = new std::string(value);

            }
            else if(key == "-DISTANCE_FUNCTION")
            {

                distanceFunction = new std::string(value);

            }
            else if(key == "-PIVOT_TYPE")
            {

                pivot_type = new std::string(value);

            }
            else if(key == "-SAMPLE_SIZE_PIVOT")
            {

                sample_size_pivot = new double(std::stod(value));

            }
            else if(key == "-SEED")
            {

                seed = new size_t(std::stoi(value));

            }
            else if(key == "-REP")
            {

                rep = new size_t(std::stoi(value));

            }
            else if(key == "-PATH_SAVE_RESULTS")
            {

                path_save_results = new std::string(value);

            }
            else
            {

                throw std::invalid_argument("Invalid key !_!");

            }

        }

        //Read seed
        if(seed == nullptr)
        {

            seed = new size_t(100);

        }

        //Path save results
        if(path_save_results == nullptr)
        {

            path_save_results = new std::string("../results/");

        }

        //Read dataset train
        if((dataset_train_cardinality != nullptr) && (dataset_train_dimensionality != nullptr))
        {

            Dataset<double>::loadNumericDataset(train, *dataset_train, *dataset_train_cardinality, *dataset_train_dimensionality);

        }
        else
        {

            if(dataset_train_separator != nullptr)
            {

                Dataset<double>::loadNumericDataset(train, *dataset_train, *dataset_train_separator);

            }
            else
            {

                Dataset<double>::loadNumericDataset(train, *dataset_train, ",");

            }

        }

        //Read distance function
        if(distanceFunction == nullptr)
        {

            df = new EuclideanDistance<BasicArrayObject<double>>();

        }
        else
        {

            if(*distanceFunction == "EUCLIDEAN")
                df = new EuclideanDistance<BasicArrayObject<double>>();
            else if(*distanceFunction == "CHEBYSHEV")
                df = new ChebyshevDistance<BasicArrayObject<double>>();
            else if(*distanceFunction == "MANHATTAN")
                df = new ManhattanDistance<BasicArrayObject<double>>();
            else
            {

                throw std::invalid_argument("Not find distance function !_!");

            }

        }

        //Read pivot
        if(pivot_type != nullptr)
        {

            if(sample_size_pivot == nullptr)
            {

                sample_size_pivot = new double(0.1);

            }

            if(*pivot_type == "RANDOM")
            {
                pvt = new RandomPivots<double>();
            }
            else if(*pivot_type == "GNAT")
            {
                pvt = new GnatPivots<double>();
            }
            else if(*pivot_type == "CONVEX")
            {
                pvt = new ConvexPivots<double>();
            }
            else if(*pivot_type == "KMEDOIDS")
            {
                pvt = new KmedoidsPivots<double>();
            }
            else if(*pivot_type == "MAXSEPARATED")
            {
                pvt = new MaxSeparatedPivots<double>();
            }
            else if(*pivot_type == "MAXVARIANCE")
            {
                pvt = new MaxVariancePivots<double>();
            }
            else if(*pivot_type == "SELECTION")
            {
                pvt = new SelectionPivots<double>();
            }
            else if(*pivot_type == "PCA")
            {
                pvt = new PCAPivots<double>();
            }
            else if(*pivot_type == "SSS")
            {
                pvt = new SSSPivots<double>();
            }
            else if(*pivot_type == "FFT")
            {
                pvt = new FFTPivots<double>();
            }
            else if(*pivot_type == "HFI")
            {
                pvt = new HFIPivots<double>();
            }
            else if(*pivot_type == "IS")
            {
                pvt = new ISPivots<double>();
            }
            else if(*pivot_type == "WDR")
            {
                pvt = new WDRPivots<double>();
            }
            else if(*pivot_type == "BPP")
            {
                pvt = new BPPPivots<double>();
            }
            else
            {
                throw std::invalid_argument("Pivot selection not find !_! " + *pivot_type);
            }

            if(*sample_size_pivot != 1.0)
                pvt->setSampleSize(*sample_size_pivot);

            pvt->setSeed(*seed);

        }

        PivotExperiments<double> expt = PivotExperiments<double>();
        expt.setDistanceFunctionName(*distanceFunction);
        expt.setDistanceFunction(df);
        expt.setOutputPath(*path_save_results);
        expt.setTrainDataset(train);
        expt.setPivotMethod(pvt);
        expt.setSeed(*seed);
        expt.setSampleSize(*sample_size_pivot);
        expt.modifySeed();
        expt.modifySampleSize();
        expt.setPivotNum({5, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500});
        expt.setSavePivot(true);
        if(*rep == 1)
            expt.runExperiment();
        else
            expt.runExperimentWithRepetitions(*rep);

        //cout << *pivot_type << endl;

    }

    return 0;

}
