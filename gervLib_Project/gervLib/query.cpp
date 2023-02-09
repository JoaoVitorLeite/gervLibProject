#include <iostream>
#include <Hermes.h>
#include <Pivots.h>
#include <Dataset.h>
#include <OmniKdTree.h>
#include <VpTree.h>
#include <PM_Tree.h>
#include <chrono>

using namespace std;

/*Keywords

    -INDEX => Nome do index utilizado
    -DATASET_TRAIN => Caminho do dataset de treino
    ** -DATASET_TRAIN_CARDINALITY => Cardinalidade do dataset de treino
    ** -DATASET_TRAIN_DIMENSIONALITY => Dimensionalidade do dataset de treino
    ** -DATASET_TRAIN_SEPARATOR => Caracter separador do dataset de treino !DEFAULT = COMMA
    -DATASET_TEST => Caminho de um dataset de teste
    ** -DATASET_TEST_CARDINALITY => Cardinalidade do dataset de teste
    ** -DATASET_TEST_DIMENSIONALITY => Dimensionalidade do dataset de teste
    ** -DATASET_TEST_SEPARATOR => Caracter separador do dataset de teste !DEFAULT = COMMA
    -DISTANCE_FUNCTION => Função de distância !DEFAULT = EUCLIDEAN DISTANCE
    -PIVOT_TYPE => Pivot
    ** -SAMPLE_SIZE_PIVOT => Tamanho da amostra para gerar os pivôs
    ** -NUM_PIVOTS => Número de pivôs !DEFAULT = Calcular dimensionalidade intrínseca([MIN = 2] <= NUM_PIVOTS <= [MAX = Dimensionality])
    ** -PIVOT_OPTIONAL => Parâmetro opcional(o selection e kmedoids não vão precisar, mas o SSS nem sempre funciona com parâmetro adicional fixo, veja nos testes unitários...)
    -SEED = Seed
    -BF => Branch factor
    -K_MAX => Valor máximo para os k-vizinhos mais próximos !DEFAULT = 50
    -NUM_QUERY => Quantidade de consulta a serem realizadas usando o dataset de treino
    -NUM_PER_LEAF => Quantidade máxima de elementos por nó folha
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
        Dataset<double>* test = new Dataset<double>();
        DistanceFunction<BasicArrayObject<double>>* df = nullptr;
        Pivot<double>* pvt = nullptr;

        //CONFERIR SE NN EXISTE ALGUM SWITCH POR AI SOBRE OS METODOS DE SELECAO DE PIVOS
        std::string names[] = {"RANDOM", "GNAT", "CONVEX", "KMEDOIDS", "MAXSEPARETED", "MAXVARIANCE", "SELECTION", "PCA", "SSS", "FFT", "HFI", "IS", "WDR"};

        std::string *index = nullptr,
                *dataset_train = nullptr,
                *dataset_train_separator = nullptr,
                *dataset_test = nullptr,
                *dataset_test_separator = nullptr,
                *distanceFunction = nullptr,
                *pivot_type = nullptr,
                *pivot_optional = nullptr,
                *path_save_results = nullptr;

        size_t *dataset_train_cardinality = nullptr,
                *dataset_train_dimensionality = nullptr,
                *dataset_test_cardinality = nullptr,
                *dataset_test_dimensionality = nullptr,
                *num_pivots = nullptr,
                *seed = nullptr,
                *k_max = nullptr,
                *num_query = nullptr,
                *num_per_leaf = nullptr,
                *bf = nullptr;

        double *sample_size_pivot = nullptr;


        for(int x = 1; x < argc; x += 2)
        {

            std::string key = argv[x];
            for(size_t x = 0; x < key.size(); x++)
                key[x] = std::toupper(key[x]);

            std::string value = argv[x+1];

            if(key == "-INDEX")
            {

                index = new std::string(value);

            }
            else if(key == "-DATASET_TRAIN")
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
            else if(key == "-DATASET_TEST")
            {

                dataset_test = new std::string(value);

            }
            else if(key == "-DATASET_TEST_CARDINALITY")
            {

                dataset_test_cardinality = new size_t(std::stoi(value));

            }
            else if(key == "-DATASET_TEST_DIMENSIONALITY")
            {

                dataset_test_dimensionality = new size_t(std::stoi(value));

            }
            else if(key == "-DATASET_TEST_SEPARATOR")
            {

                dataset_test_separator = new std::string(value);

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
            else if(key == "-NUM_PIVOTS")
            {

                num_pivots = new size_t(std::stoi(value));

            }
            else if(key == "-PIVOT_OPTIONAL")
            {

                pivot_optional = new std::string(value);

            }
            else if(key == "-SEED")
            {

                seed = new size_t(std::stoi(value));

            }
            else if(key == "-BF")
            {

                bf = new size_t(std::stoi(value));

            }
            else if(key == "-K_MAX")
            {

                k_max = new size_t(std::stoi(value));

            }
            else if(key == "-NUM_QUERY")
            {

                num_query = new size_t(std::stoi(value));

            }
            else if(key == "-PATH_SAVE_RESULTS")
            {

                path_save_results = new std::string(value);

            }
            else if(key == "-NUM_PER_LEAF")
            {

                num_per_leaf = new size_t(std::stoi(value));

            }
            else
            {

                throw std::invalid_argument("Invalid key !_!");

            }

        }

        if((index == nullptr) || (dataset_train == nullptr) || (dataset_test == nullptr) || (pivot_type == nullptr))
        {

            throw std::invalid_argument("Important arguments were not passed !_!");

        }


        //Read seed
        if(seed == nullptr)
        {

            seed = new size_t(100);

        }

        //Path save results
        if(path_save_results == nullptr)
        {

            path_save_results = new std::string("results/");

        }

        //Read K
        if(k_max == nullptr)
        {

            k_max = new size_t(50);

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



        //Read dataset test
        if((dataset_test_cardinality != nullptr) && (dataset_test_dimensionality != nullptr))
        {

            Dataset<double>::loadNumericDataset(test, *dataset_test, *dataset_test_cardinality, *dataset_test_dimensionality);

        }
        else
        {

            if(dataset_test_separator != nullptr)
            {

                Dataset<double>::loadNumericDataset(test, *dataset_test, *dataset_test_separator);

            }
            else
            {

                Dataset<double>::loadNumericDataset(test, *dataset_test, ",");

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
            else if(*pivot_type == "MAXSEPARETED")
            {
                pvt = new MaxSeparetedPivots<double>();
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
            else
            {
                throw std::invalid_argument("Pivot selection not find !_! " + *pivot_type);
            }

            if(*sample_size_pivot != 1.0)
                pvt->setSampleSize(*sample_size_pivot);
            pvt->setSeed(*seed);
            pvt->generatePivots(train, df, *num_pivots);

        }

        //Read index
        if(*index == "OMNIKDTREE")
        {

            if(num_per_leaf == nullptr)
            {

                *num_per_leaf = (size_t)std::ceil(0.1*train->getCardinality());

            }

            OmniKdTree<double>* index_query = new OmniKdTree<double>(train, df, pvt, *num_per_leaf);

            std::string fileName = *path_save_results + names[pvt->getPivotType()] + ".csv";
            std::ofstream file(fileName);

            file << *index
                 << ","
                 << (distanceFunction == nullptr ? "EUCLIDEANDISTANCE" : *distanceFunction)
                 << ","
                 << *pivot_type
                 << ","
                 << pvt->getNumberOfPivots()
                 << ","
                 << (num_query == nullptr ? test->getCardinality() : *num_query)
                 << ","
                 << *seed
                 << "\n";

            file << "k,time,count,disk\n";

            std::vector<PairResult> ans;

            for(size_t k = 5; k <= *k_max; k += 5)
            {

                size_t max_query = test->getCardinality();

                if(num_query != nullptr)
                {

                    max_query = std::min(*num_query, test->getCardinality());

                }


                for(size_t x = 0; x < max_query; x++)
                {

                    df->resetStatistics();
                    ans.clear();
                    auto start = std::chrono::steady_clock::now();
                    index_query->kNN(train, test->getInstance(x), k, ans);
                    auto end = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

                    file << k
                         << ","
                         << elapsed.count()
                         << ","
                         << df->getDistanceCount()
                         << ","
                         << index_query->getDiskAccess()
                         << "\n";

                }

            }

            file.close();

        }
        else if(*index == "VPTREE")
        {

            if(num_per_leaf == nullptr)
            {

                *num_per_leaf = (size_t)std::ceil(0.1*train->getCardinality());

            }

            VpTree<double, DistanceFunction<BasicArrayObject<double>>>* index_query = new VpTree<double, DistanceFunction<BasicArrayObject<double>>>(false, 0.0, *num_per_leaf, pvt, train, df);

            std::string fileName = *path_save_results + names[pvt->getPivotType()] + ".csv";
            std::ofstream file(fileName);

            file << *index
                 << ","
                 << (distanceFunction == nullptr ? "EUCLIDEANDISTANCE" : *distanceFunction)
                 << ","
                 << *pivot_type
                 << ","
                 << pvt->getNumberOfPivots()
                 << ","
                 << (num_query == nullptr ? test->getCardinality() : *num_query)
                 << ","
                 << *seed
                 << "\n";

            file << "k,time,count,disk\n";

            Dataset<double>* ans;

            for(size_t k = 5; k <= *k_max; k += 5)
            {

                size_t max_query = test->getCardinality();

                if(num_query != nullptr)
                {

                    max_query = std::min(*num_query, test->getCardinality());

                }


                for(size_t x = 0; x < max_query; x++)
                {

                    ans = new Dataset<double>();
                    df->resetStatistics();
                    auto start = std::chrono::steady_clock::now();
                    index_query->kNNInc(test->getFeatureVector(x), k, index_query->getRoot(), ans, df);
                    auto end = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

                    file << k
                         << ","
                         << elapsed.count()
                         << ","
                         << df->getDistanceCount()
                         << ","
                         << index_query->getLeafNodeAccess()
                         << "\n";

                    ans->clear();

                }

            }

            file.close();

        }
        else if(*index == "PMTREE")
        {

            if(num_per_leaf == nullptr)
            {

                *num_per_leaf = (size_t)std::ceil(0.1*train->getCardinality());

            }

            PM_Tree<double>* index_query = new PM_Tree<double>(train, df, pvt, *num_per_leaf, *num_pivots);

            std::string fileName = *path_save_results + names[pvt->getPivotType()] + ".csv";
            std::ofstream file(fileName);

            file << *index
                 << ","
                 << (distanceFunction == nullptr ? "EUCLIDEANDISTANCE" : *distanceFunction)
                 << ","
                 << *pivot_type
                 << ","
                 << pvt->getNumberOfPivots()
                 << ","
                 << (num_query == nullptr ? test->getCardinality() : *num_query)
                 << ","
                 << *seed
                 << "\n";

            file << "k,time,count,disk\n";

            std::vector<KnnEntry<double>> ans;

            for(size_t k = 5; k <= *k_max; k += 5)
            {

                size_t max_query = test->getCardinality();

                if(num_query != nullptr)
                {

                    max_query = std::min(*num_query, test->getCardinality());

                }


                for(size_t x = 0; x < max_query; x++)
                {

                    df->resetStatistics();
                    auto start = std::chrono::steady_clock::now();
                    index_query->kNN(test->getFeatureVector(x), k, ans);
                    auto end = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

                    file << k
                         << ","
                         << elapsed.count()
                         << ","
                         << df->getDistanceCount()
                         << ","
                         << index_query->getLeafNodeAccess()
                         << "\n";

                    ans.clear();

                }

            }

            file.close();

        }
        else
        {

            std::invalid_argument("Unrecognized index !_!");

        }

    }


    return 0;

}
