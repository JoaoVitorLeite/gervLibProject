//OMNI - NUM
    Dataset<double>* train = new Dataset<double>();
    Dataset<double>::loadNumericDataset(train, "../../gervLib/datasets/cities_norm.csv", ",");
    Dataset<double>* test = new Dataset<double>();
    Dataset<double>::loadNumericDataset(test, "../../gervLib/datasets/cities_norm.csv", ",");
    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    size_t k = 100, numPerLeaf = 55;
    OmniKdTree<double> index = OmniKdTree<double>(train, df, pvt, numPerLeaf, 2);
    std::vector<PairResult> ans;

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans.clear();
        index.kNN(train, test->getInstance(x), k, ans);

        std::vector<double> dist;
        for(size_t i = 0; i < train->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *train->instance(i)));

        sort(dist.begin(), dist.end());

        for(size_t z = 0; z < k; z++)
            if(dist[z] != ans[z].distance)
                cout << "ERRO EM: " << x << endl;

    }
    
//OMNI - TEXT
    Dataset<std::vector<char>>* train = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(train, "../../gervLib/datasets/sgb-words.csv", " ");
    Dataset<std::vector<char>>* test = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(test, "../../gervLib/datasets/sgb-words.csv", " ");
    EditDistance<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    size_t k = 100, numPerLeaf = 55;
    OmniKdTree<std::vector<char>> index = OmniKdTree<std::vector<char>>(train, df, pvt, numPerLeaf, 2);
    std::vector<PairResult> ans;

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans.clear();
        index.kNN(train, test->getInstance(x), k, ans);

        std::vector<double> dist;
        for(size_t i = 0; i < train->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *train->instance(i)));

        sort(dist.begin(), dist.end());

        std::cout << x << std::endl;
        for(size_t z = 0; z < k; z++)
            cout << "ERRO EM: " << x << endl;

    }
    
//VP - NUM
    Dataset<double>* train = new Dataset<double>();
    Dataset<double>::loadNumericDataset(train, "../../gervLib/datasets/cities_norm.csv", ",");
    Dataset<double>* test = new Dataset<double>();
    Dataset<double>::loadNumericDataset(test, "../../gervLib/datasets/cities_norm.csv", ",");
    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    size_t k = 100, numPerLeaf = 55;
    VpTree<double, DistanceFunction<BasicArrayObject<double>>> index = VpTree<double, DistanceFunction<BasicArrayObject<double>>>(false, 0.0, numPerLeaf, pvt, train, df);
    Dataset<double>* ans = new Dataset<double>();

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans->clear();
        index.kNNInc(test->getFeatureVector(x), k, index.getRoot(), ans, df);

        std::vector<double> dist;
        for(size_t i = 0; i < test->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *test->instance(i)));

        sort(dist.begin(), dist.end());

        for(size_t z = 0; z < k; z++)
            if(dist[z] != df->getDistance(ans->getFeatureVector(z), test->getFeatureVector(x)))
                cout << "ERRO EM: " << x << endl;

    }
    
//VP - TEXT
    Dataset<std::vector<char>>* train = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(train, "../../gervLib/datasets/sgb-words.csv", " ");
    Dataset<std::vector<char>>* test = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(test, "../../gervLib/datasets/sgb-words.csv", " ");
    EditDistance<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    size_t k = 100, numPerLeaf = 55;
    VpTree<std::vector<char>, DistanceFunction<BasicArrayObject<std::vector<char>>>> index = VpTree<std::vector<char>, DistanceFunction<BasicArrayObject<std::vector<char>>>>(false, 0.0, numPerLeaf, pvt, train, df);
    Dataset<std::vector<char>>* ans = new Dataset<std::vector<char>>();

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans->clear();
        index.kNNInc(test->getFeatureVector(x), k, index.getRoot(), ans, df);

        std::vector<double> dist;
        for(size_t i = 0; i < test->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *test->instance(i)));

        sort(dist.begin(), dist.end());

        for(size_t z = 0; z < k; z++)
            if(dist[z] != df->getDistance(ans->getFeatureVector(z), test->getFeatureVector(x)))
                cout << "ERRO EM: " << x << endl;

    }


//MVP - NUM
    Dataset<double>* train = new Dataset<double>();
    Dataset<double>::loadNumericDataset(train, "../../gervLib/datasets/cities_norm.csv", ",");
    Dataset<double>* test = new Dataset<double>();
    Dataset<double>::loadNumericDataset(test, "../../gervLib/datasets/cities_norm.csv", ",");
    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    size_t k = 100, numPerLeaf = 55;
    MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxVariancePivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS> index
        = MVPTree<BasicArrayObject<double>, EuclideanDistance<BasicArrayObject<double>>, MaxVariancePivots<double>, Dataset<double>, BF,PL,LC,LPN,FO,NS>(df, train);
    std::vector<KnnEntryMVP<BasicArrayObject<double>>> ans;

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans.clear();
        index.knn(test->getFeatureVector(x), k, ans);
        std::vector<double> dist;
        for(size_t i = 0; i < train->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *train->instance(i)));

        sort(dist.begin(), dist.end());

        for(size_t z = 0; z < k; z++)
            if(dist[z] != ans[z].distance)
                cout << "ERRO EM: " << x << endl;

    }
    
//MVP - TEXT

//PM - NUM
    Dataset<double>* train = new Dataset<double>();
    Dataset<double>::loadNumericDataset(train, "../../gervLib/datasets/cities_norm.csv", ",");
    Dataset<double>* test = new Dataset<double>();
    Dataset<double>::loadNumericDataset(test, "../../gervLib/datasets/cities_norm.csv", ",");
    EuclideanDistance<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    size_t k = 5, numPerLeaf = 55;
    std::vector<KnnEntry<double>> ans;
    PM_Tree<double> index = PM_Tree<double>(train, df, pvt, numPerLeaf, 2);

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans.clear();
        index.kNN(test->getFeatureVector(x), k, ans);
        std::vector<double> dist;
        for(size_t i = 0; i < test->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *test->instance(i)));

        sort(dist.begin(), dist.end());

        for(size_t z = 0; z < k; z++)
            if(dist[z] != ans[z].distance)
                cout << "ERRO EM: " << x << endl;

    }
    
//PM - TEXT
    Dataset<std::vector<char>>* train = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(train, "../../gervLib/datasets/sgb-words.csv", " ");
    Dataset<std::vector<char>>* test = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(test, "../../gervLib/datasets/sgb-words.csv", " ");
    EditDistance<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    size_t k = 5, numPerLeaf = 55;
    std::vector<KnnEntry<std::vector<char>>> ans;
    PM_Tree<std::vector<char>> index = PM_Tree<std::vector<char>>(train, df, pvt, numPerLeaf, 2);

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans.clear();
        index.kNN(test->getFeatureVector(x), k, ans);
        std::vector<double> dist;
        for(size_t i = 0; i < test->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *test->instance(i)));

        sort(dist.begin(), dist.end());

        for(size_t z = 0; z < k; z++)
            if(dist[z] != ans[z].distance)
                cout << "ERRO EM: " << x << endl;

    }


//SPB - NUM
    Dataset<std::vector<char>>* train = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(train, "../../gervLib/datasets/names.csv", " ");
    Dataset<std::vector<char>>* test = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(test, "../../gervLib/datasets/names.csv", " ");
    EditDistance<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    size_t k = 100, numPerLeaf = 55;
    SPBTree<std::vector<char>> index = SPBTree<std::vector<char>>(train, df, pvt, 2, 2);

    std::vector<KnnSPB<std::vector<char>>> ans;

    for(size_t x = 0; x < test->getCardinality(); x++)
    {

        ans.clear();
        index.knn(test->getFeatureVector(x), k, ans);
        std::vector<double> dist;
        for(size_t i = 0; i < test->getCardinality(); i++)
            dist.push_back(df->getDistance(*test->instance(x), *test->instance(i)));

        sort(dist.begin(), dist.end());

        for(size_t z = 0; z < k; z++)
            if(dist[z] != ans[z].distance)
                cout << "ERRO EM: " << x << endl;

    }














    
    
