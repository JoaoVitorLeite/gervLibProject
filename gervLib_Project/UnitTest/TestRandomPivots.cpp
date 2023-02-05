#include "TestRandomPivots.h"

void TestRandomPivots::initTestCase()
{

}


void TestRandomPivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_random_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_random_2.dat");

}


void TestRandomPivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(207);
    pvt->generatePivots(sample, df, 2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 8;
    pvtIndex[1] = 5;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(13);
    pvt->generatePivots(sample, df, 5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 2;
    pvtIndex[1] = 4;
    pvtIndex[2] = 9;
    pvtIndex[3] = 1;
    pvtIndex[4] = 0;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(10341);
    pvt->generatePivots(sample, df, 7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 0;
    pvtIndex[1] = 4;
    pvtIndex[2] = 8;
    pvtIndex[3] = 1;
    pvtIndex[4] = 2;
    pvtIndex[5] = 3;
    pvtIndex[6] = 7;


    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(35);
    pvt->generatePivots(sample, df, 2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 3;
    pvtIndex[1] = 4;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(7);
    pvt->generatePivots(sample, df, 5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 3;
    pvtIndex[1] = 5;
    pvtIndex[2] = 1;
    pvtIndex[3] = 9;
    pvtIndex[4] = 7;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(183);
    pvt->generatePivots(sample, df, 7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 1;
    pvtIndex[1] = 8;
    pvtIndex[2] = 4;
    pvtIndex[3] = 9;
    pvtIndex[4] = 3;
    pvtIndex[5] = 0;
    pvtIndex[6] = 2;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(40);
    pvt->generatePivots(sample, df, 2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 9;
    pvtIndex[1] = 6;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>(sample, 5, 1441);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 6;
    pvtIndex[1] = 9;
    pvtIndex[2] = 7;
    pvtIndex[3] = 4;
    pvtIndex[4] = 2;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(1);
    pvt->generatePivots(sample, df, 7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 1;
    pvtIndex[1] = 9;
    pvtIndex[2] = 2;
    pvtIndex[3] = 5;
    pvtIndex[4] = 7;
    pvtIndex[5] = 6;
    pvtIndex[6] = 3;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestRandomPivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(1);
    pvt->generatePivots(sample, df, 7);

    unsigned char* seri = pvt->serialize();

    RandomPivots<double>* test = new RandomPivots<double>();

    test->unserialize(seri);

    QCOMPARE(pvt->isEqual(test), true);

}


void TestRandomPivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    RandomPivots<double>* pvt = new RandomPivots<double>();
    pvt->setSeed(7);
    pvt->generatePivots(sample, df, 5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_random_1.dat");

    RandomPivots<double>* test = new RandomPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_random_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}



void TestRandomPivots::test12()
{

    Dataset<std::vector<char>>* sample = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(sample, "../../gervLib/datasets/sgb-words-30.csv", 30, 1);
    DistanceFunction<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    pvt->setSeed(1);
    pvt->generatePivots(sample, df, 7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 9;
    pvtIndex[1] = 23;
    pvtIndex[2] = 12;
    pvtIndex[3] = 25;
    pvtIndex[4] = 27;
    pvtIndex[5] = 16;
    pvtIndex[6] = 13;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}



void TestRandomPivots::test13()
{

    Dataset<std::vector<char>>* sample = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(sample, "../../gervLib/datasets/sgb-words-30.csv", 30, 1);
    DistanceFunction<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    pvt->setSeed(1);
    pvt->generatePivots(sample, df, 7);

    unsigned char* seri = pvt->serialize();

    RandomPivots<std::vector<char>>* test = new RandomPivots<std::vector<char>>();

    test->unserialize(seri);

    QCOMPARE(pvt->isEqual(test), true);

}



void TestRandomPivots::test14()
{

    Dataset<std::vector<char>>* sample = new Dataset<std::vector<char>>();
    Dataset<std::vector<char>>::loadTextDataset(sample, "../../gervLib/datasets/sgb-words-30.csv", 30, 1);
    DistanceFunction<BasicArrayObject<std::vector<char>>>* df = new EditDistance<BasicArrayObject<std::vector<char>>>();
    RandomPivots<std::vector<char>>* pvt = new RandomPivots<std::vector<char>>();
    pvt->setSeed(1);
    pvt->generatePivots(sample, df, 5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_random_2.dat");

    RandomPivots<std::vector<char>>* test = new RandomPivots<std::vector<char>>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_random_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);


}

