#include "TestMaxVariancePivots.h"


void TestMaxVariancePivots::initTestCase()
{

}


void TestMaxVariancePivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_max_variance_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_max_variance_2.dat");
    remove("../../gervLib/datasets/pivot_unit_test_max_variance_3.dat");

}


void TestMaxVariancePivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(503);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 2;
    pvtIndex[1] = 1;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxVariancePivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(19);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 6;
    pvtIndex[1] = 1;
    pvtIndex[2] = 4;
    pvtIndex[3] = 0;
    pvtIndex[4] = 3;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);


}



void TestMaxVariancePivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(524);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 1;
    pvtIndex[1] = 13;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxVariancePivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(8);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 12;
    pvtIndex[1] = 1;
    pvtIndex[2] = 0;
    pvtIndex[3] = 2;
    pvtIndex[4] = 14;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxVariancePivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(8037);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 9;
    pvtIndex[1] = 3;
    pvtIndex[2] = 12;
    pvtIndex[3] = 13;
    pvtIndex[4] = 4;
    pvtIndex[5] = 14;
    pvtIndex[6] = 8;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxVariancePivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(31);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 2;
    pvtIndex[1] = 1;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxVariancePivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(3148);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 0;
    pvtIndex[1] = 1;
    pvtIndex[2] = 12;
    pvtIndex[3] = 3;
    pvtIndex[4] = 5;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxVariancePivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(555);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 14;
    pvtIndex[1] = 2;
    pvtIndex[2] = 6;
    pvtIndex[3] = 10;
    pvtIndex[4] = 3;
    pvtIndex[5] = 7;
    pvtIndex[6] = 1;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxVariancePivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(503);
    pvt->generatePivots(sample,df,2);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_max_variance_1.dat");

    MaxVariancePivots<double>* test = new MaxVariancePivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_max_variance_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestMaxVariancePivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(8);
    pvt->generatePivots(sample,df,5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_max_variance_2.dat");

    MaxVariancePivots<double>* test = new MaxVariancePivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_max_variance_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestMaxVariancePivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxVariancePivots<double>* pvt = new MaxVariancePivots<double>();
    pvt->setSeed(555);
    pvt->generatePivots(sample,df,7);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_max_variance_3.dat");

    MaxVariancePivots<double>* test = new MaxVariancePivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_max_variance_3.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}
