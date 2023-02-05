#include "TestPCAPivots.h"


void TestPCAPivots::initTestCase()
{

}


void TestPCAPivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_sss_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_sss_2.dat");
    remove("../../gervLib/datasets/pivot_unit_test_sss_3.dat");

}


void TestPCAPivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 0;
    pvtIndex[1] = 8;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 0;
    pvtIndex[1] = 8;
    pvtIndex[2] = 7;
    pvtIndex[3] = 9;
    pvtIndex[4] = 6;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 0;
    pvtIndex[1] = 8;
    pvtIndex[2] = 7;
    pvtIndex[3] = 9;
    pvtIndex[4] = 6;
    pvtIndex[5] = 5;
    pvtIndex[6] = 4;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 0;
    pvtIndex[1] = 8;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 0;
    pvtIndex[1] = 8;
    pvtIndex[2] = 11;
    pvtIndex[3] = 14;
    pvtIndex[4] = 13;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 0;
    pvtIndex[1] = 8;
    pvtIndex[2] = 11;
    pvtIndex[3] = 14;
    pvtIndex[4] = 13;
    pvtIndex[5] = 12;
    pvtIndex[6] = 10;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 0;
    pvtIndex[1] = 9;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 0;
    pvtIndex[1] = 9;
    pvtIndex[2] = 8;
    pvtIndex[3] = 12;
    pvtIndex[4] = 13;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 0;
    pvtIndex[1] = 9;
    pvtIndex[2] = 8;
    pvtIndex[3] = 12;
    pvtIndex[4] = 13;
    pvtIndex[5] = 16;
    pvtIndex[6] = 18;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestPCAPivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    //pvt->setSeed(1105);
    pvt->generatePivots(sample,df,2);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_sss_1.dat");

    PCAPivots<double>* test = new PCAPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_sss_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestPCAPivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    //pvt->setSeed(36);
    pvt->generatePivots(sample,df,5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_sss_2.dat");

    PCAPivots<double>* test = new PCAPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_sss_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestPCAPivots::test12()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    PCAPivots<double>* pvt = new PCAPivots<double>();
    //pvt->setSeed(915);
    pvt->generatePivots(sample,df,7);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_sss_3.dat");

    PCAPivots<double>* test = new PCAPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_sss_3.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}
