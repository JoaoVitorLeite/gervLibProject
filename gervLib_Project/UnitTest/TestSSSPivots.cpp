#include "TestSSSPivots.h"


void TestSSSPivots::initTestCase()
{

}


void TestSSSPivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_sss_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_sss_2.dat");
    remove("../../gervLib/datasets/pivot_unit_test_sss_3.dat");

}


void TestSSSPivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(517);
    pvt->generatePivots(sample, df,2, 0.35);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 2;
    pvtIndex[1] = 7;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(703);
    pvt->generatePivots(sample, df,5,0.2);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 2;
    pvtIndex[1] = 5;
    pvtIndex[2] = 6;
    pvtIndex[3] = 7;
    pvtIndex[4] = 9;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(1853);
    pvt->generatePivots(sample, df,7,0.1);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 1;
    pvtIndex[1] = 2;
    pvtIndex[2] = 3;
    pvtIndex[3] = 4;
    pvtIndex[4] = 6;
    pvtIndex[5] = 7;
    pvtIndex[6] = 9;


    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(103);
    pvt->generatePivots(sample, df,2,0.35);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 1;
    pvtIndex[1] = 5;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(8);
    pvt->generatePivots(sample, df,5,0.35);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 4;
    pvtIndex[1] = 5;
    pvtIndex[2] = 7;
    pvtIndex[3] = 9;
    pvtIndex[4] = 10;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(78);
    pvt->generatePivots(sample, df,7,0.25);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 1;
    pvtIndex[1] = 2;
    pvtIndex[2] = 3;
    pvtIndex[3] = 4;
    pvtIndex[4] = 5;
    pvtIndex[5] = 8;
    pvtIndex[6] = 10;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(876);
    pvt->generatePivots(sample, df,2,0.35);

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


void TestSSSPivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(909);
    pvt->generatePivots(sample, df,5,0.35);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 3;
    pvtIndex[1] = 4;
    pvtIndex[2] = 6;
    pvtIndex[3] = 9;
    pvtIndex[4] = 12;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(99);
    pvt->generatePivots(sample, df,7,0.35);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 3;
    pvtIndex[1] = 4;
    pvtIndex[2] = 6;
    pvtIndex[3] = 8;
    pvtIndex[4] = 9;
    pvtIndex[5] = 14;
    pvtIndex[6] = 16;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSSSPivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(1105);
    pvt->generatePivots(sample,df,2,0.35);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_sss_1.dat");

    SSSPivots<double>* test = new SSSPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_sss_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestSSSPivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(36);
    pvt->generatePivots(sample,df,5,0.35);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_sss_2.dat");

    SSSPivots<double>* test = new SSSPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_sss_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestSSSPivots::test12()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SSSPivots<double>* pvt = new SSSPivots<double>();
    pvt->setSeed(915);
    pvt->generatePivots(sample,df,7,0.35);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_sss_3.dat");

    SSSPivots<double>* test = new SSSPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_sss_3.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}
