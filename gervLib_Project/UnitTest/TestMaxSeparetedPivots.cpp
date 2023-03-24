#include "TestMaxSeparetedPivots.h"


void TestMaxSeparetedPivots::initTestCase()
{

}


void TestMaxSeparetedPivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_max_separeted_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_max_separeted_2.dat");
    remove("../../gervLib/datasets/pivot_unit_test_max_separeted_3.dat");

}


void TestMaxSeparetedPivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(1627);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 6;
    pvtIndex[1] = 2;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(3508);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 2;
    pvtIndex[1] = 6;
    pvtIndex[2] = 9;
    pvtIndex[3] = 8;
    pvtIndex[4] = 4;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(41);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 6;
    pvtIndex[1] = 7;
    pvtIndex[2] = 1;
    pvtIndex[3] = 8;
    pvtIndex[4] = 9;
    pvtIndex[5] = 4;
    pvtIndex[6] = 0;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(1024);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 12;
    pvtIndex[1] = 11;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(597);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 0;
    pvtIndex[1] = 9;
    pvtIndex[2] = 7;
    pvtIndex[3] = 11;
    pvtIndex[4] = 4;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(66);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 12;
    pvtIndex[1] = 11;
    pvtIndex[2] = 9;
    pvtIndex[3] = 7;
    pvtIndex[4] = 1;
    pvtIndex[5] = 4;
    pvtIndex[6] = 14;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(1211);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 18;
    pvtIndex[1] = 14;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(2480);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 14;
    pvtIndex[1] = 12;
    pvtIndex[2] = 18;
    pvtIndex[3] = 6;
    pvtIndex[4] = 7;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(12615);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 4;
    pvtIndex[1] = 14;
    pvtIndex[2] = 18;
    pvtIndex[3] = 1;
    pvtIndex[4] = 16;
    pvtIndex[5] = 6;
    pvtIndex[6] = 9;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestMaxSeparetedPivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(1627);
    pvt->generatePivots(sample,df,2);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_max_separeted_1.dat");

    MaxSeparatedPivots<double>* test = new MaxSeparatedPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_max_separeted_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestMaxSeparetedPivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(597);
    pvt->generatePivots(sample,df,5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_max_separeted_2.dat");

    MaxSeparatedPivots<double>* test = new MaxSeparatedPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_max_separeted_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestMaxSeparetedPivots::test12()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    MaxSeparatedPivots<double>* pvt = new MaxSeparatedPivots<double>();
    pvt->setSeed(12615);
    pvt->generatePivots(sample,df,7);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_max_separeted_3.dat");

    MaxSeparatedPivots<double>* test = new MaxSeparatedPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_max_separeted_3.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}
