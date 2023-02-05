#include "TestConvexPivots.h"


void TestConvexPivots::initTestCase()
{

}


void TestConvexPivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_convex_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_convex_2.dat");
    remove("../../gervLib/datasets/pivot_unit_test_convex_3.dat");

}


void TestConvexPivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(1105);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 7;
    pvtIndex[1] = 1;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestConvexPivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(154);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 2;
    pvtIndex[1] = 6;
    pvtIndex[2] = 9;
    pvtIndex[3] = 8;
    pvtIndex[4] = 0;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);


}


void TestConvexPivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(1324);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 7;
    pvtIndex[1] = 1;
    pvtIndex[2] = 9;
    pvtIndex[3] = 0;
    pvtIndex[4] = 4;
    pvtIndex[5] = 3;
    pvtIndex[6] = 5;


    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestConvexPivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(3731);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 11;
    pvtIndex[1] = 12;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestConvexPivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(36);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 0;
    pvtIndex[1] = 9;
    pvtIndex[2] = 7;
    pvtIndex[3] = 6;
    pvtIndex[4] = 4;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestConvexPivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(387);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 0;
    pvtIndex[1] = 5;
    pvtIndex[2] = 7;
    pvtIndex[3] = 9;
    pvtIndex[4] = 1;
    pvtIndex[5] = 13;
    pvtIndex[6] = 14;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestConvexPivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(1500);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 6;
    pvtIndex[1] = 16;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestConvexPivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(795);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 12;
    pvtIndex[1] = 7;
    pvtIndex[2] = 6;
    pvtIndex[3] = 16;
    pvtIndex[4] = 9;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestConvexPivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(915);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 12;
    pvtIndex[1] = 6;
    pvtIndex[2] = 16;
    pvtIndex[3] = 7;
    pvtIndex[4] = 14;
    pvtIndex[5] = 9;
    pvtIndex[6] = 19;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}

void TestConvexPivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(1105);
    pvt->generatePivots(sample,df,2);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_convex_1.dat");

    ConvexPivots<double>* test = new ConvexPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_convex_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestConvexPivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(36);
    pvt->generatePivots(sample,df,5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_convex_2.dat");

    ConvexPivots<double>* test = new ConvexPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_convex_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestConvexPivots::test12()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    ConvexPivots<double>* pvt = new ConvexPivots<double>();
    pvt->setSeed(915);
    pvt->generatePivots(sample,df,7);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_convex_3.dat");

    ConvexPivots<double>* test = new ConvexPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_convex_3.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}
