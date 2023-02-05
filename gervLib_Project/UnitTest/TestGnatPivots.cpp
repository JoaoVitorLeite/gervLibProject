#include "TestGnatPivots.h"



void TestGnatPivots::initTestCase()
{

}



void TestGnatPivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_gnat_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_gnat_2.dat");
    remove("../../gervLib/datasets/pivot_unit_test_gnat_3.dat");

}



void TestGnatPivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(0);
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



void TestGnatPivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(17);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 7;
    pvtIndex[1] = 6;
    pvtIndex[2] = 1;
    pvtIndex[3] = 0;
    pvtIndex[4] = 9;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);


}



void TestGnatPivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(1983);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 1;
    pvtIndex[1] = 7;
    pvtIndex[2] = 3;
    pvtIndex[3] = 9;
    pvtIndex[4] = 5;
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



void TestGnatPivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(7659);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 0;
    pvtIndex[1] = 5;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}



void TestGnatPivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(946);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 1;
    pvtIndex[1] = 7;
    pvtIndex[2] = 10;
    pvtIndex[3] = 5;
    pvtIndex[4] = 9;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}



void TestGnatPivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(1000);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 1;
    pvtIndex[1] = 7;
    pvtIndex[2] = 12;
    pvtIndex[3] = 0;
    pvtIndex[4] = 5;
    pvtIndex[5] = 8;
    pvtIndex[6] = 4;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}



void TestGnatPivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(510);
    pvt->generatePivots(sample, df,2);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 18;
    pvtIndex[1] = 6;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}



void TestGnatPivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(295);
    pvt->generatePivots(sample, df,5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 7;
    pvtIndex[1] = 1;
    pvtIndex[2] = 19;
    pvtIndex[3] = 18;
    pvtIndex[4] = 9;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}



void TestGnatPivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(112);
    pvt->generatePivots(sample, df,7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 1;
    pvtIndex[1] = 18;
    pvtIndex[2] = 6;
    pvtIndex[3] = 9;
    pvtIndex[4] = 14;
    pvtIndex[5] = 16;
    pvtIndex[6] = 19;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}



void TestGnatPivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(0);
    pvt->generatePivots(sample,df,2);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_gnat_1.dat");

    GnatPivots<double>* test = new GnatPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_gnat_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}



void TestGnatPivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(946);
    pvt->generatePivots(sample,df,5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_gnat_2.dat");

    GnatPivots<double>* test = new GnatPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_gnat_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}



void TestGnatPivots::test12()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    GnatPivots<double>* pvt = new GnatPivots<double>();
    pvt->setSeed(112);
    pvt->generatePivots(sample,df,7);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_gnat_3.dat");

    GnatPivots<double>* test = new GnatPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_gnat_3.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


