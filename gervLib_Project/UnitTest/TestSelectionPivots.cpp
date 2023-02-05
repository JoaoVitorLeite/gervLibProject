#include "TestSelectionPivots.h"


void TestSelectionPivots::initTestCase()
{

}


void TestSelectionPivots::cleanupTestCase()
{

    remove("../../gervLib/datasets/pivot_unit_test_selection_1.dat");
    remove("../../gervLib/datasets/pivot_unit_test_selection_2.dat");
    remove("../../gervLib/datasets/pivot_unit_test_selection_3.dat");

}


void TestSelectionPivots::test1()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(111);
    pvt->generatePivots(sample, df, 2, 5);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 3;
    pvtIndex[1] = 5;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test2()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(48);
    pvt->generatePivots(sample, df, 5, 5);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 2;
    pvtIndex[1] = 1;
    pvtIndex[2] = 5;
    pvtIndex[3] = 0;
    pvtIndex[4] = 7;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);


}


void TestSelectionPivots::test3()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(159);
    pvt->generatePivots(sample, df, 7, 5);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 0;
    pvtIndex[1] = 8;
    pvtIndex[2] = 5;
    pvtIndex[3] = 3;
    pvtIndex[4] = 9;
    pvtIndex[5] = 4;
    pvtIndex[6] = 1;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test4()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(604);
    pvt->generatePivots(sample, df, 2, 7);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 5;
    pvtIndex[1] = 6;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test5()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(14);
    pvt->generatePivots(sample, df, 5, 7);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 7;
    pvtIndex[1] = 11;
    pvtIndex[2] = 6;
    pvtIndex[3] = 1;
    pvtIndex[4] = 14;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test6()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(223);
    pvt->generatePivots(sample, df, 7, 7);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 5;
    pvtIndex[1] = 13;
    pvtIndex[2] = 4;
    pvtIndex[3] = 8;
    pvtIndex[4] = 10;
    pvtIndex[5] = 2;
    pvtIndex[6] = 12;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test7()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(3);
    pvt->generatePivots(sample, df, 2, 10);

    size_t* pvtIndex = new size_t[2];
    pvtIndex[0] = 2;
    pvtIndex[1] = 0;

    for(size_t x = 0; x < 2; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test8()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(10376);
    pvt->generatePivots(sample, df, 5, 10);

    size_t* pvtIndex = new size_t[5];
    pvtIndex[0] = 7;
    pvtIndex[1] = 3;
    pvtIndex[2] = 2;
    pvtIndex[3] = 15;
    pvtIndex[4] = 5;

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test9()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(20);
    pvt->generatePivots(sample, df, 7, 10);

    size_t* pvtIndex = new size_t[7];
    pvtIndex[0] = 9;
    pvtIndex[1] = 1;
    pvtIndex[2] = 6;
    pvtIndex[3] = 10;
    pvtIndex[4] = 0;
    pvtIndex[5] = 13;
    pvtIndex[6] = 3;

    for(size_t x = 0; x < 7; x++)
    {

        QCOMPARE(pvt->getPivot(x)->getOID(), pvtIndex[x]);

    }

    delete (sample);
    delete (df);
    delete [] (pvtIndex);

}


void TestSelectionPivots::test10()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset1.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(111);
    pvt->generatePivots(sample,df,2,5);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_selection_1.dat");

    SelectionPivots<double>* test = new SelectionPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_selection_1.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestSelectionPivots::test11()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset2.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(14);
    pvt->generatePivots(sample,df,5,7);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_selection_2.dat");

    SelectionPivots<double>* test = new SelectionPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_selection_2.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}


void TestSelectionPivots::test12()
{

    Dataset<double>* sample = new Dataset<double>();
    Dataset<double>::loadNumericDataset(sample, "../../gervLib/datasets/open2/Dataset3.csv", ",");
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    SelectionPivots<double>* pvt = new SelectionPivots<double>();
    pvt->setSeed(20);
    pvt->generatePivots(sample,df,7,10);
    pvt->saveToFile("../../gervLib/datasets/pivot_unit_test_selection_3.dat");

    SelectionPivots<double>* test = new SelectionPivots<double>();
    test->loadFromFile("../../gervLib/datasets/pivot_unit_test_selection_3.dat");

    QCOMPARE(pvt->isEqual(test), true);

    delete (sample);
    delete (df);

}
