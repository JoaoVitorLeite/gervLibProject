#include "TestDataset.h"



void TestDataset::initTestCase()
{



}



void TestDataset::cleanupTestCase()
{

    remove("../../gervLib/datasets/dataset_unit_test_1.dat");
    remove("../../gervLib/datasets/dataset_unit_test_2.dat");

}



void TestDataset::test1()
{

    Dataset<double> data;
    Dataset<double>::loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset1.csv", 10, 2);
    Dataset<double> data2;
    Dataset<double>::loadNumericDataset(&data2, "../../gervLib/datasets/open2/Dataset1.csv", ",");

    QCOMPARE(data.getCardinality(), (size_t)10);
    QCOMPARE(data.getDimensionality(), (size_t)2);
    QCOMPARE(data2.getCardinality(), (size_t)10);
    QCOMPARE(data2.getDimensionality(), (size_t)2);
    QCOMPARE(data.isEqual(&data2), true);

    Dataset<double> data3 = Dataset<double>();
    data3.setData(data.getElements());

    QCOMPARE(data2.isEqual(&data3),true);

}



void TestDataset::test2()
{

    Dataset<double> data;
    Dataset<double>::loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset2.csv", 15, 3);
    Dataset<double> data2;
    Dataset<double>::loadNumericDataset(&data2, "../../gervLib/datasets/open2/Dataset2.csv", ",");

    QCOMPARE(data.getCardinality(), (size_t)15);
    QCOMPARE(data.getDimensionality(), (size_t)3);
    QCOMPARE(data2.getCardinality(), (size_t)15);
    QCOMPARE(data2.getDimensionality(), (size_t)3);
    QCOMPARE(data.isEqual(&data2), true);

    Dataset<double> data3 = Dataset<double>();
    data3.setData(data.getElements());

    QCOMPARE(data2.isEqual(&data3),true);


}



void TestDataset::test3()
{

    Dataset<double> data;
    Dataset<double>::loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset3.csv", 20, 4);
    Dataset<double> data2;
    Dataset<double>::loadNumericDataset(&data2, "../../gervLib/datasets/open2/Dataset3.csv", ",");

    QCOMPARE(data.getCardinality(), (size_t)20);
    QCOMPARE(data.getDimensionality(), (size_t)4);
    QCOMPARE(data2.getCardinality(), (size_t)20);
    QCOMPARE(data2.getDimensionality(), (size_t)4);
    QCOMPARE(data.isEqual(&data2), true);

    Dataset<double> data3 = Dataset<double>();
    data3.setData(data.getElements());

    QCOMPARE(data2.isEqual(&data3),true);

}



void TestDataset::test4()
{

    Dataset<double> data;
    Dataset<double>::loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset3.csv", 20, 4);

    size_t validPos = 2, invalidPos = 200;

    BasicArrayObject<double>* aux = data.getInstance(validPos);

    QCOMPARE(aux->isEqual(data.getInstance(validPos)), true);

    QVERIFY_EXCEPTION_THROWN(data.instance(invalidPos), std::invalid_argument);

}



void TestDataset::test5()
{

    Dataset<double> data;
    Dataset<double>::loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset3.csv", 20, 4);

    std::vector<BasicArrayObject<double>*> sp = data.sample(5, false, 87);

    for(size_t x = 0; x < 5; x++)
    {

        QCOMPARE(data.contains(*sp[x]), true);

    }

}



void TestDataset::test6()
{

    Dataset<double> data;
    Dataset<double>::loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset3.csv", 20, 4);
    Dataset<double> data2 = data;
    data.shuffle(76);

    for(size_t x = 0; x < data.getCardinality(); x++)
    {

        QCOMPARE(data2.contains(*data.instance(x)),true);

    }

}



void TestDataset::test7()
{

    Dataset<double> data;
    data.loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset1.csv", " ");
    data.setSeed(7);

    unsigned char* seri = data.serialize();

    Dataset<double> data2 = Dataset<double>();

    data2.unserialize(seri);

    QCOMPARE(data.isEqual(&data2), true);

}




void TestDataset::test8()
{

    Dataset<std::vector<char>> data;
    data.loadTextDataset(&data, "../../gervLib/datasets/names.csv", 2472, 2);
    data.setSeed(7);

    unsigned char* seri = data.serialize();

    Dataset<std::vector<char>> data2 = Dataset<std::vector<char>>();

    data2.unserialize(seri);

    QCOMPARE(data.isEqual(&data2), true);

}



void TestDataset::test9()
{

    Dataset<double> data;
    data.loadNumericDataset(&data, "../../gervLib/datasets/open1/Dataset3.csv", " ");
    data.setSeed(7);

    data.saveToFile("../../gervLib/datasets/dataset_unit_test_1.dat");

    Dataset<double> data2 = Dataset<double>();

    data2.loadFromFile("../../gervLib/datasets/dataset_unit_test_1.dat");

    QCOMPARE(data.isEqual(&data2), true);

}



void TestDataset::test10()
{

    Dataset<std::vector<char>> data;
    data.loadTextDataset(&data, "../../gervLib/datasets/names.csv", 2472, 2);
    data.setSeed(7);

    data.saveToFile("../../gervLib/datasets/dataset_unit_test_2.dat");

    Dataset<std::vector<char>> data2 = Dataset<std::vector<char>>();

    data2.loadFromFile("../../gervLib/datasets/dataset_unit_test_2.dat");

    QCOMPARE(data.isEqual(&data2), true);


}



void TestDataset::test11()
{

    Dataset<std::vector<char>> data;
    data.loadTextDataset(&data, "../../gervLib/datasets/names.csv", 2472, 2);

    Dataset<std::vector<char>> data2;
    data2.loadTextDataset(&data2, "../../gervLib/datasets/names.csv", " ");

    QCOMPARE(data.isEqual(&data2), true);

}














































