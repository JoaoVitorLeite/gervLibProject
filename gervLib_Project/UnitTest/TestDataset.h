#ifndef TESTDATASET_H
#define TESTDATASET_H

#include <QtTest>
#include <Dataset.h>
#include <cstdio>


class TestDataset: public QObject
{
    Q_OBJECT

private slots:

    void initTestCase();
    void cleanupTestCase();
    //void init();
    //void cleanup();
    void test1();
    void test2();
    void test3();
    void test4();
    void test5();
    void test6();
    void test7();
    void test8();
    void test9();
    void test10();
    void test11();


};

#endif // TESTDATASET_H
