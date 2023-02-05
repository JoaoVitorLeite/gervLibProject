#include <iostream>
#include <QTest>
#include <TestDataset.h>
#include <TestRandomPivots.h>
#include <TestGnatPivots.h>
#include <TestConvexPivots.h>
#include <TestMaxSeparetedPivots.h>
#include <TestMaxVariancePivots.h>
#include <TestSelectionPivots.h>
#include <TestSSSPivots.h>
#include <TestPCAPivots.h>
#include <TestKmedoidsPivots.h>

int main(int argc, char *argv[])
{

    int status = 0;

//    auto runTest = [&status, argc, argv](QObject* obj) {
//        status |= QTest::qExec(obj, argc, argv);
//    };

//    // run suite
//    // runTest(new TestObject);

    status |= QTest::qExec(new TestDataset, argc, argv);
    status |= QTest::qExec(new TestRandomPivots, argc, argv);
    status |= QTest::qExec(new TestGnatPivots, argc, argv);
    status |= QTest::qExec(new TestConvexPivots, argc, argv);
    status |= QTest::qExec(new TestMaxSeparetedPivots, argc, argv);
    status |= QTest::qExec(new TestMaxVariancePivots, argc, argv);
    status |= QTest::qExec(new TestSelectionPivots, argc, argv);
    status |= QTest::qExec(new TestSSSPivots, argc, argv);
    status |= QTest::qExec(new TestPCAPivots, argc, argv);
    status |= QTest::qExec(new TestKmedoidsPivots, argc, argv);

    return status;
}
