QT += testlib
QT -= gui

CONFIG += qt console warn_on depend_includepath testcase
CONFIG -= app_bundle

TEMPLATE = app

DESTDIR = $$PWD/bin
OBJECTS_DIR = $$PWD/bin
MOC_DIR = $$PWD/bin
RCC_DIR = $$PWD/bin
UI_DIR = $$PWD/bin

INCLUDEPATH += \
            ../gervLib/dataset/include/ \
            ../gervLib/Distance/ \
            ../gervLib/Distance/util/include/ \
            ../gervLib/hermes/util/include/ \
            ../gervLib/utils/include/ \
            ../gervLib/pivots/include/ \
            ../gervLib/eigen/ \

SOURCES +=  tst_test_gervlib.cpp \
    TestConvexPivots.cpp \
    TestDataset.cpp \
    TestGnatPivots.cpp \
    TestKmedoidsPivots.cpp \
    TestMaxSeparetedPivots.cpp \
    TestMaxVariancePivots.cpp \
    TestPCAPivots.cpp \
    TestRandomPivots.cpp \
    TestSSSPivots.cpp \
    TestSelectionPivots.cpp \
    main.cpp


HEADERS += \
    TestConvexPivots.h \
    TestDataset.h \
    TestGnatPivots.h \
    TestKmedoidsPivots.h \
    TestMaxSeparetedPivots.h \
    TestMaxVariancePivots.h \
    TestPCAPivots.h \
    TestRandomPivots.h \
    TestSSSPivots.h \
    TestSelectionPivots.h






HEADERS += \
        ../gervLib/utils/include/Util.h \

SOURCES += \
        ../gervLib/utils/include/Util.cpp \
