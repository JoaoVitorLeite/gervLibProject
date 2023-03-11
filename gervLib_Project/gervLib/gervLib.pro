QT -= gui

CONFIG += c++17 console
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

DESTDIR = $$PWD/bin
OBJECTS_DIR = $$PWD/bin
MOC_DIR = $$PWD/bin
RCC_DIR = $$PWD/bin
UI_DIR = $$PWD/bin

INCLUDEPATH += \
            Distance/ \
            Distance/util/include/ \
            utils/include/ \
            dataset/include/ \
            pivots/include/ \
            eigen/ \
            VpTree/Node/Bucket/ \
            VpTree/Node/ \
            VpTree/ \
            pm_tree/ \
            omni/ \
            kdtree/ \
            mvptree/ \


HEADERS += \
        VpTree/Node/Bucket/Bucket.h \
        VpTree/Node/Bucket/Pair.h \
        VpTree/Node/DirectorNode.h \
        VpTree/Node/LeafNodeVPTree.h \
        VpTree/Node/Node.h \
        VpTree/QueueItem.h \
        VpTree/VpTree.h \
        kdtree/DirectoryNodeKdTree.h \
        kdtree/KdTree.h \
        kdtree/LeafNodeKdTree.h \
        kdtree/NodeKdTree.h \
    #memory/MemoryManager.h \
        mvptree/datapoint.h \
        mvptree/mvpnode.h \
        mvptree/mvptree.h \
        omni/OmniKdTree.h \
        pivots/include/Cluster.h \
        pivots/include/ConvexPivots.h \
        pivots/include/FFTPivots.h \
        pivots/include/GnatPivots.h \
        dataset/include/Dataset.h \
        Distance/BasicArrayObject.h \
        pivots/include/HFIPivots.h \
        pivots/include/ISPivots.h \
        pivots/include/Kmedoids.h \
        pivots/include/KmedoidsPivots.h \
        pivots/include/MaxSeparetedPivots.h \
        pivots/include/MaxVariancePivots.h \
        pivots/include/PCAPivots.h \
        pivots/include/Pivot.h \
        pivots/include/Pivots.h \
        pivots/include/RandomPivots.h \
        pivots/include/SSSPivots.h \
        pivots/include/SelectionPivots.h \
        pivots/include/WDRPivots.h \
        pm_tree/PM_Tree.h \
        utils/include/Util.h \




SOURCES += \
        main.cpp \
        utils/include/Util.cpp \
#        query.cpp




TRANSLATIONS += \
    gervLib_en_US.ts

CONFIG += lrelease
CONFIG += embed_translations

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    hermes/gervLib.cflags \
    hermes/gervLib.cxxflags

