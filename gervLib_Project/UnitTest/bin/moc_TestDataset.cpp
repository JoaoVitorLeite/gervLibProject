/****************************************************************************
** Meta object code from reading C++ file 'TestDataset.h'
**
** Created by: The Qt Meta Object Compiler version 68 (Qt 6.1.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../TestDataset.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TestDataset.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 68
#error "This file was generated using the moc from 6.1.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_TestDataset_t {
    const uint offsetsAndSize[30];
    char stringdata0[110];
};
#define QT_MOC_LITERAL(ofs, len) \
    uint(offsetof(qt_meta_stringdata_TestDataset_t, stringdata0) + ofs), len 
static const qt_meta_stringdata_TestDataset_t qt_meta_stringdata_TestDataset = {
    {
QT_MOC_LITERAL(0, 11), // "TestDataset"
QT_MOC_LITERAL(12, 12), // "initTestCase"
QT_MOC_LITERAL(25, 0), // ""
QT_MOC_LITERAL(26, 15), // "cleanupTestCase"
QT_MOC_LITERAL(42, 5), // "test1"
QT_MOC_LITERAL(48, 5), // "test2"
QT_MOC_LITERAL(54, 5), // "test3"
QT_MOC_LITERAL(60, 5), // "test4"
QT_MOC_LITERAL(66, 5), // "test5"
QT_MOC_LITERAL(72, 5), // "test6"
QT_MOC_LITERAL(78, 5), // "test7"
QT_MOC_LITERAL(84, 5), // "test8"
QT_MOC_LITERAL(90, 5), // "test9"
QT_MOC_LITERAL(96, 6), // "test10"
QT_MOC_LITERAL(103, 6) // "test11"

    },
    "TestDataset\0initTestCase\0\0cleanupTestCase\0"
    "test1\0test2\0test3\0test4\0test5\0test6\0"
    "test7\0test8\0test9\0test10\0test11"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_TestDataset[] = {

 // content:
       9,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags, initial metatype offsets
       1,    0,   92,    2, 0x08,    0 /* Private */,
       3,    0,   93,    2, 0x08,    1 /* Private */,
       4,    0,   94,    2, 0x08,    2 /* Private */,
       5,    0,   95,    2, 0x08,    3 /* Private */,
       6,    0,   96,    2, 0x08,    4 /* Private */,
       7,    0,   97,    2, 0x08,    5 /* Private */,
       8,    0,   98,    2, 0x08,    6 /* Private */,
       9,    0,   99,    2, 0x08,    7 /* Private */,
      10,    0,  100,    2, 0x08,    8 /* Private */,
      11,    0,  101,    2, 0x08,    9 /* Private */,
      12,    0,  102,    2, 0x08,   10 /* Private */,
      13,    0,  103,    2, 0x08,   11 /* Private */,
      14,    0,  104,    2, 0x08,   12 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void TestDataset::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<TestDataset *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->initTestCase(); break;
        case 1: _t->cleanupTestCase(); break;
        case 2: _t->test1(); break;
        case 3: _t->test2(); break;
        case 4: _t->test3(); break;
        case 5: _t->test4(); break;
        case 6: _t->test5(); break;
        case 7: _t->test6(); break;
        case 8: _t->test7(); break;
        case 9: _t->test8(); break;
        case 10: _t->test9(); break;
        case 11: _t->test10(); break;
        case 12: _t->test11(); break;
        default: ;
        }
    }
    (void)_a;
}

const QMetaObject TestDataset::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_meta_stringdata_TestDataset.offsetsAndSize,
    qt_meta_data_TestDataset,
    qt_static_metacall,
    nullptr,
qt_incomplete_metaTypeArray<qt_meta_stringdata_TestDataset_t

, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>


>,
    nullptr
} };


const QMetaObject *TestDataset::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *TestDataset::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_TestDataset.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int TestDataset::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<QMetaType *>(_a[0]) = QMetaType();
        _id -= 13;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
