///*********************************************************************
//*                                                                    *
//* Copyright (c)                                                      *
//* ZJU-DBL, Zhejiang University, Hangzhou, China                      *
//*                                                                    *
//* All Rights Reserved.                                               *
//*                                                                    *
//* Permission to use, copy, and distribute this software and its      *
//* documentation for NON-COMMERCIAL purposes and without fee is       *
//* hereby granted provided  that this copyright notice appears in     *
//* all copies.                                                        *
//*                                                                    *
//* THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANTIES ABOUT THE        *
//* SUITABILITY OF THE SOFTWARE, EITHER EXPRESS OR IMPLIED, INCLUDING  *
//* BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY,      *
//* FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHOR  *
//* SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A      *
//* RESULT OF USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS    *
//* DERIVATIVES.                                                       *
//*                                                                    *
//*********************************************************************/

//#include <iostream>
//#include <fstream>
//#include <math.h>
//#include <time.h>
//#include <vector>
//#include <sys/stat.h>]

//using namespace std;

//#include <cstring>

//#include "spb_tree/blockfile/blk_file.h"
//#include "spb_tree/blockfile/cache.h"
//#include "spb_tree/gadget/gadget.h"
//#include "spb_tree/heap/binheap.h"
//#include "spb_tree/heap/heap.h"
////#include "./rand/rand.h"
//#include "BasicArrayObject.h"
//#include "spb_tree/RAF.h"
//#include "spb_tree/pb-tree.h"
//#include <string>
//#include "spb_tree/spb_config.h"
//double compdists = 0;
//int func;
////double IOread = 0;
////double IOwrite = 0;
////double cc = 0;
////int MAXINT, BITS;
////double EPS, MAXDIST;
////int bolcksize;
////int dim, num_obj, func;

////Note that the id of data in the file should begin with 0
//BasicArrayObject<double> ** readobjects(char* filename, int num_obj, int dim)
//{
//    FILE *fp;

//    int record_count = 0;
//    BasicArrayObject<double> ** objset = new BasicArrayObject<double>*[num_obj];
//    for (int i = 0;i<num_obj;i++)
//    {
//        objset[i] = new BasicArrayObject<double>();
//    }

//    if ((fp = fopen(filename, "r")) == NULL)
//    {
//        for (int i = 0;i<num_obj;i++)
//            delete objset[i];
//        error("Could not open the data file", TRUE);
//    }
//    else
//    {
//        fscanf(fp, "%d %d %d", &dim, &num_obj, &func);
//        while (record_count<num_obj)
//        {
//            float d;
//            objset[record_count]->setOID(record_count);
////            objset[record_count]->resize(dim);
//            //objset[record_count]->data = new float[objset[record_count]->size];
//            for (int i = 0; i < dim; i++)
//            {
//                fscanf(fp, "%f", &d);
//                objset[record_count]->set(i, d);
//            }
//            fscanf(fp, "\n");
//            record_count++;
//        }
//    }

//    return objset;
//}

//BasicArrayObject<double> * readobjects2(char* filename)
//{
//    FILE *fp;
//    int record_count = 0;

//    if ((fp = fopen(filename, "r")) == NULL)
//    {
//        error("Could not open the data file", TRUE);
//        return NULL;
//    }
//    else
//    {
//        fscanf(fp, "%d %d %d", &dim, &num_obj, &func);

//        BasicArrayObject<double> * objset = new BasicArrayObject<double>[num_obj];
//        while (record_count < num_obj)
//        {
//            float d;
//            objset[record_count].setOID(record_count);
//            //objset[record_count].resize(dim);
//            //objset[record_count].data = new float[dim];
//            for (int i = 0; i < dim; i++)
//            {
//                fscanf(fp, "%f", &d);
//                objset[record_count].set(i,d);
//            }
//            fscanf(fp, "\n");
//            record_count++;
//        }
//        return objset;
//    }

//}

//clock_t bulidIndex(char* index_name, char * filename, char*pname,int n_p)
//{
//    /*******
//    index_name -- the path of the index to be stored
//    filename -- the path of the dataset
//    pname -- the path of the pivot file
//    n_p -- the number of pivots
//    *******/

//    //pivot selection
//    PB_Tree<double> * pb = new PB_Tree<double>();

//    //parameter settings
//    pb->df = new EuclideanDistance<BasicArrayObject<double>>();
//    pb->eps = EPS;
//    pb->num_cand = 40; // the number of candidates for selecting pivots (if we read the pivots from pivot file, it is useless here)
//    pb->num_piv = n_p;
//    pb->bits = BITS;
//    int keysize;
//    BasicArrayObject<double> * os = readobjects2(filename); // the function needed to be rewirte to read differnt formats of the dataset

//    pb->readptable(pname);

////    IOread = IOwrite = compdists = 0;
//    clock_t begin;

//    begin = clock();
//    pb->H_bulkload(os, num_obj);
//    keysize = os[0].keysize;
//    delete[] os;
//    cout << "keysize:" << keysize << endl;
//    pb->bplus = new B_Tree();
//    pb->bplus->init(index_name, bolcksize, NULL, keysize);
//    pb->bplus->bulkload("bulkload.txt");

//    //build RAF
//    pb->bplus->load_root();
//    B_Node * node = pb->bplus->root_ptr;
//    B_Node * temp;

//    while (node->level != 0)
//    {
//        node = node->entries[0]->get_son();
//    }
//    int * obj_order = new int[num_obj];
//    int k = 0;
//    for (int i = 0;i<node->num_entries;i++)
//    {
//        obj_order[k] = node->entries[i]->son;

//        k++;
//    }
//    temp = node->get_right_sibling();
//    //delete node;
//    while (temp != NULL)
//    {

//        node = temp;
//        for (int i = 0;i<node->num_entries;i++)
//        {
//            obj_order[k] = node->entries[i]->son;
//            k++;
//        }
//        temp = node->get_right_sibling();
//        delete node;
//    }

//    pb->draf = new RAF<double>();
//    pb->draf->num_obj = num_obj;
//    pb->draf->init(index_name, bolcksize, NULL);

//    BasicArrayObject<double> ** objS = readobjects(filename, num_obj, dim);
//    int * result = pb->draf->buid_from_array(objS, obj_order);

//    //delete object sets
//    for (int i = 0;i<num_obj;i++)
//        delete objS[i];
//    delete[] obj_order;

//    node = pb->bplus->root_ptr;


//    while (node->level != 0)
//    {
//        node = node->entries[0]->get_son();
//    }

//    k = 0;
//    for (int i = 0;i<node->num_entries;i++)
//    {
//        node->entries[i]->ptr = result[k];
//        k++;
//    }
//    char * buffer = new char[bolcksize];
//    node->write_to_buffer(buffer);

//    pb->bplus->file->write_block(buffer, node->block);
//    delete [] buffer;

//    temp = node->get_right_sibling();
//    //delete node;
//    while (temp != NULL)
//    {
//        buffer = new char[bolcksize];
//        node = temp;
//        for (int i = 0;i<node->num_entries;i++)
//        {
//            node->entries[i]->ptr = result[k];
//            k++;
//        }
//        node->write_to_buffer(buffer);

//        pb->bplus->file->write_block(buffer, node->block);

//        delete [] buffer;
//        temp = node->get_right_sibling();
//        delete node;
//    }
//    delete[] result;

//    pb->bplus->close();
//    return begin;
//}

//#define ARRAYLENGTH 100



//int * build_MBR(B_Node * node, PB_Tree<double> *pb)
//{
//    int* minp = new int[pb->num_piv];
//    int* maxp = new int[pb->num_piv];
//    for (int i = 0;i<pb->num_piv;i++)
//    {
//        minp[i] = MAXINT;
//        maxp[i] = 0;
//    }

//    if (node->level == 0)
//    {
//        for (int i = 0;i<node->num_entries;i++)
//        {
//            int* t = pb->R_Zconvert(node->entries[i]->key);
//            for (int j = 0;j<pb->num_piv;j++)
//            {
//                if (t[j]<minp[j])
//                    minp[j] = t[j];
//                if (t[j]>maxp[j])
//                    maxp[j] = t[j];
//            }
//            delete[] t;
//        }

//    }
//    else
//    {
//        for (int i = 0;i<node->num_entries;i++)
//        {
//            int* t = build_MBR(node->entries[i]->get_son(), pb);
//            node->entries[i]->del_son();
//            int * t1 = new int[pb->num_piv];
//            for (int j = 0;j<pb->num_piv;j++)
//            {
//                t1[j] = t[pb->num_piv + j];
//            }
//            for (int j = 0;j<pb->num_piv;j++)
//            {
//                if (t[j]<minp[j])
//                    minp[j] = t[j];
//                if (t1[j]>maxp[j])
//                    maxp[j] = t1[j];
//            }

//            node->entries[i]->min = pb->Zconvert(t);
//            node->entries[i]->max = pb->Zconvert(t1);

//            delete[]  t;
//            delete[] t1;
//        }

//        char * buffer = new char[bolcksize];
//        node->write_to_buffer(buffer);

//        pb->bplus->file->write_block(buffer, node->block);
//        delete [] buffer;


//    }
//    int * temp = new int[2 * pb->num_piv];
//    for (int i = 0;i<pb->num_piv;i++)
//    {
//        temp[i] = minp[i];
//        temp[pb->num_piv + i] = maxp[i];
//    }

//    delete[] minp;
//    delete[] maxp;
//    return temp;
//}

//int * H_build_MBR(B_Node * node, PB_Tree<double> *pb)
//{
//    int* minp = new int[pb->num_piv];
//    int* maxp = new int[pb->num_piv];
//    for (int i = 0;i<pb->num_piv;i++)
//    {
//        minp[i] = MAXINT;
//        maxp[i] = 0;
//    }

//    if (node->level == 0)
//    {
//        for (int i = 0;i<node->num_entries;i++)
//        {

//            unsigned * key = new unsigned[pb->num_piv];
//            unsigned * t = new unsigned[pb->num_piv];
//            for (int j = 0;j<pb->bplus->keysize;j++)
//            {
//                t[j] = 0;
//                key[pb->num_piv + j - pb->bplus->keysize] = node->entries[i]->key[j];
//            }
//            for (int j = pb->bplus->keysize;j<pb->num_piv;j++)
//            {
//                t[j] = 0;
//                key[j - pb->bplus->keysize] = 0;
//            }
//            pb->R_Hconvert(t, key, pb->num_piv);

//            for (int j = 0;j<pb->num_piv;j++)
//            {
//                if (t[j]<minp[j])
//                    minp[j] = t[j];
//                if (t[j]>maxp[j])
//                    maxp[j] = t[j];
//            }
//            delete[] t;
//            delete[] key;
//        }

//    }
//    else
//    {
//        for (int i = 0;i<node->num_entries;i++)
//        {
//            int* t = H_build_MBR(node->entries[i]->get_son(), pb);
//            node->entries[i]->del_son();
//            int * t1 = new int[pb->num_piv];
//            for (int j = 0;j<pb->num_piv;j++)
//            {
//                t1[j] = t[pb->num_piv + j];
//            }
//            for (int j = 0;j<pb->num_piv;j++)
//            {
//                if (t[j]<minp[j])
//                    minp[j] = t[j];
//                if (t1[j]>maxp[j])
//                    maxp[j] = t1[j];
//            }

//            unsigned * mi = new unsigned[pb->num_piv];
//            unsigned * ma = new unsigned[pb->num_piv];
//            pb->Hconvert(mi, (unsigned*)t, pb->num_piv);
//            pb->Hconvert(ma, (unsigned*)t1, pb->num_piv);
//            node->entries[i]->min = new unsigned[pb->bplus->keysize];
//            node->entries[i]->max = new unsigned[pb->bplus->keysize];
//            for (int j = 0;j<pb->bplus->keysize;j++)
//            {
//                node->entries[i]->min[j] = mi[j + pb->num_piv - pb->bplus->keysize];
//                node->entries[i]->max[j] = ma[j + pb->num_piv - pb->bplus->keysize];
//            }

//            delete[] mi;
//            delete[] ma;
//            delete[]  t;
//            delete[] t1;
//        }

//        char * buffer = new char[bolcksize];
//        node->write_to_buffer(buffer);

//        pb->bplus->file->write_block(buffer, node->block);
//        delete [] buffer;

//    }

//    int * temp = new int[2 * pb->num_piv];
//    for (int i = 0;i<pb->num_piv;i++)
//    {
//        temp[i] = minp[i];
//        temp[pb->num_piv + i] = maxp[i];
//    }

//    delete[] minp;
//    delete[] maxp;
//    return temp;
//}

//int main(int argc, char** argv)
//{

//    //******************************build the index***********
//    clock_t begin, buildEnd, queryEnd;
//    double buildComp, queryComp;
//    struct stat sdata1;
//    struct stat sdata2;

//    int buffer_size = 32;

//    char * datafile = "../datasets/LA.txt";//argv[1];// the path of input data file
//    char *pname = "../datasets/pivot_all_LA.txt";//argv[2];// the path of input pivots
//    char * indexfile = "../datasets/spb-index-file";//argv[3];// the path to store the built index
//    FILE * f = fopen("../datasets/spb-cost.txt", "w");//fopen(argv[4], "w"); // the path to store the building and query cost
//    MAXDIST = 25000;//atof(argv[5]);// the maximum distance for the input dataset
//    EPS = MAXDIST / 1000;
//    MAXINT = (MAXDIST / EPS);
//    BITS = ((int)log2(MAXINT) + 1); //  the bits to represent space filling curve values

//    double radius[7];
//    int kvalues[] = {1, 5, 10, 20, 50, 100};

//    if (string(datafile).find("LA") != -1) {
//        double r[] = {473, 692, 989, 1409, 1875, 2314, 3096 };
//        memcpy(radius, r, sizeof(r));
//    }
//    else if (string(datafile).find("integer") != -1) {
//        double r[] = {2321,2733, 3229,3843, 4614, 5613, 7090 };
//        memcpy(radius, r, sizeof(r));
//    }
//    else if (string(datafile).find("sf3.txt") != -1) {
//        double r[] = { 100, 200, 300, 400, 500, 600, 700 };
//        memcpy(radius, r, sizeof(r));
//    }
//    else if (string(datafile).find("mpeg") != -1) {

//         double r[] = {3838, 4092, 4399, 4773, 5241, 5904, 7104};
//        memcpy(radius, r, sizeof(r));
//    }

//    int pn = 5;//atoi(argv[6]);// the number of pivots
//    bolcksize = 4096;//atoi(argv[7]);	// the page size
//    char * querydata = "../datasets/LA_query.txt";//argv[8];// the path of input query data
//    fprintf(f, "pivotnum: %d\n", pn);
//    //compdists = 0;
//    IOread = IOwrite = 0;
//    begin = bulidIndex(indexfile, datafile, pname,pn);

//    PB_Tree<double> * pb = new PB_Tree<double>();

//    Cache* c = new Cache(buffer_size, bolcksize);
//    pb->c = c;

//    pb->df = new EuclideanDistance<BasicArrayObject<double>>();
//    pb->num_piv = pn;
//    pb->readptable(pname);
//    pb->eps = EPS;
//    pb->bits = BITS;

//    pb->bplus = new B_Tree();
//    pb->bplus->init_restore(indexfile, NULL);
//    pb->bplus->load_root();
//    H_build_MBR(pb->bplus->root_ptr, pb);
//    buildEnd = clock() - begin;
//    buildComp = pb->df->getDistanceCount();
//    fprintf(f, "building...\n");
//    fprintf(f, "finished... %f build time\n", (double)buildEnd / CLOCKS_PER_SEC);
//    fprintf(f, "finished... %f distances computed\n", buildComp);
//    fprintf(f, "finished... %f IO times\n", IOread + IOwrite);

//    char * bfile = new char[strlen(indexfile) + 2];
//    strcpy(bfile, indexfile);
//    strcat(bfile, ".b");
//    char * raffile = new char[strlen(indexfile) + 4];
//    strcpy(raffile, indexfile);
//    strcat(raffile, ".raf");
//    stat(bfile, &sdata1);
//    stat(raffile, &sdata2);
//    fprintf(f, "saved... %lli bytes\n", (long long)(sdata1.st_size + sdata2.st_size));
//    fflush(f);
//    //************end of build index*************************

//    //************************ similarity searh***********************
//    fprintf(f, "\nquerying...\n");

//    pb->draf = new RAF<double>();
//    pb->draf->init_restore(indexfile, c);


//    BasicArrayObject<double> * q = new BasicArrayObject<double>();
//    ifstream in;
//    double io = 0;
//    double dists = 0;
//    int qcount = 100; // the number of quereis
//    q->resize(dim);   // the dimension of the query object
//    //q->data = new float[q->size];
//    double rad;
//    cout << "start knnSearching......" << endl;
//    for (int k = 0; k < 6; k++) {
//        in.open(querydata);
//        begin = clock();
//        IOread = IOwrite = 0;
//        dists = 0;
//        rad = 0;
//        double pf = 0;
//        for (int j = 0;j < qcount; j++)
//        {

//            c->clear();
//            compdists = 0;

//            double dIn;
//            for (int i = 0;i < q->getSize();i++)
//            {
//                in >> dIn;
//                q->Set(i, dIn);
//            }

//            rad += pb->BFkNN(q, kvalues[k]); //kNN query function

//            pf += c->page_faults;
//            dists += compdists;
//        }
//        queryEnd = clock() - begin;
//        queryComp = dists;
//        fprintf(f, "k: %d\n", kvalues[k]);
//        fprintf(f, "finished... %f query time\n", (double)queryEnd / CLOCKS_PER_SEC / qcount);
//        fprintf(f, "finished... %f distances computed\n", queryComp / qcount);
//        fprintf(f, "finished... %f IO times\n", IOread / qcount);
//        fprintf(f, "finished... %f radius\n", rad / qcount);
//        in.close();
//    }
//    cout << "start rangeSearching......" << endl;
//    for (int k = 0; k < 7; ++k) {
//        in.open(querydata);
//        begin = clock();
//        IOread = IOwrite = 0;
//        dists = 0;
//        rad = 0;
//        double pf = 0;
//        for (int j = 0;j < qcount; j++)
//        {

//            c->clear();
//            compdists = 0;

//            double dIn;
//            for (int i = 0;i < q->GetSize();i++)
//            {
//                in >> dIn;
//                q->set(i, dIn);
//            }

//            rad += pb->ImprovedRQ_new(q, radius[k]);  // range query function

//            pf += c->page_faults;
//            dists += compdists;
//        }
//        queryEnd = clock() - begin;
//        queryComp = dists;
//        fprintf(f, "r: %f\n", radius[k]);
//        fprintf(f, "finished... %f query time\n", (double)queryEnd / CLOCKS_PER_SEC / qcount);
//        fprintf(f, "finished... %f distances computed\n", queryComp / qcount);
//        fprintf(f, "finished... %f IO times\n", IOread / qcount);

//        fprintf(f, "finished... %f objs\n", rad / qcount);
//        in.close();
//    }
//        delete q;
//        q = NULL;
//        pb = NULL;
//        fprintf(f, "\n");

//    //***********************end of similarity search*******************
//}

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <sys/stat.h>

using namespace std;

#include <cstring>

#include "spb_tree/blockfile/blk_file.h"
#include "spb_tree/blockfile/cache.h"
#include "spb_tree/gadget/gadget.h"
#include "spb_tree/heap/binheap.h"
#include "spb_tree/heap/heap.h"
//#include "./rand/rand.h"
#include "BasicArrayObject.h"
#include "spb_tree/RAF.h"
#include "spb_tree/pb-tree.h"
#include <string>
#include "spb_tree/spb_config.h"
double compdists = 0;
int func;
//double IOread = 0;
//double IOwrite = 0;
//double cc = 0;
//int MAXINT, BITS;
//double EPS, MAXDIST;
//int bolcksize;
//int dim, num_obj, func;

template <class type>
clock_t buildIndex(char* index_name, char * filename, std::string sep, DistanceFunction<BasicArrayObject<type>>* df, char* pname, int n_p, Pivot<type>* pvt = nullptr)
{
    /*******
    index_name -- the path of the index to be stored
    filename -- the path of the dataset
    pname -- the path of the pivot file
    n_p -- the number of pivots
    *******/

    //pivot selection
    PB_Tree<type> * pb = new PB_Tree<type>();

    //parameter settings
    pb->pvtMethod = pvt;
    pb->df = df;
    pb->eps = EPS;
    pb->num_cand = 40; // the number of candidates for selecting pivots (if we read the pivots from pivot file, it is useless here)
    pb->num_piv = n_p;
    pb->bits = BITS;
    int keysize;
    Dataset<type>* data = new Dataset<type>();

    if constexpr (std::is_same<type, std::vector<char>>::value)
    {

        Dataset<type>::loadTextDataset(data, filename, sep);

    }
    else
    {

        Dataset<type>::loadNumericDataset(data, filename, sep);

    }

    dim = data->getDimensionality();
    num_obj = data->getCardinality();
    BasicArrayObject<type>* os = data->getInstances(); // the function needed to be rewirte to read differnt formats of the dataset

    pb->readptable(pname, data);

//    IOread = IOwrite = compdists = 0;
    clock_t begin;

    begin = clock();
    pb->H_bulkload(os, num_obj);
    keysize = os[0].keysize;
    delete[] os;
    cout << "keysize:" << keysize << endl;
    pb->bplus = new B_Tree();
    pb->bplus->init(index_name, bolcksize, NULL, keysize);
    pb->bplus->bulkload((char*)"bulkload.txt");

    //build RAF
    pb->bplus->load_root();
    B_Node * node = pb->bplus->root_ptr;
    B_Node * temp;

    while (node->level != 0)
    {
        node = node->entries[0]->get_son();
    }
    int * obj_order = new int[num_obj];
    int k = 0;
    for (int i = 0;i<node->num_entries;i++)
    {
        obj_order[k] = node->entries[i]->son;

        k++;
    }
    temp = node->get_right_sibling();
    //delete node;
    while (temp != NULL)
    {

        node = temp;
        for (int i = 0;i<node->num_entries;i++)
        {
            obj_order[k] = node->entries[i]->son;
            k++;
        }
        temp = node->get_right_sibling();
        delete node;
    }

    pb->draf = new RAF<type>();
    pb->draf->num_obj = num_obj;
    pb->draf->init(index_name, bolcksize, NULL);

    BasicArrayObject<type> ** objS = data->getInstancesArray();
    int * result = pb->draf->buid_from_array(objS, obj_order);

    //delete object sets
    for (int i = 0;i<num_obj;i++)
        delete objS[i];
    delete[] obj_order;

    node = pb->bplus->root_ptr;


    while (node->level != 0)
    {
        node = node->entries[0]->get_son();
    }

    k = 0;
    for (int i = 0;i<node->num_entries;i++)
    {
        node->entries[i]->ptr = result[k];
        k++;
    }
    char * buffer = new char[bolcksize];
    node->write_to_buffer(buffer);

    pb->bplus->file->write_block(buffer, node->block);
    delete [] buffer;

    temp = node->get_right_sibling();
    //delete node;
    while (temp != NULL)
    {
        buffer = new char[bolcksize];
        node = temp;
        for (int i = 0;i<node->num_entries;i++)
        {
            node->entries[i]->ptr = result[k];
            k++;
        }
        node->write_to_buffer(buffer);

        pb->bplus->file->write_block(buffer, node->block);

        delete [] buffer;
        temp = node->get_right_sibling();
        delete node;
    }
    delete[] result;

    pb->bplus->close();

    delete data;

    return begin;

}

template clock_t buildIndex<double>(char*, char*, std::string , DistanceFunction<BasicArrayObject<double>>* , char* , int , Pivot<double>*);
template clock_t buildIndex<std::vector<char>>(char*, char*, std::string , DistanceFunction<BasicArrayObject<std::vector<char>>>* , char* , int , Pivot<std::vector<char>>*);

template <class type>
int * build_MBR(B_Node * node, PB_Tree<type> *pb)
{
    int* minp = new int[pb->num_piv];
    int* maxp = new int[pb->num_piv];
    for (int i = 0;i<pb->num_piv;i++)
    {
        minp[i] = MAXINT;
        maxp[i] = 0;
    }

    if (node->level == 0)
    {
        for (int i = 0;i<node->num_entries;i++)
        {
            int* t = pb->R_Zconvert(node->entries[i]->key);
            for (int j = 0;j<pb->num_piv;j++)
            {
                if (t[j]<minp[j])
                    minp[j] = t[j];
                if (t[j]>maxp[j])
                    maxp[j] = t[j];
            }
            delete[] t;
        }

    }
    else
    {
        for (int i = 0;i<node->num_entries;i++)
        {
            int* t = build_MBR(node->entries[i]->get_son(), pb);
            node->entries[i]->del_son();
            int * t1 = new int[pb->num_piv];
            for (int j = 0;j<pb->num_piv;j++)
            {
                t1[j] = t[pb->num_piv + j];
            }
            for (int j = 0;j<pb->num_piv;j++)
            {
                if (t[j]<minp[j])
                    minp[j] = t[j];
                if (t1[j]>maxp[j])
                    maxp[j] = t1[j];
            }

            node->entries[i]->min = pb->Zconvert(t);
            node->entries[i]->max = pb->Zconvert(t1);

            delete[]  t;
            delete[] t1;
        }

        char * buffer = new char[bolcksize];
        node->write_to_buffer(buffer);

        pb->bplus->file->write_block(buffer, node->block);
        delete [] buffer;


    }
    int * temp = new int[2 * pb->num_piv];
    for (int i = 0;i<pb->num_piv;i++)
    {
        temp[i] = minp[i];
        temp[pb->num_piv + i] = maxp[i];
    }

    delete[] minp;
    delete[] maxp;
    return temp;
}

template <class type>
int * H_build_MBR(B_Node * node, PB_Tree<type> *pb)
{
    int* minp = new int[pb->num_piv];
    int* maxp = new int[pb->num_piv];
    for (int i = 0;i<pb->num_piv;i++)
    {
        minp[i] = MAXINT;
        maxp[i] = 0;
    }

    if (node->level == 0)
    {
        for (int i = 0;i<node->num_entries;i++)
        {

            unsigned * key = new unsigned[pb->num_piv];
            unsigned * t = new unsigned[pb->num_piv];
            for (int j = 0;j<pb->bplus->keysize;j++)
            {
                t[j] = 0;
                key[pb->num_piv + j - pb->bplus->keysize] = node->entries[i]->key[j];
            }
            for (int j = pb->bplus->keysize;j<pb->num_piv;j++)
            {
                t[j] = 0;
                key[j - pb->bplus->keysize] = 0;
            }
            pb->R_Hconvert(t, key, pb->num_piv);

            for (int j = 0;j<pb->num_piv;j++)
            {
                if (t[j]<minp[j])
                    minp[j] = t[j];
                if (t[j]>maxp[j])
                    maxp[j] = t[j];
            }
            delete[] t;
            delete[] key;
        }

    }
    else
    {
        for (int i = 0;i<node->num_entries;i++)
        {
            int* t = H_build_MBR(node->entries[i]->get_son(), pb);
            node->entries[i]->del_son();
            int * t1 = new int[pb->num_piv];
            for (int j = 0;j<pb->num_piv;j++)
            {
                t1[j] = t[pb->num_piv + j];
            }
            for (int j = 0;j<pb->num_piv;j++)
            {
                if (t[j]<minp[j])
                    minp[j] = t[j];
                if (t1[j]>maxp[j])
                    maxp[j] = t1[j];
            }

            unsigned * mi = new unsigned[pb->num_piv];
            unsigned * ma = new unsigned[pb->num_piv];
            pb->Hconvert(mi, (unsigned*)t, pb->num_piv);
            pb->Hconvert(ma, (unsigned*)t1, pb->num_piv);
            node->entries[i]->min = new unsigned[pb->bplus->keysize];
            node->entries[i]->max = new unsigned[pb->bplus->keysize];
            for (int j = 0;j<pb->bplus->keysize;j++)
            {
                node->entries[i]->min[j] = mi[j + pb->num_piv - pb->bplus->keysize];
                node->entries[i]->max[j] = ma[j + pb->num_piv - pb->bplus->keysize];
            }

            delete[] mi;
            delete[] ma;
            delete[]  t;
            delete[] t1;
        }

        char * buffer = new char[bolcksize];
        node->write_to_buffer(buffer);

        pb->bplus->file->write_block(buffer, node->block);
        delete [] buffer;

    }

    int * temp = new int[2 * pb->num_piv];
    for (int i = 0;i<pb->num_piv;i++)
    {
        temp[i] = minp[i];
        temp[pb->num_piv + i] = maxp[i];
    }

    delete[] minp;
    delete[] maxp;
    return temp;
}

int main(int argc, char** argv)
{

    //******************************build the index***********
    clock_t begin, buildEnd, queryEnd;
    double buildComp, queryComp;
    struct stat sdata1;
    struct stat sdata2;

    int buffer_size = 32;

    char * datafile = (char*)"../datasets/LA.txt";//argv[1];// the path of input data file
    char *pname = (char*)"../datasets/pivot_all_LA.txt";//argv[2];// the path of input pivots
    char * indexfile = (char*)"../datasets/spb-index-file";//argv[3];// the path to store the built index
    FILE * f = fopen("../datasets/spb-cost.txt", "w");//fopen(argv[4], "w"); // the path to store the building and query cost
    MAXDIST = 25000;//atof(argv[5]);// the maximum distance for the input dataset
    EPS = MAXDIST / 1000;
    MAXINT = (MAXDIST / EPS);
    BITS = ((int)log2(MAXINT) + 1); //  the bits to represent space filling curve values

    double radius[7];
    int kvalues[] = {1, 5, 10, 20, 50, 100};

    if (string(datafile).find("LA") != -1) {
        double r[] = {473, 692, 989, 1409, 1875, 2314, 3096 };
        memcpy(radius, r, sizeof(r));
    }
    else if (string(datafile).find("integer") != -1) {
        double r[] = {2321,2733, 3229,3843, 4614, 5613, 7090 };
        memcpy(radius, r, sizeof(r));
    }
    else if (string(datafile).find("sf3.txt") != -1) {
        double r[] = { 100, 200, 300, 400, 500, 600, 700 };
        memcpy(radius, r, sizeof(r));
    }
    else if (string(datafile).find("mpeg") != -1) {
        double r[] = {3838, 4092, 4399, 4773, 5241, 5904, 7104};
        memcpy(radius, r, sizeof(r));
    }

    int pn = 5;//atoi(argv[6]);// the number of pivots
    bolcksize = 4096;//atoi(argv[7]);	// the page size
    char * querydata = (char*)"../datasets/LA_query.txt";//argv[8];// the path of input query data
    fprintf(f, "pivotnum: %d\n", pn);
    //compdists = 0;
    IOread = IOwrite = 0;
    DistanceFunction<BasicArrayObject<double>>* df = new EuclideanDistance<BasicArrayObject<double>>();
    begin = buildIndex<double>(indexfile, datafile, " ", df, pname, pn);

    PB_Tree<double> * pb = new PB_Tree<double>();

    Cache* c = new Cache(buffer_size, bolcksize);
    pb->c = c;

    pb->df = new EuclideanDistance<BasicArrayObject<double>>();
    pb->num_piv = pn;
    pb->readptable(pname);
    pb->eps = EPS;
    pb->bits = BITS;

    pb->bplus = new B_Tree();
    pb->bplus->init_restore(indexfile, NULL);
    pb->bplus->load_root();
    H_build_MBR(pb->bplus->root_ptr, pb);
    buildEnd = clock() - begin;
    buildComp = pb->df->getDistanceCount();
    fprintf(f, "building...\n");
    fprintf(f, "finished... %f build time\n", (double)buildEnd / CLOCKS_PER_SEC);
    fprintf(f, "finished... %f distances computed\n", buildComp);
    fprintf(f, "finished... %f IO times\n", IOread + IOwrite);

    char * bfile = new char[strlen(indexfile) + 2];
    strcpy(bfile, indexfile);
    strcat(bfile, ".b");
    char * raffile = new char[strlen(indexfile) + 4];
    strcpy(raffile, indexfile);
    strcat(raffile, ".raf");
    stat(bfile, &sdata1);
    stat(raffile, &sdata2);
    fprintf(f, "saved... %lli bytes\n", (long long)(sdata1.st_size + sdata2.st_size));
    fflush(f);
    //************end of build index*************************

    //************************ similarity searh***********************
    fprintf(f, "\nquerying...\n");

    pb->draf = new RAF<double>();
    pb->draf->init_restore(indexfile, c);

    BasicArrayObject<double> * q = new BasicArrayObject<double>(-1, dim);
    q->resize(dim);   // the dimension of the query object
    ifstream in;
    double io = 0;
    double dists = 0;
    int qcount = 100; // the number of quereis
    //q->data = new float[q->size];
    double rad;
    cout << "start knnSearching......" << endl;
    for (int k = 0; k < 6; k++) {
        in.open(querydata);
        begin = clock();
        IOread = IOwrite = 0;
        dists = 0;
        rad = 0;
        double pf = 0;

        for (int j = 0;j < qcount; j++)
        {

            c->clear();
            compdists = 0;

            double dIn;
            for (int i = 0;i < q->getSize();i++)
            {
                in >> dIn;
                q->set(i, dIn);
            }

            rad += pb->BFkNN(q, kvalues[k]); //kNN query function

            pf += c->page_faults;
            dists += df->getDistanceCount();
        }
        queryEnd = clock() - begin;
        queryComp = dists;
        fprintf(f, "k: %d\n", kvalues[k]);
        fprintf(f, "finished... %f query time\n", (double)queryEnd / CLOCKS_PER_SEC / qcount);
        fprintf(f, "finished... %f distances computed\n", queryComp / qcount);
        fprintf(f, "finished... %f IO times\n", IOread / qcount);
        fprintf(f, "finished... %f radius\n", rad / qcount);
        in.close();
    }
    cout << "start rangeSearching......" << endl;
    for (int k = 0; k < 7; ++k) {
        in.open(querydata);
        begin = clock();
        IOread = IOwrite = 0;
        dists = 0;
        rad = 0;
        double pf = 0;
        for (int j = 0;j < qcount; j++)
        {

            c->clear();
            compdists = 0;

            double dIn;
            for (int i = 0;i < q->GetSize();i++)
            {
                in >> dIn;
                q->set(i, dIn);
            }

            rad += pb->ImprovedRQ_new(q, radius[k]);  // range query function

            pf += c->page_faults;
            dists += df->getDistanceCount();
        }
        queryEnd = clock() - begin;
        queryComp = dists;
        fprintf(f, "r: %f\n", radius[k]);
        fprintf(f, "finished... %f query time\n", (double)queryEnd / CLOCKS_PER_SEC / qcount);
        fprintf(f, "finished... %f distances computed\n", queryComp / qcount);
        fprintf(f, "finished... %f IO times\n", IOread / qcount);

        fprintf(f, "finished... %f objs\n", rad / qcount);
        in.close();
    }
    delete q;
    q = NULL;
    pb = NULL;
    fprintf(f, "\n");

    //***********************end of similarity search*******************

    return 0;

}
