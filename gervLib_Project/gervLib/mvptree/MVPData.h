#ifndef MVPDATA_H
#define MVPDATA_H

#include <Dataset.h>

template <class type>
struct vp_t{

    BasicArrayObject<type>* key;

    vp_t(BasicArrayObject<type>* _key)
    {

        key = _key;

    }

};

template <class type>
struct datapoint_t
{

    BasicArrayObject<type>* key;
    double* dists;

    datapoint_t(BasicArrayObject<type>* _key, size_t pl)
    {

        key = _key;
        dists = new double[pl];

        for(size_t i = 0; i < pl; i++)
        {

            dists[i] = -1.0;

        }

    }

};


#endif // MVPDATA_H
