#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP

#include <EuclideanDistance.h>
#include <ManhattanDistance.h>
#include <ChebyshevDistance.h>
#include <CanberraDistance.h>
#include <BasicArrayObject.h>

template <class FeatureVector>
class Evaluator{

private:
    uint16_t types;
    uint32_t ndf;

public:
    static const u_int16_t EUCLIDEAN = 1;
    static const u_int16_t CITYBLOCK = 2;
    static const u_int16_t CHEBYSHEV = 3;
    static const u_int16_t JEFFREY = 4;
    static const u_int16_t CANBERRA = 5;
    static const u_int16_t BRAYCURTIS = 6;
    static const u_int16_t QUISQUARE = 7;

public:
    /**
    * Constructor.
    *
    * @param types The distance function number.
    */
    Evaluator(uint16_t types){
        resetStatistics();
        setType(types);
    }

    /**
    * Destructor.
    *
    */
    Evaluator(){
        resetStatistics();
        setType(0);
    }

    /**
    * Increments the statistics by one.
    *
    */
    void updateStatistics(){

        ndf++;
    }


    /**
    * Sets the statistics value to 0.
    *
    */
    void resetStatistics(){

        ndf = 0;
    }


    /**
    * Returns the statistics value.
    *
    * @return The number of distance function calculations.
    */
    u_int32_t getStatistics(){

        return ndf;
    }


    /**
    * Sets the distance function
    *
    * @param distanceFunction The number of distance function.
    */
    void setType(uint16_t distanceFunction){

        this->types = distanceFunction;
    }

    void setDistanceFunction(uint16_t distanceFunction){

        setType(distanceFunction);
    }


    /**
    * Returns the distance function type.
    *
    * @return The distance function type.
    */
    uint16_t getType(){

        return types;
    }

    uint16_t distanceType(){

        return getType();
    }


    /**
    * Calculates the similarity between two feature vectors.
    *
    * @return The distance between the two compared objects.
    */
    double getDistance(FeatureVector &obj1, FeatureVector &obj2){

        return GetDistance(&obj1, &obj2);
    }


    /**
    * Calculates the distance between two feature vectors.
    *
    * @return The distance between the two compared objects.
    */
    double GetDistance(FeatureVector *obj1, FeatureVector *obj2){

        // Update stats
        updateStatistics();

        double answer = 0.0;

        if (getType() == Evaluator::EUCLIDEAN){
            EuclideanDistance<FeatureVector> d;
            answer = d.getDistance(*obj1, *obj2);
        }
        if (getType() == Evaluator::CITYBLOCK){
            ManhattanDistance<FeatureVector>  d;
            answer = d.getDistance(*obj1, *obj2);
        }
        if (getType() == Evaluator::CHEBYSHEV){
            ChebyshevDistance<FeatureVector>  d;
            answer = d.getDistance(*obj1, *obj2);
        }
        return answer;
    }
};

typedef Evaluator<FeatureVector> WrapperDF;

#endif // EVALUATOR_HPP
