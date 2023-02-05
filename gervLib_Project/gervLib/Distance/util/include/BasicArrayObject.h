#ifndef BASICARRAYOBJECT_H
#define BASICARRAYOBJECT_H

#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <sstream>

/**
* For illustration, consider the feature vector as follows:
* +-----+------+------------------+
* | OID | Size | Vector Data []   |
* +-----+------+------------------+
*
* @brief This class implements a generic feature vector
* @arg DType The data type stored by each position of the feature vector
*/
template <class DType>
class BasicArrayObject{


    private:
        //The internal variable that really stores the data
        std::vector<DType> data;
        //The OID that identifies the feature vector
        size_t OID;
        //A previous directive that can allow store and retrieve
        //the feature vector from BLOB or FILE
        unsigned char *serialized;

        /**
        * Gets the equivalent byte vector of the object.
        * @return The equivalent byte vector of the object.
        */
        unsigned char *serializeNumeric(){

            if (serialized == NULL){
                serialized = new unsigned char[getSerializedSize()];
                size_t size = getSize();
                memcpy(serialized, &OID, sizeof(size_t));
                memcpy(serialized + sizeof(size_t), &size, sizeof(size_t));
                for (size_t x = 0; x < size; x++){
                    memcpy(serialized + sizeof(size_t) + sizeof(size_t)+(sizeof(DType)*x), &data[x], sizeof(DType));
                }
            }
            return serialized;
        }

        unsigned char* serializeVecChar()
        {

            if (serialized == NULL)
            {

                size_t tamanhoTotal = 2*sizeof(size_t), total = 0;

                for(size_t x = 0; x < getSize(); x++)
                {

                    tamanhoTotal += sizeof(size_t) + data[x].size()*sizeof(char);

                }

                serialized = new unsigned char[tamanhoTotal];
                size_t size = getSize();

                memcpy(serialized + total, &OID, sizeof(size_t));
                total += sizeof(size_t);

                memcpy(serialized + total, &size, sizeof(size_t));
                total += sizeof(size_t);


                for (size_t x = 0; x < getSize(); x++)
                {

                    size_t sizeVec = data[x].size();
                    memcpy(serialized + total, &sizeVec, sizeof(size_t));
                    total += sizeof(size_t);

                    memcpy(serialized + total, data[x].data(), data[x].size()*sizeof(char));
                    total += data[x].size()*sizeof(char);

                }

            }

            return serialized;

        }

        /**
        * Transform a byte vector into an object.
        * @param dataIn The byte vector.
        * @param dataSize The byte vector size.
        */
        void unserializeNumeric(unsigned char *dataIn, size_t dataSize = 0){

            DType *d = 0;
            size_t size_vector;

            // This is the reverse of Serialize(). So the steps are similar.
            // Remember, the format of the serizalized object is
            // +----------+-----------+--------+
            // | OID | Size | Vector Data []   |
            // +----------+-----------+--------+

            memcpy(&OID, dataIn, sizeof(size_t));

            if (dataSize != 0) {
                size_vector = (dataSize - sizeof(size_t) - sizeof(size_t)) / sizeof(DType);
            } else {
                memcpy(&size_vector, dataIn + sizeof(size_t), sizeof(size_t));
            }

            d = new DType[size_vector];

            for (size_t x = 0; x < size_vector; x++){
                memcpy(&d[x], dataIn + sizeof(size_t) + sizeof(size_t) + (sizeof(DType)*x), sizeof(DType));
            }

            data.clear();

            for (size_t x = 0; x < size_vector; x++){
                data.push_back(d[x]);
            }

            // Since we have changed the object contents, we must invalidate the old
            // serialized version if it exists. In fact we, may copy the given serialized
            // version of tbe new object to the buffer but we don't want to spend memory.
            if (serialized != NULL){
                delete [] serialized;
                serialized = NULL;
            }//end if
            delete[] d;
        }

        void unserializeVecChar(unsigned char *dataIn, size_t dataSize = 0)
        {

            size_t total = 0;
            size_t vecSize = OID = -1;

            memcpy(&OID, dataIn + total, sizeof(size_t));
            total += sizeof(size_t);

            memcpy(&vecSize, dataIn + total,sizeof(size_t));
            total += sizeof(size_t);

            data = std::vector<std::vector<char>>(vecSize);

            for(size_t x = 0; x < vecSize; x++)
            {

                size_t size;
                memcpy(&size, dataIn + total, sizeof(size_t));
                total += sizeof(size_t);

                char* aux = new char[size];

                memcpy(aux, dataIn + total, size*sizeof(char));
                total += size*sizeof(char);

                std::vector<char> vecChar(aux,aux+size);
                data[x] = vecChar;

            }

        }


    public:

        /**
        * Constructor Method.
        * Sets data and size to empty and 0, respectively.
        */
        BasicArrayObject(){
            data.clear();
            serialized = NULL;
        }

        /**
        * Constructor Method.
        * Sets the values of the vector to current.
        */
        BasicArrayObject(const size_t OID, const std::vector<DType> &data){

            this->OID = OID;
            for(size_t x = 0; x < data.size(); x++)
                set(data[x]);
            serialized = NULL;
        }


        BasicArrayObject(const size_t OID, size_t size){

            this->OID = OID;
            data = std::vector<DType>();
            data.reserve(size);
            serialized = NULL;

        }

        /**
        * Destructor.
        */
        ~BasicArrayObject(){

            data.clear();
            if (serialized != NULL)
                delete[] serialized;
        }

        /**
        * @deprecated
        * @copydoc setOID(size_t OID).
        */
        void SetOID(size_t OID){

            setOID(OID);
        }


        /**
        * Sets the feature vector OID.
        * @param OID The OID of the feature vector.
        * Caution: The OID is not checked as unique.
        */
        void setOID(size_t OID){
            this->OID = OID;
        }

        /**
        * Gets the feature vector OID.
        * @return The feature vector OID.
        */
        size_t getOID() const{
            return OID;
        }

        /**
        * @brief This method is deprecated, but still usefull for arboretum
        * template functions.
        * @deprecated This method is deprecated. Use getOID() instead.
        * @copydoc getOID().
        */
        size_t GetOID() const{
            return getOID();
        }

        /**
        * Re-sizes a Basic Array Object.
        * All previous stored values are cleaned according to previous size.
        * @param size The new Basic Array size.
        */
        void resize(size_t size){

            DType aux2;
            std::vector<DType> aux;
            for (size_t x = 0; x < std::min((size_t) size, data.size()); x++){
                aux.push_back(data[x]);
            }
            data.clear();
            aux.clear();
            data = aux;
            for (size_t x = std::min((size_t) size, data.size()); x < size; x++){
                data.push_back(aux2);
            }
        }

        /**
        * Re-sizes a Basic Array Object.
        * All previous stored values are cleaned according to previous size.
        * New elements receives 'value' as the default value.
        * @param size The new Basic Array size.
        * @param value The default value to be set.
        */
        void resize(size_t size, DType value){

            std::vector<DType> aux;
            for (size_t x = 0; x < std::min((size_t) size, data.size()); x++){
                aux.push_back(data[x]);
            }
            data.clear();
            aux.clear();
            data = aux;
            for (size_t x = std::min((size_t) size, data.size()); x < size; x++){
                data.push_back(value);
            }
        }


        /**
        * @copydoc set(DType value).
        */
        void Set(DType value){
            set(value);
        }

        /**
        * Set the next (mandatory) position of the feature vector.
        * @param value The value to be pushed.
        */
        void add(DType value){
            set(value);
        }

        /**
        * @deprecated
        * @copydoc add(DType value).
        */
        void set(DType value){
            data.push_back(value);
        }

        /**
        * @deprecated
        * @copydoc set(size_t pos, DType value).
        */
        void Set(size_t pos, DType value){

            set(pos, value);
        }

        /**
        * Sets a specific value in a specific position.
        * @todo Throw a out-of-bounds exception.
        * @param pos The position of the insertion.
        * @param value The value to be pushed.
        */
        void set(size_t pos, DType value){

            if (pos == data.size()){
                data.push_back(value);
            } else {
                if (pos > data.size()){
                    DType aux;
                    for (int x = 0; x < (pos-data.size()); x++)
                        data.push_back(aux);
                    data[pos] = value;
                } else {
                    data[pos] = value;
                }
            }
        }

        /**
        * Gets the entire stored data.
        * @return The entire stored data.
        */
        std::vector<DType> getData(){

            return data;
        }

        /**
        * Overloaded operator allowing modifications.
        * @param idx The index to be queried.
        */
        DType& operator[] (size_t idx) {

            return data[idx];
        }

        /**
        * Overloaded operator that not allows modifications.
        * @param idx The index to be queried.
        */
        const DType& operator[] (size_t idx)  const{

            return data[idx];
        }

        /**
        * Get a value in a specific position.
        * @todo Throw a out-of-bounds exception.
        * @param idx The position value to be retrieved.
        * @return The value of the position idx.
        */
        DType *get(size_t idx){

            return (&data[idx]);
        }

        /**
        * @brief This method is deprecated, but still necessary for arboretum
        * template functions.
        * @deprecated This method is deprecated. Use get(size_t idx) instead.
        * @copydoc get(size_t idx).
        */
        DType *Get(size_t idx){

            return get(idx);
        }

        /**
        * Gets the number of elements in the feature vector.
        * @return The number of elements of the feature vector.
        */
        size_t getSize(){

            return data.size();
        }

        size_t size(){

            return getSize();
        }

        /**
        * @brief This method is deprecated, but still usefull for arboretum
        * template functions.
        * @deprecated This method is deprecated. Use getSize() instead.
        * @copydoc getSize().
        */
        size_t GetSize(){

            return getSize();
        }

        /**
        * @deprecated
        * @copydoc getObject().
        */
        BasicArrayObject<DType> GetObject(){

            return getObject();
        }

        /**
        * Return the instance of the current Basic Array Object.
        * @return The current instance of Basic Array Object.
        */
        BasicArrayObject<DType> getObject(){

            return this;
        }

        /**
        * Gets an instantied copy of the object.
        * @return A copy of the object.
        */
        BasicArrayObject<DType> *clone(){

            return new BasicArrayObject<DType>(getOID(), getData());
        }

        /**
        * @brief This method is deprecated, but still usefull for arboretum
        * template functions.
        * @deprecated This method is deprecated. Use clone() instead.
        * @copydoc clone().
        */
        BasicArrayObject<DType> *Clone(){

            return clone();
        }

        /**
        * Check if the obj is equal to the current object.
        * @param obj The object to be compared.
        * @return True if the objects are equal, else otherwise.
        */
        bool isEqual(BasicArrayObject<DType> *obj){

            if ((getOID() != obj->GetOID()) || (getSize() != obj->GetSize()))
                return false;

            for (int x = 0; x < getSize(); x++)
                if (*(this->get(x)) != *((DType*) (obj->Get(x))))
                    return false;

            return true;
        }

        /**
        * @brief This method is deprecated, but still usefull for arboretum
        * template functions.
        * @deprecated This method is deprecated. Use isEqual(stObject *obj) instead.
        * @copydoc isEqual(stObject *obj).
        */
        bool IsEqual(BasicArrayObject<DType> *obj){

            return isEqual(obj);
        }

        /**
        * Gets the size of the byte vector.
        * @return The size of the bytes vector.
        */
        size_t getSerializedSize()
        {

            if constexpr (std::is_same<DType, std::vector<char>>::value)
            {

                size_t total = 2*sizeof(size_t);

                for(size_t x = 0; x < getSize(); x++)
                {

                    total += sizeof(size_t) + data[x].size()*sizeof(char);

                }

                return total;

            }
            else
            {

                return (sizeof(size_t) + sizeof(size_t) +  (sizeof(DType) * data.size()));

            }

        }

        /**
        * @brief This method is deprecated, but still usefull for arboretum
        * template functions.
        * @deprecated This method is deprecated. Use getSerializedSize() instead.
        * @copydoc getSerializedSize().
        */
        size_t GetSerializedSize()
        {

            return getSerializedSize();

        }

        unsigned char *serialize(){

            if(serialized == NULL){

                if constexpr (std::is_same<DType, std::vector<char>>::value)
                {
                    serialized = serializeVecChar();
                }
                else
                {
                    serialized = serializeNumeric();
                }

            }

            return serialized;


        }

        /**
        * @brief This method is deprecated, but still usefull for arboretum
        * template functions.
        * @deprecated This method is deprecated. Use serialize() instead.
        * @copydoc serialize().
        */
        unsigned char *Serialize(){

            return serialize();
        }

        /**
        * Gets the equivalent byte vector of the object.
        * @return The equivalent string byte vector of the object.
        */
        std::string serializeToString(){

            std::string answer;

            serialize();
            for (int x = 0; x < getSerializedSize(); x++){
                answer += serialized[x];
            }

            return answer;
        }

        void unserialize(unsigned char *dataIn, size_t dataSize = 0)
        {

            if constexpr (std::is_same<DType, std::vector<char>>::value)
            {
                unserializeVecChar(dataIn, dataSize);
            }
            else
            {
                unserializeNumeric(dataIn, dataSize);
            }

        }

        /**
        * @brief This method is deprecated, but still usefull for arboretum
        * template functions.
        * @deprecated This method is deprecated.
        * Use unserialize(const unsigned char *dataIn, size_t dataSize) instead.
        * @copydoc unserialize(const unsigned char *dataIn, size_t dataSize).
        */
        void Unserialize(unsigned char *dataIn, size_t dataSize = 0){

            unserialize(dataIn, dataSize);
        }

        /**
        * @copydoc unserialize(const unsigned char *dataIn, size_t dataSize).
        */
        void unserializeFromString(std::string dataIn){

            size_t size_vector;
            DType *d;

            // This is the reverse of Serialize(). So the steps are similar.
            // Remember, the format of the serizalized object is
            // +----------+-----------+--------+
            // | OID | Size | Vector Data []   |
            // +----------+-----------+--------+

            memcpy(&OID, dataIn.c_str(), sizeof(size_t));
            memcpy(&size_vector, dataIn.c_str() + sizeof(size_t), sizeof(size_t));
            d = new DType[size_vector];

            for (size_t x = 0; x < size_vector; x++){
                memcpy(&d[x], dataIn.c_str() + sizeof(size_t) + sizeof(size_t) + (sizeof(DType)*x), sizeof(DType));
            }

            data.clear();

            for (size_t x = 0; x < size_vector; x++){
                data.push_back(d[x]);
            }

            // Since we have changed the object contents, we must invalidate the old
            // serialized version if it exists. In fact we, may copy the given serialized
            // version of tbe new object to the buffer but we don't want to waste memory.
            if (serialized != NULL){
                delete [] serialized;
                serialized = NULL;
            }//end if
        }

        static std::string base64Chars(){

            return "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-";
        }



        static bool isBase64(char c){

            return (isalnum(c) || (c == '+') || (c == '-'));
        }


        static std::string toBase64(std::string input){

            unsigned char * bytes_to_encode = (unsigned char *) input.c_str();
            int in_len = input.size();

            std::string ret;
            int i = 0;
            int j = 0;
            unsigned char char_array_3[3];
            unsigned char char_array_4[4];

            while (in_len--) {
                char_array_3[i++] = *(bytes_to_encode++);
                if (i == 3) {
                    char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
                    char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
                    char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
                    char_array_4[3] = char_array_3[2] & 0x3f;

                    for(i = 0; (i <4) ; i++)
                        ret += base64Chars()[char_array_4[i]];
                    i = 0;
                }
            }

            if (i)
            {
                for(j = i; j < 3; j++)
                    char_array_3[j] = '\0';

                char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
                char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
                char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
                char_array_4[3] = char_array_3[2] & 0x3f;

                for (j = 0; (j < i + 1); j++)
                    ret += base64Chars()[char_array_4[j]];

                while((i++ < 3))
                    ret += '=';

            }

            return ret;
        }

        static std::string fromHexaDecimal(std::string input){

            static const char* const lut = "0123456789ABCDEF";
            size_t len = input.length();
            //if (len & 1) throw std::invalid_argument("odd length");

            std::string output;
            output.reserve(len / 2);
            for (size_t i = 0; i < len; i += 2)
            {
                char a = (char) input[i];
                const char* p = std::lower_bound(lut, lut + 16, a);
                //if (*p != a) throw std::invalid_argument("not a hex digit");

                char b = input[i + 1];
                const char* q = std::lower_bound(lut, lut + 16, b);
                //if (*q != b) throw std::invalid_argument("not a hex digit");

                output.push_back(((p - lut) << 4) | (q - lut));
            }

            return output;
        }

        static std::string fromBase64(std::string input){

            size_t in_len = input.size();
            size_t i = 0;
            size_t j = 0;
            int in_ = 0;
            unsigned char char_array_4[4], char_array_3[3];
            std::string ret;

            while (in_len-- && ( input[in_] != '=') && isBase64(input[in_])) {
                char_array_4[i++] = input[in_]; in_++;
                if (i ==4) {
                    for (i = 0; i <4; i++)
                        char_array_4[i] = static_cast<unsigned char>(base64Chars().find(char_array_4[i]));

                    char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
                    char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
                    char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

                    for (i = 0; (i < 3); i++)
                        ret += char_array_3[i];
                    i = 0;
                }
            }

            if (i) {
                for (j = i; j <4; j++)
                    char_array_4[j] = 0;

                for (j = 0; j <4; j++)
                    char_array_4[j] = static_cast<unsigned char>(base64Chars().find(char_array_4[j]));

                char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
                char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
                char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

                for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
            }

            return ret;
        }

        std::string toString(std::string separator)
        {

            std::stringstream str;

            if constexpr (std::is_same<DType, std::vector<char>>::value)
            {

                str << "[";

                for(size_t x = 0; x < getSize(); x++)
                {

                    size_t vecSize = data[x].size();

                    str << "[";

                    for(size_t y = 0; y < vecSize; y++)
                    {

                        str << data[x][y];

                        if(y != (vecSize-1))
                        {

                            str << separator;

                        }
                        else
                        {

                            str << "]";

                        }

                    }

                    str << "]";

                    if(x != (getSize()-1))
                    {

                        str << ", ";

                    }

                }

                str << "]";

            }
            else
            {

                str << "[";

                for(size_t i = 0; i < size(); i++)
                {

                    str << std::to_string(*get(i));

                    if(i != (size()-1))
                    {

                        str << separator;

                    }
                    else
                    {

                        str << "]";

                    }

                }

            }

            return str.str();

        }

        std::string toString()
        {

            return toString(" ");

        }

        std::string toStringWithOID(std::string separator=",")
        {

            return std::to_string(this->OID) + " -> " + toString(separator);

        }        

};


typedef BasicArrayObject<double> FeatureVector;
typedef std::vector<FeatureVector> FeatureVectorList;

#endif
