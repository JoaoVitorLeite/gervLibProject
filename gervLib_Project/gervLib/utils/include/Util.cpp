#include "Util.h"


size_t* uniqueRandomNumber(size_t start, size_t end, size_t nNumber, size_t seed)
{

   if(nNumber > (end-start+1)) throw std::invalid_argument("Too many numbers !_!");
   std::unordered_set<size_t> set;
   srand(seed);

   while (set.size() < nNumber)
   {

       size_t randomNumber = rand()%end + start;
       if(set.find(randomNumber) == set.end()) set.insert(randomNumber);

   }

   size_t* ans = new size_t[nNumber];

   size_t x = 0;
   for(size_t i : set)
       ans[x++] = i;

   set.clear();

   return ans;

}


size_t* randomNumber(size_t start, size_t end, size_t nNumber, size_t seed)
{

   size_t* ans = new size_t[nNumber];

   srand(seed);

   for(size_t x = 0; x < nNumber; x++)
       ans[x] = rand()%end + start;

   return ans;

}


size_t* shuffleIndex(size_t start, size_t end, size_t seed)
{

   srand(seed);

   size_t* ans = new size_t[end-start+1];
   size_t y, temp;

   for(size_t x = 0; x < (end-start+1); x++) ans[x] = x;

   for(size_t x = end - 1; x > 0; x--)
   {

      y = rand() % end + start;

      temp = ans[x];
      ans[x] = ans[y];
      ans[y] = temp;

   }

   return ans;

}



std::string selectFolderName(std::string directory, std::string startwith)
{

   std::vector<std::string> foldersName;
   std::string ans;
   DIR *dir;
   struct dirent *ent;
   if((dir = opendir(directory.c_str())) != NULL)
   {

       while ((ent = readdir(dir)) != NULL)
       {

           if(ent->d_type == DT_DIR)
           {

               std::string name = ent->d_name;

               if((name.rfind(startwith, 0) == 0) && (name != ".") && (name != ".."))
               {

                   foldersName.push_back(ent->d_name);

               }

           }

       }

       closedir(dir);

   }
   else
   {

       std::cout << std::filesystem::current_path() << std::endl;
       throw std::invalid_argument("Could not open directory");

   }


   if(foldersName.size() > 0)
   {

       std::sort(foldersName.begin(), foldersName.end(), [](const std::string & a, const std::string & b)
       {

           std::string number_a = std::regex_replace(a, std::regex("[^0-9]*"), std::string("$1")), number_b = std::regex_replace(b, std::regex("[^0-9]*"), std::string("$1"));

           return std::stoi(number_a) < std::stoi(number_b);

       });

       int num = std::stoi(std::regex_replace(foldersName.back(), std::regex("[^0-9]*"), std::string("$1"))) + 1;
       ans = startwith + "_" + std::to_string(num);

   }
   else
   {

       ans = startwith + "_1";

   }

   return ans;

}


