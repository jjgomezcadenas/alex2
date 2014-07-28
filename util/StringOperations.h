#ifndef STRINGSOP_H
#define STRINGSOP_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <numeric>

#include <iterator>


namespace alex {

  std::string MergeStrings(std::string s1,std::string s2);
  std::string PathFromStrings(std::string path,std::string name);
  std::vector<std::string> SplitString(std::string str);

//--------------------------------------------------------------------
  template< typename T>
  std::string SerializeVectorInList(std::vector<T> v)
//--------------------------------------------------------------------
  {
    std::ostringstream s;
    s << std::endl;
    for (auto elem : v)
      s<< elem << std::endl << std::ends;

    return s.str();
  }

} // namespace 

#endif 
