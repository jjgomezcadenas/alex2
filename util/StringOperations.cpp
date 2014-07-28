
#include <alex/StringOperations.h>
using namespace std;

namespace alex {

//--------------------------------------------------------------------
  string MergeStrings(string s1,string s2)
//--------------------------------------------------------------------
  {
    ostringstream s;
    s << s1 << s2 ;
    return s.str();
  }
//--------------------------------------------------------------------
  string PathFromStrings(string path,string name)
//--------------------------------------------------------------------
  {
    ostringstream s;
    s << path << "/" << name ;
    return s.str();
  }

//--------------------------------------------------------------------
  std::vector<std::string> SplitString(string str)
//--------------------------------------------------------------------
  {
    std::istringstream buf(str);
    std::istream_iterator<std::string> beg(buf), end;

    std::vector<std::string> tokens(beg, end); // done!

    // for(auto& s: tokens)
    //     std::cout << '"' << s << '"' << '\n';
    return tokens;
  }


} // namespace 

