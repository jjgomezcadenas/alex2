
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



} // namespace 

