
#ifndef ALGOAEx2_
#define ALGOAEx2_
// Generated by AlexConf: do not edit
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>
#include <TH1F.h>
#include <TH2F.h>
#include <alex/IAlgorithm.h>
namespace alex {
  class AEx2: public IAlgorithm {
  public:
    AEx2();
    ~AEx2(){}
    bool Init() ;
    bool Execute() ;
    bool End() ;
    std::string  Name() const {return fName;}
    void SetName(std::string name) {fName = name;}
  private:
    std::string fName;
    std::string fFileName;
    int fNEvents;
    double fBz;
    std::vector<double> fBfield;
    std::vector<int> fBVector;
    TH1F* fH1_Y;
    TH2F* fH2_YZ;
  };
}
#endif 