#include "LogUtil.h"

using std::string; 
using std::cout; 
using std::endl; 
using std::ostream;

namespace alex {

  void InitLogger()
  {
    log4cpp::Appender *appender1 = new log4cpp::OstreamAppender("console", &std::cout);
    appender1->setLayout(new log4cpp::BasicLayout());
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog.addAppender(appender1);
    klog.setPriority(log4cpp::Priority::DEBUG);
  }

  void SetDebugLevel(string debugLevel) 
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG 
          << "In DebugSvc::Configure SetLogger to " << debugLevel;

    if(debugLevel == "FATAL")
      klog.setPriority(log4cpp::Priority::FATAL);
    else if (debugLevel == "ERROR")
      klog.setPriority(log4cpp::Priority::ERROR);
   else if (debugLevel == "WARN")
      klog.setPriority(log4cpp::Priority::WARN);
    else if (debugLevel == "NOTICE")
      klog.setPriority(log4cpp::Priority::NOTICE);
    else if (debugLevel == "INFO")
      klog.setPriority(log4cpp::Priority::INFO);
    else
      klog.setPriority(log4cpp::Priority::DEBUG);
  }

}
