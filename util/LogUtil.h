#ifndef LOGUTIL_
#define LOGUTIL_
/*
 LOGUTIL: Utilities for Logging

 JJGC, July, 2014.
*/


#include "log4cpp/Category.hh"
#include "log4cpp/Appender.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/OstreamAppender.hh"
#include "log4cpp/Layout.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/Priority.hh"
#include <string>

namespace alex {
	void InitLogger();
	void SetDebugLevel(std::string debugLevel);

}
#endif