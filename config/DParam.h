#ifndef DPARAM_H
#define DPARAM_H

#include <string>
#include <alex/IData.h>
#include <cstdlib>


//Describes a parameter.
//   <Param>
    //   <name>dataPath</name>
    //   <dataType>string</dataType>
    //   <value>/Users/jjgomezcadenas/Development/NEXT/DATA</value>
    // </Param> 

namespace alex {

	class DParam : public IData
	{
		public:
				DParam() {};
				virtual ~DParam(){};
		
				std::string Name() const {return fName;}
				std::string DataType() const {return fDataType;}
				double GetValueAsDouble() const {return atof(fValue.c_str());}
				int GetValueAsInt() const {return atoi(fValue.c_str());}
				std::string Value() const {return fValue;}
				std::string Serialize() const ;
	
				void SetData(std::string name, std::string dataType, std::string value)
				{fName = name; fDataType=dataType; fValue=value;}
		
		private:
				std::string fName;
				std::string fDataType;
				std::string fValue;
	};
}
#endif //IDATA_H