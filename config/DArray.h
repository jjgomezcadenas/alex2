#ifndef DARRAY_H
#define DARRAY_H

#include <cstdlib>
#include <string>
#include <alex/IData.h>


//Describes an array.
//   <Array>
    //   <name>P</name>
    //   <dataType>double</dataType>
    //   <dim>3</dim>
    //   <value>0.1 1.0 3</value>
    // </Array>

namespace alex {

	class DArray : public IData
	{
		public:
			DArray() {};
			virtual ~DArray(){};
		
			std::string Name() const {return fName;}
			std::string DataType() const {return fDataType;}
			int Dim() const {return fDim;}
			std::string Value() const {return fValue;}
			std::vector<std::string>  Split() ;
			std::string Serialize() const ;
	
			void SetData(std::string name, std::string dataType, 
									std::string dim, std::string value)
			{fName = name; fDim = atoi(dim.c_str()); fDataType=dataType; fValue=value;}
		
		private:
			std::string fName;
			std::string fDataType;
			std::string fValue;
			int fDim;
	};
}
#endif 