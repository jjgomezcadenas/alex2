#ifndef IALGORITHM_H
#define IALGORITHM_H

#include <string>

//Interface for Alex algorithms. 

namespace alex {

	class IAlgorithm {
	public:
	virtual ~IAlgorithm() {};
	virtual bool Init() = 0;
	virtual bool Execute() = 0;
	virtual bool End() = 0;
	virtual std::string  Name() const = 0;
	virtual void SetName(std::string) = 0;
	};
}
#endif //IALGORITHM_H