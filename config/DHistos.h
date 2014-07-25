#ifndef DHISTOS_H
#define DHISTOS_H

#include <cstdlib>
#include <string>
#include <alex/IData.h>


//Describes h1 and h2 histos.
//   <H1D>
    //     <name>fHY</name>
    //     <title>Y distribution</title>
    //     <nbinsx>50</nbinsx>
    //     <xlow>-10.0</xlow>
    //     <xup>10.0</xup>
    //   </H1D>

//   <H2D>
    //     <name>fHYZ</name>
    //     <title>Y vs Z distribution</title>
    //     <nbinsx>10</nbinsx>
    //     <nbinsy>10</nbinsy>
    //     <xlow>-10.0</xlow>
    //     <xup>10.0</xup>
    //     <ylow>-10.</ylow>
    //     <yup>10.0</yup>
    //   </H2D>
    // </ToyAnalysis2>

namespace alex {

	class DHisto 
	{
		public:
			DHisto() {};
			virtual ~DHisto(){};
		
			std::string Name() const {return fName;}
			std::string Title() const {return fTitle;}

			void SetData(std::string name, std::string title)
			{fName = name; fTitle = title;}


		private:
			std::string fName;
			std::string fTitle;
	};
	class DH1: public IData, public DHisto
	{
		public:
			DH1() {};
			virtual ~DH1(){};
		
			int Nbinsx() const {return fNbinsx;}
			double Xlow() const {return fXlow;}
			double Xup() const {return fXup;}

			std::string Serialize() const ;
	
			void SetData(std::string name, std::string title, 
						std::string nbinsx, std::string xlow,std::string xup )
				{
					DHisto::SetData(name,title); 
					fNbinsx=atoi(nbinsx.c_str()); 
					fXlow=atof(xlow.c_str()); fXup= atof(xup.c_str());
				}
		
		private:
			int fNbinsx;
			double fXlow;
			double fXup;
	};
	class DH2 : public DH1
	{
		public:
			DH2() {};
			virtual ~DH2(){};
		
			int Nbinsy() const {return fNbinsy;}
			double Ylow() const {return fYlow;}
			double Yup() const {return fYup;}

			std::string Serialize() const ;
	
			void SetData(std::string name, std::string title, 
						std::string nbinsx, std::string xlow,std::string xup,
						std::string nbinsy, std::string ylow,std::string yup)

				{
					DH1::SetData(name, title, nbinsx, xlow, xup);
					fNbinsy=atoi(nbinsy.c_str()); 
					fYlow=atof(ylow.c_str()); fYup= atof(yup.c_str());
				}
		
		private:
			int fNbinsy;
			double fYlow;
			double fYup;
	};
}
#endif 