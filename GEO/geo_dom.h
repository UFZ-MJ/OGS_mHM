/**************************************************************************
GeoLib - Object: Domain
Task: 
Programing:
09/2003 OK Implementation
09/2005 CC GeoLib2
**************************************************************************/

#ifndef gs_dom_INC

#define gs_dom_INC
  /* Schutz gegen mehrfaches Einfuegen */

// C++ STL
#include <string>
#include <vector>
#include <cstdio>

/* Class */
/*---------------------------------------------------------------*/

class CGLDomain
{
  private:
    std::string name;
    long Insert(CGLDomain *);
    std::vector<CGLDomain*> GetVector(void);
  public:
    double x_min,x_max;
    double y_min,y_max;
    double z_min,z_max;
    // constructor
    CGLDomain(void);
    // destructor
    ~CGLDomain(void);
    int Read(char *,FILE *);
    CGLDomain* Get(std::string);
};
extern std::vector<CGLDomain*> domain_vector;
extern int GEOReadDomain(char*,int,FILE*);
#endif
