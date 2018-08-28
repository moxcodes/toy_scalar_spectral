#include <vector>
#include <boost/math/special_functions/legendre.hpp>
#include "scalarFunction.hpp"



class scalarHistory
{
public:
  int spectralPoints=10;
  std::vector<scalarFunction::scalarFunction> data;
 
}
