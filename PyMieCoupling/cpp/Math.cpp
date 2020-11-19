
#include <vector>
#include <complex>


typedef std::vector<std::complex<double>> iVec;
typedef std::complex<double> complex128;


std::vector<double> linespace(double start, double end, int N)
{
std::vector<double> vector  = std::vector<double>(N);
double delta = (end-start)/N;

for (int i = 0; i < N; i++)
  {
    vector[i] = delta*i;
  }
  return vector;
}



double Sum(std::vector<double> vector)
{
   int N = vector.size();
   double sum = 0.;
   for (int i = 0; i < N; i++)
   {
     sum += vector[i];
   }
   return sum;
}



complex128 Sum(std::vector<complex128> vector)
{
   int N = vector.size();
   complex128 sum = 0.;
   for (int i = 0; i < N; i++)
   {
     sum += vector[i];
   }
   return sum;
}



std::tuple<std::vector<double> , std::vector<double>> Arrange(double start, double stop)
{
  std::vector<double> Vec0;
  std::vector<double> Vec1;
  for (double i = start; i < stop; i++)
  {
    Vec0.push_back(i);
    Vec1.push_back( ( 2 * (i + 1) + 1) / ( (i + 1) * (i + 2) ) ) ;
  }
  return std::make_tuple(Vec0, Vec1);
}






void print(std::vector <double> const &a) {
   std::cout << "The vector elements are : ";
   for(long unsigned int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
}



void print(std::vector <complex128> const &a) {
   std::cout << "The vector elements are : ";
   for(long unsigned int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
}












// -
