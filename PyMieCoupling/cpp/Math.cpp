
#include <vector>
#include <complex>
#include <tuple>

typedef std::vector<std::complex<double>> iVec;
typedef std::complex<double> complex128;


std::vector<double> linespace(const double start,
                              const double end,
                              const long unsigned int N)
{
std::vector<double> vector = std::vector<double>(N);

const double delta = (end-start)/N;

for (long unsigned int i = 0; i < N; i++)
  {
    vector[i] = delta*i;
  }
  return vector;
}



complex128 Sum(complex128* vector, const long unsigned int N)
{

   complex128 sum = 0.;
   for (long unsigned int i = 0; i < N; i++)
   {
     sum += vector[i];
   }
   return sum;
}




double Sum(const std::vector<double>* vector)
{
   const long unsigned int N = vector->size();
   double sum = 0.;
   for (long unsigned int i = 0; i < N; i++)
   {
     sum += (*vector)[i];
   }
   return sum;
}






complex128 Sum(const std::vector<complex128>* vector)
{
   const long unsigned int N = vector->size();
   complex128 sum = 0.;
   for (long unsigned int i = 0; i < N; i++)
   {
     sum += (*vector)[i];
   }
   return sum;
}




void print(std::vector<double>* a) {
   std::cout << "The vector elements are : ";
   for(long unsigned int i=0; i < a->size(); i++)
   std::cout << a->at(i) << ' ';
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



void print(std::vector <complex128>* a) {
   std::cout << "The vector elements are : ";
   for(long unsigned int i=0; i < a->size(); i++){
   std::cout << a->at(i) << ' ';}
}



std::pair<std::vector<double>* , std::vector<double>*> Arrange(const double start,
                                                              const double stop)
{
  std::vector<double> *Vec0 = new std::vector<double>;
  std::vector<double> *Vec1 = new std::vector<double>;
  for (double i = start; i < stop; i++)
  {
    Vec0->push_back(i);
    Vec1->push_back( ( 2 * (i) + 1) / ( (i) * (i + 1) ) ) ;
  }
  return std::make_pair(Vec0, Vec1);
}














// -
