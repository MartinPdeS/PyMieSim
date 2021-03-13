#include <iostream>
#include <complex>
#include <vector>
#include <tuple>

typedef std::vector<double> Vec;
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;

#define PI (double)3.14159265358979323846264338

std::vector<double>
linespace(const double start,
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



complex128
Sum(complex128* vector, const long unsigned int N)
{
   complex128 sum = 0.;
   for (long unsigned int i = 0; i < N; i++)
   {
     sum += vector[i];
   }
   return sum;
}




double
Sum(const std::vector<double>* vector)
{
   const long unsigned int N = vector->size();
   double sum = 0.;
   for (long unsigned int i = 0; i < N; i++)
   {
     sum += (*vector)[i];
   }
   return sum;
}


template <class T>
T Sum(const std::vector<T>* vector)
{
   const long unsigned int N = vector->size();
   T sum = 0.;
   for (long unsigned int i = 0; i < N; i++)
   {
     sum += (*vector)[i];
   }
   return sum;
}


template <class T>
void printVector(std::vector <T> const &a) {
   std::cout << "The vector elements are : ";
   for(long unsigned int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
}



std::pair<std::vector<double>* , std::vector<double>*>
Arrange(const double start,
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


template<typename T>
Vec
Linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);

  return linspaced;
}


std::tuple<Vec, Vec>
Meshgrid(int num)
{
  Vec phi = Linspace(-PI/2, PI/2, num);
  Vec theta = Linspace(-PI, PI, num);

  Vec Phi, Theta;

  for(int i=0; i < num; ++i)
  {
    Phi.insert( end(Phi), begin(phi), end(phi) );
    for(int j=0; j < num; ++j)
    {
      Theta.push_back(theta[i]);
    }
  }


  return std::make_tuple(Phi, Theta);
}



std::tuple<Vec, Vec, Vec>
FullMesh(int num)
{
  Vec Phi, Theta;
  std::tie(Phi, Theta) = Meshgrid(num);
  Vec R(num*num, 10);

  return std::make_tuple(R, Phi, Theta);
}



// -
