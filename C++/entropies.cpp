#include <Rcpp.h>
using namespace Rcpp;

#include <math.h>
#include <cstdlib>
#include <functional>
#include <algorithm>


double calcMean(NumericVector x) {
  double mean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
  
  return mean; 
}


double calcSD(NumericVector x)
{
  double mean = calcMean(x);
  int n = x.size();
  
  double sum_diff = 0.0; 
  for(unsigned int i = 0; i < n; ++i)
  {
    sum_diff += pow(x[i] - mean, 2);
  }
  
  double sq_sd = (double)1/(n-1)*sum_diff;
  
  return sqrt(sq_sd);
}


// [[Rcpp::export]]
double sampleEn(NumericVector x, int m, double r){
  int N = x.size();
  
  if(N<2)
  {
    return NAN; 
  }
  
  // Code adapted from: 
  // http://blog.schochastics.net/post/sample-entropy-with-rcpp/
  
  int cm = 0, cm_1 = 0;
  double tol = 0.0;
  
  double sd = calcSD(x);
  
  // tolerance
  tol = sd * r;
  
  for (unsigned int i = 0; i < N - m; i++) {
    for (unsigned int j = i + 1; j < N - m; j++) {      
      bool eq = true;
      
      // Chebyshev distance criteria
      for (unsigned int k = 0; k < m; k++) {
        if (abs(x[i+k] - x[j+k]) > tol) {
          eq = false;
          break;
        }
      }
      if (eq) cm++;
      
      // check for length m+1
      int k = m;
      if (eq && abs(x[i+k] - x[j+k]) <= tol)
        cm_1++;
    }
  }
  
  if (cm > 0 && cm_1 > 0)
    return log((double)cm / (double)cm_1);
  else
    return 0.0; 
}


// [[Rcpp::export]]
double shannonEnDiscrete(NumericVector x){
  double entropy=0;
  int n = x.size();
  
  if(n<2)
  {
    return NAN; 
  }
  
  // Code adapted from: 
  // https://stackoverflow.com/questions/20965960/shannon-entropy
  
  // For discrete data (few different doubles)
  std::map<double, int> counts;
  for (unsigned int i = 0; i < n; i++) {
    counts[x[i]]++; 
  }
  
  typename std::map<double, int>::iterator it;
  it = counts.begin();
  while(it != counts.end()){
    double p_x = (double)it->second/n; 
    if (p_x>0) entropy-=p_x*log2(p_x);
    it++;
  }
  return entropy;
}
