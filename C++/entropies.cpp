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


std::map<int, double> prob(NumericVector x, bool discretize){
  int n = x.size();
  
  
  // Discretize: Scale and tranform to int to avoid the use of doubles as keys.
  if(discretize)
  {
    std::transform(x.begin(), x.end(), x.begin(), [](double d) { return (int)(d*100);});
  }
  
  std::map<int, int> counts; 
  for (unsigned int i = 0; i < n; i++)
  {
    counts[x[i]]++;
  }
  
  std::map<int, double> probs; 
  
  std::map<int, int>::iterator it;
  it = counts.begin();
  while(it != counts.end()){
    int count = it->second;
    probs[it->first] = (double)count/n;
    it++; 
  }
  
  return probs; 
}


// [[Rcpp::export]]
// Shannon Entropy: 
// http://people.math.harvard.edu/~ctm/home/text/others/shannon/entropy/entropy.pdf
double shannonEnDiscrete(NumericVector x){
  double entropy=0;
  int n = x.size();
  
  if(n<2)
  {
    return NAN; 
  }
  
  // Code adapted from: 
  // https://stackoverflow.com/questions/20965960/shannon-entropy
  
  std::map<int, double> probs = prob(x, true);
  
  typename std::map<int, double>::iterator it;
  // Calculate Shannon Entropy
  it = probs.begin();
  while(it != probs.end()){
    double p_x = it->second;
    if (p_x>0) entropy-=p_x*log2(p_x);
    it++;
  }
  return entropy;
}