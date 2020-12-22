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


std::map<int, double> prob(NumericVector x){
  int n = x.size();
  
  std::map<int, int> counts; 
  for (unsigned int i = 0; i < n; i++)
  {
    counts[x[i]]++;
  }
  
  std::map<int, double> probs; 
  
  typename std::map<int, int>::iterator it;
  it = counts.begin();
  while(it != counts.end()){
    int count = it->second;
    probs[it->first] = (double)count/n;
    it++; 
  }
  
  return probs; 
}


// [[Rcpp::export]]
// Sample Entropy:
// https://journals.physiology.org/doi/full/10.1152/ajpheart.2000.278.6.h2039
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
  
  std::map<int, double> probs = prob(x);
  
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


std::vector<int> biDerivative(NumericVector x)
{
  // Could also be done with bit-wise operators for potential speed up
  int n = x.size();

  auto XOR = [](int x_1, int x_2)->int {int d_i = ((x_1+x_2)==1) ? 1 : 0;
                                        return d_i;};
  
  std::vector<int> d;
  for(int i = 0; i<(n-1); i++)
  {
    int d_i = XOR(x[i], x[i+1]);
    d.push_back(d_i);
  }
  
  return d; 
}


 // [[Rcpp::export]]
 // Binary Entropy: https://arxiv.org/ftp/arxiv/papers/1305/1305.0954.pdf
double biEn(NumericVector x, bool tresBin){
  
  int n = x.size();
  double sc_f = 1; // Scaling factor initialization
  std::function<double(double, int)> en_k; // Entropy of the k-th derivative
  
  if(!tresBin)
  {
    // scaling factor
    sc_f = (1/(pow(2,n-1)-1));
  
    // biEntropy term of k-th derivation
    en_k = [](double p_1, int k)->double {return (-p_1*log2(p_1)-(1-p_1)*log2(1-p_1))*pow(2,k);};
  }
  else
  {
    // scaling factor
    double norm = 0.0;
    for(unsigned int k=0;k<=(n-2);k++)
    {
      norm += log2(k+2);
    }
    sc_f = 1/norm;
    
    // Tres biEntropy term of k-th derivation
    en_k = [](double p_1, int k)->double {return (-p_1*log2(p_1)-(1-p_1)*log2(1-p_1))*log2(k+2);};
  }

  double biEn = 0.0; 
  NumericVector dk = x; 
  for(unsigned int k = 0; k<=(n-2);k++)
  {
    std::map<int, double> probs = prob(dk);
    double p_1 = probs[1];
    
    if((p_1==0.0) || (p_1==1.0))
    {
      break;
    }
    
    // add entropy term of k-th derivation
    biEn += en_k(p_1, k);
    
    dk = biDerivative(dk);
  }

  // scale Binary Entropy
  return sc_f*biEn; 
}

// [[Rcpp::export]]

double msEn(NumericVector x, int nScales)
{
  
  
  
}

/*** R
x <- c(1,0,0,1)

biEntropy <- biEn(x, tresBin=T)
*/