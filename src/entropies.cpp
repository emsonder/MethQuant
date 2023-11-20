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
  if(!discretize)
  {
    std::transform(x.begin(), x.end(), x.begin(), [](double d) { return (int)(d*100000);});
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

//' Shannon Entropy:
//'
//' Rcpp Implementation of a Shannon Entropy
//' http://people.math.harvard.edu/~ctm/home/text/others/shannon/entropy/entropy.pdf
//'
//' @name conditionalEntropy
//' @param x Numeric vector
//' @param y Numeric vector
//' @return conditional entropy H(X|Y)
//'
//' @export
// [[Rcpp::export]]
double shannonEnDiscrete(NumericVector x, bool normalize, bool discretize){
  double entropy=0;
  int n = x.size();

  if(n<2)
  {
    return NAN;
  }

  // Code adapted from:
  // https://stackoverflow.com/questions/20965960/shannon-entropy
  std::map<int, double> probs = prob(x, discretize);

  typename std::map<int, double>::iterator it;
  // Calculate Shannon Entropy
  it = probs.begin();
  while(it != probs.end()){
    double p_x = it->second;
    if (p_x>0) entropy-=p_x*log(p_x);
    it++;
  }

  if(normalize)
  {
    entropy=entropy/log(probs.size());
  }

  return entropy;
}

//' Sample Entropy
//'
//' Rcpp Implementation of a Sample Entropy (// https://journals.physiology.org/doi/full/10.1152/ajpheart.2000.278.6.h2039)
//' Implementation adapted from: http://blog.schochastics.net/post/sample-entropy-with-rcpp/
//'
//' @name conditionalEntropy
//' @param x Numeric vector
//' @param y Numeric vector
//' @return conditional entropy H(X|Y)
//'
//' @export
// [[Rcpp::export]]
double sampleEn(NumericVector x, int m, double r){
  int N = x.size();

  if(N<2)
  {
    return R_NaN;
  }

  // Code adapted from:
  // http://blog.schochastics.net/post/sample-entropy-with-rcpp/

  int cm = 0, cm_1 = 0;
  double tol = 0.0;

  double sd = calcSD(x);

  // tolerance
  // tol = sd * r;
  tol = r;

  for (unsigned int i = 0; i < N - m; i++) {
    for (unsigned int j = i + 1; j < N-m; j++) {
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
  else if(cm>0)
    return R_PosInf; // or NAN
  else
    return R_NaN; // This case cannot happen!
}
