#include <Rcpp.h>
using namespace Rcpp;

#include <math.h>
#include <cstdlib>
#include <functional>
#include <algorithm>

// [[Rcpp::export]]
double mhl(NumericVector x){
  
  int r = x.size();
  int L = x.size()-1;
  
  int norm_factor = (r*(r+1)/2); 
  double met_sub_frac = 0.0;
  
  for (unsigned int l = 0; l <= L; l++) {
      
      int n_met_cpgs = (l+1);
      int n_met_sub = 0; 
      
      // check all substrings of length l 
      for(unsigned int i=0; i<(r-l); i++) {
        
          NumericVector substring = x[Rcpp::Range(i, i+l)];
          
          int n_met = std::accumulate(substring.begin(), 
                                      substring.end(), 0);
          if(n_met==n_met_cpgs)
          {
            n_met_sub +=1;
          }
      }
      
      met_sub_frac += ((double)n_met_cpgs*(double)n_met_sub/(double)(r-l));
  }
  
  double mhl=(double)met_sub_frac/norm_factor;
  
  return mhl;
}
