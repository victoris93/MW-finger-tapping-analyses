#include <Rcpp.h>
using namespace Rcpp;


/*
 * This function measures approximate entropy (Pincus + Kalman, 1997, PNAS).
 * We count only exact matches, so d(x(i),y(i)) is either 0 or 1 (0 for exact matches),
 * 1 otherwise. Therefore, the Cm count the number of exact matches of sub-sequences of
 * length m (including the matching sequence itself).
 * 
 * Output is: ApEn(0), ..., ApEn(maxM)
 */

// [[Rcpp::export]]
NumericVector apen_int(IntegerVector x, int maxM){
  NumericVector Phi(maxM); // from 1, ..., maxM+1
  NumericVector ApEn(maxM+1); // from 0, ..., maxM
  int N=x.size();
  IntegerVector C(N);
  
  for(int i=0; i<N; i++){
    if(IntegerVector::is_na(x[i])){
      Rprintf("WARN\n");
    }
  }
  
  int d;
  for(int m=1; m<=maxM+1; m++){ 
    Phi[m-1]=0.0;
    for(int i=0; i<N-m+1; i++){
      C[i]=0;
      for(int j=0; j<N-m+1; j++){
        
        // check if subseq of length m starting at i is identical to that at j
        d=0;
        for(int p=0; p<m; p++){
          if(x[i+p]!=x[j+p]){
            d=1; break;
          }
        }
        C[i] += int(d==0);
      }
      Phi[m-1] += log(float(C[i])/float(N-m+1));  // this is base-e log
    }
    Phi[m-1] *= 1./(N-m+1);
  }
  ApEn[0] = -Phi[0];
  for(int i=1; i<=maxM; i++ ){
    ApEn[i] = Phi[i-1] - Phi[i];  // ApEn(1) = Phi(1)-Phi(2)
  }
  //Rprintf("%i\n", maxM);
  //Rprintf("%f", log(2));
  return ApEn;
}