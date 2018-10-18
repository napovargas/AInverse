#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <time.h>
#include <algorithm>
#include <fstream>
#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
void AinverseCpp(arma::ivec ivecs, arma::ivec ivecd, arma::sword ina){
  std::ofstream iun22;
  double      dz          = 0.0;
  //double      zero        = 1e-12;
  //double      deta        = 0.0;
  double      den         = 0.0;
  double      rden        = 0.0;
  double      r25         = 0.0;
  double      r50         = 0.0;
  double      cx          = 0.0;
  double      summ        = 0.0;
  //arma::uword iun22       = 22;
  //arma::uword icnta       = 0;
  arma::uword maxrow      = 10;
  arma::uword maxan       = 1e5;
  arma::sword k           = 0;
  arma::sword maxgh       = -1;
  arma::uword maxg1       = maxgh + 1;
  arma::uword maxent      = maxan*maxrow;
  arma::sword i, ii       = 0;
  arma::sword isir        = 0;
  arma::sword idam        = 0;
  arma::sword j, jj       = 0;
  arma::ivec  isub        = zeros<ivec>(maxent);
  arma::ivec  jsub        = zeros<ivec>(maxent);
  arma::vec   U           = zeros(maxan);
  arma::vec   V           = zeros(maxan);
  arma::vec   AIDIAG      = zeros(maxan);
  arma::vec   COEF        = zeros(maxan);
  List        Out;
  //ivecs.replace(0, -1);
  //ivecd.replace(0, -1);
  iun22.open ("AInverse.txt", ios::trunc);
  /* Main loop */
  //deta                    = dz;
  //icnta                   = 0;
  k                       = -1;
  for(i = maxg1; i < ina; i++){
    isir                  = ivecs(i) - 1;
    idam                  = ivecd(i) - 1;
    den                   = 1.0;
    if (isir > maxgh) {
      den                 = den - .25 * U(isir);
    }
    if (idam > maxgh) {
      den                 = den - .25 * U(idam);
    }
    rden                  = 1.0 / den;
    r25                   = .25 * rden;
    r50                   = -.50 * rden;
    U(i)                  += den;
    V(i)                  = sqrt(den);
    AIDIAG(i)             += rden;
    /*
     * Both parents known
     */
    if (isir > maxgh && idam > maxgh) {
      k++;
      isub(k)             = i;
      jsub(k)             = isir;
      COEF(k)             = r50;
      k++;
      isub(k)             = i;
      jsub(k)             = idam;
      COEF(k)             = r50;
      AIDIAG(isir)        += r25;
      AIDIAG(idam)        += r25;
      if (isir < idam) {
        ii = idam;
        idam = isir;
        isir = ii;
      }
      k++;
      isub(k)             = isir;
      jsub(k)             = idam;
      COEF(k)             = r25;
      if (isir == idam) {
        COEF(k)           += r25;
      }
    }
    /*
     * Uknown dam
     */
    else if (isir > maxgh && idam == maxgh) {
      k++;
      isub(k)             = i;
      jsub(k)             = isir;
      COEF(k)             = r50;
      AIDIAG(isir)        += r25;
    }
    /*
     * Unknown sire
     */
    else if (isir == maxgh && idam > maxgh) {
      k++;
      isub(k)             = i;
      jsub(k)             = idam;
      COEF(k)             = r50;
      AIDIAG(idam)        += r25;
    }
    if (i == ina) {
      goto statement_3;
    }
    /*
     * Inner loop
     */
    for(j = i + 1; j < ina; j++){
      isir                = ivecs(j) - 1;
      idam                = ivecd(j) - 1;
      summ                = dz;
      if (isir >= i) {
        summ              = .5e0 * V(isir);
      }
      if (idam >= i) {
        summ              += .5e0 * V(idam);
      }
      V(j)                = summ;
      U(j)                += summ * summ;
    }
    /*
     * End main loop
     */
  }
  statement_3:
    goto statement_4400;
  statement_4400:
    for(i = 0; i < ina; i++){
      //iun22.open ("Ainv.txt", ios::trunc);
      iun22 << i + 1 << " " << i + 1 << " " << AIDIAG(i) << std::endl;
      //iun22.close();
    }
    for(i = 0; i <= k; i++){
      ii                = isub(i);
      jj                = jsub(i);
      cx                = COEF(i);
      //iun22.open ("Ainv.txt", ios::trunc);
      iun22 << ii + 1 << " " << jj + 1 << " " << cx << std::endl;
      //iun22.close();
    }
  iun22.close();
  Out["ii"] = isub.head(k + 1);
  Out["jj"] = jsub.head(k + 1);
  Out["cx"] = COEF.head(k + 1);
  Out["D"]  = AIDIAG.head(ina);
 // return(Out);
}