#include <Rcpp.h>
#include <math.h>
#include <string>
using namespace Rcpp;

double estimate_beta(NumericVector group, NumericVector atrisk0, NumericVector atrisk1, NumericVector status) {
  
  double beta, beta_old = 1000;
  double beta_low = -2;
  double beta_up = 4;
  double eps = 0.0001;
  double score, score_deriv, score_low, score_up;
  double w, w_deriv;
  double S0, S1;
  int n = group.size();
  int k=0;
  
  // Find starting value using bisection method
  // score at lower limit
  score_low = 0;
  score_up = 0;
  
  for(int i = 0; i < n; i++) {
    if (status[i] == 1) {
      // score at lower limit
      S0 = (atrisk0[i] + exp(beta_low) * atrisk1[i]);
      S1 = exp(beta_low) * atrisk1[i];
      w = group[i] - S1 / S0;
      score_low += w;
      // score at upper limit
      S0 = (atrisk0[i] + exp(beta_up) * atrisk1[i]);
      S1 = exp(beta_up) * atrisk1[i];
      w = group[i] - S1 / S0;
      score_up += w;
    }
  }
  // Stopping criteria for initial boudaries
  // Uses that the score function decreases when beta increases
  if (score_low < 0) {
    beta = beta_low;
    k = 5;
  }
  if (score_up > 0) {
    beta = beta_up;
    k = 5;
  }

  while ((k < 4)) {
    k++;
    // score at middle limit
    score = 0;
    beta = (beta_low + beta_up) / 2;
    for(int i = 0; i < n; i++) {
      if (status[i] == 1) {
        S0 = (atrisk0[i] + exp(beta) * atrisk1[i]);
        S1 = exp(beta) * atrisk1[i];
        w = group[i] - S1 / S0;
        score += w;
      }
    }
    //printf("score:%f %f %f, beta: %f %f %f\n", score_low, score, score_up, beta_low, beta, beta_up);
    
    if (score > 0) {
      beta_low = beta;
      score_low = score;
    } else {
      beta_up = beta;
      score_up = score;
    }
  }
  
  
  // Newton iteration
  while ((fabs(beta - beta_old) > eps) && (k < 15)) {
    k++;
    score = 0;
    score_deriv = 0;
    
    for(int i = 0; i < n; i++) {
      if (status[i] == 1) {
        S0 = (atrisk0[i] + exp(beta) * atrisk1[i]);
        S1 = exp(beta) * atrisk1[i];
        w = group[i] - S1 / S0;
        w_deriv = - (S1 * S0 - S1 * S1) / (S0 * S0);
        score += w;
        score_deriv += w_deriv;
      }
    }
    
    beta_old = beta;
    beta = beta - score / score_deriv;
//    printf("%f\n", beta);
  }
  
  return(beta);
}


// [[Rcpp::export]]
List analyze_lwyy(NumericMatrix x_in, int n_id, int n_total) {
  
  NumericMatrix id_mat(n_id, 5); // id group start stop Bi
  std::fill(id_mat.begin(), id_mat.end(), 0);
  
  int i_id;
  int atrisk0 = 0;
  int atrisk1 = 0;
  double eventtime;
  double A = 0; // Used to calculate robust variance estimator
  double B = 0; // Used to calculate robust variance estimator
  double beta_est;
  // Initialize the matrix x   
  NumericMatrix x(x_in.nrow(), 10);
  x(_, 0) = x_in(_, 0);
  x(_, 1) = x_in(_, 7);
  x(_, 2) = x_in(_, 2);
  x(_, 3) = x_in(_, 3);
  x(_, 4) = x_in(_, 6);
  
  int col_id = 0;
  int col_group = 1;
  int col_start = 2; 
  int col_stop = 3;
  int col_status = 4;
  int col_atrisk0 = 5;
  int col_atrisk1 = 6;
  int col_S0 = 7;
  int col_S1 = 8;
  int col_W = 9;
  // END Initialize the matrix x   
  
  // Initialize id, group, start time, and stop time in matrix id_mat
  for (int i = 0; i < x.nrow(); i++) { // loop over the event times
    i_id = x(i, 0) - 1;
    
    if (id_mat(i_id, 0) == 0) {
      id_mat(i_id, 0) = x(i, col_id);
      id_mat(i_id, 1) = x(i, col_group);
    }
    
    // Grab start time of individual (column 3, i.e. col index 2) 
    if (id_mat(i_id, 2) == -1)
      id_mat(i_id, 2) = x(i, col_start);
    else {
      if (x(i, col_start) < id_mat(i_id, 2)) 
        id_mat(i_id, 2) = x(i, col_start);
    }
    
    // Grab stop time of individual (column 4 of id_mat, i.e. col index 3) 
    if (id_mat(i_id, 3) == -1)
      id_mat(i_id, 3) = x(i, col_stop);
    else {
      if (x(i, col_stop) > id_mat(i_id, 3)) 
        id_mat(i_id, 3) = x(i, col_stop);
    }
  }  // END Initialize id, group, start time, and stop time in matrix id_mat
  
  // Calculate the number of at-risk subjects for every eventtime
  for (int i = 0; i < x.nrow(); i++) {
    eventtime = x(i, col_stop);
    atrisk0 = atrisk1 = 0;
    // At-risk number calculation
    for (int j = 0; j < id_mat.nrow(); j++) {
      if ((id_mat(j, 2) < eventtime) && ((id_mat(j, 3) >= eventtime))) {
        if (id_mat(j, 1) == 0)
          atrisk0++;
        if (id_mat(j, 1) == 1)
          atrisk1++;
      }
    }
    // Calculations: S0, S1=S2, W
    x(i, col_atrisk0) = atrisk0;
    x(i, col_atrisk1) = atrisk1;
  } // END at-risk calculation
  
  
   beta_est = estimate_beta(x(_, col_group), x(_, col_atrisk0), x(_, col_atrisk1), x(_, col_status));
  
  // Calculate S0, S1=S2, W, A, and Bi
  for (int i = 0; i < x.nrow(); i++) {
    eventtime = x(i, col_stop);
    // S0, S1=S2, W
    atrisk0 = x(i, col_atrisk0);
    atrisk1 = x(i, col_atrisk1);
    x(i, col_S0) = atrisk0 + exp(beta_est) * atrisk1;
    x(i, col_S1) = exp(beta_est) * atrisk1;
    x(i, col_W) = x(i, col_group) - x(i, col_S1) / x(i, col_S0);
    
    A += x(i, col_status) * x(i, col_S1) / x(i, col_S0) * (1 - x(i, col_S1) / x(i, col_S0));
    
    // Calculation of Bi
    if (x(i, col_status) == 1) {
      
      id_mat(x(i, col_id)-1, 4) += x(i, col_W);
      
      for (int j = 0; j < id_mat.nrow(); j++) {
        if ((id_mat(j, 2) < eventtime) && ((id_mat(j, 3) >= eventtime))) {
          id_mat(j, 4) -= (id_mat(j, 1) - x(i, col_S1) / x(i, col_S0)) * exp(beta_est * id_mat(j, 1)) / x(i, col_S0);
        }
      }
    } // END Calculation of Bi
  }  // END for-loop
  
  
  // Calcualte sum_i B_i^2
  for (int j = 0; j < id_mat.nrow(); j++) {
    B += id_mat(j, 4) * id_mat(j, 4); 
  }
  // Finalize calculation of A and B
  B = B / n_total;
  A = A / n_total;
  
  // Prepare and return output variable
  List out;
  colnames(id_mat) = CharacterVector::create("id", "group", "study_start", "study_stop", "Bi");
  //colnames(x) = CharacterVector::create("id", "group", "study_start", "study_stop", "Bi");
  // id group study_start study_stop status atrisk0 atrisk1 S0 S1 W
  // out["x"] = x;
  out["id_mat"] = id_mat;
  out["A"] = A;
  out["B"] = B;
  out["robustSE"] = sqrt( B / (A*A) / n_total);
  out["information"] = n_total * A * A / B ;
  out["beta_est"] = beta_est;
  return(out);
}