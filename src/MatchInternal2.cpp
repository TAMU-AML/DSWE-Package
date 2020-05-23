/*
// MIT License
//   
//   Copyright (c) 2020 Nitesh Kumar, Abhinav Prakash, and Yu Ding
//   
//   Permission is hereby granted, free of charge, to any person obtaining a copy
//   of this software and associated documentation files (the "Software"), to deal
//   in the Software without restriction, including without limitation the rights
//   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//   copies of the Software, and to permit persons to whom the Software is
//   furnished to do so, subject to the following conditions:
//     
//     The above copyright notice and this permission notice shall be included in all
//     copies or substantial portions of the Software.
//   
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//     SOFTWARE.
*/

// function to match data set
// match.cov : main matching function given data sets
// ref : data set considered as a reference
// obj : data set whose each observation searches for a matching observation in ref
// thres : threshold against which similarity is obtained between ref and obj data sets
// circ_pos : position of circular variables (such as wind direction etc), if supplied by user

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec matchcov(arma::mat& ref , arma::mat& obj, arma::rowvec& thres, arma::rowvec& circ_pos, int flag)
{
  // Variables to store no of rows
  int row_ref = ref.n_rows;
  int row_obj = obj.n_rows;

  // Vector to store matched obj data set index
  arma::vec match(row_ref);
  match.fill(0.0);

  // Vector to store index of obj data set
  arma::vec index = arma::linspace(1, row_obj, row_obj);

  // Looping through each element of reference set
  for(int i = 0; i < row_ref; i++){

    // Calculating ratio between obj set and each ref observation
    arma::rowvec ref_i = ref.row(i);
    arma::mat score = arma::abs(obj.each_row() - ref_i);

    // If circular variable supplied, extracting corresponding data sets
    if(flag > 0){

       // Circular variable as vector in ref and matrix in obj
       arma::vec circref = ref_i(arma::conv_to<arma::uvec>::from(circ_pos - 1));
       arma::mat circdata = obj.cols(arma::conv_to<arma::uvec>::from(circ_pos - 1));
       int cols = circdata.n_cols;

       // Loop to update circular variable score
       for(int j = 0; j < cols; j++){

         // Updating score of circular variable
         arma::uvec id_dcor = arma::find((circdata.col(j) - circref(j)) > 180);
         arma::vec score_vec = score.col(j);
         arma::vec circdata_vec = circdata.col(j);
         score_vec(id_dcor) = (360 - (circdata_vec(id_dcor) - circref(j))) / circref(j);
         score.col(j) = score_vec;

          }
       }

    // Checking calculated score against threshold
    arma::umat decision(obj.n_rows, obj.n_cols);
    int obj_col = obj.n_cols;

    for(int k = 0; k < obj_col; k++) {

      arma::uvec des = score.col(k) < thres(k);
      decision.col(k) = des;

    }

    // Finding index of obj data set whose score <= threshold
    arma::mat id = arma::conv_to<arma::mat>::from(decision);
    arma::vec id_sum = sum(id, 1);
    arma::uvec id_index = find(id_sum == obj.n_cols);

    // Vector to score adjusted score
    arma::mat score_adjusted;

    // vector to store maxscore out of matched index
    arma::vec maxscore;

    // Variables and vectors to store min score out of max score and index
    double min_score;
    arma::uvec min_id;
    arma::uvec id_min;

    // Vector to store unmatched observations index
    arma::uvec unmatched_id;

    // filtering incase of multiple match
    int id_num = id_index.n_elem;
    if(id_num > 0){

      // Calculating adjusted score and filtering most similar
      score_adjusted = score.each_row() / thres;
      maxscore = arma::max(score_adjusted.rows(id_index), 1);
      min_score = arma::min(maxscore);
      min_id = find(maxscore == min_score);
      id_min = id_index(min_id(0));

      // Storing the index of matched observation in obj
      match(i) = index(arma::conv_to<int>::from(id_min));

      // Index of unmatched observation in obj
      unmatched_id = find(index != match(i));

      // Removing matched observation from obj
      obj = obj.rows(unmatched_id);

      // Removing matched index from index vector created before
      index = index.elem(unmatched_id);

     }

  }

  return match;

}
