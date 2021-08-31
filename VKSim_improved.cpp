#include <Rcpp.h>
using namespace Rcpp;
#include <tuple>


double calcDBE(NumericVector x){
  return x[0] - (x[1]/2) + 1;
}



std::tuple<NumericVector, bool> reactMolecule(NumericVector x, NumericVector y, double min_dbe) {
  std::size_t n = x.size();
  NumericVector out(n); 
  int temp = -1;
  for (std::size_t i = 0; i < n; ++i) {
    temp = x[i] + y[i];
    if(temp < 0) break;
    out[i] = temp;
  }
  
  double H_C = out[1]/out[0];
  double dbe = calcDBE(out);
  
  bool unchanged = ((temp < 0) | (H_C < 0.2) | (H_C > 3.1) | (dbe < min_dbe)); 
  if(unchanged){out = x;}
  
  return std::forward_as_tuple(out, unchanged);
}

// [[Rcpp::export]]
List generateReaction(List start, NumericVector reaction, int max_reaction, double min_dbe){
  List curr(start);
  List outputs;
  std::size_t sn = start.size();
  //Get element labels
  NumericVector first = start[0];
  std::vector<std::string> names = first.attr("names");
  Rcpp::Rcout << "Reaction Vector: "<< reaction << " Max " << max_reaction << " iterations" <<"\n";
  for(int i = 0; i < max_reaction; ++i){
    int counter = 0;
    List out(sn); 
    for (std::size_t j = 0; j < sn; ++j) {
      auto raw = reactMolecule(curr[j], reaction, min_dbe);
      NumericVector first = std::get<0>(raw);
      int unchanged = std::get<1>(raw);
      counter = counter + unchanged;
      first.attr("names") = names;
      out[j] = first;
    }
    Rcpp::Rcout << "Unreacted/Molecules: "<< counter<< "/" << sn << "\n";
    outputs.push_back(out);
    curr = out;
    if(counter == static_cast<int>(sn)){
      Rcpp::Rcout << "Breaking on iteration: "<< i + 1<< "\n";
      break;
    }
  }
  return outputs;
}


// [[Rcpp::export]]
List generateReactions(List start, List reactions, double min_dbe) {
  List curr(start);
  std::size_t rn = reactions.size();
  List outputs(rn);
  
  for (std::size_t i = 0; i < rn; ++i) {
    List current_rxn = reactions[i];
    NumericVector reaction = current_rxn[0];
    int max_reaction = current_rxn[1];
    auto temp = generateReaction(curr, reaction, max_reaction, min_dbe);
    outputs[i] = temp;
    curr = temp[temp.size() - 1];
  }
  
  return outputs;
}

// [[Rcpp::export]]
String createReaction(CharacterVector name_vec, std::vector<int> num_vec){
  String pos_formula("+");
  String neg_formula("-");
  std::size_t n = name_vec.size();
  for(std::size_t i = 0; i < n; ++i){
    if(num_vec[i] != 0) {
      String element = name_vec[i];
      element += std::to_string(abs(num_vec[i]));
      if(num_vec[i] < 0) {
        neg_formula.push_back(element);
      } else {
        pos_formula.push_back(element);
      }
    }
    else
    {continue;}
  }
  String complete("");
  if(pos_formula != "+"){complete.push_back(pos_formula);}
  if(neg_formula != "-"){complete.push_back(neg_formula);}
  
  return complete;
}
