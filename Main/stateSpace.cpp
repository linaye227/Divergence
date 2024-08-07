/*******************************************************************************
 *                                                                             *
 * Cosmos:(C)oncept et (O)utils (S)tatistique pour les (Mo)deles               *
 * (S)tochastiques                                                             *
 *                                                                             *
 * Copyright (C) 2009-2012 LSV & LACL                                          *
 * Authors: Paolo Ballarini Benoît Barbot & Hilal Djafri                       *
 * Website: http://www.lsv.ens-cachan.fr/Software/cosmos                       *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify        *
 * it under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation; either version 3 of the License, or           *
 * (at your option) any later version.                                         *
 *                                                                             *
 * This program is distributed in the hope that it will be useful,             *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of              *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
 * GNU General Public License for more details.                                *
 *                                                                             *
 * You should have received a copy of the GNU General Public License along     *
 * with this program; if not, write to the Free Software Foundation, Inc.,     *
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                 *
 * file stateSpace.cpp created by Benoit Barbot on 01/12/11.                   *
 *******************************************************************************
 */

#include "stateSpace.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <deque>
#include <fstream>
#include <iostream>
#include <set>
// #include <bits/stdc++.h>
#include "sparse_io.hpp"

namespace boostmat = boost::numeric::ublas;
// using namespace boost::numeric::ublas;

stateSpace::stateSpace() {
  nbState = 0;
  nbTrans = 0;
  maxRate = 0;
}

int stateSpace::findHash(const vector<int> *vect) const {
  const auto it = S.find(vect);
  if (it != S.end())
    return (it->second); // found
  else
    return (-1);
}

double stateSpace::getMu(int state) const { return (*muvect)[state]; }

void stateSpace::add_state(vector<int> v) {
  vector<int> *v2 = new vector<int>(v);
  S[v2] = (int)nbState;
  // for(int i = 0  ; i< v2->size(); i++)cerr << (*v2)[i]<< ",";
  // cerr << endl;
  findstate->push_back(*v2);
  nbState++;
  if ((nbState % 100000) == 0)
    cerr << "Number of states :" << nbState << endl;
}
set <vector <int>> Aprime;
set <vector <int>> AlmostLoose1;
map <vector<int>,int> mp;
int cmax = 2;
int sno = 2;
int pmax = 5;
void stateSpace::exploreStateSpace() {
  // apply a Dijkstra algorithm on the product of the SPN an the LHA to produce
  // the state space. The list of state is store in the hash table S and
  // the transition list is stored in transitionList.
  N.reset();
  cerr << "Exploring state space" << endl;
  stack<vector<int>, vector<vector<int>>> toBeExplore;
  auto init = N.getState().getVector();
  cerr << "Today\n";
  cerr << init.size() << endl;
  for (auto i : init) cerr << i <<" ";
  cerr << endl;
  A.reset(N.getState());
  mp[init]++;
  pmax = init[init.size()-2] + 1; // Equal to the maximum possible clients in the system + 1 to reduce the size of the marking
  cmax = init[init.size()-3]; // Equal to the number of different types of clients available
  init.push_back(A.CurrentLocation);
  toBeExplore.push(init);
  findstate = new vector<vector<int>>(0);
  add_state(init);
  while (!toBeExplore.empty()) {
    vector<int> place = toBeExplore.top();
    toBeExplore.pop();
    vector<int> currentstate = place;
    int lhaloc = place.back();
    A.CurrentLocation = lhaloc;
    place.pop_back();
    N.setState(place);
    const auto savMark = N.getState();
    for (size_t t = 0; t < N.tr; t++){
      for (let b : N.Transition[t].bindingList) {
        N.setState(savMark);
        A.CurrentLocation = lhaloc;
        if (N.IsEnabled(t, b)) {
          N.fire(t, b, 0.0);
          int SE = A.synchroniseWith(t, N.getState(), b);
          if (SE > -1) {
            string lab = N.Transition[t].label;
            vector<int> marking = N.getState().getVector(); // new marking
            if (lab[0]=='I'){
              // If it is a input channel transition, no change has to be made in the marking
            }
            else if (lab[0]=='D'){
              // Client Departure Transition
              string num1 = "";
              for (int i=4;i<lab.size();i++){
                if (lab[i]!='_') num1+=lab[i];
                else break;
              }
              int numb1 = stoi(num1); // Channel number from which client has departed
              if (numb1==0){

              }
              else {
                for (int i=2+(numb1-1)*((cmax+2)*pmax);i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1;i++){
                  marking[i] = marking[i+1]; // altering the marking (making it such that everything in the channel is shifted back by 1 because the client in the first position has left)
                }
                marking[2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1] = 0;
                for (int i=2+(numb1-1)*((cmax+2)*pmax) + cmax*pmax;i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax + pmax - 1;i++){
                  if (marking[i]==0 && marking[i+1]==1 ) {
                    marking[i] = 1;
                    marking[i+1]  =0; // Chamging the markings corresponding to the "In Place" in the petri net
                    break;
                  }
                }
              }
            }
            else {
              string num1 = "";
              int numb1=0,numb2;
              int undercount = 0;
              for (int i=3;i<lab.size();i++){
                if (lab[i]!='_') {
                  num1+=lab[i];
                }
                else {
                  undercount++;
                  if (undercount>1) break;
                  numb1 = stoi(num1);
                  num1 = "";
                  i++;
                }
              }
              numb2 = stoi(num1); // Channel to which the client arrives
              if (numb1==0){
                // If the Client came from Input channel, do nothing (marking is already correctly built)
              }
              else if (numb1==numb2){
                for (int i=2+(numb1-1)*((cmax+2)*pmax);i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1;i++){
                  marking[i] = marking[i+1]; // make the marking shift back by 1
                }
                marking[2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1] = 0;
                for (int i=2+(numb1-1)*((cmax+2)*pmax) + cmax*pmax;i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax + pmax - 1;i++){
                  if (marking[i]==0 && marking[i+1]==1 ) {
                    marking[i] = 1;
                    marking[i+1]  =0; // The "In Place" is also shifted back by 1
                    break;
                  }
                }
              }
              else {
                for (int i=2+(numb1-1)*((cmax+2)*pmax);i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1;i++){
                  marking[i] = marking[i+1];
                }
                marking[2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1] = 0;
                for (int i=2+(numb1-1)*((cmax+2)*pmax) + cmax*pmax;i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax + pmax - 1;i++){
                  if (marking[i]==0 && marking[i+1]==1 ) {
                    marking[i] = 1;
                    marking[i+1]  =0; // This case is also similar to the previous cases
                    break;
                  }
                }
              }
            }
            if (mp[marking]==0){
              mp[marking]++;
              nbTrans++;
              if (marking[0]-1<=marking[marking.size()-1] && marking[0]!=1){     
                // These are the states having number of clients <= the minimum number of clients required for the given theta and epsylon values
                marking.push_back(A.CurrentLocation);
                auto its = S.find(&marking); // checking if the marking already exists
                if (its == S.end()) {
                  toBeExplore.push(marking);
                  add_state(marking); 
                }
              }
              else if (marking[0]==1){
                // Target Set
                marking.push_back(A.CurrentLocation);
                auto its = S.find(&marking);
                if (its == S.end()) {
                  Aprime.insert(marking);
                  add_state(marking);
                }
              }
              else if (marking[0]-1 >=marking[marking.size()-2]){
                // If the total number of clients is equal to the maximum number of clients corresponding to the given theta and epsylon, no more clients must be added (f1)
                marking.push_back(A.CurrentLocation);
                auto its = S.find(&marking);
                if (its == S.end()) {
                  AlmostLoose1.insert(marking);
                  add_state(marking);
                }
              }
              else{
                marking.push_back(A.CurrentLocation);
                auto its = S.find(&marking); 
                if (its == S.end()) {
                  toBeExplore.push(marking);
                  add_state(marking); 
                }
              }
            }
          }
        }
      }
    }
  }
  cerr << nbState << " states found" << endl
       << nbTrans << " transitions found" << endl;
  cerr << "Size of AlmostLoose1: " << AlmostLoose1.size() << endl;
  cerr << "Size of Aprime " << Aprime.size() << endl;
}

void stateSpace::computeDecisiveness(double threshold) {
  // apply a Dijkstra algorithm on the product of the SPN an the LHA to produce
  // the state space. The list of state is store in the hash table S and
  // the transition list is stored in transitionList.

  N.reset();
  cerr << "Exploring state space" << endl;

  stack<vector<int>, vector<vector<int>>> toBeExplore;
  auto init = N.getState().getVector();

  A.reset(N.getState());
  init.push_back(A.CurrentLocation);
  toBeExplore.push(init);
  findstate = new vector<vector<int>>(0);
  add_state(init);

  vector<size_t> immTrans;
  for (size_t t = 0; t < N.tr; t++)
    if (N.Transition[t].DistTypeIndex == IMMEDIATE ||
        N.Transition[t].DistTypeIndex == DETERMINISTIC)
      immTrans.push_back(t);

  while (!toBeExplore.empty()) {
    vector<int> place = toBeExplore.top();
    toBeExplore.pop();

    vector<int> currentstate = place;

    int lhaloc = place.back();
    A.CurrentLocation = lhaloc;
    place.pop_back();
    N.setState(place);

    bool existImm = false;

    for (size_t t = 0; t < immTrans.size() && !existImm; t++)
      for (let b : N.Transition[immTrans[t]].bindingList)
        if (N.IsEnabled(immTrans[t], b)) {
          existImm = true;
          break;
        }
    const auto savMark = N.getState();

    for (size_t t = 0; t < N.tr; t++)
      for (let b : N.Transition[t].bindingList) {

        N.setState(savMark);
        A.CurrentLocation = lhaloc;

        if (N.IsEnabled(t, b) &&
            (!existImm || N.Transition[t].DistTypeIndex == IMMEDIATE ||
             N.Transition[t].DistTypeIndex == DETERMINISTIC)) {

          N.fire(t, b, 0.0);
          // N.symmetrize();

          int SE = A.synchroniseWith(t, N.getState(), b);

          if (SE > -1) {
            // cerr << "transition:" << *it << endl;
            vector<int> marking = N.getState().getVector();

            nbTrans++;
            marking.push_back(A.CurrentLocation);
            // vector<double> Param = N.GetDistParameters(*it);
            // transitionsList.push( make_pair(make_pair(currentstate,
            // marking),Param[0] ));
          
            auto its = S.find(&marking);
            if (its == S.end()) {

              /*
              //N.Marking.printHeader(cerr);
               cerr << "state:"<< nbState << " -> ";
                  for( let i : marking)cerr << i << " ";
                  //N.Marking.print(cerr, 0.0);
               cerr << endl;
              */

              toBeExplore.push(marking);
              add_state(marking);
            } // else { cerr << " -> " << its->second  << endl;}
          }
          // N.unfire(t,b);
        }
      }
  }
  cerr << nbState << " states found" << endl
       << nbTrans << " transitions found" << endl;
}
set <int> rows;
void stateSpace::buildTransitionMatrix() {
  cerr << "Building transition matrix" << endl;

  // transform the transition list into a sparse transition probability matrix (the probability matrix that i need to find the reachability probability)
  transitionsMatrix = new boost::numeric::ublas::compressed_matrix<double>(
      nbState, nbState, nbTrans);
  auto &mat = *transitionsMatrix;

  cerr << "Exploring graph" << endl;
  for (size_t i = 0; i < nbState; i++) {
    // iterating in the generated state space
    vector<int> place = (*findstate)[i]; // ith state
    auto it1 = AlmostLoose1.find(place);
    auto it2 = Aprime.find(place);
    if (it1!=AlmostLoose1.end() || it2!=Aprime.end()) {
      mat(i,i) = 1;
      continue;
    }
    int lhaloc = place.back();
    place.pop_back();
    pmax = place[place.size()-2]  + 1;
    cmax = place[place.size()-3];
    N.setState(place);
    mat(i,i) = 0.0;
    for (size_t t = 0; t < N.tr; t++) {
      // Loop over binding here
      abstractBinding b;
      if (N.IsEnabled(t, b)) {
        A.CurrentLocation = lhaloc;
        N.fire(t, b, 0.0);
        vector<int> marking = N.getState().getVector();
        int SE = A.synchroniseWith(t, N.getState(), b);
        string lab = N.Transition[t].label;
        if (SE > -1) {
          N.unfire(t, b);
          if (lab[0]=='I'){
            
          }
          else if (lab[0]=='D'){
            string num1 = "";
            for (int i=4;i<lab.size();i++){
              if (lab[i]!='_') num1+=lab[i];
              else break;
            }
            int numb1 = stoi(num1);
            if (numb1==0){

            }
            else {
              for (int i=2+(numb1-1)*((cmax+2)*pmax);i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1;i++){
                marking[i] = marking[i+1];
              }
              marking[2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1] = 0;
              for (int i=2+(numb1-1)*((cmax+2)*pmax) + cmax*pmax;i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax + pmax - 1;i++){
                if (marking[i]==0 && marking[i+1]==1 ) {
                  marking[i] = 1;
                  marking[i+1]  =0;
                  break;
                }
              }
            }
          }
          else {
            string num1 = "";
            int numb1=0,numb2;
            int undercount = 0;
            for (int i=3;i<lab.size();i++){
              if (lab[i]!='_') {
                num1+=lab[i];
              }
              else {
                undercount++;
                if (undercount>1) break;
                numb1 = stoi(num1);
                num1 = "";
                i++;
              }
            }
            numb2 = stoi(num1);
            if (numb1==0){

            }
            else if (numb1==numb2){
              for (int i=2+(numb1-1)*((cmax+2)*pmax);i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1;i++){
                marking[i] = marking[i+1];
              }
              marking[2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1] = 0;
              for (int i=2+(numb1-1)*((cmax+2)*pmax) + cmax*pmax;i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax + pmax - 1;i++){
                if (marking[i]==0 && marking[i+1]==1 ) {
                  marking[i] = 1;
                  marking[i+1]  =0;
                  break;
                }
              }
            }
            else {
              for (int i=2+(numb1-1)*((cmax+2)*pmax);i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1;i++){
                marking[i] = marking[i+1];
              }
              marking[2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax-1] = 0;
              for (int i=2+(numb1-1)*((cmax+2)*pmax) + cmax*pmax;i<2+(numb1-1)*((cmax+2)*pmax)+(cmax)*pmax + pmax - 1;i++){
                if (marking[i]==0 && marking[i+1]==1 ) {
                  marking[i] = 1;
                  marking[i+1]  =0;
                  break;
                }
              }
            }
          }
          marking.push_back(A.CurrentLocation);
          // if (find(AlmostLoose1.begin(),AlmostLoose1.end(),marking)!=AlmostLoose1.end() || find(Aprime.begin(),Aprime.end(),marking)!=Aprime.end()) rows.insert(i);
          N.GetDistParameters(t, b);
          // cerr << N.ParamDistr[0] << ", ";
          int j = findHash(&marking);
          mat(i, j) = N.ParamDistr[0];
          /*cerr << "->state:";
           for (auto it2=marking.begin(); it2!= marking.end() ; it2++) {
           cerr << *it2 << ":";
           }
           cerr << N.ParamDistr[0] << endl;*/
        } else
          N.unfire(t, b);
      }
    }
  }
  // for (size_t i = 0; i < transitionsMatrix->size1(); i++) {
  //   long double sum1 = 0;
  //   if (rows.find(i) != rows.end()) {
  //     for (size_t j = 0; j < transitionsMatrix->size2(); j++) (*transitionsMatrix)(i, j)= 0;
  //     (*transitionsMatrix)(i, i) = 1;
  //   }
  //   else {
  //     for (size_t j = 0; j < transitionsMatrix->size2(); j++) sum1 += (*transitionsMatrix)(i, j);
  //     for (size_t j = 0; j < transitionsMatrix->size2(); j++) (*transitionsMatrix)(i, j)/=sum1;
  //   }
  // }
  /* Quick fix to redo */

  /*
   cerr << "uniformize to 1" << endl;
   for (it1_t it1 = mat.begin1(); it1 != mat.end1(); it1++)
   {
   double sum = 0.0;
   for (it2_t it2 = it1.begin(); it2 != it1.end(); it2++){
   //cerr << "non null:" << it2.index1() << ":" << it2.index2() << endl;
   if(it2.index1()!= it2.index2())sum += *it2;
   }
   for (it2_t it2 = it1.begin(); it2 != it1.end(); it2++){
   //cerr << "non null:" << it2.index1() << ":" << it2.index2() << endl;
   *it2 /= sum;
   }
   mat(it1.index1(),it1.index1())= 1.0;
   }
   */

  cerr << "Adding self loop" << endl;
  // Add self loop to ensure that mat is a probability matrix.
  // If the model is a CTMC the value on diagonal are wrong.

  for (auto it1 = mat.begin1(); it1 != mat.end1(); it1++) {
    double sum = 1.0;
    for (auto it2 = it1.begin(); it2 != it1.end(); it2++) {
      // cerr << "non null:" << it2.index1() << ":" << it2.index2() << endl;
      if (it2.index1() != it2.index2())
        sum -= *it2;
    }
    mat(it1.index1(), it1.index1()) = sum;
  }

  cerr << " copying" << endl;

  finalVector = new boost::numeric::ublas::vector<double>(nbState);
  for (hash_state::iterator it = S.begin(); it != S.end(); it++) {
    A.CurrentLocation = it->first->back();
    if (A.isFinal()) {
      (*finalVector)(it->second) = 1.0;
      // cerr << "final:" << it->second << endl;
    } else {
      (*finalVector)(it->second) = 0.0;
    }
  }
}

/*double stateSpace::maxRate(){
 double t = 0.0;
 for(boost::numeric::ublas::compressed_matrix<double>::iterator1 it
 = transitionsMatrix->begin1(); it!= transitionsMatrix->end1(); it++){
 for(boost::numeric::ublas::compressed_matrix<double>::iterator2 it2
 = it.begin(); it2!= it.end(); it2++){
 if(it.index1() != it2.index2())t = max(t,*it2);
 };
 };
 return(t);
 }*/

double stateSpace::uniformizeMatrix() {
  // First Compute infinitesimal generator
  // replace all value on the diagonal by opposite of the sum

  double lambda = 0.0;
  for (auto it1 = transitionsMatrix->begin1(); it1 != transitionsMatrix->end1();
       it1++) {
    double sum = 0.0;
    for (auto it2 = it1.begin(); it2 != it1.end(); it2++) {
      // cerr << "non null:" << it2.index1() << ":" << it2.index2() << endl;
      if (it2.index1() != it2.index2())
        sum += *it2;
    }
    lambda = max(lambda, sum);
    (*transitionsMatrix)(it1.index1(), it1.index1()) = -sum;
  }
  // Divide each coefficient of the matrix by lambda
  // and add 1 on the diagonal

  for (auto it1 = transitionsMatrix->begin1(); it1 != transitionsMatrix->end1();
       it1++) {
    for (auto it2 = it1.begin(); it2 != it1.end(); it2++) {
      // cerr << "non null:" << it2.index1() << ":" << it2.index2() << endl;
      *it2 /= lambda;
      if (it2.index1() == it2.index2())
        *it2 += 1.0;
    }
  }
  maxRate = lambda;
  return lambda;
}

void stateSpace::printP() {
  cerr << "Probability transition matrix:" << endl;
  for (size_t i = 0; i < transitionsMatrix->size1(); i++) {
    for (size_t j = 0; j < transitionsMatrix->size2(); j++)
      cerr << (*transitionsMatrix)(i, j) << "\t";
    cerr << endl;
  }
  cerr << endl << "Final Vector" << endl;
  for (size_t i = 0; i < finalVector->size(); i++) {
    cerr << (*finalVector)(i);
    cerr << endl;
  }
}

void stateSpace::outputMat() {
  cerr << "Exporting the transition matrix" << endl;

  fstream outputFile;
  outputFile.open("matrixFile", fstream::out);

  outputFile << boostmat::io::sparse(*transitionsMatrix);
  outputFile << *finalVector << endl;

  for (hash_state::iterator it = S.begin(); it != S.end(); it++) {
    outputFile << "(";
    vector<int> vect = *(*it).first;
    for (size_t i = 0; i < N.Msimpletab.size(); i++) {
      if (i > 0)
        outputFile << ",";
      outputFile << vect[N.Msimpletab[i]];
    };

    /*for(int i =0; i< vect.size()-1; i++){
     if(i>0)outputFile << ",";
     outputFile << vect[i];
     }*/
    outputFile << ")=";
    outputFile << (*it).second << endl;
  }

  outputFile.close();
}

void stateSpace::outputPrism() {
  cerr << "Exporting the model for Prism" << endl;

  fstream outputFile;
  outputFile.open("prismStates.sta", fstream::out);

  outputFile << "(";
  for (size_t i = 0; i < N.Place.size(); i++) {
    outputFile << N.Place[i].label.substr(1,N.Place[i].label.size()-1);
    outputFile << ",";
  };
  outputFile << "automata)" << endl;

  for (size_t it = 0; it < findstate->size(); it++) {
    outputFile << it << ":(";
    vector<int> vect = (*findstate)[it];
    for (size_t i = 0; i < N.Place.size(); i++) {
      outputFile << vect[i];
      outputFile << ",";
    };
    outputFile << vect[vect.size() - 1];

    outputFile << ")" << endl;
  }

  outputFile.close();

  fstream outputMatrixFile;
  outputMatrixFile.open("prismMatrix.tra", fstream::out);
  long long transnum =0;
  for (auto it1 = transitionsMatrix->begin1(); it1 != transitionsMatrix->end1();
       it1++) {
    for (auto it2 = it1.begin(); it2 != it1.end(); it2++) {
      if (*it2 >= (1e-16)) {
        transnum+=1;
      }
    }
  }
  outputMatrixFile << nbState << " " << transnum << endl;
  for (auto it1 = transitionsMatrix->begin1(); it1 != transitionsMatrix->end1();
       it1++) {
    for (auto it2 = it1.begin(); it2 != it1.end(); it2++) {
      if (*it2 >= (1e-16)) {
        outputMatrixFile << it2.index1() << " " << it2.index2() << " " << (*it2)
                         << endl;
      }
    }
  }
  outputMatrixFile.close();

  fstream outputProperty;
  outputProperty.open("prismProperty.pctl", fstream::out);
  outputProperty << "P=? [ (true) U  (Count=1"; // Change according to the target set
  outputProperty << ")]";
  outputProperty.close();

  fstream outputLabel;
  outputLabel.open("prismLabel.lbl", fstream::out);
  outputLabel << "0='init' 1='deadlock'\n 1: 0"; // change according to the initial state.
  outputLabel.close();
}

void stateSpace::launchPrism(string prismPath) {
  cerr << "Starting Prism" << endl;
  string cmd = prismPath +
               " -gs -maxiters 1000000000 -ctmc -importtrans prismMatrix.tra "
               "-importstates prismStates.sta -importlabels prismLabel.lbl -v "
               "-cuddmaxmem 10000000 prismProperty.pctl > prismOutput";
  if (0 != system(cmd.c_str())) {
    cerr << "Fail to launch prism" << endl;
    exit(EXIT_FAILURE);
  }
  cerr << "Prism finish" << endl;
}

void stateSpace::importPrism() {
  cerr << "Importing Prism result" << endl;
  string line;
  size_t poseq = 1;
  string pos;
  string prob;
  ifstream myfile("prismOutput");
  if (myfile.is_open()) {
    do {
      getline(myfile, line);
    } while (myfile.good() &&
             line != "Probabilities (non-zero only) for all states:");

    muvect = new vector<double>(nbState, 0.0);
    // int n=0;
    while (myfile.good() && poseq > 0) {
      getline(myfile, line);
      // cerr << line << endl;
      poseq = line.find("=");

      if (poseq != string::npos) {
        // cerr << line << endl;
        size_t si = 1 + line.find("(", 0);
        pos = line.substr(si, poseq - 1 - si);
        // cerr << "pos:" << pos << endl;
        prob = line.substr(poseq + 1, line.size());

        vector<int> vect;
        size_t it = 0;
        // cerr << "v:";
        while (it < pos.length()) {
          size_t it2 = pos.find(",", it);
          if (it2 == string::npos)
            it2 = pos.length();
          // cerr << "test:" << it<< ":" << it2 << endl;
          vect.push_back(atoi((pos.substr(it, it2 - it)).c_str()));
          // cerr << atoi((pos.substr(it,it2-it)).c_str() ) << ",";
          it = it2 + 1;
        }
        // cerr << endl;

        int state = findHash(&vect);
        // cerr << "state" << state << ":";
        (*muvect)[state] = atof(prob.c_str());
        // cerr << atof(prob.c_str());
        // muvect->push_back(atof(prob.c_str()));
        // S[new vector<int>(vect)] = n;
        // n++;
      }
    }
    myfile.close();
    // nbState = n;
  }
}

void stateSpace::outputTmpLumpingFun() {
  cerr << "Exporting the temporary lumping function" << endl;
  fstream outputlumptmp;

  outputlumptmp.open("lumpingfunTmp.cpp", fstream::out);
  outputlumptmp << "#include \"markingImpl.hpp\"" << endl << endl;

  for (size_t i = 0; i < N.Msimpletab.size(); i++) {
    int j = N.Msimpletab[i];
    outputlumptmp << "const int reducePL_"
                  << N.Place[j].label.substr(1, N.Place[j].label.length() - 1)
                  << " = " << i << ";" << endl;
  };

  outputlumptmp << "void SPN::print_state(const vector<int> &vect){" << endl;
  for (size_t i = 0; i < N.Msimpletab.size(); i++) {
    int j = N.Msimpletab[i];
    outputlumptmp << "\tcerr << \""
                  << N.Place[j].label.substr(1, N.Place[j].label.length() - 1)
                  << " = \" << vect[reducePL_"
                  << N.Place[j].label.substr(1, N.Place[j].label.length() - 1)
                  << "] << endl;" << endl;
  };
  outputlumptmp << "}" << endl;

  outputlumptmp
      << "bool SPN::precondition(const abstractMarking &Marking){return true;}"
      << endl;

  outputlumptmp << endl
                << "void SPN::lumpingFun(const abstractMarking "
                   "&Marking,vector<int> &vect){"
                << endl;
  for (size_t i = 0; i < N.Msimpletab.size(); i++) {
    int j = N.Msimpletab[i];
    outputlumptmp << "\tvect[reducePL_"
                  << N.Place[j].label.substr(1, N.Place[j].label.length() - 1)
                  << "] = Marking.P->_PL_"
                  << N.Place[j].label.substr(1, N.Place[j].label.length() - 1)
                  << "; //To Complete" << endl;
  };
  outputlumptmp << "}" << endl;
  outputlumptmp.close();
}

void stateSpace::outputVect() {
  cerr << "Exporting the probability vector" << endl;

  fstream outputFile;
  outputFile.open("muFile", fstream::out);
  outputFile.precision(15);

  outputFile << "[" << muvect->size() << "](";
  for (size_t i = 0; i < muvect->size(); i++) {
    if (i > 0)
      outputFile << ",";
    outputFile << (*muvect)[i];
  }
  outputFile << ")" << endl;

  for (hash_state::iterator it = S.begin(); it != S.end(); it++) {
    outputFile << "(";
    vector<int> vect = *(*it).first;
    for (size_t i = 0; i < N.Msimpletab.size(); i++) {
      if (i > 0)
        outputFile << ",";
      outputFile << vect[N.Msimpletab[i]];
    };

    /*for(int i =0; i< vect.size()-1; i++){
     if(i>0)outputFile << ",";
     outputFile << vect[i];
     }*/
    outputFile << ")=";
    outputFile << (*it).second << endl;
  }

  outputFile.close();
}

double stateSpace::returnPrismResult() { return (*muvect)[0]; }

void stateSpace::inputVect() {
  cerr << "Start reading muFile" << endl;
  ifstream inputFile("muFile", fstream::in);

  if (!inputFile.good()) {
    cerr << "Fail to open muFile at" << endl;
    exit(EXIT_FAILURE);
  }

  boostmat::vector<double> v1;
  inputFile >> v1;
  nbState = v1.size();
  muvect = new vector<double>(nbState);
  for (size_t i = 0; i < nbState; i++) {
    (*muvect)[i] = v1(i);
  }

  string line;
  size_t poseq;
  string pos;
  string stateid;
  while (inputFile.good()) {
    getline(inputFile, line);
    poseq = line.find("=");

    if (poseq != string::npos) {
      pos = line.substr(1, poseq - 2);
      stateid = line.substr(poseq + 1, line.size());

      vector<int> vect;
      size_t it = 0;
      while (it < pos.length()) {
        size_t it2 = pos.find(",", it);
        if (it2 == string::npos)
          it2 = pos.length();

        vect.push_back(atoi((pos.substr(it, it2 - it)).c_str()));
        it = it2 + 1;
      }

      S[new vector<int>(vect)] = (int)atoi(stateid.c_str());
    }
  }

  inputFile.close();
  if (S.empty()) {
    cerr << "muFile empty" << endl;
    exit(EXIT_FAILURE);
  }
  cerr << "Finished reading muFile with " << nbState << " states" << endl;
}

void stateSpace::inputMat() {
  fstream inputFile("matrixFile", fstream::in);

  if (!inputFile.good()) {
    cerr << "Fail to open matrixFile" << endl;
    exit(EXIT_FAILURE);
    return;
  }

  /*boostmat::matrix<double> m1;
   inputFile >> m1;
   nbState = m1.size1();*/
  boostmat::compressed_matrix<double, boostmat::row_major> m;
  inputFile >> boostmat::io::sparse(m);

  /*for (unsigned i = 0; i < nbState; ++ i)
   for (unsigned j = 0; j < nbState; ++ j)
   if(m1 (i,j) != 0.)  m (i, j) = m1 (i,j);*/
  transitionsMatrix = new boostmat::compressed_matrix<double>(m);

  // cerr << *transitionsMatrix << endl;

  finalVector = new boostmat::vector<double>();
  inputFile >> (*finalVector);

  string line;
  size_t poseq;
  string pos;
  string stateid;
  while (inputFile.good()) {
    getline(inputFile, line);
    // cerr << line << endl;
    poseq = line.find("=");

    if (poseq != string::npos) {
      pos = line.substr(1, poseq - 2);
      stateid = line.substr(poseq + 1, line.size());

      vector<int> vect;
      size_t it = 0;
      while (it < pos.length()) {
        size_t it2 = pos.find(",", it);
        if (it2 == string::npos)
          it2 = pos.length();
        vect.push_back(atoi((pos.substr(it, it2 - it)).c_str()));
        it = it2 + 1;
      }

      S[new vector<int>(vect)] = (int)atoi(stateid.c_str());
    }
  }

  inputFile.close();

  cerr << "DTMC size:" << finalVector->size() << endl;
}
