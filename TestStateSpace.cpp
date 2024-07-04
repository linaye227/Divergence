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
 #include <array>
 #include <boost/numeric/ublas/io.hpp>
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/matrix_expression.hpp>
 #include <boost/numeric/ublas/matrix_sparse.hpp>
 #include <deque>
 #include <fstream>
 #include <iomanip>
 #include <iostream>
 #include <queue>
 #include <set>
 #include <stack>
 
 #include "parameters.hpp"
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
 int cmax = 2;
 int pmax = 48;
 int sno = 2;
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
 
 void stateSpace::exploreStateSpace() {
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
           N.symmetrize();
 
           int SE = A.synchroniseWith(t, N.getState(), b);
 
           if (SE > -1) {
             // cerr << "transition:" << *it << endl;
             vector<int> marking = N.getState().getVector();
 
             nbTrans++;
             marking.push_back(A.CurrentLocation);
             // vector<double> Param = N.GetDistParameters(*it);
             // transitionsList.push( make_pair(make_pair(currentstate,
             // marking),Param[0] ));
             string lab = N.Transition[t].abel
            if (lab[0]=='I'){
              
            }
            else if (lab[0]=='D'){
              string num1 = "";
              for (int i=4;i<lab.size();i++){
                if (lab[i]!='_') num1+=lab[i];
                else break;
              }
              int numb1 = stoi(num1);
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
              // cerr << "HERE\n";
              numb2 = stoi(num1);
              if (numb1==0){
                // for (auto i : marking) cerr << i << " ";
                // cerr << endl;
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
 
 class Precise_float {
   static const int nb_bucket = 1024;
 
 public:
   Precise_float() {
     value.fill(0.0);
     mine = nb_bucket - 1;
     maxe = 0;
   };
   Precise_float(double d) {
     Precise_float();
     *this += d;
   };
 
   Precise_float &operator=(double d) {
     value.fill(0.0);
     mine = nb_bucket - 1;
     maxe = 0;
     *this += d;
     return *this;
   }
 
   Precise_float &operator+=(double d) {
     if (d == 0.0)
       return *this;
     int exp = 0;
     double di = frexp(d, &exp);
     double di2 = this->value[exp + nb_bucket / 2] + di;
     if (abs(di2) >= 0.5 && abs(di2) < 1.0) {
       this->value[exp + nb_bucket / 2] = di2;
       this->mine = std::max(0, std::min(mine, exp + (nb_bucket / 2)));
       this->maxe =
           std::min(nb_bucket - 1, std::max(maxe, exp + (nb_bucket / 2)));
     } else {
       this->value[exp + nb_bucket / 2] = 0;
       double d2 = ldexp(di2, exp);
       *this += d2;
     }
 
     return *this;
   };
 
   Precise_float &operator+=(const Precise_float &rhs) {
     for (int i = rhs.mine; i <= rhs.maxe; i++) {
       if (rhs.value[i] != 0) {
         *this += ldexp(rhs.value[i], i - nb_bucket / 2);
       }
     }
     return *this;
   }
 
   friend Precise_float operator+(const Precise_float &lhs,
                                  const Precise_float &rhs) {
     Precise_float pf = lhs;
     pf += rhs;
     return pf;
   }
 
   double to_double() {
     double acc = 0.0;
     for (int i = mine; i <= maxe; i++) {
       acc += ldexp(value[i], i - nb_bucket / 2);
     }
     return acc;
   }
 
   friend std::ostream &operator<<(std::ostream &os, const Precise_float &pf) {
     os << "[";
     for (int i = pf.mine; i <= pf.maxe; i++) {
       if (pf.value[i] != 0.0)
         os << ",2^" << (i - nb_bucket / 2) << "=" << pf.value[i];
     }
     os << "]";
     return os;
   }
 
 private:
   std::array<double, nb_bucket> value;
   int mine, maxe;
 };
 
 struct sys_state {
   const vector<int> t;
   sys_state(const vector<int> &v) : t(v){};
   bool operator==(const sys_state &other) const {
     return t == other.t;
     /*assert(t.size() == other.t.size());
     for (unsigned int i = 0; i < t.size(); i++)
       if (t[i] != other.t[i])
         return false;
     return true;*/
   }
 };
 std::ostream &operator<<(std::ostream &os, const sys_state &st) {
   os << "[";
   for (int i = 0; i < st.t.size(); i++) {
     if (i > 0)
       os << "; ";
     os << st.t[i];
   }
   os << "]";
   return os;
 }
 
 namespace std {
 template <> struct hash<sys_state> {
   std::size_t operator()(const sys_state &k) const {
     using std::hash;
     using std::size_t;
 
     // Compute individual hash values for first,
     // second and third and combine them using XOR
     // and bit shifting:
     size_t res = hash<int>()(k.t[0]);
     for (unsigned int i = 1; i < k.t.size(); i++)
       res = (res ^ (hash<int>()(k.t[i]) << 1)) >> 1;
 
     return res;
   }
 };
 } // namespace std
 
 template <typename TF> class traj_weight {
 public:
   const sys_state states;
   TF weight;
   TF likelihood;
   int position;
   traj_weight(sys_state states, TF weight, TF likelihood, int pos)
       : states(states), weight(weight), likelihood(likelihood), position(pos) {}
 
   bool is_larger(traj_weight &rhs) { return this->weight > rhs.weight; }
 
   /*friend void swap(traj_weight *lhs, traj_weight *rhs) {
     std::swap(lhs, rhs);
     int pi = lhs->position;
     lhs->position = rhs->position;
     rhs->position = pi;
   }*/
 };
 template <typename TF>
 std::ostream &operator<<(std::ostream &os, const traj_weight<TF> &tw) {
   os << "states: " << tw.states << ", weight: " << tw.weight
      << ", likelyhood: " << tw.likelihood << ", position: " << tw.position;
   return os;
 }
 
 template <typename TF> class trajs_heap {
   vector<traj_weight<TF> *> heap;
   unordered_map<sys_state, traj_weight<TF> *> hash;
 
 public:
   trajs_heap() {}
   tuple<sys_state, TF, TF> pop() {
     traj_weight<TF> &traj = *heap[0];
     const TF pf = traj.weight;
     const TF likelihood = traj.likelihood;
     traj.weight = 0.0;
     traj.likelihood = 0.0;
     swap(0, heap.size() - 1);
     traj.position = -1;
     heap.pop_back();
     siftDown(0);
     return make_tuple(traj.states, pf, likelihood);
   }
   void add(const sys_state &states, const TF &weight, const TF &likelihood) {
     // cerr << "add to heap" << weight;
     auto tw = hash.find(states);
     if (tw != hash.end()) {
       traj_weight<TF> *traj = (tw->second);
       traj->likelihood =
           (traj->weight * traj->likelihood + weight * likelihood) /
           (traj->weight + weight);
       traj->weight += weight;
 
       if (traj->position < 0) {
         traj->position = heap.size();
         heap.push_back(traj);
       }
       siftUp(traj->position);
     } else {
       traj_weight<TF> *tw =
           new traj_weight<TF>(states, weight, likelihood, heap.size());
       hash.emplace(states, tw);
       heap.push_back(tw);
       siftUp(heap.size() - 1);
     }
   }
   bool empty() { return heap.size() == 0; }
   long size() { return heap.size(); }
 
 private:
   size_t getLeftChildIndex(size_t nodeIndex) const { return 2 * nodeIndex + 1; }
   size_t getRightChildIndex(size_t nodeIndex) const {
     return 2 * nodeIndex + 2;
   }
   size_t getParentIndex(size_t nodeIndex) const { return (nodeIndex - 1) / 2; }
 
   void swap(size_t a, size_t b) {
     auto t = heap[a];
     heap[a] = heap[b];
     heap[b] = t;
     int pi = heap[a]->position;
     heap[a]->position = heap[b]->position;
     heap[b]->position = pi;
   }
 
   void siftUp(size_t i) {
     size_t parentIndex;
 
     if (i != 0) {
       parentIndex = getParentIndex(i);
 
       if (heap[i]->is_larger(*heap[parentIndex])) {
         swap(i, parentIndex);
         siftUp(parentIndex);
       }
     }
   }
 
   void siftDown(size_t i) {
     size_t leftChildIndex, rightChildIndex, minIndex;
     leftChildIndex = getLeftChildIndex(i);
     rightChildIndex = getRightChildIndex(i);
     if (rightChildIndex >= heap.size()) {
       if (leftChildIndex >= heap.size())
         return;
       else
         minIndex = leftChildIndex;
     } else {
 
       if (heap[leftChildIndex]->is_larger(*heap[rightChildIndex]))
         minIndex = leftChildIndex;
       else
         minIndex = rightChildIndex;
     }
 
     if (heap[minIndex]->is_larger(*heap[i])) {
       swap(minIndex, i);
       siftDown(minIndex);
     }
   }
   friend ostream &operator<<(ostream &os, const trajs_heap<TF> &th) {
     os << "Heap contents:" << endl;
     for (auto tw : th.heap) {
       os << *tw << endl;
     }
     // os << "Hash contents:" << endl;
     // for (const auto &h : th.hash) {
     // os << h.first << ": " << *(h.second) << endl;
     //}
     return os;
   }
 };
 
 bool comp_traj(const pair<vector<int>, double> &lhs,
                const pair<vector<int>, double> &rhs) {
   return get<1>(lhs) < get<1>(rhs);
 }
 
 void stateSpace::computeDecisiveness(double threshold) {
   REHandling lumpingFun;
   // apply a Dijkstra algorithm on the product of the SPN an the LHA to produce
   // the state space.
   N.reset();
   cerr << "Exploring state space for decisiveness" << endl;
 
   trajs_heap<double> toBeExplore;
   auto init = N.getState().getVector();
 
   A.reset(N.getState());
   init.push_back(A.CurrentLocation);
   toBeExplore.add(init, 1.0, 1.0);
   // cerr << toBeExplore << endl << endl;
   findstate = new vector<vector<int>>(0);
 
   vector<size_t> immTrans;
   for (size_t t = 0; t < N.tr; t++)
     if (N.Transition[t].DistTypeIndex == IMMEDIATE ||
         N.Transition[t].DistTypeIndex == DETERMINISTIC)
       immTrans.push_back(t);
 
   Precise_float lower_bound_dec;
   Precise_float upper_bound_dec;
   Precise_float mid_value;
 
   while (!toBeExplore.empty() &&
          (1.0 - (lower_bound_dec.to_double() + upper_bound_dec.to_double()) >
           threshold)) {
     if (nbTrans % 100 == 0) {
       cerr << "\033[A\033[2K" << nbTrans
            << " transition taken, toExplore:" << toBeExplore.size()
            << " deciveness_int :[" << upper_bound_dec.to_double() << ","
            << 1.0 - lower_bound_dec.to_double()
            << "] width: " << std::setprecision(16)
            << (1.0 -
                (lower_bound_dec.to_double() + upper_bound_dec.to_double()))
            << endl;
     }
 
     // pop_heap(toBeExplore.begin(), toBeExplore.end(), comp_traj);
     auto state_prob = toBeExplore.pop();
     // cerr << toBeExplore << endl << endl;
     //   toBeExplore.pop_back();
 
     vector<int> place = get<0>(state_prob).t;
 
     vector<int> currentstate = place;
 
     int lhaloc = place.back();
     A.CurrentLocation = lhaloc;
     place.pop_back();
     N.setState(place);
 
     double total_weight = 0.0;
     double total_original_weight = 0.0;
     bool is_IS_enabled = lumpingFun.precondition(N.Marking);
     double mu_c =
         (is_IS_enabled ? lumpingFun.direct_computation(N.Marking) : 1.0);
     for (size_t t = 0; t < N.tr; t++)
       for (let b : N.Transition[t].bindingList)
         if (N.IsEnabled(t, b)) {
           N.GetDistParameters(t, b);
           total_original_weight += N.ParamDistr[0];
           N.fire(t, b, 0.0);
           double mu =
               (is_IS_enabled ? lumpingFun.direct_computation(N.Marking) : mu_c);
           N.unfire(t, b);
           total_weight += N.ParamDistr[0] * mu / mu_c;
         }
 
     bool existImm = false;
     for (size_t t = 0; t < immTrans.size() && !existImm; t++)
       for (let b : N.Transition[immTrans[t]].bindingList)
         if (N.IsEnabled(immTrans[t], b)) {
           existImm = true;
           break;
         }
     const auto savMark = N.getState();
     bool deadlock = true;
     if (P.verbose > 5) {
       cerr << "state " << get<1>(state_prob) << endl;
     }
     for (size_t t = 0; t < N.tr; t++)
       for (let b : N.Transition[t].bindingList) {
         N.setState(savMark);
         A.CurrentLocation = lhaloc;
 
         if (N.IsEnabled(t, b) &&
             (!existImm || N.Transition[t].DistTypeIndex == IMMEDIATE ||
              N.Transition[t].DistTypeIndex == DETERMINISTIC)) {
           deadlock = false;
 
           N.GetDistParameters(t, b);
           double p_orig = N.ParamDistr[0] / total_original_weight;
           N.fire(t, b, 0.0);
           double mu =
               (is_IS_enabled ? lumpingFun.direct_computation(N.Marking) : mu_c);
           double p_bias = (mu / mu_c) * p_orig;
           double likelyhood = get<2>(state_prob) * mu_c / mu;
           double prob = get<1>(state_prob) * p_bias;
 
           int SE = A.synchroniseWith(t, N.getState(), b);
 
           if (P.verbose > 5) {
             cerr << "from mu_c " << mu_c << " " << t << ":"
                  << N.Transition[t].label;
             b.print();
             // N.Marking.print(cerr, 0.0);
             cerr << "-> " << mu << " with original prob: " << p_orig
                  << " bias_prob: " << p_bias << " likelyhood: " << likelyhood
                  << " deciveness_int :[" << upper_bound_dec.to_double() << ","
                  << 1.0 - lower_bound_dec.to_double() << "] -> " << endl;
           }
 
           if (SE > -1 && likelyhood > 0.0) {
             // cerr << "transition:" << N.Transition[t].label << endl;
 
             if (A.isFinal()) {
               upper_bound_dec += prob;
               mid_value += prob * likelyhood;
               // cerr << "accepting" << endl;
             } else {
 
               vector<int> marking = N.getState().getVector();
 
               nbTrans++;
               marking.push_back(A.CurrentLocation);
 
               toBeExplore.add(marking, prob, likelyhood);
             }
           } else {
             // cerr << "transition:" << N.Transition[t].label << endl;
             lower_bound_dec += prob;
           }
         }
       }
     // Go to sink state
 
     if (P.verbose > 5)
       cerr << "-> to sink " << mu_c << " prob: "
            << (total_original_weight - total_weight) / total_original_weight
            << endl;
     if (deadlock) {
       lower_bound_dec += get<1>(state_prob); //* get<2>(state_prob);
     } else {
       if (total_original_weight >= total_weight) {
         lower_bound_dec += get<1>(state_prob) *
                            (total_original_weight - total_weight) /
                            total_original_weight;
       } else {
         N.Marking.print(cerr, 0.0);
         cerr << "Property not guaranteed: " << endl;
         cerr << "Property not guaranteed: " << endl;
         exit(1);
       }
     }
   }
   cerr << "Probability :" << nbTrans << " transition taken, decisiveness_int :["
        << mid_value.to_double() - threshold / 2.0 << ","
        << mid_value.to_double() + threshold * 2.0 << "]" << endl;
   //  << toBeExplore << endl;
 }
 
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
    N.setState(place);

    /*cerr << "state:";
     for (auto it2=currentstate.begin(); it2!= currentstate.end() ; it2++) {
     cerr << *it2 << ":";
     }
     cerr << endl;*/
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
            // cerr << "HERE\n";
            numb2 = stoi(num1);
            if (numb1==0){
              // for (auto i : marking) cerr << i << " ";
              // cerr << endl;
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
     outputFile << N.Place[i].label;
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
   outputMatrixFile << nbState << " " << nbTrans << endl;
 
   for (auto it1 = transitionsMatrix->begin1(); it1 != transitionsMatrix->end1();
        it1++) {
     for (auto it2 = it1.begin(); it2 != it1.end(); it2++) {
       if (*it2 >= (1e-16)) {
         outputMatrixFile << it2.index1() << " " << it2.index2() << " " << *it2
                          << endl;
       }
     }
   }
   outputMatrixFile.close();
 
   fstream outputProperty;
   outputProperty.open("prismProperty.ctl", fstream::out);
   outputProperty << "P=? [ (true) U (";
   bool first = true;
   int lhaloc = A.CurrentLocation;
   for (int it = 0; it < (int)A.NbLoc; ++it) {
     A.CurrentLocation = it;
     if (A.isFinal()) {
       if (first) {
         first = false;
       } else {
         outputProperty << "|";
       }
       outputProperty << "(automata = " << it << ")";
     }
   }
   A.CurrentLocation = lhaloc;
   outputProperty << ")]";
   outputProperty.close();
 
   fstream outputLabel;
   outputLabel.open("prismLabel.lbl", fstream::out);
   outputLabel << "0='init' 1='deadlock'\n0: 0";
   outputLabel.close();
 }
 
 void stateSpace::launchPrism(string prismPath) {
   cerr << "Starting Prism" << endl;
   string cmd = prismPath +
                " -gs -maxiters 1000000000 -ctmc -importtrans prismMatrix.tra "
                "-importstates prismStates.sta -importlabels prismLabel.lbl -v "
                "-cuddmaxmem 10000000 prismProperty.ctl > prismOutput";
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
 
   outputlumptmp << "bool SPN::precondition(const abstractMarking "
                    "&Marking){return true;}"
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
     cerr << "Fail to open muFile, load empty values" << endl;
     return;
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
 