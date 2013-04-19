
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "hloc.hpp"
#include "./iterator_python_sequence.hpp"
#include <algorithm>
#include <functional>
#include <numeric>

using namespace triqs::python_tools;
using namespace std;
using namespace boost;

// I define a more tolerant comparison between vectors for the AllBlocs map.
// This is really not so nice, using a key which is a vector<double> is
// probably a bad idea
struct lt {
  bool operator()(vector<double> v1, vector<double> v2) const {
    int i=0;
    for(vector<double>::iterator it = v1.begin(); it != v1.end(); ++it, ++i) {
      if(*it < (v2[i] - 1e-8)) {
        return true;
      } else if(v2[i] < (*it - 1e-8)) {
        return false;
      }
    }
    return false;
  }
};

typedef Hloc::Bloc Bloc;
typedef Hloc::Operator Operator;
typedef Hloc::REAL_OR_COMPLEX REAL_OR_COMPLEX;
typedef Hloc::Operator::BlocMatrixElement BlocMatrixElement;

// ******************************************************************

template<typename KEY, typename VAL>
inline VAL & myfind(map<KEY,VAL> & MAP, const KEY &  s) {
  typename map<KEY, VAL>::iterator p;
  p= MAP.find(s);
  assert (p != MAP.end());
  return p->second;
}

template<typename KEY, typename VAL>
inline const VAL & myfind(const map<KEY,VAL> & MAP, const KEY &  s) {
  typename map<KEY, VAL>::const_iterator p;
  p= MAP.find(s);
  assert (p != MAP.end());
  return p->second;
}

template<typename T>
bool operator < (const vector<T> & A, const vector<T> & B) { 
  assert(A.size() == B.size());
  for (uint i=0; i<A.size(); ++i) { if (abs(A[i] -B[i])>1.e-10) return (A[i] < B[i]); }
  return false; // they are equal, not <
};

// ******************************************************************

Bloc::Bloc(int dim_, int num_):dim(dim_), num(num_),
			       H_(dim,0),deltaH_(dim,0),
			       H(&H_[0]), deltaH(&deltaH_[0]) {}

Bloc::Bloc(const Bloc & B):
  dim(B.dim), num(B.num),
  H_(B.H_),deltaH_(B.deltaH_),
  H(&H_[0]), deltaH(&deltaH_[0]) {}

Bloc& Bloc::operator=(const Bloc&X) {
  dim = X.dim;
  num = X.num;
  H_ = X.H_;
  deltaH_ = X.deltaH_;
  H = &H_[0]; deltaH = &deltaH_[0];
  return *this;
}

ostream & operator<< (ostream & out, const Bloc & B) {
  out<< "number :"<<B.num<< ", dimension : "<<B.dim<<endl;
  out<< "Eigenvalues : ";
  std::copy(B.H_.begin(),B.H_.end(),ostream_iterator<double>(out," "));
  out<<endl;
  return out;
}

// ******************************************************************
int Operator::_number =0;

Operator::Operator(string name_, StatisticType stat,
		   const vector<vector<Hloc::REAL_OR_COMPLEX> > & MatrixElements):
  name(name_), Number(_number++),Statistic(stat),
  BlocMatrixElements(MatrixElements),
  BlocCorrespondance(MatrixElements.size()),
  transpose(NULL) {}

Operator::Operator(const Operator & Op):
  name(Op.name),Number(Op.Number),Statistic(Op.Statistic),
  BlocMatrixElements(Op.BlocMatrixElements),
  BlocCorrespondance(Op.BlocCorrespondance),
  transpose(NULL) {}

ostream & operator<< (ostream & out, const Operator & Op) {
  out<<"Operator "<< Op.name << " number : "<<Op.Number<<endl<<" statistic : "<<(Op.Statistic==Hloc::Operator::Fermionic ? "Fermionic" : "Bosonic")<<endl;
  for (uint num =0; num<Op.BlocCorrespondance.size(); ++num) { 
    const Bloc *Bp = Op.BlocCorrespondance[num];
    out<< "  "<< num<< " --> "; 
    if (Bp)  {
      out<<Bp->num;
      //int n1 = Op.BlocMatrixElements[num].size()/Bp->dim;// what a mess !
      /*
      out <<"  with matrix element :  "<<endl
	  <<BlocMatrixElement(n1,Op.BlocCorrespondance[num], 
			      &Op.BlocMatrixElements[num][0]).M<<endl;
      */


    }
    else out<<"None"<<endl; 
  }
  return out;
}

// ******************************************************************

// I put here data structure useful only for this file. Hence in a special name space.
namespace Hloc_construction { 

  // The Fundamental (canonical) operators bosonic or fermionic
  struct OpFundamental { 
    int number; 
    static int NF,NB;      // number of fermions, bosons
    static int MaxNbosons; // maximum number of bosons (to truncate the Hilbert Space)
    bool isdagger; 
    bool isfermion() { return (number <NF); }
  };
  int OpFundamental::NF, OpFundamental::NB,OpFundamental::MaxNbosons;

  /*
    A pure state of the local Hilbert space
    F : for the fermionic Fock space : An integer n in [0,2^Nc], where Nc is the number of C operators
    is a basis state (using usual binary decomposition)
    B : for the bosons.
  */
  struct PureState { 
    int F;
    vector<int> B; int NBmax;
    // Order to use them as keys in maps.
    bool operator < (const PureState & X) const { 
      if (F != X.F) return (F<X.F);
      return (B<X.B);
    }

    // action of a fundamental operator : returns the sign : -1, 1 or 0 if 0 vector is 0.
    int apply(OpFundamental C) { 
      //cout<<"apply "<<"  "<<F << " "<<C.number <<endl;
      int sign = 0;
      if (C.isfermion()) { 
	int NC= 1<<C.number;
	if (C.isdagger) // F & NC : bit number C.number is at 1
	  { if (F & NC) return 0; else F ^= NC;}
	else
	  { if (!(F & NC)) return 0; else F ^= NC;}
	// compute signe
	for (int i =0, k = F; i<C.number; ++i, k>>=1) { if (k & 1) sign++;}
      }
      else { // operator is bosonic
	if (C.isdagger)
	  { if (B[C.number] < NBmax) B[C.number]++; else  return 0;}
	else
	  { if (B[C.number] >0)  B[C.number]--; else return 0;}
      }
      //cout<<"apply "<< F << " "<<C.number << " IS "<< (sign%2 == 0 ? 1 : -1)<<endl;
      return (sign%2 == 0 ? 1 : -1);
    }

    /// Given a symmetry defined by : number of the Cdagger -> character
    std::complex<double> ActWithSymmetry(const map<int,std::complex<double>> & S) const { 
      std::complex<double> res=1;
      for (map<int,std::complex<double>>::const_iterator p=S.begin(); p!=S.end(); ++p) { 
	int NC = 1<<p->first;
	// F & NC : bit number C.number is at 1
	if (F & NC) res *= p->second;
      }
      return res;
    }

    // Just to print the state clearly
    friend std::ostream & operator << (std::ostream & out, const PureState & p) {
      for (int i = OpFundamental::NB-1; i >= 0; --i) out << p.B[i] << " ";
      for (int i = OpFundamental::NF-1; i >= 0; --i) out << (p.F>>i & 1) << " ";
      return out;
    }

  };
  
  //--------------------------------------
  // An "iterator" that sweeps over all PureStates
  class PureStateGenerator { 
    PureState ps;
    bool fini;
  public : 
    PureStateGenerator(int NmaxBosons){
      OpFundamental::MaxNbosons = NmaxBosons;
      ps.F = 0; ps.B.resize(OpFundamental::NB,0);
      fini=false;
    }
    //------------------
    const PureState & operator*(){ return (ps);}
    //------------------
    const PureState * operator->(){ return (&ps);}
    //------------------
    inline bool atEnd(){ return fini;}
    //------------------
    PureStateGenerator& operator++() {

      int NF = OpFundamental::NF, NB = OpFundamental::NB;
      ++ps.F;
      if (ps.F == 1<<NF) {
        ps.F = 0; int i = 0;
        while ( (i < NB) && (ps.B[i] == OpFundamental::MaxNbosons) ) { ps.B[i] = 0; ++i; }
        if (i == NB) { fini = true; return *this; }
        else { ++ps.B[i]; }
      }
      return *this;

    }
  };
  //-----------------------------------------------
  
  /* 
     A state of the Fock space is a sum_i lambda_i |i> where |i> is a pure state
     so it is simply a map i-> lambda_i
  */
  typedef map<PureState, REAL_OR_COMPLEX> State;
  typedef map<PureState, REAL_OR_COMPLEX>::const_iterator StateIterator;
  
  /*
    Now I need to define the operators and to apply them to the states.
    This class is different from the Operator class, which stores
    only the minimal amount of information, to optimize the code.
  */
  struct FullOperator { 
    
    // structure of the operator is [ alpha * c1*c2*c3*c4], where ci are OpFundamental
    vector< std::pair<REAL_OR_COMPLEX, vector<OpFundamental > > > mystruct;
    // Images of a bloc by the Operator
    vector<int> BlocCorrespondance;
    Hloc::Operator::StatisticType Statistic;

    //----------------------------------------------------------

    FullOperator(python::list L){
      // statistic : 0, undefined, 1, Fermion, 2 boson
      int stat=0;
      // transcribe the python structure
      for (IteratorOnPythonListOf2Tuples< REAL_OR_COMPLEX, python::list > p(L); !p.atEnd(); ++p) {
	int newstat = 2 - python::len(p->x2) %2; 
	if (stat !=0) { if (stat !=newstat) TRIQS_RUNTIME_ERROR<<"An Operator must be Fermionic or Bosonic";}
	else stat = newstat;
	vector<OpFundamental> res;
	for (IteratorOnPythonList<int> C(p->x2); !C.atEnd(); ++C) { 
	  OpFundamental OP; 
	  OP.number = abs(*C) - 1; // start from 1 in the python to store the sign !
	  OP.isdagger = (*C<0); // dagger have negative number. Cf python procedure.
	  res.push_back(OP);
	}
	mystruct.push_back(make_pair(p->x1,res));
      }
      Statistic = (stat==1 ? Operator::Fermionic : Operator::Bosonic); // empty is one i.e. bosonic.
    }

    //-----------------------------------------------------------

    /// Application of the Operator on a Pure State  : we don't need a general state
    void operator() (const PureState & psi, State & res) const { 
      res.clear();
      for (vector< std::pair<REAL_OR_COMPLEX, vector<OpFundamental> > >::const_iterator term = mystruct.begin(); 
	   term != mystruct.end(); ++term) { 
	PureState etat(psi);
	int sign = 1;
	// iterate on the monomial
	for (vector<OpFundamental>::const_iterator C = term->second.begin(); sign && (C !=term->second.end()); ++C) 
	  sign *= etat.apply(*C);
	if (sign !=0) res[etat] += sign*term->first;
      }
    }
  }; // end FullOperator
  
  //***********************************************************************
  // this struct is a temporary one where all the calculation are actually done.
  // afterwards, all vector will be copied into the Hloc ones, and hence
  // compacted in memory (better for cache).
  
  struct mydata {
    
    // BlocContents [B->num] is the list of the Pure Fock state contained in B. This vector defined their order.
    vector< vector<PureState> > BlocContents;
    
    // number of the pure Fock state  ---> (Bloc number, number in the bloc)
    map<PureState, std::pair<int, int> > MapStateBlocs;
    
    // is a state s entirely in Bloc B ?
    inline bool IsStateInBloc(const State & s, const Bloc * B) {
      for (StateIterator p=s.begin(); p!=s.end(); ++p)
      	if (MapStateBlocs[p->first].first != B->num) return false;
      return true;
    }
    
    // Data constructed for Hloc : 
    vector<Bloc> BlocList;
    map<string,Operator> OperatorMap,OperatorMapTranspose;
    double E_GS; int dimtot,MaxDimBlock;
    
    // temporary data 
    vector < triqs::arrays::matrix<REAL_OR_COMPLEX> *> Pinv, P;
    map<string,FullOperator> FullOperatorMap;

    //---------------------------

    ~mydata() { 
      assert(Pinv.size()==P.size());
      for (uint i=0; i<Pinv.size(); ++i) {delete Pinv[i]; delete P[i];}
    }
     
    //---------------------------------
    
    mydata(int NF, int NB, 
	   python::dict AllOperatorDict, 
	   python::dict QuantumNumbersList, 
	   python::list Symmetries, 
	   python::object SelectQN,
	   int Nmaxbosons)
    { 
      OpFundamental::NF = NF;OpFundamental::NB = NB;
      for (IteratorOnPythonDict < string , python::list > p(AllOperatorDict); !p.atEnd(); ++p)
	FullOperatorMap.insert(make_pair(p->key,FullOperator(p->val)));
      
      // reinit the Operator::_number
      Operator::_number=0;

      // I construct the bloc by sorting the pure state by their quantum numbers.
      // AllBlocs a map (string key computed from quantum numbers) -> a list of PureStates
      // with these quantum numbers, after truncation by SelectQN
      map<vector<double> , vector<PureState>, lt > AllBlocs; 
      
      // transcribe the symmetries into C++
      vector< map<int,std::complex<double>> > SymChar;
      for (IteratorOnPythonList<python::dict> p(Symmetries); !p.atEnd(); ++p) {
	SymChar.push_back(map<int,std::complex<double>>());
	for (IteratorOnPythonDict<int,std::complex<double>> it(*p); !it.atEnd(); ++it) 
	  SymChar.back()[it->key] = it->val;
      }

      // We iterate on all PureStates
      for (PureStateGenerator PS(Nmaxbosons); !PS.atEnd(); ++PS) { 
	python::list allqns_py; vector<double> allqns;
	stringstream fs; fs.setf(ios::fixed,ios::floatfield); fs.precision(3);
	State s2;
	// iterate on QN.	
	for (IteratorOnPythonDict<string,python::object> p(QuantumNumbersList); !p.atEnd(); ++p) {
	  myfind(FullOperatorMap,p->key)(*PS,s2);
	  if (s2.size()>1) TRIQS_RUNTIME_ERROR<<"Hloc : "<<p->key<<" is not a quantum number or  it does not leave pure state pure";
	  std::complex<double> val = (s2.size()==0 ? 0 : s2.begin()->second);
	  allqns_py.append(val);
	  allqns.push_back(real(val));allqns.push_back(imag(val));
	}
	// iterate on QN given by .Symmetry action
	for (uint u =0; u<SymChar.size(); ++u) { 
	  std::complex<double> val = PS->ActWithSymmetry(SymChar[u]);
	  allqns_py.append(val);
	  allqns.push_back(real(val));allqns.push_back(imag(val));
	}
	//copy(allqns.begin(), allqns.end(),ostream_iterator<double>(cout, " "));cout<<endl;
	// I filter out the pure state where SelectQN  is false if it exists.
	//if  ( (SelectQN ==python::PythonNone()) || python::extract<bool>(SelectQN(allqns_py))) 
	if  ( (SelectQN.is_none()) || python::extract<bool>(SelectQN(allqns_py))) 
	  AllBlocs[allqns].push_back(*PS);
      }
          
      // I now have the blocks, so I can build BlocList and MapStateBlocs
      // Build the blocs BlocList :  bloc.number -> bloc *
      int NBlocs =AllBlocs.size();
      P.resize(NBlocs,NULL); Pinv.resize(NBlocs,NULL);
      int num=0; // number of the bloc
      for (map< vector<double> , vector< PureState >, lt >::const_iterator p=AllBlocs.begin(); p!=AllBlocs.end(); ++p, ++num) {
	Bloc B(p->second.size(),num);    
	assert(B.dim>0);
	BlocList.push_back(B);
	BlocContents.push_back(p->second);
	for (int i=0; i<B.dim; ++i)
	  MapStateBlocs[p->second[i]] = std::make_pair(num,i);
      } 
    
      // Now diagonalize the H
      for (vector<Bloc>::iterator B = BlocList.begin(); B!=BlocList.end(); ++B) {
        triqs::arrays::matrix<REAL_OR_COMPLEX> Hmat(B->dim,B->dim); Hmat() =0;
	const vector<PureState> & BContent (BlocContents[B->num]);  // A basis of the blocs as pure states.
	FullOperator & Hop(myfind(FullOperatorMap,string("Hamiltonian")));    // the operator
	for (uint j=0;j<BContent.size(); ++j) {  // for all vector of the basis
	  State S2;
	  Hop(BContent[j],S2);                  // action of Operator
	  if (S2.size()==0) continue; // operator gives 0. no matrixelement computed (more exactly an empty one)
	  if (!IsStateInBloc(S2,&(*B)))  TRIQS_RUNTIME_ERROR<<"Hamiltonian is not diagonal in the blocks";
	  for (State::const_iterator p=S2.begin(); p != S2.end(); ++p)  { 
	    Hmat(MapStateBlocs[p->first].second, int(j) ) = p->second;
	  }
	}
	// call lapack to diagonalize
	triqs::arrays::vector<double> ev(B->dim);
        triqs::arrays::matrix<REAL_OR_COMPLEX> temp(Hmat);
        std::tie(ev,Hmat) = triqs::arrays::linalg::eigenelements(temp());
	// ev is the list of eigenvalues.
	// prepare to sort them
	vector< std::pair<double,int> > tmp(B->dim);
	for (int i =0; i<B->dim; ++i)  { 
	  tmp[i] = std::make_pair(ev(i),i);
	  // cout<<" PERM"<<i <<"  "<< tmp[i].first<< "  "<<tmp[i].second<<endl;
	}

	std::sort(tmp.begin(),tmp.end());
	for (int i=0; i<B->dim; i++) { B->H_[i]= tmp[i].first; B->deltaH_[i]=B->H[i] - B->H[0];}
      
	// I need to modify the P matrix
	//Hmat = Hmat.transpose();
	triqs::arrays::matrix<REAL_OR_COMPLEX> Mtmp(Hmat);
	for (int i=0; i<B->dim; i++)
	  for (int j=0; j<B->dim; j++)
	    {
	      //   cout<<" PERM"<<i <<"  "<< tmp[i].first<< "  "<<tmp[i].second<<endl;
	      Mtmp ( tmp[i].second, tmp[j].second) = Hmat(i,j);
	    }

	Pinv[B->num] = new triqs::arrays::matrix<REAL_OR_COMPLEX>(Mtmp);
	P[B->num]    = new triqs::arrays::matrix<REAL_OR_COMPLEX>(Mtmp);
	*(P[B->num]) = inverse(*(P[B->num]));

      }
      // end diagonalization of Hamiltonian
 
      // Compute the elements matrices of all operators
      for (map<string,FullOperator>::iterator OP= FullOperatorMap.begin(); OP != FullOperatorMap.end(); ++OP) {
      
	// a vector of MatrixElements stored in 1d vectors
	vector<vector<REAL_OR_COMPLEX> > AllMatrixElements;
      
	for (vector<Bloc>::const_iterator B = BlocList.begin(); B != BlocList.end(); ++B) {
	  //Compute the matrix element of an operator OP 
	  //between bloc B and its image Bp in the Fock basis.
	  ///It checks that B is mapped to exactly one bloc.
	  const Bloc * Bp=NULL;
	  for (uint j=0;j<uint(B->dim); ++j) {  // for all vector of the basis
	    State S2;
	    OP->second(BlocContents[B->num][j],S2);                  // action of Operator
	    if (S2.size()==0) { continue;} // operator gives 0. no matrix element computed (more exactly an empty one)
	    if (Bp==NULL)  { // first determination of Bp
	      Bp= &BlocList[MapStateBlocs[S2.begin()->first].first]; 
	      AllMatrixElements .push_back(vector<REAL_OR_COMPLEX> (B->dim *Bp->dim,0));
	    }
	    if (!IsStateInBloc(S2,Bp))  TRIQS_RUNTIME_ERROR<<"Operator "<<OP->first<<" does not connect one bloc to one bloc";
	    
	    for (State::const_iterator p=S2.begin(); p != S2.end(); ++p)  {
	      SmallMatrix<REAL_OR_COMPLEX,ByLines> SM(Bp->dim, B->dim,AllMatrixElements.back());
	      SM(MapStateBlocs[p->first].second, int(j) ) = p->second;
	    }
	  }
	  if (Bp==NULL)  { 
	    AllMatrixElements .push_back(vector<REAL_OR_COMPLEX> (0));
	    OP->second.BlocCorrespondance.push_back(-1);
	  }
	  else { 
	    OP->second.BlocCorrespondance.push_back(Bp->num);
	    SmallMatrix<REAL_OR_COMPLEX,ByLines> SM  (Bp->dim, B->dim,AllMatrixElements .back());
            triqs::arrays::matrix<REAL_OR_COMPLEX> res(Bp->dim, B->dim); res() = 0.0;
            for (int i=0; i<Bp->dim; i++)
              for (int j=0; j<B->dim; j++)
                for (int k=0; k<Bp->dim; k++)
                  for (int l=0; l<B->dim; l++)
                    res(i,j) += (*(Pinv[Bp->num]))(i,k) * SM(k,l) * (*(P[B->num]))(l,j);
            for (int i=0; i<Bp->dim; i++)
              for (int j=0; j<B->dim; j++) SM(i,j) = res(i,j);
	  }
	}
	OperatorMap.insert(make_pair(OP->first,Operator(OP->first,OP->second.Statistic,AllMatrixElements )));
      } // loop on operators.
    
      // Construct the transpose of the operators
      for (map<string,Operator>::iterator OP = OperatorMap.begin(); OP != OperatorMap.end(); ++OP) {

        vector<vector<REAL_OR_COMPLEX> > AllMatrixElementsT(BlocList.size());

	for (vector<Bloc>::const_iterator B = BlocList.begin(); B != BlocList.end(); ++B) {

          int n = myfind(FullOperatorMap,OP->first).BlocCorrespondance[B->num];
          if (n != -1) {
            const Bloc * Bp = &BlocList[n];
	    AllMatrixElementsT[Bp->num].resize(B->dim *Bp->dim,0);
            SmallMatrix<REAL_OR_COMPLEX,ByLines> SM (B->dim, Bp->dim, AllMatrixElementsT[Bp->num]);
            SmallMatrix<REAL_OR_COMPLEX,ByLines> SM1 (Bp->dim, B->dim, OP->second.BlocMatrixElements[B->num]);
            for (int i=0; i<B->dim; i++)
              for (int j=0; j<Bp->dim; j++) SM(i,j) = SM1(j,i);
          }
        }
	OperatorMapTranspose.insert(make_pair(OP->first,Operator(OP->first,OP->second.Statistic,AllMatrixElementsT)));
      }

      E_GS = BlocList[0].H[0];
      dimtot = 0; 
      vector<int> vn;
      for (Hloc::BlocIterator B(BlocList); !B.atEnd(); ++B) { 
	E_GS = min(E_GS, B->H[0]);
	dimtot += B->dim;
	vn.push_back(B->dim);
      }
      MaxDimBlock = *max_element(vn.begin(), vn.end());

    }
  }; //struct
}; //namespace

using namespace Hloc_construction;

// ******************************************************************

Hloc::Hloc(int NF, int NB, 
	   python::dict AllOperatorDict, 
	   python::dict  QuantumNumbersList, 
	   python::list Symmetries, 
	   python::object SelectQN,
	   int Nmaxbosons):
  datatmp(new Hloc_construction::mydata(NF,NB,AllOperatorDict,QuantumNumbersList,Symmetries,SelectQN,Nmaxbosons)),
  BlocList(datatmp->BlocList),
  OperatorMap(datatmp->OperatorMap),  
  OperatorMapTranspose(datatmp->OperatorMapTranspose),  
  E_GS(datatmp->E_GS),
  NBlocks(BlocList.size()),
  NopC(NF),
  DimHilbertSpace(datatmp->dimtot),
  MaxDimBlock(datatmp->MaxDimBlock)
{ 
  // Now reconnect all the bloc correspondance of operators.
  for ( map<string,Operator>::iterator p = OperatorMap.begin(); p != OperatorMap.end(); ++p) {
    FullOperator & fullop(myfind(datatmp->FullOperatorMap,p->second.name));
    for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B) { 
      int n = fullop.BlocCorrespondance[B->num];
      p->second.BlocCorrespondance[B->num] = (n==-1 ? NULL : &BlocList[n]);
    }
    // link the tranposed operator
    p->second.transpose = &(myfind(OperatorMapTranspose,p->first));
  }

 // Now reconnect all the bloc correspondance of the TRANSPOSED operators.
  for ( map<string,Operator>::iterator p = OperatorMapTranspose.begin(); p != OperatorMapTranspose.end(); ++p) {
    FullOperator & fullop(myfind(datatmp->FullOperatorMap,p->second.name));
    vector <int> BlocCorresInverse(fullop.BlocCorrespondance.size(),-1);
    // cout<<"INVERINTB"<< p->second.name<<endl;
    //std::copy(fullop.BlocCorrespondance.begin(),fullop.BlocCorrespondance.end(),ostream_iterator<int>(cout," "));
    // inverse the bloc correspondance table of the block of Operator fullop
    for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B) {
       int Bpnum = fullop.BlocCorrespondance[B->num];
       if (Bpnum!=-1)  { 
	 assert(BlocCorresInverse[Bpnum] ==-1);
	 BlocCorresInverse[Bpnum] = B->num;
       }
    }
    for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B) { 
      int n = BlocCorresInverse[B->num];
      p->second.BlocCorrespondance[B->num] = (n==-1 ? NULL : &BlocList[n]);
    }
    // link the tranposed operator
    p->second.transpose = &(myfind(OperatorMap,p->first));
  }
  delete datatmp;

  // just to check that the iterator (or BlocList) is really in the order of B->num
  int u =0;
  for (BlocIterator B = BlocBegin(); !B.atEnd(); ++B, ++u) assert(B->num==u);

}
// ******************************************************************

Hloc::Hloc(const Hloc & H1):
  BlocList(H1.BlocList),
  OperatorMap(H1.OperatorMap),  
  OperatorMapTranspose(H1.OperatorMapTranspose),  
  E_GS(H1.E_GS),
  NBlocks(BlocList.size()),
  NopC(H1.NopC),
  DimHilbertSpace(H1.DimHilbertSpace),
  MaxDimBlock(H1.MaxDimBlock)
{ 
  assert(0); // not tested
  // Now reconnect all the bloc correspondance of operators.
  for ( map<string,Operator>::iterator p = OperatorMap.begin(); p != OperatorMap.end(); ++p) {
    for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B) { 
      const Bloc * Bp = myfind(H1.OperatorMap,p->first).BlocCorrespondance[B->num];
      p->second.BlocCorrespondance[B->num] = (Bp == NULL ? NULL : &BlocList[Bp->num]);
    }
    // link the tranposed operator
    p->second.transpose = &(myfind(OperatorMapTranspose,p->first));
  }

  // same for the transposed
  // Now reconnect all the bloc correspondance of operators.
  for ( map<string,Operator>::iterator p = OperatorMapTranspose.begin(); p != OperatorMapTranspose.end(); ++p) {
    for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B) { 
      const Bloc * Bp = myfind(H1.OperatorMapTranspose,p->first).BlocCorrespondance[B->num];
      p->second.BlocCorrespondance[B->num] = (Bp == NULL ? NULL : &BlocList[Bp->num]);
    }
    // link the tranposed operator
    p->second.transpose = &(myfind(OperatorMap,p->first));
  }
 
}

//-------------------------------------------------

ostream & operator<< (ostream & out, const Hloc &hloc) {
  // Final print
  vector<int> dims;
  for (Hloc::BlocIterator B = hloc.BlocBegin(); B != hloc.BlocEnd(); ++B) dims.push_back(B->dim);
  //std::transform(hloc.BlocList.begin(), hloc.BlocList.end(),std::back_inserter(dims),bind(&Bloc::dim, boost::lambda::_1));
  std::sort(dims.begin(),dims.end());
  int sum = std::accumulate(dims.begin(),dims.end(),0);
  out << "All blocs dimensions :";
  std::copy(dims.begin(),dims.end(),ostream_iterator<int>(out,","));
  out<<endl;
  out << "Check total dimension : "<<sum<<endl<<"Diagonalisation done "<<endl;

  // Now print all the blocs....
  out << "Details of the blocs"<<endl;
  std::copy(hloc.BlocList.begin(), hloc.BlocList.end(),ostream_iterator<Bloc>(out,"---------\n"));

 // Now print all the operators
  out << "Details of the operators"<<endl;
  for ( map<string,Operator>::const_iterator p = hloc.OperatorMap.begin(); p != hloc.OperatorMap.end(); ++p) {
    out<<p->second<<endl<<"--------------"<<endl;
    out<<p->second.Transpose()<<endl<<"--------------"<<endl;
  }
  return out;
}

// ******************************************

int Hloc::MaxOp_DimAllMatrixElements() const { 
  int res =0;
  for (map<string,Operator>::const_iterator p = OperatorMap.begin(); p!= OperatorMap.end(); ++p) { 
    int s =0;
    for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B)
      s += B->dim * p->second[B].Btarget->dim;
    res+=s;
  }
  return res;
}

// ********************************************************************

double Hloc::PartitionFunction(double Beta) const { 

  
  // Compute the partition function 
  double Z = 0;
  for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B)
    for (int i =0;  i<B->dim; ++i)
      Z += exp( -Beta * ( B->H[i] - E_GS )); 
  return Z;
  
}

// ********************************************************************


void Hloc::LocalGreenFunction(const Operator & Op1, const Operator & Op2, triqs::gf::gf_view<triqs::gf::imfreq> & G) const {
  
  const double Beta(G.domain().beta);

  // Compute the partition function 
  double Z = PartitionFunction(Beta);
  
  // sum over blocs A, B
  for (BlocIterator A = BlocBegin(); A != BlocEnd(); ++A)
    for (BlocIterator B = BlocBegin(); B != BlocEnd(); ++B)
      if ( (Op1[A].Btarget ==  B) && ((Op2[B].Btarget == A) )) {
	// sum over all states in the blocs
	for (int ai=0; ai<A->dim; ++ai)
	  for (int bi=0; bi<B->dim; ++bi) {
	    double Ea = A->H[ai] - E_GS, Eb = B->H[bi] - E_GS;
	    double fact = Op1[A].M(bi,ai) * Op2[B].M(ai,bi) * (exp(-Beta*Eb) + exp(-Beta*Ea ))/Z ;
	    cout<<fact<<endl;
	    // for omega : G(omega) = fact/(omega + Ea - Eb);
	  }
      }
}

// ********************************************************************

const Operator & Hloc::operator[] (const string & s) const { 
  map<string,Operator>::const_iterator p = OperatorMap.find(s);
  if (p==OperatorMap.end()) TRIQS_RUNTIME_ERROR<<"Operator "<<s<<" not found";
  return p->second;
}



