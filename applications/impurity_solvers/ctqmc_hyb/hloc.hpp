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
#ifndef HLOC_H
#define HLOC_H

#include <map>
#include <set>
#include <vector>
#include <boost/python.hpp>
#include <triqs/arrays.hpp>
#include <triqs/arrays/linalg/det_and_inverse.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>
#include <triqs/gf/imfreq.hpp>
#include "small_matrix.hpp"
#include "ctqmc_utils.hpp"

// internal class for construction. See cpp file
namespace Hloc_construction { struct mydata; }

/**
   Blocks and local operators matrix elements for a local 
   Hamiltonian (with a small Hilbert space).
   Contains : 
    - the blocks (class Hloc::Bloc) with their dimensions, a unique number
      and the eigenvalues in ascending order.
    - a list of local operators, with the property that they connect 
      one bloc to one bloc, and their matrix elements.
*/
class Hloc { 
public : 
  /**
     The type of the matrices. Change this to use complex matrices
  */
  typedef double REAL_OR_COMPLEX;

  /**
     @param Hamiltonian : The Hamiltonian as python operator
     @param QuantumNumbers : a python dict name -> operator as python operator
     @param Symmetries : a list of dict :Cnumber -> Character as prepared by Operator.Transcribe_to_C_OpList
     @param SelectQN : truncation of the local Hilbert space. Default is None, i.e. lambda qn : true
  */
  Hloc(int NF, int NB, 
       boost::python::dict AllOperatorDict, 
       boost::python::dict QuantumNumbersList, 
       boost::python::list Symmetries, 
       boost::python::object SelectQN,
       int Nmaxbosons);

  /** Copy construction
   */
  Hloc (const Hloc & H1);

  /** 
      Blocks of the diagonalised local Hamiltonian.
  */
  class Bloc {
  public:
    /// Dimension of the bloc
    int dim;
    /// Number (unique) of the bloc
    int num;
    Bloc(int dim_, int num_);
    Bloc(const Bloc & B);
    Bloc& operator=(const Bloc&);
  protected:
    std::vector<double> H_;
    std::vector<double> deltaH_;
    friend class Hloc;
    friend struct Hloc_construction::mydata;
  public:
    /// Eigenvalues of the Hamiltonian in the bloc H[0] to H[dim-1]
    const double * restrict H;
    /// Eigenvalues of the Hamiltonian in the bloc shifted by the minimum of H
    const double * restrict deltaH;
    friend std::ostream & operator<< (std::ostream & out, const Bloc & B);
  };

 //****************************************************************

  /**
     Iterator on the Blocs
     Usage : (!! <b> A bit different from STL iterators </b>!!)
     for (Hloc::BlocIterator B(HLOC_OBJECT); !B.atEnd(); ++B)
       B->num is the number, etc....
   */
  typedef ctqmc_utils::myConstIteratorVector<Bloc> BlocIterator;


  //****************************************************************

  /** 
      An Operator in the local Hilbert Space connecting 
      one bloc to *exactly* ONE other bloc.
      Op[B].Btarget is the image of B by Op.
      Op[B].M is the corresponding matrix element.
  */
  class Operator { 
  public:
    /// Type of Statistic
    enum StatisticType {Fermionic,Bosonic};

    // to gather the results of  Op[B]
    class BlocMatrixElement {
    public :
      const Bloc * Btarget; 
      const SmallMatrix<REAL_OR_COMPLEX,ByLines> M;

      BlocMatrixElement():
        Btarget(NULL),M(0,0,NULL)  {}

      BlocMatrixElement(int Bdim, const Bloc * Btarget_, const REAL_OR_COMPLEX * M_ ) : 
	Btarget(Btarget_),M(Btarget->dim, Bdim, const_cast<REAL_OR_COMPLEX *>(M_)) { assert(Btarget);}

      BlocMatrixElement(const BlocMatrixElement & B) : Btarget(B.Btarget), M(B.M) {}
    };


    // BlocCorrespondance constructed as a vector of NULL. Filled later by Hloc.
    Operator(std::string name_, StatisticType stat, const std::vector<std::vector<REAL_OR_COMPLEX> > & MatrixElements);

    Operator(const Operator & Op);

    /// name of the Operator
    const std::string name;

    /// A unique number for each operator
    const int Number; 

    /// Is the operator bosonic or fermionic
    const StatisticType Statistic; 

    /** 
	Matrix elements of the operator. 
	Op[B].Btarget is the image of B by Op, maybe NULL.
	Op[B].M is the corresponding matrix element.
    */
    inline const BlocMatrixElement operator[] (const Bloc * B) const { 
      assert(B); 
      return (BlocCorrespondance[B->num] ? 
	      BlocMatrixElement(B->dim,BlocCorrespondance[B->num],&BlocMatrixElements[B->num][0])
	      :BlocMatrixElement() );
    }

    friend std::ostream & operator<< (std::ostream & out, const Operator & Op);
 
    const Operator & Transpose() const {return *transpose; }
 
  protected:
    std::vector<std::vector<REAL_OR_COMPLEX> > BlocMatrixElements; 
    std::vector<const Bloc *> BlocCorrespondance;
    friend class Hloc;
    friend struct Hloc_construction::mydata;
    friend class BlocMatrixElement;
    Operator * transpose; // where is my tranpose operator ?
 private:
    Operator & operator=(const Operator & Op);// private, so forbidden
    static int _number;
  };

  

  //****************************************************************

  /**
     DiagonalOperator : D[B] is a vector<VALTYPE> : longeur B->dim
  */
  class DiagonalOperator {
  public:
    // We need this ref to Hloc in the inner product of TraceSlice
    const Hloc & H;
  private:
    typedef std::vector<std::vector<REAL_OR_COMPLEX> > VV;
    VV data;
    static inline VV construct_data(const Hloc & h) { 
      VV res;
      for (BlocIterator B = h.BlocBegin(); B != h.BlocEnd(); ++B)
	res.push_back(std::vector<REAL_OR_COMPLEX>(B->dim,0));
      return res;
    }
  public:
    DiagonalOperator(const Hloc & h_):H(h_),data(construct_data(h_)) {}
    ///
    inline std::vector<REAL_OR_COMPLEX> & operator[] (const Bloc * B) { return data[B->num]; }
    inline const std::vector<REAL_OR_COMPLEX> & operator[] (const Bloc * B) const { return data[B->num]; }
  };

  //****************************************************************

  typedef std::map<std::string,Operator>::const_iterator OperatorIterator;
  inline OperatorIterator OperatorIteratorBegin() const { return OperatorMap.begin();}
  inline OperatorIterator OperatorIteratorEnd() const { return OperatorMap.end();}
 
  //****************************************************************

private:
  Hloc_construction::mydata * datatmp;
protected: 
  std::vector<Bloc> BlocList;
  std::map<std::string,Operator> OperatorMap,OperatorMapTranspose;
public:

  // do we guarantee that the Bloc will come in order of their ->num ?? ok ?
  inline BlocIterator BlocBegin() const { return BlocIterator(this->BlocList);} 
  inline BlocIterator BlocEnd() const { return BlocIterator(this->BlocList,true);}

  /// Returns the Operator of name s
  const Operator & operator[] (const std::string & s) const;
  
  /// Ground state energy
  const double E_GS;
  
  /**  Number of blocks */
  const int NBlocks;
  
  /**  Number C operators */
  const int NopC;
  
  /**  Number of operators */
  int N_Operators() const { return OperatorMap.size();}
 
  /** Sum of the dimension of the blocs */
  const int DimHilbertSpace;

  /**  Maximum Dimension of the blocks */
  const int MaxDimBlock;

  /// Max_OP sum _B B->dim * Op[B]->dim;
  int MaxOp_DimAllMatrixElements() const;
  
  /// Returns the local Green function
  void LocalGreenFunction(const Operator & OP1, const Operator & Op2, triqs::gf::gf_view<triqs::gf::imfreq> & G) const;

  /// Compute the paritition function at inverse temperature beta
  double PartitionFunction(double Beta) const; 

  friend std::ostream & operator<< (std::ostream & out, const Hloc &hloc);
  inline std::string print(){std::stringstream fs; fs<<*this<<std::endl;return fs.str();}
};

#endif
