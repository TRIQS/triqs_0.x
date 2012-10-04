
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

#include "GF_Bloc_Base.hpp"
#include <sys/stat.h>
#include <list>
#include <boost/python/stl_iterator.hpp>

using python::extract;
using python::object;
using python::dict;

template<typename T> 
inline T & _check_sp(boost::shared_ptr<T> const & p) { assert(p); return *p;}

template <typename DataType> 
GF_Bloc_Base<DataType>::GF_Bloc_Base (python::object IndicesL_,
				      python::object IndicesR_,
				      PyObject * Data,
				      boost::shared_ptr<MeshGF> Mesh,
				      boost::shared_ptr<TailGF> Tail):
  IndicesL(IndicesL_),
  IndicesR(IndicesR_),
  mesh_ptr(Mesh),
  tail_ptr(Tail),
  mesh(_check_sp(mesh_ptr)),
  tail(_check_sp(tail_ptr)),
  Beta(mesh.Beta),
  Statistic(mesh.Statistic),
  N1(int(PySequence_Size(IndicesL.ptr()))),// int for 64bits
  N2(int(PySequence_Size(IndicesR.ptr()))),// int for 64bits
  data_ptr( Data==Py_None ? 
	    new PyArray<DataType,3>(N1,N2,mesh.index_max -  mesh.index_min +1  ,FortranOrder):
	    new PyArray<DataType,3>(Data,"",SameOrder,AssertView)),
  data(*data_ptr),
  data_const(*data_ptr)
{
  if (Data!=Py_None) assert (data_ptr->as_PyObjectPtrBorrowedRef()  == Data);
  //data.dumpStructureInformation(cout);

  //if (!is_fortran_ordering(data)) TRIQS_RUNTIME_ERROR<<"Oops... Data is not Fortran order ! Internal error";
 
  // check the size !
  if (data.shape() !=TinyVector<int,3>(N1,N2,mesh.index_max -  mesh.index_min +1)) {
    cout << data.shape() << " versus " << TinyVector<int,3>(N1,N2,mesh.index_max -  mesh.index_min +1) << endl;
    TRIQS_RUNTIME_ERROR<<"GF_Bloc_Base : Data array provided is of incorrect dimension !!";
    }
  data.reindexSelf(TinyVector<int,3>(1,1,mesh.index_min));
}
 
//------------------------------------------------------------------------------------------------

template <typename DataType>
GF_Bloc_Base<DataType>::GF_Bloc_Base(const GF_Bloc_Base & Gin):
  IndicesL(Gin.IndicesL),
  IndicesR(Gin.IndicesR),
  mesh_ptr(Gin.mesh_ptr),
  tail_ptr(Gin.tail_ptr),
  mesh(_check_sp(mesh_ptr)),
  tail(_check_sp(tail_ptr)),
  Beta(Gin.Beta),
  Statistic(Gin.Statistic),
  N1(Gin.N1),
  N2(Gin.N2),
  data_ptr(new PyArray<DataType,3>(Gin.data)),
  data(*data_ptr),
  data_const(*data_ptr)
{ 
  data.reindexSelf(TinyVector<int,3>(1,1,mesh.index_min));
}

//------------------------------------------------------------------------------------------------

template <typename DataType>
GF_Bloc_Base<DataType>::~GF_Bloc_Base() { delete data_ptr;}

//------------------------------------------------------------------------------------------------

template <typename DataType>
void GF_Bloc_Base<DataType>::operator= (const GF_Bloc_Base<DataType> & Gin)
{
  check_have_same_structure(*this,Gin);
  //IndicesList = Gin.IndicesList; Indices are NOT copied
  (Array<DataType,3>)data = (Array<DataType,3>)Gin.data;
  tail = Gin.tail;
}

//------------------------------------------------------------------------------------------------

inline vector<string>  my_indices(const python::object & IND)
{
  vector<string> v;
  v.push_back(""); // to start at 1
  python::stl_input_iterator<python::object> begin(IND), end;
  std::list<python::object>  l(begin,end);
  
  for (std::list<python::object>::const_iterator p = l.begin(); p !=l.end(); ++p) { 
    python::str s(*p);
    s =  s.replace(" ","").replace("(","").replace(")","").replace("'","").replace(",","-");
    v.push_back( PyString_AsString(s.ptr()));
  }
  return v;
}

//-------------------------------------------------------------

inline double real(double x){ return x;}
inline double imag(double x){ return 0;}
template <typename DataType> inline DataType myset(double r2, double r3);
template<> inline double myset<double>(double r2, double r3) { return r2;}
template<> inline COMPLEX myset<COMPLEX> (double r2, double r3) {return (r2 + I *r3);}

template <typename DataType>
void GF_Bloc_Base<DataType>::save(string file, bool accumulate) const
{

  vector<string> indices_namesL = my_indices(IndicesL),indices_namesR = my_indices(IndicesR);
  int step = 1;
  bool NewStyleSave = true;
  
  if ( (N1*N2>1) && NewStyleSave) system("mkdir -p " + file);
  for (int n1=1; n1<=N1;n1++) {
    for (int n2=1; n2<=N2;n2++) {
      stringstream fs; 
      if ( N1*N2>1) fs<<file<<(NewStyleSave ? "/" : "_") <<indices_namesL[n1]<<"_"<<indices_namesR[n2]<<".dat";
      else fs<<file<<".dat";
      ofstream f(fs.str().c_str(), (accumulate ? ios::out|ios::app : ios::out));
      f.setf(ios::fixed,ios::floatfield);f.precision(PRECISION_OUTPUT);f<<endl<<endl;
      for (int i=mesh.index_min;i<=mesh.index_max;i+=step)  {
        switch (mesh.typeGF) { 
          case Imaginary_Time : 
          case Imaginary_Legendre :
            f<<real(mesh[i])<<"  "<<data(n1,n2,i) <<endl;
            break;
          case Imaginary_Frequency : 
            f<<imag(mesh[i])<<"  "<<real(data(n1,n2,i))<<"   "<<imag(data(n1,n2,i)) <<endl;
            break;
          default : 
            f<<real(mesh[i])<<"  "<<real(data(n1,n2,i))<<"   "<<imag(data(n1,n2,i)) <<endl;
            break;
        }      
      }
    }
  }
  tail.save(file,accumulate);

}


//-------------------------------------------------------------

template <typename DataType>
void GF_Bloc_Base<DataType>::load(string filename)
{
  vector<string> indices_namesL = my_indices(IndicesL),indices_namesR = my_indices(IndicesR);
  tail.load(filename);
  bool NewStyleSave = true;

  Array<DataType,1> G_complete(Range(mesh.index_min,mesh.index_max));
  double r1,r2,r3,r_check=0;
  for (int n1=1; n1<=N1;n1++) 
    for (int n2=1; n2<=N2;n2++)
      {
        G_complete = 0;
          {
            stringstream fs;
            if (N1*N2>1) fs<<filename<<(NewStyleSave ? "/" : "_")<<indices_namesL[n1]<<"_"<<indices_namesR[n2]<<".dat";
            else fs<<filename<<".dat";
            ifstream f(fs.str().c_str());
            if (f.fail()) 
              { 
                TRIQS_RUNTIME_ERROR<<"Error in GFBloc::load : file " << fs.str() << "Not found : Putting 0";
              }
            else
              {
                for (int i=mesh.index_min;i<=mesh.index_max;i++)
                  {
                    if (f.fail()) { REPORT <<"NOT enough data : completing with 0"<<endl; break;}
		    if ((mesh.typeGF == Imaginary_Time)||(mesh.typeGF == Imaginary_Legendre)) {
                      f >> r1 >> r2;
                      r_check = real(mesh[i]);
                      r3=0;
                    } else {
                      f >> r1 >> r2 >> r3;
                      if (mesh.typeGF == Imaginary_Frequency) r_check = imag(mesh[i]);
                      if (mesh.typeGF == Real_Frequency) r_check = real(mesh[i]);
                    }
		    if (abs(r1 - r_check) >1e-6) 
		      {
			TRIQS_RUNTIME_ERROR<<"Error in GFB_array::load : mesh index incorrect ";
			G_complete = 0;
			break;
		      }
		    else
		      G_complete(i) = myset<DataType>(r2,r3); //r2 + I *r3;
		    
                  }
              }
          }
        data(n1,n2,ALL) = G_complete(Range(mesh.index_min,mesh.index_max));
      }
}

//-------------------------------------------------------------------------

template <typename DataType>
void GF_Bloc_Base<DataType>::MPI_reduce_sum_onsite() { 
  myMPI_reduce_sum_on_site((Array<DataType,3>&)data);
  // WHAT ABOUT THE TAIL ?
}

//-------------------------------------------------------------------------

template <typename DataType>
void GF_Bloc_Base<DataType>::MPI_bcast() { 
  myMPI_bcast((Array<DataType,3>&)data);
  // BCAST THE TAIL ?
}

//-------------------------------------

template class GF_Bloc_Base<COMPLEX>;
template class GF_Bloc_Base<double>;

 
sum (expr,  x_ << domain )
sum (expr,  x_ | domain )
sum (expr,  x_ |= domain )
sum (expr,  x_ & y_ |= domain )

auto i_ = first_index(A);
auto k_ = first_index(B);

// HELP ! : i_== k_ !!!!!

auto i_= placeholder<0> & A.index_domain(1);

auto i_= A.index<0>(1); // obscur ?

auto j_ = second_index(A);

sum (A(i_,j_) * B(j_,1), j_ = A.index_range(1) );
sum (A(i_,j_) * B(j_,1), j_ |= range(1,n) );
sum (A(i_,j_) * B(j_,1), j_ <<= range(1,n) );
sum (A(i_,j_) * B(j_,1), j_ << range(1,n) );

sum (A(i_,j_) * B(j_,1), j_);
sum (A(i_,j_) * B(j_,i_), i_ & j_);


// GOOD 
sum (A(i_,j_) * B(j_,1), j_ = range(1,n) );
sum (A(i_,j_) * B(j_,1), i_ = A.index_range(0), j_ = A.index_range(1) );

// is the same : ok, 3 --> range(3,4)
sum (A(i_,j_) * B(j_,1), j_ = 3 ) == eval (A(i_,j_) * B(j_,1), j_ = 3 );  

sum( expr, i_= D) ----> triqs::lazy::sum_impl<D_type>::invoke(i_>> expr, D);

integrate( g(k0,om_), om_ = segment[a,b]); 

//options of the sum ?
sum(expr, i_=D, Riemann( ));

auto om_ = omega_index(G);
sum(G(k0,om_), om_);






