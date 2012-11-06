
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

#include "TailGF.hpp"

TailGF::TailGF(int OrderMinMIN_, int OrderMaxMAX_,
	       python::list IndicesL_,  python::list IndicesR_ ) :
  IndicesL(IndicesL_),IndicesR(IndicesR_),
  N1(int(PySequence_Size(IndicesL_.ptr()))),// int for 64bits
  N2(int(PySequence_Size(IndicesR_.ptr()))),// int for 64bits
  OrderMinMIN(OrderMinMIN_),OrderMaxMAX(OrderMaxMAX_),
  M(OrderMaxMAX - OrderMinMIN +1,N1,N2,COrder),
  OrderMaxArray(N1,N2,COrder)
 {
   M.reindexSelf(TinyVector<int,3>(OrderMinMIN,1,1));
   OrderMaxArray = OrderMaxMAX;
 } 

TailGF::TailGF(int OrderMinMIN_, int OrderMaxMAX_, python::object Arr, int OrderMax_, 
	       python::list IndicesL_,  python::list IndicesR_ ) :
  IndicesL(IndicesL_),IndicesR(IndicesR_),
  N1(int(PySequence_Size(IndicesL_.ptr()))),// int for 64bits
  N2(int(PySequence_Size(IndicesR_.ptr()))),// int for 64bits
  OrderMinMIN(OrderMinMIN_),OrderMaxMAX(OrderMaxMAX_),
  M(Arr.ptr()),
  OrderMaxArray(N1,N2,COrder)
 {
   if (M.shape() != TinyVector<int,3>(OrderMaxMAX - OrderMinMIN +1,N1,N2)) TRIQS_RUNTIME_ERROR<<"Reconstruction : size mismatch";
   M.reindexSelf(TinyVector<int,3>(OrderMinMIN,1,1));
   OrderMaxArray = OrderMax_;
 } 

inline python::object aux1(python::object o1) {
  python::extract<int> R(o1);
  return ( R.check() ? python::slice(R(),R()+1,1) : o1);
}

TailGF::TailGF(const TailGF & t, python::object slL, python::object slR):
  IndicesL(t.IndicesL[slL]),IndicesR(t.IndicesR[slR]),
  N1(int(PySequence_Size(IndicesL.ptr()))),// int for 64bits
  N2(int(PySequence_Size(IndicesR.ptr()))),// int for 64bits
  OrderMinMIN(t.OrderMinMIN),OrderMaxMAX(t.OrderMaxMAX),
  M(t.M.as_BoostObject() [python::make_tuple( python::slice(), aux1(slL),aux1(slR))]),
  OrderMaxArray(t.OrderMaxArray.as_BoostObject() [python::make_tuple(aux1(slL),aux1(slR))]) //VIEW
{
   // will never happen
   if (M.extent(1)!=N1) TRIQS_RUNTIME_ERROR<<"Tail GF constructor :"<<M.extent(1) <<"  "<<N1; 
   if (M.extent(2)!=N2) TRIQS_RUNTIME_ERROR<<"Tail GF constructor :"<<M.extent(2) <<"  "<<N2; 
   M.reindexSelf(TinyVector<int,3>(OrderMinMIN,1,1));
 } 


TailGF::TailGF(const TailGF & t) : 
  IndicesL(t.IndicesL),IndicesR(t.IndicesR),
  N1(t.N1),N2(t.N2),
  OrderMinMIN(t.OrderMinMIN),OrderMaxMAX(t.OrderMaxMAX),
  M(t.M.as_PyObjectPtrBorrowedRef()), //VIEW
  OrderMaxArray(t.OrderMaxArray.as_PyObjectPtrBorrowedRef()) //VIEW
{
   M.reindexSelf(TinyVector<int,3>(OrderMinMIN,1,1));
}

Array<COMPLEX,2> TailGF::eval(COMPLEX omega) const {
  const int omin(OrderMin()), omax(OrderMax());
  Array<COMPLEX,2> Res(N1,N2);
  COMPLEX z= pow(omega , -omin);
  Res = 0;
  for (int u =omin; u<= omax;u++,z /=omega ) 
    Res += z* (*this)[u];
  return Res;
}

python::object TailGF::__getitem__(int i) const {
  if (!has_coef(i)) TRIQS_RUNTIME_ERROR<<"TailGF:: operator[] :: order is incorrect";
  python::object a2 = M.as_BoostObject() [python::make_tuple( i-OrderMinMIN, python::slice(),python::slice())];
  python::object cls = python::import("pytriqs.Base.GF_Local.ArrayViewWithIndexConverter").attr("ArrayViewWithIndexConverter");
  return cls(a2, IndicesL, IndicesR);
}

void TailGF::__setitem__(int i,const PyArray<COMPLEX,2> & val) {
  if (!has_coef(i))  TRIQS_RUNTIME_ERROR<<"TailGF:: operator[] :: order is incorrect";
  if (val.shape() !=TinyVector<int,2>(N1, N2)) TRIQS_RUNTIME_ERROR<<"Tail : __setitem__ : size mismatch";
  M(i,ALL,ALL) = val;
}

python::object TailGF::AllCoefs() const {
  python::dict R;
  const int omin(OrderMin()), omax(OrderMax());
  for (int r=omin; r<=omax; r++) R[r]= this->__getitem__(r);
  return R;
}
python::object TailGF::__repr__() const {
  //  if (Nnd>1) return object(string("TailGF class"));
  stringstream fs;
  const int omin(OrderMin()), omax(OrderMax());
  for (int r=omin; r<=omax; r++){
    if ((N1==1) && (N2==1) ) fs<< M(r,1,1); else fs<< M(r,ALL,ALL);
    if (r<0) fs<<" Om^"<<abs(r)<<"  ";
    if (r>0) fs<<" /Om^"<<r<<"  ";
    if (r!=omax) fs<<" + ";
  }
  return python::object(string(fs.str()));
}

python::object TailGF::__call__(COMPLEX omega) const {
  PyArray<COMPLEX,2> Res(N1,N2,COrder);
  Res = eval(omega);
  python::object cls = python::import("pytriqs.Base.GF_Local.ArrayViewWithIndexConverter").attr("ArrayViewWithIndexConverter");
  return cls(Res, IndicesL,IndicesR);
}

void TailGF::from_L_T_R(const PyArray<COMPLEX,2> & L, const TailGF & T, const PyArray<COMPLEX,2> & R) {
  const int omin(T.OrderMin()), omax(T.OrderMax());
  //  Lnp * G2 pp * R pn
  if (L.shape() != TinyVector<int,2>(N1,T.N1)) TRIQS_RUNTIME_ERROR<<"The left hand site matrix has incorrect dimentions";
  if (R.shape() != TinyVector<int,2>(T.N2,N2)) TRIQS_RUNTIME_ERROR<<"The right hand site matrix has incorrect dimentions";

  Array<COMPLEX,2> tmp(T.N1,N2, fortranArray); // tmp = TR
  for (int r=omin; r<=omax; r++){
    matmul_lapack ( T.M(r,ALL,ALL), R, tmp ) ;  //T * R -> tmp
    Array<COMPLEX,2> RES(M(r,ALL,ALL)); // view
    matmul_lapack ( L,tmp, RES ) ; // L * tmp -> this
  }
  cleanM_outside(omin,omax);
}

TailGF TailGF::copy() {
  TailGF t(OrderMinMIN, OrderMaxMAX, IndicesL, IndicesR);
  t.M = M;
  t.OrderMaxArray = OrderMaxArray;
  return t;
}
 
void TailGF::invert() {
  const int omin(OrderMin()), omax(OrderMax());
  if (N1!=N2) TRIQS_RUNTIME_ERROR<<"Inversion can only be done for square matrix !";
  int new_OrderMin = - omin;
  if (new_OrderMin <OrderMinMIN) TRIQS_RUNTIME_ERROR<<" I can not inverse with the OrderMinMIN and OrderMaxMAX provided";
  int new_OrderMax = min(OrderMaxMAX,omax - omin + new_OrderMin);

  Array<COMPLEX,3> newM(Range(OrderMinMIN, OrderMaxMAX),Range(0,N1-1),Range(0,N2-1)); //new_OrderMax - new_OrderMin +1,N1,N2);
  newM=0;

  Array<COMPLEX,2> tmp(N1,N2,fortranArray),S(N1,N2,fortranArray);
  Array<COMPLEX,2> M0( newM(new_OrderMin,ALL,ALL));
  M0 = M(omin,ALL,ALL); 
  if (!Inverse_Matrix(M0)) TRIQS_RUNTIME_ERROR<<"TailGF::invert : This tail is 0  !!!";

  // b_n = - a_0^{-1} * sum_{p=0}^{n-1} b_p a_{n-p} for n>0
  // b_0 = a_0^{-1}
  for (int n = 1; n<=new_OrderMax-new_OrderMin; ++n) {
    S = 0;
    for (int p=0; p<n;p++) {
      matmul_lapack(newM(new_OrderMin + p,ALL,ALL),M(omin+ n-p ,ALL,ALL),tmp);
      S +=tmp;
    }
    matmul_lapack(M0,S,tmp);
    newM(n+new_OrderMin,ALL,ALL) = -tmp;
  }
  M= newM;
  OrderMaxArray = new_OrderMax;
}


void TailGF::save(string file, bool accumulate) const {
  const int omin(OrderMin()), omax(OrderMax());
  stringstream fs; 
  fs<<file<<(( ( N1*N2>1) ) ? "/" : ".") <<"Moments.dat";
  ofstream f(fs.str().c_str(), (accumulate ? ios::out|ios::app : ios::out));
  f.setf(ios::fixed,ios::floatfield);f.precision(PRECISION_OUTPUT);
  f<<OrderMin()<<endl<<OrderMax()<<endl;
  for (int r=omin; r<=omax; r++) f << (*this)[r];
  f.close();
}
 
void TailGF::load(string filename) {

  bool NewStyleSave = true;
  stringstream fs1; 
  fs1<<filename<<((NewStyleSave && ( N1*N2>1) ) ? "/" : ".") <<"Moments.dat";
  ifstream mom; mom.open(fs1.str().c_str());
  if (mom.fail()) { 
    REPORT <<"CAN NOT FIND THE MOMENTS : GF_Bloc_w::load : Putting 0 !"<<endl;M=0;
  }
  else {
    int omin,omax;
    mom>>omin>>omax;
    if (omin<OrderMinMIN) TRIQS_RUNTIME_ERROR<<"TailGF::load : OrderMin too small";
    if (omax>OrderMaxMAX) TRIQS_RUNTIME_ERROR<<"TailGF::load : OrderMax too large";
    zero();
    Array <COMPLEX,2> MM(N1,N2,fortranArray);
    for (int i=omin; i<=omax;i++) {mom>>MM; (*this)[i] = MM;}
    OrderMaxArray = omax;
  }

}

python::object TailGF::__reduce_to_dict__() const {
  // construct a new dictionnary and fill it
  python::dict d; 
  d["OrderMinMIN"] = OrderMinMIN;
  d["OrderMax"] = OrderMax();
  d["OrderMaxMAX"] = OrderMaxMAX;
  // correction ticket #94
  d["IndicesL"] = python::tuple(IndicesL).attr("__repr__")();
  d["IndicesR"] = python::tuple(IndicesR).attr("__repr__")();
  //d["IndicesL"] = python::import("numpy").attr("array")(IndicesL);
  //d["IndicesR"] = python::import("numpy").attr("array")(IndicesR);
  d["array"] = python::object(M);
  return d;
}

#include <boost/python/exec.hpp>
// correction ticket #94
inline python::list obj_or_string_to_list(python::object x) { 
 python::extract<std::string> a(x); 
 bool is_string = a.check();
 python::object r = (is_string ? python::eval(python::str(x)) : x); // if x is a string, eval it
 return python::list(r);
}

python::object TailGF::__factory_from_dict__(const python::object & dic){
  python::dict d(dic);
  int omin= python::extract<int>(d["OrderMinMIN"]);
  int omax= python::extract<int>(d["OrderMaxMAX"]);
  int omax2= python::extract<int>(d["OrderMax"]);
  python::list indL, indR;
  if (d.has_key("Indices")) {
    indL = obj_or_string_to_list( d["Indices"]);
    //indL = python::list(d["Indices"]);
    indR = indL;
  } else {
    indL = obj_or_string_to_list(d["IndicesL"]);
    indR = obj_or_string_to_list(d["IndicesR"]);
    //indL = python::list(d["IndicesL"]);
    //indR = python::list(d["IndicesR"]);
  }
  python::object arr(d["array"]);
  PyArray<COMPLEX,3> ARR(arr);
  if (ARR.extent(0) !=( omax - omin +1)) TRIQS_RUNTIME_ERROR<<"Reconstruction : Inconsistent data for TailGF";
  return python::object(TailGF(omin,omax,arr,omax2,indL,indR));
}


python::tuple TailGF::__reduce__() const { 
  return python::make_tuple( python::import("pytriqs.Base.Utility.myUtils").attr("call_factory_from_dict"),
			     python::make_tuple(python::object(*this).attr("__class__"),__reduce_to_dict__()));

}

bool TailGF::operator==(const TailGF & other) const {
  const int omin(OrderMin()), omax(OrderMax());
  const int tomin(other.OrderMin()), tomax(other.OrderMax());
  return (omin == tomin) && (omax == tomax) &&
    (N1 ==other.N1) &&(N2 ==other.N2) &&  (max(abs(M(Range(omin,omax),ALL,ALL) - other.M(Range(omin,omax),ALL,ALL)))) < ZERO; 
}  


TailGF & TailGF::operator *= (const TailGF & t) {

    const int omin(OrderMin()), omax(OrderMax());
    const int new_ordermin(omin+t.OrderMin());
    const int new_ordermax(min(omax+t.OrderMax(),OrderMaxMAX));

    if (N1!=t.N1 || N2!=t.N2) TRIQS_RUNTIME_ERROR<<"Multiplication is valid only for similar tail shapes !";
    if (new_ordermin < OrderMinMIN) TRIQS_RUNTIME_ERROR<<"The multiplication makes the new tail have a too small OrderMin";

    Array<COMPLEX,3> newM(Range(OrderMinMIN, OrderMaxMAX), Range(1,N1), Range(1,N2));
    Array<COMPLEX,2> tmp(N1, N2, fortranArray);

    for (int n = new_ordermin; n <= new_ordermax; ++n) {
      const int kmin(max(0,n-t.OrderMax()-omin));
      const int kmax(min(omax-omin,n-t.OrderMin()-omin));
    //  cout << "n, kmin, kmax" << n << " " << kmin << " " << kmax << endl;
      for (int k = kmin; k <= kmax; ++k) {
        matmul_lapack(M(omin+k,ALL,ALL), t.M(n-omin-k,ALL,ALL), tmp);
        newM(n,ALL,ALL) += tmp;
      }
    }

    OrderMaxArray = new_ordermax;
    M = newM;
    
    return *this;
}  

TailGF TailGF::transpose() const{ 
 // return a new view with transpose array, as numpy
 boost::python::object Mt= M.as_BoostObject().attr("transpose")(boost::python::make_tuple(0,2,1));
 return  TailGF(OrderMinMIN, OrderMaxMAX, Mt , this->OrderMax(),IndicesL,IndicesR );
}

TailGF TailGF::conjugate(bool is_matsubara_expansion) const {
 // return a new view with conjugate array, as numpy
 boost::python::object Mt= M.as_BoostObject().attr("conjugate")();
 TailGF r = TailGF(OrderMinMIN, OrderMaxMAX, Mt , this->OrderMax(),IndicesL,IndicesR );
 // In the Matsubara frequency case
 if (is_matsubara_expansion) 
  for (int n= OrderMinMIN; n<=OrderMaxMAX; ++n) {
    r.M(n,ALL,ALL) *= ( ((n%2 ==1) || (n%2 == -1)) ? -1 : 1);
  }
 return r;
}

void TailGF::check(const TailGF & t) {
 if (t.N1 !=N1) TRIQS_RUNTIME_ERROR<<"TailGF::operator =  N1 mismatch : t.N1 ="<<t.N1 <<" this.N1 =  "<<N1 ;
 if (t.N2 !=N2) TRIQS_RUNTIME_ERROR<<"TailGF::operator =  N2 mismatch : t.N2 ="<<t.N2 <<" this.N2 =  "<<N2;
 const int tomin(t.OrderMin()), tomax(t.OrderMax());
 if ((tomin < OrderMinMIN) ||  (tomin > OrderMaxMAX) || 
   (tomax < OrderMinMIN) ||  (tomax > OrderMaxMAX) )
  TRIQS_RUNTIME_ERROR
   <<"TailGF::operator = : OrderMin/Max mismatch\n"
   << " When the tailGF is a view, the size of the RHS must match the size of this\n"
   << " t Min/max : "<< tomin << "  "<< tomax << " this :  "<< OrderMinMIN<< "  "<<OrderMaxMAX<<"\n";
}

void TailGF::cleanM_outside(int wmin, int wmax) {
 if (M.lbound(0) < wmin) M(Range(M.lbound(0),wmin-1),ALL,ALL) = 0;
 if (M.ubound(0) > wmax) M(Range(wmax+1,M.ubound(0)),ALL,ALL) = 0;
}

TailGF & TailGF::operator = (const TailGF & t) {   
 check(t);
 const int tomin(t.OrderMin()), tomax(t.OrderMax());
 //  M = 0; // wrong if M a (partial) view of t.M !
 M(Range(tomin,tomax),ALL,ALL) = t.M(Range(tomin,tomax),ALL,ALL);
 cleanM_outside(tomin,tomax);
 OrderMaxArray=tomax;
 return *this;
}

TailGF & TailGF::operator += (const TailGF & t) {
 check(t);
 OrderMaxArray=min(OrderMax(),t.OrderMax());
 M(ALL,ALL,ALL) += t.M(ALL,ALL,ALL);
 return *this;
}

TailGF & TailGF::operator -= (const TailGF & t) {
 check(t);
 OrderMaxArray=min(OrderMax(),t.OrderMax());
 M(ALL,ALL,ALL) -= t.M(ALL,ALL,ALL);
 return *this;
}


