
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

#include "my_mpi.hpp"
#include <fstream>
#include <complex>

using namespace std;

std::ostream *REPORT_ptr_ =NULL;

void REPORT_FINALIZE() {}
//{if (IS_MASTER_NODE) {REPORT<<"Ending at : "<<endl; cout<<"okok"<<endl; system("date"); cout<<"finifini"<<endl; } }

 
std::ostream & REPORT_PTR() {
  if (!REPORT_ptr_) 
    {
      cout.setf(ios::fixed,ios::floatfield);cout.precision(PRECISION_OUTPUT);
      REPORT_ptr_ =  (IS_MASTER_NODE ?  &(cout) :  new ofstream("/dev/null",ios::out));
      //if (IS_MASTER_NODE) { 
      //(*REPORT_ptr_)<<"Starting on "<<myMPI_number_proc()<< "Nodes at : "<<endl;
      //system("date");
      //}
    }
  return(*REPORT_ptr_); }



#ifdef HAVE_MPI


/************************************************************
   C++ bindings 
*********************************************************/

#ifdef MPI_BINDINGS_CXX

template <> const MPI::Datatype  mpi_datatype<bool>::value = MPI::LOGICAL;
template <> const MPI::Datatype  mpi_datatype<int>::value = MPI::INT;
template <> const MPI::Datatype  mpi_datatype<double>::value = MPI::DOUBLE;
template <> const MPI::Datatype  mpi_datatype<float>::value = MPI::FLOAT;
template <> const MPI::Datatype  mpi_datatype<complex<double> >::value = MPI::DOUBLE_COMPLEX;

void myMPI_send(const string & A, int dest)
{ 
  int size = A.length()+1; // take the 0 at the end of the string....
  myMPI_send(size,dest);
  MPI::COMM_WORLD.Send((char *) A.c_str(),size, MPI::CHAR, dest, 99);
} 

void myMPI_recv(string & A,int from)
{
  MPI::Status status; int size;
  myMPI_recv(size,from);
  char *Value= new char[size];  
  MPI::COMM_WORLD.Recv(Value,size, MPI::CHAR, from, 99,status);
  A.assign(Value);
  delete[] Value;
}

#else

template <> const MPI_Datatype  mpi_datatype<bool>::value = MPI_LOGICAL;
template <> const MPI_Datatype  mpi_datatype<int>::value = MPI_INT;
template <> const MPI_Datatype  mpi_datatype<double>::value = MPI_DOUBLE;
template <> const MPI_Datatype  mpi_datatype<float>::value = MPI_FLOAT;
template <> const MPI_Datatype  mpi_datatype<complex<double> >::value = MPI_DOUBLE_COMPLEX;

void myMPI_send(const string & A, int dest)
{ 
  //  cout<<"Node "<<myMPI_myid()<<" sending "<< A<<"  "<<A.length()<<endl;
  int size = A.length()+1; // take the 0 at the end of the string....
  myMPI_send(size,dest);
  MPI_Send((char *) A.c_str(),size, MPI_CHAR, dest, 99,MPI_COMM_WORLD);
} 

void myMPI_recv(string & A,int from)
{
  MPI_Status status; int size;
  myMPI_recv(size,from);
  char *Value= new char[size];  
  MPI_Recv(Value,size, MPI_CHAR, from, 99,MPI_COMM_WORLD,&status);
  A.assign(Value);
  //  cout<<"Node "<<myMPI_myid()<<" received "<< A<<"  "<<A.length()<<endl;
  delete[] Value;
}

#endif

#endif
