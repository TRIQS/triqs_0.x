
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

#ifndef my_MPI_H
#define my_MPI_H

#include <triqs/utility/report_stream.hpp>
#include <cstring>


inline int system(std::string s){return(system(s.c_str()));}

/**
 * \defgroup mympi Interface for MPI calls.
 *          - Can be used on a MPI or non-MPI machine as well
 *          - Extended for Blitz Array (if -DUSE_BLITZ)
 *          - Simpler calls than the C calls.
 @{
*/

#define FATAL(s) { std::stringstream fs_for_fatal; fs_for_fatal<<s; throw fs_for_fatal.str();}
#define FIN flush<<endl
#define ERROR(s) { std::clog<<"ERROR : "<<s<<endl; return;}
#define FLAG(x)     cout<<x<<endl
//#define myassert(x,message)  if (!(x)) FATAL(message);
inline void myassert(bool x,std::string message)  {if (!(x)) FATAL(message);}

const int PRECISION_OUTPUT=15;

void REPORT_FINALIZE();
std::ostream & REPORT_PTR();

/** 
    REPORT << blabla 
 */
#define REPORT REPORT_PTR()

/// Id of the master
const int myMPI_master=0;

/// Test : are we on master node ?
#define is_master_node  (myMPI_myid()==myMPI_master)
#define REPORT_WARNING(s){if is_master_node std::clog<<"WARNING : "<<s<<endl;}

/// Beginning of section : master only
#define BEGIN_MASTER_ONLY  if is_master_node {

/// End of "master only" section
#define END_MASTER_ONLY  }

#ifndef HAVE_MPI

/// \internal
struct MPI_Request {};
const MPI_Request MPI_REQUEST_NULL=MPI_Request();


					
///   Bcast 
template <class T>   void  myMPI_bcast(T & x, int s=1) {}			

/// Reduce with sum
template <class T>   void  myMPI_reduce_sum(T & a, T & b, int s=1) {b=a;}

/// Reduce with sum, result in a
template <class T>   void  myMPI_reduce_sum_on_site(T & a, int s=1) {}

/// A in from --> B in dest
template <class T>   void  myMPI_exchange(T & A,int from,T & B,int dest){B=A;}

/// Send A to TO
template <class T>   void  myMPI_send(const T & A,int to){}

/// Receive A from FROM
template <class T>   void  myMPI_recv(T & A,int from){}

/// Non blocking receive
template <class T>   void  myMPI_irecv(T & A,int from,MPI_Request & request){}

/// Init MPI
inline void myMPI_init(int argc, char *argv[] ){}

/// Number of processors
inline int myMPI_number_proc(){ return (1);}

/// What is the number of the processor ?
inline int myMPI_myid(){ return (myMPI_master);}

/// End of MPI
inline void myMPI_finalize(){REPORT_FINALIZE();}

/// Simply BARRIER
inline void myMPI_barrier() {}

/// To slice an interval over the nodes : min bound
inline int myMPI_slice_inf(int imin,int imax) {return(imin);}

/// To slice an interval over the nodes : max bound
inline int myMPI_slice_sup(int imin,int imax){return (imax);}

///
template<class T> void myMPI_scatterv(const T & A, T & B) {B=A;}

#else 


/************************************************************
   C++ bindings 
*********************************************************/

#ifdef MPI_BINDINGS_CXX

#include "mpicxx.h"

/**
   Type of MPI data.
 */
template <class T> struct mpi_datatype { static const MPI::Datatype value;};

// we have to redefine the FATAL to end properly in MPI.
//#undef FATAL
//#define FATAL(s) { if (is_master_node) std::clog<<"FATAL : "<<s<<endl; myMPI_barrier(); myMPI_finalize(); exit(1);}

inline void myMPI_init(int argc, char *argv[] ){MPI::Init(argc,argv);}
inline int myMPI_myid(){return MPI::COMM_WORLD.Get_rank();} 
inline int myMPI_number_proc(){ return MPI::COMM_WORLD.Get_size();}
inline void myMPI_finalize(){ MPI::COMM_WORLD.Barrier(); REPORT_FINALIZE();MPI::Finalize();}
inline void myMPI_barrier(){ MPI::COMM_WORLD.Barrier();}


inline int myMPI_slice_inf(int imin,int imax){
  int j=(imax - imin + 1)/myMPI_number_proc();
  int i= imax - imin + 1 -myMPI_number_proc()*j;
  return ( myMPI_myid()<=i-1 ?  imin + myMPI_myid()*(j+1) :  imin + myMPI_myid()*j + i);
}


inline int myMPI_slice_sup(int imin,int imax){
  int j=(imax - imin + 1)/myMPI_number_proc();
  int i= imax - imin + 1 -myMPI_number_proc()*j;
  return ( myMPI_myid()<=i-1 ?  imin + (myMPI_myid()+1)*(j+1) -1 : imin + (myMPI_myid()+1)*j  + i - 1 );
}


template<class T>
inline void myMPI_bcast(T & x,int s=1)			
{
  myMPI_barrier();
  MPI::COMM_WORLD.Bcast(& x, s, mpi_datatype<T>::value, myMPI_master);
  myMPI_barrier();
}

  
template<class T>
inline void  myMPI_reduce_sum(T & a, T & b,int s=1)				
{
  myMPI_barrier();
  MPI::COMM_WORLD.Reduce (&a,&b,s, mpi_datatype<T>::value, MPI::SUM, myMPI_master);
  myMPI_barrier();
}
 
template<class T>
inline void  myMPI_reduce_sum_on_site(T & a, int s=1)				
 {
   T  b; 
   myMPI_barrier();  
   MPI::COMM_WORLD.Reduce (&a,&b,s, mpi_datatype<T>::value, MPI::SUM, myMPI_master);
   myMPI_barrier();
   a=b;
 }


inline void myMPI_bcast(string & s)
{ int i = s.length()+1; myMPI_bcast(i); char *cc= new char[i];strcpy(cc,s.c_str());
  myMPI_barrier();
  MPI::COMM_WORLD.Bcast(cc, i, MPI::CHAR, myMPI_master);
  myMPI_barrier();
  s.assign(cc);delete cc;}

template<class T>
void myMPI_send(const T & A, int dest)
{ 
  MPI::COMM_WORLD.Send(const_cast<T*>(&A),1, mpi_datatype<T>::value, dest, 99);
}

template<class T>
void myMPI_recv(T & A,int from)
{
  MPI::Status status; 
  MPI::COMM_WORLD.Recv(&A,1, mpi_datatype<T>::value, from, 99,status);
}


// Special case for the strings.

void myMPI_send(const string & A, int dest);

void myMPI_recv(string & A,int from);

/*
template<class T>
void myMPI_irecv(T & A,int from,MPI_Request & request)
{
  MPI_Irecv(&A,1, mpi_datatype<T>::value, from, 99,MPI_COMM_WORLD,&request);
}
*/

#else

/************
    SIMPLE C bindings.
****************/


#include "mpi.h"

/**
   Type of MPI data.
 */
template <class T> struct mpi_datatype { static const MPI_Datatype value;};

// we have to redefine the FATAL to end properly in MPI.
//#undef FATAL
//#define FATAL(s) { if (is_master_node) std::clog<<"FATAL : "<<s<<endl; myMPI_barrier(); myMPI_finalize(); exit(1);}

inline void myMPI_init(int argc, char *argv[] ){MPI_Init(&argc,&argv);}
inline int myMPI_myid(){int num; MPI_Comm_rank(MPI_COMM_WORLD,&num);return(num);} 
inline int myMPI_number_proc(){ int num; MPI_Comm_size(MPI_COMM_WORLD,&num); return(num);}
inline void myMPI_finalize(){ MPI_Barrier(MPI_COMM_WORLD); REPORT_FINALIZE();MPI_Finalize();}
inline void myMPI_barrier(){MPI_Barrier(MPI_COMM_WORLD);}


inline int myMPI_slice_inf(int imin,int imax){
  int j=(imax - imin + 1)/myMPI_number_proc();
  int i= imax - imin + 1 -myMPI_number_proc()*j;
  return ( myMPI_myid()<=i-1 ?  imin + myMPI_myid()*(j+1) :  imin + myMPI_myid()*j + i);
}


inline int myMPI_slice_sup(int imin,int imax){
  int j=(imax - imin + 1)/myMPI_number_proc();
  int i= imax - imin + 1 -myMPI_number_proc()*j;
  return ( myMPI_myid()<=i-1 ?  imin + (myMPI_myid()+1)*(j+1) -1 : imin + (myMPI_myid()+1)*j  + i - 1 );
}


template<class T>
inline void myMPI_bcast(T & x,int s=1)			
{
  MPI_Barrier( MPI_COMM_WORLD);
  MPI_Bcast(& x, s, mpi_datatype<T>::value, myMPI_master, MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD);
}

  
template<class T>
inline void  myMPI_reduce_sum(T & a, T & b,int s=1)				
{
  MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Reduce (&a,&b,s, mpi_datatype<T>::value, MPI_SUM, myMPI_master, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
}
 
template<class T>
inline void  myMPI_reduce_sum_on_site(T & a, int s=1)				
 {
   T  b; 
   MPI_Barrier(MPI_COMM_WORLD);  
   MPI_Reduce (&a,&b,s, mpi_datatype<T>::value, MPI_SUM, myMPI_master, MPI_COMM_WORLD); 
   MPI_Barrier(MPI_COMM_WORLD);
   a=b;
 }


inline void myMPI_bcast(std::string & s)
{ int i = s.length()+1; myMPI_bcast(i); char *cc= new char[i];strcpy(cc,s.c_str());
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(cc, i, MPI_CHAR, myMPI_master, MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD);
  s.assign(cc);delete cc;}

template<class T>
void myMPI_send(const T & A, int dest)
{ 
  MPI_Send(const_cast<T*>(&A),1, mpi_datatype<T>::value, dest, 99,MPI_COMM_WORLD);
}

template<class T>
void myMPI_recv(T & A,int from)
{
  MPI_Status status; 
  MPI_Recv(&A,1, mpi_datatype<T>::value, from, 99,MPI_COMM_WORLD,&status);
}


// Special case for the strings.

void myMPI_send(const std::string & A, int dest);

void myMPI_recv(std::string & A,int from);

template<class T>
void myMPI_irecv(T & A,int from,MPI_Request & request)
{
  MPI_Irecv(&A,1, mpi_datatype<T>::value, from, 99,MPI_COMM_WORLD,&request);
}

#endif

#endif

/* @} */

#endif
