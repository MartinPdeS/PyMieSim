#pragma once
#include<malloc.h>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

namespace py = pybind11 ;
using namespace pybind11::literals ;
typedef std::complex<float> complex_f ;
typedef std::complex<double> complex_d ;
typedef unsigned int uint ;


template <size_t alignment= 64>
void* _64B_aligned_malloc( size_t size );









template<class Type,uint8_t Dim,class IndexType=uint>
class Multi_array
{

};


template<class Type, class IndexType>
class Multi_array<Type,1,IndexType>
{
	public :
    /* Default constructor */
    Multi_array();

	/* constructor default stride */
	Multi_array
	(
		IndexType n_i , // Number of elements in i
		void* (*alloc_func)(size_t size) = &_64B_aligned_malloc, // Custom allocation function
		void (*free_func)(void* ptr) = &free
	);

	/* usual constructor */
	Multi_array
	(
		IndexType n_i , // Number of elements in i
		size_t stride_i , // The number of Bytes of one element
		void* (*alloc_func)(size_t size) = &_64B_aligned_malloc, // Custom allocation function
		void (*free_func)(void* ptr) = &free
	);

	/* Constructing from an existing pointer */
	Multi_array ( Type* prt, IndexType n_i , size_t stride_i = sizeof(Type) );

	/* Constructing from a 1D Numpy array */
	static Multi_array numpy_copy	( py::array_t<Type,py::array::c_style>& np_array );
	static Multi_array numpy_share	( py::array_t<Type,py::array::c_style>& np_array );

	/* Copy constructor */
	Multi_array( Multi_array& Mom );
	Multi_array( const Multi_array& Mom );

	/* Move constructor */
	// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2027.html#Move_Semantics
	// https://en.cppreference.com/w/cpp/language/move_constructor
	Multi_array( Multi_array&& Mom );

	/* Destructor */
	~Multi_array();

    /* Assignment operators (not implemented)*/
	/*
		Copy assignment operator
			https://en.cppreference.com/w/cpp/language/copy_assignment
	*/

	/*
		Move assignment operator
			https://en.cppreference.com/w/cpp/language/move_assignment
	*/

    /* Indexing operators */
	Type& operator()( IndexType i ); /* Returns a reference to an element */
	Type& operator[]( IndexType i ); /* Returns a reference to an element */

	/*
		Same behavior for const Multi_array
		See : https://en.cppreference.com/w/cpp/language/member_functions#const-_and_volatile-qualified_member_functions
	*/
	Type& operator()( IndexType i ) const ; /* Returns a reference to an element */
	Type& operator[]( IndexType i ) const ; /* Returns a reference to an element */


	Type* get_ptr(){ return ptr; } ;
	Type* get_ptr()const{ return ptr; }  ;

	/* Copy to numpy methods */

	/*
		move_py :
		Same meaning as move constructor.
		"steal" the resources held by the current object and gives them to a numpy array
		and leave the current object in valid but empty state.
		Python inherits all memory responsibilities
	*/
	py::array_t<Type, py::array::c_style> move_py();
	py::array_t<Type, py::array::c_style> move_py(IndexType n_i); // subset move : numpy array has shape {n_i}

	// 	copy_py : Same meaning as copy constructor. Makes a copy of the current object in a  numpy array format and gives it to python (with all responsibilities)
	py::array_t<Type, py::array::c_style> copy_py();
	py::array_t<Type, py::array::c_style> copy_py(IndexType n_i); // subset copy : partial copy. numpy array has shape {n_i}

	//	share_py : Doesn't make a copy and give access to current memmory to python (with no responsibilities)
	py::array_t<Type, py::array::c_style> share_py();
	py::array_t<Type, py::array::c_style> share_py(IndexType n_i); // subset share

	// 	affect_py : Doesn't make a copy and give access to current memmory to python (with all responsibilities)
	py::array_t<Type, py::array::c_style> affect_py();
	py::array_t<Type, py::array::c_style> affect_py(IndexType n_i); // subset affect


	IndexType get_n_i(){return n_i;};
	IndexType get_n_i()const{return n_i;}  ;

	size_t get_stride_i(){return stride_i;};
	size_t get_stride_i()const{return stride_i;};

	uint64_t get_alloc_memory_size(){return n_i*stride_i;}; /*This does not take strides into account...*/

	// alloc_ptr get_alloc_func(){return alloc_func;};
	// free_ptr get_free_func(){return free_func};

	private :
	void* (*alloc_func)(size_t size) ;
	void (*free_func)(void* ptr) ;

	Type* ptr ;
	IndexType n_i ;
	size_t stride_i ;

    char* displace( IndexType n_Bytes );
	char* displace( IndexType n_Bytes ) const;

    // void check_overflow();

	void _free_func();

};
