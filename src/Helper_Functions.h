#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H
#include "Var.h"

template<typename T>
inline void copy_omp(const std::vector<T> &in, std::vector<T> &out)
{
	#pragma omp parallel for
	for(size_t i = 0; i < in.size(); ++i)
		out[i] = in[i];
		
	return;
}

  /**
   *  @brief Copies the range [first,last) into result.
   *  @ingroup mutating_algorithms
   *  @param  __first  An input iterator.
   *  @param  __last   An input iterator.
   *  @param  __result An output iterator.
   *  @return   result + (first - last)
   *
   *  This inline function will boil down to a call to @c memmove whenever
   *  possible.  Failing that, if random access iterators are passed, then the
   *  loop count will be known (and therefore a candidate for compiler
   *  optimizations such as unrolling).  Result may not be contained within
   *  [first,last); the copy_backward function should be used instead.
   *
   *  Note that the end of the output range is permitted to be contained
   *  within [first,last).
  */
template<typename iterator>
inline void copy_omp(iterator const& __first, iterator const& __last, iterator const& __result)
{
	#pragma omp parallel
	{   
		int tid = omp_get_thread_num();
		auto chunksize = (__last - __first) / omp_get_num_threads();
		iterator begin = __first + chunksize * tid;
		iterator begin_new = __result + chunksize * tid;
		iterator end = (tid == omp_get_num_threads() -1) ? __last : begin + chunksize;
		std::copy(begin, end, begin_new);
	}
}

template<typename T>
inline void copy_n_omp(const std::vector<T> &in, const size_t n, std::vector<T> &out)
{
	#pragma omp parallel for
	for(size_t i = 0; i < n; ++i)
		out[i] = in[i];

	return;
}

  /**
   *  @brief Copies the range [first,first+n) into [result,result+n).
   *  @ingroup mutating_algorithms
   *  @param  __first  An input iterator.
   *  @param  __n      The number of elements to copy.
   *  @param  __result An output iterator.
   *  @return  result+n.
   *
   *  This inline function will boil down to a call to @c memmove whenever
   *  possible.  Failing that, if random access iterators are passed, then the
   *  loop count will be known (and therefore a candidate for compiler
   *  optimizations such as unrolling).
  */
template<typename iterator>
inline void copy_n_omp(iterator const& __first, size_t const& __n, iterator const& __result)
{
	#pragma omp parallel
	{   
		int tid = omp_get_thread_num();
		size_t chunksize = __n / omp_get_num_threads();
		iterator begin = __first + chunksize * tid;
		iterator begin_new = __result + chunksize * tid;
		size_t nt = (tid == omp_get_num_threads() -1) ? chunksize + __n % omp_get_num_threads() : chunksize;
		std::copy_n(begin, nt, begin_new);
	}
}


  /**
   *  @brief Fills the range [first,last) with copies of value.
   *  @ingroup mutating_algorithms
   *  @param  __first  A forward iterator.
   *  @param  __last   A forward iterator.
   *  @param  __value  A reference-to-const of arbitrary type.
   *  @return   Nothing.
   *
   *  This function fills a range with copies of the same value.  For char
   *  types filling contiguous areas of memory, this becomes an inline call
   *  to @c memset or @c wmemset.
  */
template<typename T, typename iterator>
inline void fill_omp(iterator const& __first, iterator const& __last, T const& __value)
{
	#pragma omp parallel
	{   
		auto tid = omp_get_thread_num();
		auto chunksize = (__last - __first) / omp_get_num_threads();
		auto begin = __first + chunksize * tid;
		auto end = (tid == omp_get_num_threads() -1) ? __last : begin + chunksize;
		std::fill(begin, end, __value);
	}
}

#endif