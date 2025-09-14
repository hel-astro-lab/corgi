#pragma once

#include <type_traits>
#include <array>
#include <utility>
#include <concepts>



/*! Internal meta-programming tools
 *
 * See especially
 *  - https://codereview.stackexchange.com/questions/107877/simple-multi-dimensional-array-class-in-c11
 *  - https://github.com/maddouri/hyper_array/
 *
 */

namespace corgi {

template <std::size_t D, typename... I>
concept indices_for = (sizeof...(I) == D) and (std::integral<I> and ...);

  namespace internals {

/// building block of a neat trick for checking multiple types against a given trait
template <bool...>
struct bool_pack
{};

/// neat trick for checking multiple types against a given trait
/// https://codereview.stackexchange.com/a/107903/86688
template <bool... bs>
using are_all_true = std::is_same<bool_pack<true, bs...>,
                                  bool_pack<bs..., true>>;

/// checks that all the template arguments are integral types
/// @note `T&` where `std::is_integral<T>::value==true` is considered integral
/// by removing any reference then using `std::is_integral`
template <typename... Ts>
using are_integral = are_all_true<
    std::is_integral<
        typename std::remove_reference<Ts>::type
    >::value...
>;

/// compile-time sum
template <typename T>
constexpr T ct_plus(const T x, const T y) { return x + y; }

/// compile-time product
template <typename T>
constexpr T ct_prod(const T x, const T y) { return x * y; }

/// compile-time equivalent to `std::accumulate()`
template <
    typename    T,  ///< result type
    std::size_t N,  ///< length of the array
    typename    O   ///< type of the binary operation
>
constexpr
T ct_accumulate(const ::std::array<T, N>& arr,  ///< accumulate from this array
                const size_t first,             ///< starting from this position
                const size_t length,            ///< accumulate this number of elements
                const T      initial_value,      ///< let this be the accumulator's initial value
                const O&     op                 ///< use this binary operation
               )
{
    // https://stackoverflow.com/a/33158265/865719
    //return (first < (first + length))
    //     ? op(arr[first],
    //          ct_accumulate(arr,
    //                        first + 1,
    //                        length - 1,
    //                        initial_value,
    //                        op))
    //     : initial_value;

    // TODO this one is allowed in >c++14
    T ret = initial_value;
    for(size_t i=0; i<length; i++){
        ret = op(ret, arr[first + i]);
    }
    return ret;
}


/// compile-time equivalent to `std::inner_product()`
template <
    typename T,      ///< the result type
    typename T_1,    ///< first array's type
    size_t   N_1,    ///< length of the first array
    typename T_2,    ///< second array's type
    size_t   N_2,    ///< length of the second array
    typename O_SUM,  ///< summation operation's type
    typename O_PROD  ///< multiplication operation's type
>
constexpr
T ct_inner_product(const ::std::array<T_1, N_1>& arr_1,  ///< calc the inner product of this array
                   const size_t  first_1,        ///< from this position
                   const ::std::array<T_2, N_2>& arr_2,  ///< with this array
                   const size_t  first_2,        ///< from this position
                   const size_t  length,         ///< using this many elements from both arrays
                   const T       initial_value,   ///< let this be the summation's initial value
                   const O_SUM&  op_sum,         ///< use this as the summation operator
                   const O_PROD& op_prod         ///< use this as the multiplication operator
                  )
{
    // same logic as `ct_accumulate()`
    //return (first_1 < (first_1 + length))
    //     ? op_sum(op_prod(arr_1[first_1],
    //                      arr_2[first_2]),
    //              ct_inner_product(arr_1, first_1 + 1,
    //                               arr_2, first_2 + 1,
    //                               length - 1,
    //                               initial_value,
    //                               op_sum, op_prod))
    //     : initial_value;
      
    // TODO this one is allowed in >c++14
    T ret = initial_value;
    for(size_t i=0; i<length; i++){
        ret = op_sum(ret, op_prod(arr_1[first_1 + i], arr_2[first_2 + i]));
    }
    return ret;
}



/*! computes the index coefficients assuming row-major order
 *
 *
 *  what we compute:
 *        \f[
 *            \begin{cases}
 *            C_i = \prod_{j=i+1}^{n-1} L_j
 *            \\
 *            \begin{cases}
 *                i   &\in [0, \text{Dimensions - 1}] \\
 *                C_i &: \text{\_coeffs[i]}           \\
 *                L_j &: \text{\_lengths[j]}
 *            \end{cases}
 *            \end{cases}
 *        \f]
 */
//template <typename size_type, std::size_t D>
//::std::array<size_type, D>
//compute_index_coeffs(const ::std::array<size_type, D>& dimension_lengths) noexcept
//{
//    ::std::array<size_type, D> coeffs;
//    for (size_type i = 0; i < D; ++i)
//    {
//        coeffs[i] = ct_accumulate(dimension_lengths,
//                                  i + 1,
//                                  D - i - 1,
//                                  static_cast<size_type>(1),
//                                  ct_prod<size_type>);
//    }
//    return coeffs;
//}

  } // end of namespace internals
} // end of corgi
