#pragma once

#include <type_traits>
#include <array>



/*! Internal meta-programming tools
 *
 * See especially
 *  - https://codereview.stackexchange.com/questions/107877/simple-multi-dimensional-array-class-in-c11
 *  - https://github.com/maddouri/hyper_array/
 *
 */

namespace corgi {
  namespace internals {


/// shorthand for the enable_if syntax
/// @see http://en.cppreference.com/w/cpp/types/enable_if#Helper_types
//
// TODO: move to c++-14 where this is in std::
template <bool b, typename T=void>
using enable_if_t = typename std::enable_if<b, T>::type;

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
                const T      initialValue,      ///< let this be the accumulator's initial value
                const O&     op                 ///< use this binary operation
               )
{
    // https://stackoverflow.com/a/33158265/865719
    return (first < (first + length))
         ? op(arr[first],
              ct_accumulate(arr,
                            first + 1,
                            length - 1,
                            initialValue,
                            op))
         : initialValue;
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
                   const T       initialValue,   ///< let this be the summation's initial value
                   const O_SUM&  op_sum,         ///< use this as the summation operator
                   const O_PROD& op_prod         ///< use this as the multiplication operator
                  )
{
    // same logic as `ct_accumulate()`
    return (first_1 < (first_1 + length))
         ? op_sum(op_prod(arr_1[first_1],
                          arr_2[first_2]),
                  ct_inner_product(arr_1, first_1 + 1,
                                   arr_2, first_2 + 1,
                                   length - 1,
                                   initialValue,
                                   op_sum, op_prod))
         : initialValue;
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
//compute_index_coeffs(const ::std::array<size_type, D>& dimensionLengths) noexcept
//{
//    ::std::array<size_type, D> coeffs;
//    for (size_type i = 0; i < D; ++i)
//    {
//        coeffs[i] = ct_accumulate(dimensionLengths,
//                                  i + 1,
//                                  D - i - 1,
//                                  static_cast<size_type>(1),
//                                  ct_prod<size_type>);
//    }
//    return coeffs;
//}


//--------------------------------------------------
// tyepdef shortcut for checking indicies
//
// TODO: does not work; complains about <anonymous> type?
//
// template<std::size_t Dim, typename... Args>
// using check_index_length_t = 
//   typename  corgi::internals::enable_if_t<(sizeof...(Args) == Dim) &&
//             corgi::internals::are_integral<Args...>::value, void>;

template<std::size_t Dim, typename... Args>
using check_index_length_t = 
  typename enable_if_t<(sizeof...(Args) == Dim), void>::type;


//--------------------------------------------------
// N-length tuple of type T, i.e., tuple_of<3, int> = tuple<int, int, int>
// see: 
//  - https://stackoverflow.com/questions/38885406/produce-stdtuple-of-same-type-in-compile-time-given-its-length-by-a-template-a

template <size_t I,typename T> 
struct tuple_n{
    template< typename...Args> using type = typename tuple_n<I-1, T>::template type<T, Args...>;
};

template <typename T> 
struct tuple_n<0, T> {
    template<typename...Args> using type = std::tuple<Args...>;   
};
template <size_t I,typename T>  using tuple_of = typename tuple_n<I,T>::template type<>;




//--------------------------------------------------
// implementations of std::apply (from c++-17)
// See:
//  - https://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer
//
// ver1
/*
template<typename Function, typename Tuple, size_t ... I>
auto apply(Function f, Tuple t, std::index_sequence<I ...>)
{
     return f(std::get<I>(t) ...);
}

template<typename Function, typename Tuple>
auto apply(Function f, Tuple t)
{
    static constexpr auto size = std::tuple_size<Tuple>::value;
    return apply(f, t, std::make_index_sequence<size>{});
}
*/

//--------------------------------------------------
// ver2
  
template<std::size_t...Is>
auto index_over(std::index_sequence<Is...>){
  return [](auto&&f)->decltype(auto){
    return decltype(f)(f)( std::integral_constant<std::size_t, Is>{}... );
  };
}
template<std::size_t N>
auto index_upto(std::integral_constant<std::size_t, N> ={}){
  return index_over( std::make_index_sequence<N>{} );
}

//template< class T >
//constexpr std::size_t tuple_size_v = tuple_size<T>::value;
//typename std::tuple_size<std::remove_reference_t<Tuple>::value>()

template<class T>
constexpr auto tuple_size_v = std::tuple_size<T>::value;
template<class F, class Tuple>
decltype(auto) apply( F&& f, Tuple&& tup ) {
  auto indexer = index_upto<
    tuple_size_v<std::remove_reference_t<Tuple>>
  >();
  return indexer(
    [&](auto...Is)->decltype(auto) {
      return std::forward<F>(f)(
        std::get<Is>(std::forward<Tuple>(tup))...
      );
    }
  );
}


// C++-14 index_sequence (maybe?); toggle to get C++-11 compatibility
// used together with tuple unpacking in corgi::Node::id()
//template <size_t... Is>
//struct index_sequence;
// See also: 
//  - https://stackoverflow.com/questions/17424477/implementation-c14-make-integer-sequence/17426611#17426611
//  - http://aherrmann.github.io/programming/2016/02/28/unpacking-tuples-in-cpp14/

//--------------------------------------------------
// array into tuple conversion
//
// See:
//  - https://stackoverflow.com/questions/37029886/how-to-construct-a-tuple-from-an-array
//
//  try1
//template <class... Formats, size_t N, size_t... Is>
//std::tuple<Formats...> as_tuple(std::array<char*, N> const& arr,
//                                std::index_sequence<Is...>)
//{
//    return std::make_tuple(Formats{arr[Is]}...);
//}
//
//template <class... Formats, size_t N,
//          class = enable_if_t<(N == sizeof...(Formats))>>
//std::tuple<Formats...> as_tuple(std::array<char*, N> const& arr)
//{
//    return as_tuple<Formats...>(arr, std::make_index_sequence<N>{});
//}

// try2
// See:
//  - https://stackoverflow.com/questions/41207774/how-do-i-create-a-tuple-of-n-ts-from-an-array-of-t
template<std::size_t... I, typename U>
constexpr auto into_tuple(const U &arr, std::index_sequence<I...>) {
    return std::make_tuple(arr[I]...);
}

template<typename T, std::size_t N>
constexpr auto into_tuple(const T (&arr)[N]) {
    return into_tuple(arr, std::make_index_sequence<N>{});
}

template<typename T, std::size_t N>
constexpr auto into_tuple(const std::array<T, N> &arr) {
    return into_tuple(arr, std::make_index_sequence<N>{});
}

//--------------------------------------------------
// tuple into array
// See:
//  - https://stackoverflow.com/questions/10604794/convert-stdtuple-to-stdarray-c11
//template<typename First, typename... Rem>
//std::array<First, 1+sizeof...(Rem)>
//into_array(const std::tuple<First, Rem...>& t) {
//  std::array<First, 1+sizeof...(Rem)> arr;
//  ArrayFiller<First, decltype(t), 1+sizeof...(Rem)>::into_array(t, arr);
//  return arr;
//}

// Convert tuple into a array implementation
template<typename T, std::size_t N, typename Tuple,  std::size_t... I>
constexpr decltype(auto) into_array_impl(const Tuple& a, std::index_sequence<I...>)
{
  return std::array<T,N>{std::get<I>(a)...};
}

// Convert tuple into a array
template<typename Head, typename... T>
constexpr decltype(auto) into_array(const std::tuple<Head, T...>& a)
{
  using Tuple = std::tuple<Head, T...>;
  constexpr auto N = sizeof...(T) + 1;
  return into_array_impl<Head, N, Tuple>(a, std::make_index_sequence<N>());
}






  } // end of namespace internals
} // end of corgi
