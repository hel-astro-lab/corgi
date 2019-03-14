#pragma once

#include <functional>


// ref: https://stackoverflow.com/questions/40700181/sum-the-components-of-a-tuple-up-by-using-stdget-stdtuple-size-stdtuple


// namespace for utility code:
namespace utility {

  template<std::size_t...Is>
  auto index_over( std::index_sequence<Is...> ) {
    return [](auto&&f)->decltype(auto){
      return decltype(f)(f)( std::integral_constant<std::size_t,Is>{}... );
    };
  }
  template<std::size_t N>
  auto index_upto() {
    return index_over( std::make_index_sequence<N>{} );
  }

}



// namespace for semantic-equivalent replacements of `std` code:
namespace notstd {
  template<class F, class Tuple>
  decltype(auto) apply( F&& f, Tuple&& tuple ) {
    using dTuple = std::decay_t<Tuple>;
    auto index = utility::index_upto< std::tuple_size<dTuple>{} >();
    return index( [&](auto...Is)->decltype(auto){
      auto target=std::ref(f);
      return target( std::get<Is>( std::forward<Tuple>(tuple) )... );
    } ); 
  }
}


