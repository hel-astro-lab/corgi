#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

namespace corgi { namespace internals {

//! Efficiently remove an element from a vector without
//! preserving order. If the element is not the last element
//! in the vector, transfer the last element into its position
//! using a move if possible.
//! Regardless, we then shrink the size of the vector deleting
//! the element at the end, which will either be destructed or
//! the element we were deleting.
//! @note: Effectively invalidates the current iterator.
template<class value_t>
bool unstable_remove(
    typename std::vector<value_t>& container,
    typename std::vector<value_t>::iterator it
    )
{
    // Leave in-situ if we are already the tail element.
    auto last_element = container.end() - 1;
    if (it != last_element) {
        // overwrite this element with what is in the last,
        // which should have the same effect as deleting this.
        *it = std::move(*last_element);
    }
    // release the last cell of the vector, because it should
    // now either be destructed or contain the value we were
    // deleting.
    container.pop_back();
}

} } // ns corgi::internals
