/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef _MISC_HASH_HPP_
#define _MISC_HASH_HPP_

#include <cstdint>
#include <memory>
#include <type_traits>
#include "ska/bytell_hash_map.hpp"

/*template <typename T, int TWidth = sizeof(T)>
struct HashHelpers;

template <typename T>
struct HashHelpers<T, 4> {
   static uint32_t fibonacci_muliplier() { return uint32_t(0x9e3779b9); }

   static uint32_t rotate_left(uint32_t x, int n) {
      return (x << n) | (x >> (32 - n));
   }
};

template <typename T>
struct HashHelpers<T, 8> {
   static uint64_t fibonacci_muliplier() {
      return uint64_t(0x9e3779b97f4a7c15);
   }

   static uint64_t rotate_left(uint64_t x, int n) {
      return (x << n) | (x >> (64 - n));
   }
};

template <typename T, typename U = typename std::make_unsigned<T>::type>
struct Hasher;

// only add specialization for unsigned result types
template <typename T>
struct Hasher<T, T> {
   T state;

   Hasher(T init = 0) : state(init) {}

   template <typename U,
             typename std::enable_if<std::is_integral<U>::value, int>::type = 0>
   void addValue(U val) {
      state = (HashHelpers<T>::rotate_left(state, 5) ^ T(val)) *
              HashHelpers<T>::fibonacci_muliplier();
   }

   T getHash() const { return state; }
};*/

template <typename K, typename V, typename H = std::hash<K>,
          typename E = std::equal_to<K>>
using HashMap =
    ska::bytell_hash_map<K, V, H, E, std::allocator<std::pair<K, V>>>;

template <typename T, typename H = std::hash<T>, typename E = std::equal_to<T>>
using HashSet = ska::bytell_hash_set<T, H, E, std::allocator<T>>;

#endif
