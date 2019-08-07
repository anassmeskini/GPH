#ifndef HASH_HPP_
#define HASH_HPP_

#include "ska/bytell_hash_map.hpp"
#include <cstdint>
#include <memory>
#include <type_traits>

template<typename K,
         typename V,
         typename H = std::hash<K>,
         typename E = std::equal_to<K>>
using HashMap =
  ska::bytell_hash_map<K, V, H, E, std::allocator<std::pair<K, V>>>;

template<typename T, typename H = std::hash<T>, typename E = std::equal_to<T>>
using HashSet = ska::bytell_hash_set<T, H, E, std::allocator<T>>;

#endif
