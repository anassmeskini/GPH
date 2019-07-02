#ifndef BITSET_HPP
#define BITSET_HPP

#ifndef BOOST_FOUND
#include <vector>
using bitset = std::vector<bool>;
#else
#include "boost/dynamic_bitset.hpp"
using bitset = boost::dynamic_bitset;
#endif

#endif
