#pragma once
#include <folly/AtomicHashMap.h>
#include <utility>

using namespace folly;
using namespace std;

class Trie {
 private:
  AtomicHashMap<int64_t, int32_t> tree;

 public:
  explicit Trie(size_t nmax) : tree((size_t)nmax) {}

  void insert(const int64_t key_hash, const int32_t rank) {
    tree.insert({key_hash, rank});
  }

  pair<int32_t, bool> try_insert(const int64_t key_hash, const int32_t rank) {
    const auto [rank_find, success] = tree.insert({key_hash, rank});
    return make_pair(rank_find->second, success);
  }

  int32_t find(const int64_t key_hash) const {
    auto search = tree.find(key_hash);
    if (search != tree.end()) {
      return search->second;
    } else {
      return -1;
    }
  }

  void clear() { tree.clear(); }

  int size() { return tree.size(); }
};