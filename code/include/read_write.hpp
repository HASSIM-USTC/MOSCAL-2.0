#include "deom.hpp"
#include "mr_level_1.hpp"

void write_ddos_key(DEOM_DATA *d, const char *fname_ddos, const char *fname_keys) {
  ofstream fstream_ddos(fname_ddos);
  ofstream fstream_keys(fname_keys);
  IOFormat HeavyFmt(FullPrecision);

  for (int i = 0; i < d->nddo; i++) {
    if (is_valid(d->ddos[i], d->ferr)) {
      fstream_ddos << d->ddos[i].format(HeavyFmt) << '\n';
      fstream_keys << d->keys[i] << '\n';
    }
  }
  fstream_ddos.close();
  fstream_keys.close();
}