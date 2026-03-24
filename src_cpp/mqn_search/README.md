## MQN C++ Encoder

This directory contains the C++ implementation used to compute MQN descriptors for the preprocessing pipeline.

Main files:

- [main.cpp](main.cpp): file-based executable
- [main.cpp_stdio](main.cpp_stdio): stdio-oriented variant
- [reader.cpp](reader.cpp)
- [reader.h](reader.h)

The file-based binary is used by [1_add_MQN2file.sh](../../scripts/prepare_data/1_add_MQN2file.sh).

Expected behavior:

- input: one SMILES per line
- output: semicolon-separated records containing the original SMILES, the MQN sum bin, and the MQN vector

The repo currently expects a compiled executable at `build/bin/mqn_new`.

There is no checked-in build system in the repository at the moment, but the binary has been built successfully with:

```bash
mkdir -p build/bin
cd src_cpp/mqn_search
g++ -std=c++17 -Wall -Wextra -O2 *.cpp -o ../../build/bin/mqn_new
```

The prep pipeline and benchmark helper scripts expect the resulting executable at `build/bin/mqn_new`.
