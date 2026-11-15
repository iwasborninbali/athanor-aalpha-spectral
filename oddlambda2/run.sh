#!/bin/sh
# Replay from this directory: sh run.sh
set -eu
cd "$(dirname "$0")"
python3 - <<'PY'
import networkx as nx
atlas = nx.graph_atlas_g()
for n in (5, 7):
    with open('atlas%d.g6' % n, 'wb') as out:
        for g in atlas:
            if len(g) == n:
                out.write(nx.to_graph6_bytes(g, header=False))
PY
cc -O3 -std=c11 -o enumerate9 enumerate9.c
c++ -O3 -std=c++17 -o exact_scan exact_scan.cpp
./enumerate9 atlas7.g6 graphs9.g6 2> enumerate9.replay.log
{ ./exact_scan 5 atlas5.g6; ./exact_scan 7 atlas7.g6; ./exact_scan 9 graphs9.g6; } > exact_scan.replay.log
python3 sympy_check.py > sympy_check.replay.log
diff -u enumerate9.log enumerate9.replay.log
diff -u exact_scan.log exact_scan.replay.log
diff -u sympy_check.log sympy_check.replay.log
printf '%s  %s\n' '7d9f2373e60e51f1170f5c9e93bcd76cea0491245eb95d00685bb98ba775aa9e' graphs9.g6 | sha256sum -c -
sha256sum -c SHA256SUMS
