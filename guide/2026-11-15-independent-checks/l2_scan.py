"""Guide's floating-point corroboration of oddlambda2: max over all graphs of odd order n of the second largest
eigenvalue of A_alpha equals (n-3)/2 + alpha, attained uniquely by K_1 join 2K_m for alpha > 0 (numpy, margin 1e-9)."""
import sys
import numpy as np
n = int(sys.argv[1]); alphas = [0.0, 0.05, 0.125, 0.25, 1/3, 0.4, 0.49]
def decode(line):
    data = [c - 63 for c in line]; m = data[0]; bits = []
    for x in data[1:]: bits.extend((x >> k) & 1 for k in range(5, -1, -1))
    A = np.zeros((m, m)); i = 0
    for v in range(1, m):
        for u in range(v):
            if bits[i]: A[u, v] = A[v, u] = 1
            i += 1
    return A
lines = [l.strip() for l in sys.stdin.buffer if l.strip()]
A = np.array([decode(l) for l in lines]); D = np.zeros_like(A); idx = np.arange(n); D[:, idx, idx] = A.sum(axis=2)
for a in alphas:
    l2 = np.linalg.eigvalsh(a * D + (1 - a) * A)[:, -2]
    target = (n - 3) / 2 + a
    eq = np.where(np.abs(l2 - target) < 1e-9)[0]
    print(f"n={n} alpha={a:.4f} max={l2.max():.9f} target={target:.9f} violations={(l2 > target + 1e-9).sum()} equality={len(eq)} {[lines[i].decode() for i in eq[:3]]}")
