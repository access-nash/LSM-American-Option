import numpy as np

# parameters (from Longstaff Schwarz paper)
S0 = 36.0
K = 40.0
T = 1.0
r = 0.06
sigma = 0.2
M = 50          # number of time steps
I = 100000      # number of paths
D = 3           # number of basis functions
seed = 42

########################################################################################################################
# the section below uses algo from Yves Hilpisch book

# American put payoff h_t(S)
def payoff(S):
    return np.maximum(K - S, 0.0)

# MC
def simulate_paths(S0, r, sigma, T, M, I, seed):
    rng = np.random.default_rng(seed)
    dt = T / M
    S = np.zeros((M + 1, I))
    S[0] = S0
    for t in range(1, M + 1):
        z = rng.standard_normal(I)
        S[t] = S[t - 1] * np.exp((r - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * z)
    return S

    # b_d(S) = S^d, d = 0..D-1 # this is th eversion of basis function in Yves Hilpsich book
def basis(S, D):
    return np.column_stack([S ** d for d in range(D)])


#LSM
dt = T / M
disc = np.exp(-r * dt)

S = simulate_paths(S0, r, sigma, T, M, I, seed)

V = payoff(S[M]) # terminal value

# backward induction for T > t > 0
for t in range(M - 1, 0, -1):
    Y = disc * V                                   # discounted continuation value Y_{t,i}
    B = basis(S[t], D)
    alpha, *_ = np.linalg.lstsq(B, Y, rcond=None)
    C_hat = B @ alpha
    h = payoff(S[t])
    V = np.where(h > C_hat, h, Y)                  # exercise vs continue


# LSM estimator at t = 0
V0_i = disc * V
V0_LSM = V0_i.mean()
se = V0_i.std(ddof=1) / np.sqrt(I)

print(f"LSM American put value: {V0_LSM:.4f} (s.e. {se:.4f})")
########################################################################################################################
# The section below here uses algo from original LSM paper - 2 main differences compared to book algo -
# i) only in the money payoffs are used in LSM paper ii) basis function in original paper uses some Laguerre polynomials

def simulate_paths_antithetic(S0, r, sigma, T, M, I, seed):
    rng = np.random.default_rng(seed)
    dt = T / M
    S = np.zeros((M + 1, I))
    S[0] = S0
    for t in range(1, M + 1):
        z = rng.standard_normal(I // 2)
        z = np.concatenate([z, -z])
        S[t] = S[t - 1] * np.exp((r - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * z)
    return S


def basis_laguerre(S, K):
    X = S / K
    w = np.exp(-X / 2)
    return np.column_stack([np.ones_like(X), w, w * (1 - X), w * (1 - 2 * X + X ** 2 / 2)])


# LSM
dt = T / M
disc = np.exp(-r * dt)

S = simulate_paths_antithetic(S0, r, sigma, T, M, I, seed)

V = payoff(S[M])

for t in range(M - 1, 0, -1):
    Y = disc * V
    h = payoff(S[t])

    itm = h > 0
    V = Y.copy()
    if itm.sum() > 0:
        B = basis_laguerre(S[t][itm], K)
        alpha, *_ = np.linalg.lstsq(B, Y[itm], rcond=None)
        C_hat = B @ alpha
        exercise = h[itm] > C_hat
        idx = np.where(itm)[0][exercise]
        V[idx] = h[idx]

V0_i = disc * V
V0_LSM = V0_i.mean()
se = V0_i.std(ddof=1) / np.sqrt(I)

print(f"LSM (paper) American put: {V0_LSM:.4f} (s.e. {se:.4f})")

# checking whether basis function makes any difference compared to book
def basis_power(x, D):
    return np.column_stack([x ** d for d in range(D)])

V = payoff(S[M])

for t in range(M - 1, 0, -1):
    Y = disc * V
    h = payoff(S[t])
    itm = h > 0

    V = Y.copy()
    if itm.sum() > 0:
        x = S[t][itm]
        B = basis_power(x, D)
        alpha, *_ = np.linalg.lstsq(B, Y[itm], rcond=None)
        C_hat = B @ alpha
        exercise = h[itm] > C_hat
        idx = np.where(itm)[0][exercise]
        V[idx] = h[idx]

V0_i = disc * V
V0_LSM = V0_i.mean()
se = V0_i.std(ddof=1) / np.sqrt(I)

print(f"LSM (paper) American put with Yves basis function: {V0_LSM:.4f} (s.e. {se:.4f})")