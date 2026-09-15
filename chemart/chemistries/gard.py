"""GARD (Graded Autocatalysis Replication Domain), basic lipid-world model. Catalog id: gard.

For one assembly, n_i is the count of lipid type i inside it (species A_i)
and rho_i the buffered outside concentration (species L_i):

    dn_i/dt = k_f rho_i (1 - N/N_max) + sum_j beta_ij n_j rho_i (1 - N/N_max) - k_b n_i

Each term is one mass-action reaction. The crowding factor (1 - N/N_max),
N = sum_i n_i, multiplies every joining reaction; it is recorded on those
rates as `crowding_capacity` because it is not mass action.
"""

from chemart.helpers.explicit import network


def generate(p, rng):
    G = p.N_G
    beta = rng.lognormal(mean=p.beta_mu, sigma=p.beta_sigma, size=(G, G))
    L = [f"L{i + 1}" for i in range(G)]
    A = [f"A{i + 1}" for i in range(G)]

    def join(k):
        return {"law": "mass-action", "k": float(k), "crowding_capacity": float(p.N_max)}

    reactions = []
    for i in range(G):
        reactions += [(f"{L[i]} -> {A[i]}", join(p.k_f)), (f"{A[i]} -> {L[i]}", p.k_b)]
    for i in range(G):
        for j in range(G):
            reactions.append((f"{L[i]} + {A[j]} -> {A[i]} + {A[j]}", join(beta[i, j])))
            if p.catalysed_leaving:
                reactions.append((f"{A[i]} + {A[j]} -> {L[i]} + {A[j]}", float(p.k_b * beta[i, j])))
    return network(
        reactions,
        species=L + A,
        initial_state={s: p.rho for s in L},
        extras={
            "buffered": L,
            "compartments": {
                "assembly": {
                    "members": A,
                    "N_max": p.N_max,
                    "crowding": "joining rates are multiplied by (1 - N/N_max), N = sum of A_i counts",
                    "fission": "when N reaches N_max the assembly splits into two halves; each daughter inherits a sample of the composition",
                }
            },
        },
    )
