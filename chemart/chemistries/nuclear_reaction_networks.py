"""Nuclear reaction networks (book 20.2). Catalog id: nuclear-reaction-networks.

Topology only: reaction rates depend on temperature and density through
cross sections (e.g. REACLIB), which the book does not give. Baryon number,
charge (fully ionised nuclei) and electron lepton number are exact
conservation laws and are attached to every network.
"""

from chemart.helpers.explicit import network

#: species -> (baryon number, charge, electron lepton number)
PARTICLES = {
    "n": (1, 0, 0), "1H": (1, 1, 0), "2H": (2, 1, 0), "3H": (3, 1, 0),
    "3He": (3, 2, 0), "4He": (4, 2, 0), "7Li": (7, 3, 0), "7Be": (7, 4, 0),
    "12C": (12, 6, 0), "13C": (13, 6, 0), "13N": (13, 7, 0), "14N": (14, 7, 0),
    "15N": (15, 7, 0), "15O": (15, 8, 0),
    "e+": (0, 1, -1), "e-": (0, -1, 1), "nu_e": (0, 0, 1), "anti_nu_e": (0, 0, -1),
    "gamma": (0, 0, 0),
}

NETWORKS = {
    # Book eq. 20.4 (Bethe-Weizsaecker cycle).
    "cno-cycle": [
        "12C + 1H -> 13N + gamma",
        "13N -> 13C + e+ + nu_e",
        "13C + 1H -> 14N + gamma",
        "14N + 1H -> 15O + gamma",
        "15O -> 15N + e+ + nu_e",
        "15N + 1H -> 12C + 4He",
    ],
    # Book eq. 20.3 is the first step; the rest is the textbook pp-I branch.
    "pp-chain": [
        "1H + 1H -> 2H + e+ + nu_e",
        "2H + 1H -> 3He + gamma",
        "3He + 3He -> 4He + 2 1H",
    ],
    # The 12 key reactions of big-bang nucleosynthesis (book fig. 20.4).
    "big-bang-nucleosynthesis": [
        "n -> 1H + e- + anti_nu_e",
        "1H + n -> 2H + gamma",
        "2H + 1H -> 3He + gamma",
        "2H + 2H -> 3He + n",
        "2H + 2H -> 3H + 1H",
        "3H + 2H -> 4He + n",
        "3He + n -> 3H + 1H",
        "3He + 2H -> 4He + 1H",
        "3He + 4He -> 7Be + gamma",
        "3H + 4He -> 7Li + gamma",
        "7Be + n -> 7Li + 1H",
        "7Li + 1H -> 2 4He",
    ],
}

#: book 20.2: the CNO structure is open, with 1H flowing in and 4He, e+, nu_e out
BOUNDARY = {"cno-cycle": {"inflow": ["1H"], "outflow": ["4He", "e+", "nu_e"]}}


def generate(p, rng):
    net = network([(text, None) for text in NETWORKS[p.network]])
    present = [s.id for s in net.species]
    conservation = []
    for index, name in enumerate(("baryon number", "charge", "electron lepton number")):
        vector = {s: PARTICLES[s][index] for s in present if PARTICLES[s][index]}
        if vector:
            conservation.append({"name": name, "vector": vector})
    net.extras = {"conservation": conservation}
    if p.network in BOUNDARY:
        net.extras["boundary"] = BOUNDARY[p.network]
    return net
