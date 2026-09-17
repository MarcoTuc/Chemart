"""Adleman's DNA computation of a Hamiltonian path (Science 266:1021, 1994).

Catalog id: dna-hpp. Book section 19.3.1.

The in-silico reaction network of the published wet experiment. Molecules are
oligonucleotides, all sequences stored 5' -> 3':

- ``O<i>``    the code word of vertex i, a random 20-mer (fig. 2 gives O2, O3, O4);
- ``Obar<i>`` its Watson-Crick complement, which serves as the splint of vertex i;
- ``E<i>-<j>`` the edge oligonucleotide O_{i->j}: the 3' 10-mer of O_i followed by
  the 5' 10-mer of O_j, except that it carries all of O_i when i = v_in and all
  of O_j when j = v_out (so that O_{v_in} and Obar_{v_out} can prime the PCR);
- ``P<v0>-<v1>-...`` a path duplex: the ligated sense strand of the edge
  oligonucleotides of the walk, held together by the splints of its interior
  vertices (a two-vertex path is the bare edge oligonucleotide).

Step 1 of Adleman's algorithm is the ligation reaction, and it is the only step
that makes molecules, so it is the reaction set:

    P(v0..vk) + Obar_vk + E_vk-v  ->  P(v0..vk,v)

The splint Obar_vk bridges the junction: its 5' half pairs with the 3' 10-mer of
the growing strand and its 3' half with the 5' 10-mer of the incoming edge
oligonucleotide, so ligation extends the path by one vertex. Nucleotides are
conserved (``extras["conservation"]``). The walk set of a graph with a cycle is
infinite, so the closure is cut off at ``max_path_vertices`` vertex slots per
path and the network is reported as truncated.

Steps 2-5 select, they do not react, so they are reported in
``extras["analysis"]``: PCR with primers O_{v_in} and Obar_{v_out} (step 2), the
agarose gel band at 20 * n bp = 140 bp for the seven-vertex graph (step 3),
successive affinity purification on Obar_i beads for every intermediate vertex
(step 4) and the graduated-PCR readout (step 5), which prints a path as bands at
20 (k+1) bp for a vertex sitting in slot k.
"""

from __future__ import annotations

from chemart.helpers import params as params_helper
from chemart.kinetics import AVOGADRO
from chemart.network import Network, Reaction, Species

BASES = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

#: Figure 1: the directed graph Adleman solved, 7 vertices and 14 edges.
ADLEMAN_GRAPH = [[0, 1], [0, 3], [0, 6], [1, 2], [1, 3], [2, 1], [2, 3],
                 [3, 2], [3, 4], [4, 1], [4, 5], [5, 1], [5, 2], [5, 6]]
#: Figure 1 legend: with v_in = 0 and v_out = 6 this is the unique Hamiltonian path.
ADLEMAN_PATH = [0, 1, 2, 3, 4, 5, 6]
#: Figure 2: the three vertex code words the paper prints (5' -> 3').
PUBLISHED_SEQUENCES = {
    "2": "TATCGGATCGGTATATCCGA",
    "3": "GCTATTCGAGCTTAAAGCTA",
    "4": "GGCTAGGTACCAGCATGCTT",
}
#: Figure 2: the two edge oligonucleotides and the splint it prints. Obar3 is
#: printed 3' -> 5' in the paper so that it aligns under the sense strand; it is
#: the base-by-base complement of O3.
PUBLISHED_EDGES = {"2-3": "GTATATCCGAGCTATTCGAG", "3-4": "CTTAAAGCTAGGCTAGGTAC"}
PUBLISHED_SPLINT_3_PRINTED = "CGATAAGCTCGAATTTCGAT"


# ---------------------------------------------------------------------------
# Sequences
# ---------------------------------------------------------------------------
def complement(seq: str) -> str:
    """Base-by-base Watson-Crick complement, i.e. the opposite strand 3' -> 5'."""
    return "".join(COMPLEMENT[c] for c in seq)


def revcomp(seq: str) -> str:
    """The opposite strand written 5' -> 3'."""
    return complement(seq)[::-1]


def edge_oligo(seq_i: str, seq_j: str, i: int, j: int, v_in: int, v_out: int) -> str:
    """O_{i->j}: the 3' half of O_i then the 5' half of O_j (whole code words at the ends)."""
    half = len(seq_i) // 2
    left = seq_i if i == v_in else seq_i[half:]
    right = seq_j if j == v_out else seq_j[:half]
    return left + right


def graduated_pcr(walk, oligo_length: int = 20,
                  vertices: int | None = None) -> dict[str, list[int]]:
    """Band lengths per lane: lane i holds one band of 20 (k+1) bp per slot k of vertex i.

    Graduated PCR runs one reaction per vertex with O_{v_in} as one primer and
    Obar_i as the other, so a vertex that the path enters twice prints two bands
    and a vertex it never enters prints none (the empty lane the paper writes `x`).
    There is one lane per vertex other than the start; `vertices` gives their
    number, and defaults to the vertices the walk itself visits.
    """
    lanes = range(vertices) if vertices else sorted(set(walk))
    return {str(v): [oligo_length * (k + 1) for k, x in enumerate(walk) if x == v]
            for v in lanes if v != walk[0]}


# ---------------------------------------------------------------------------
# Species ids
# ---------------------------------------------------------------------------
def vertex_id(i: int) -> str:
    return f"O{i}"


def splint_id(i: int) -> str:
    return f"Obar{i}"


def edge_id(i: int, j: int) -> str:
    return f"E{i}-{j}"


def path_id(walk) -> str:
    return "P" + "-".join(str(v) for v in walk)


def molecule_id(walk) -> str:
    """A two-vertex path is the bare edge oligonucleotide."""
    return edge_id(*walk) if len(walk) == 2 else path_id(walk)


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
def _check(p, n: int) -> list[tuple[int, int]]:
    if p.oligo_length % 2:
        raise ValueError(
            f"oligo_length must be even (it is split into two halves), got {p.oligo_length}"
        )
    edges = params_helper.edges("graph", p.graph)
    outside = sorted({v for e in edges for v in e if not 0 <= v < n})
    if outside:
        raise ValueError(
            f"graph uses vertices {outside}, which are not in 0..{n - 1}; "
            f"raise vertices (now {n}) or renumber the graph"
        )
    for name in ("v_in", "v_out"):
        v = getattr(p, name)
        if not 0 <= v < n:
            raise ValueError(f"{name}={v} is not a vertex of the graph (0..{n - 1})")
    if p.v_in == p.v_out:
        raise ValueError(f"v_in and v_out must be different vertices, both are {p.v_in}")
    return list(dict.fromkeys(edges))


def _sequences(p, n: int, rng) -> list[str]:
    """The vertex code words: the given ones, the rest random as in the paper."""
    given: dict[int, str] = {}
    for key, value in p.sequences.items():
        if not (isinstance(key, str) and key.isdigit() and 0 <= int(key) < n):
            raise ValueError(
                f"sequences keys must be vertex ids '0'..'{n - 1}' as strings, got {key!r}"
            )
        ok = (isinstance(value, str) and len(value) == p.oligo_length
              and set(value.upper()) <= set(BASES))
        if not ok:
            raise ValueError(
                f"sequences[{key!r}] must be a string of {p.oligo_length} bases over "
                f"ACGT (oligo_length), got {value!r}"
            )
        given[int(key)] = value.upper()
    return [given.get(i) or "".join(rng.choice(list(BASES), p.oligo_length))
            for i in range(n)]


# ---------------------------------------------------------------------------
def generate(p, rng) -> Network:
    n = p.vertices
    edges = _check(p, n)
    code = _sequences(p, n, rng)
    L = p.oligo_length

    oligo = {(i, j): edge_oligo(code[i], code[j], i, j, p.v_in, p.v_out) for i, j in edges}
    out_edges: dict[int, list[int]] = {}
    for i, j in edges:
        out_edges.setdefault(i, []).append(j)

    # Step 1: every walk of at most max_path_vertices vertex slots, in order of length.
    walks: list[list[int]] = [[i, j] for i, j in edges]
    frontier = walks
    while frontier:
        grown: list[list[int]] = []
        for walk in frontier:
            if len(walk) >= p.max_path_vertices:
                continue
            grown += [walk + [k] for k in out_edges.get(walk[-1], [])]
        walks += grown
        frontier = grown
    truncated = any(len(w) == p.max_path_vertices and out_edges.get(w[-1]) for w in walks)

    sense = {molecule_id(w): "".join(oligo[(a, b)] for a, b in zip(w, w[1:])) for w in walks}
    reactions = [
        Reaction.of([molecule_id(w), splint_id(w[-1]), edge_id(w[-1], k)],
                    [molecule_id(w + [k])])
        for w in walks if len(w) < p.max_path_vertices
        for k in out_edges.get(w[-1], [])
    ]

    species = [Species(vertex_id(i), code[i]) for i in range(n)]
    species += [Species(splint_id(i), revcomp(code[i])) for i in range(n)]
    species += [Species(edge_id(i, j), oligo[(i, j)]) for i, j in edges]
    species += [Species(path_id(w), sense[path_id(w)]) for w in walks if len(w) > 2]
    # A path duplex carries its sense strand plus one splint per interior vertex.
    nucleotides = [L] * (2 * n) + [len(oligo[e]) for e in edges]
    nucleotides += [len(sense[path_id(w)]) + L * (len(w) - 2) for w in walks if len(w) > 2]

    initial = {edge_id(i, j): float(p.pmol) for i, j in edges}
    initial.update({splint_id(i): float(p.pmol)
                    for i in range(n) if i not in (p.v_in, p.v_out)})

    return Network(
        species=species,
        reactions=reactions,
        status="truncated" if truncated else "complete",
        initial_state=initial,
        extras={
            "graph": {"vertices": n, "edges": [list(e) for e in edges],
                      "v_in": p.v_in, "v_out": p.v_out},
            "encoding": {
                "oligo_length": L,
                "scheme": (
                    "vertex i -> random code word O_i; edge i->j -> the 3' half of "
                    "O_i followed by the 5' half of O_j (all of O_i when i = v_in, "
                    "all of O_j when j = v_out); the complement Obar_i splints the "
                    "junction at vertex i"
                ),
                "figure_2": {
                    "O2": PUBLISHED_SEQUENCES["2"], "O3": PUBLISHED_SEQUENCES["3"],
                    "O4": PUBLISHED_SEQUENCES["4"],
                    "O2->3": PUBLISHED_EDGES["2-3"], "O3->4": PUBLISHED_EDGES["3-4"],
                    "O3bar_printed_3_to_5": PUBLISHED_SPLINT_3_PRINTED,
                },
            },
            "conservation": [{"name": "nucleotides", "vector": nucleotides}],
            "analysis": _analysis(p, n, walks, sense),
        },
    )


def _analysis(p, n: int, walks, sense) -> dict:
    """Adleman's five steps: step 1 built the molecules, steps 2-5 select."""
    L = p.oligo_length
    from_start = [w for w in walks if w[0] == p.v_in]
    step2 = [w for w in from_start if w[-1] == p.v_out]
    band = L * n
    step3 = [w for w in step2 if len(sense[molecule_id(w)]) == band]
    order = [i for i in range(n) if i not in (p.v_in, p.v_out)]
    kept, kept_after = step3, []
    for i in order:
        kept = [w for w in kept if i in w]
        kept_after.append(len(kept))
    step5 = kept

    return {
        "algorithm": (
            "1 ligate random paths; 2 keep the paths from v_in to v_out; "
            "3 keep the paths with exactly n vertices; 4 keep the paths that "
            "enter every vertex; 5 read out what is left"
        ),
        "steps": [
            {"step": 1, "operation": "ligation (T4 DNA ligase)",
             "detail": "edge oligonucleotides joined by the vertex splints",
             "molecules": len(walks), "paths_from_v_in": len(from_start)},
            {"step": 2, "operation": "PCR",
             "detail": f"primers {vertex_id(p.v_in)} and {splint_id(p.v_out)}",
             "kept": len(step2)},
            {"step": 3, "operation": "agarose gel",
             "detail": f"excise the {band} bp band = {L}-mer code words x {n} vertices",
             "band_bp": band, "kept": len(step3)},
            {"step": 4, "operation": "affinity purification on magnetic beads",
             "detail": "successive capture on biotinylated "
                       + ", ".join(splint_id(i) for i in order),
             "order": order, "kept_after": kept_after, "kept": len(step5)},
            {"step": 5, "operation": "PCR and graduated PCR readout",
             "detail": "one band per vertex slot at 20 (k+1) bp", "kept": len(step5)},
        ],
        "hamiltonian_paths": [list(w) for w in step5],
        "solution_species": [molecule_id(w) for w in step5],
        "gel_survivors": [list(w) for w in step3],
        "graduated_pcr": {molecule_id(w): graduated_pcr(w, L, n) for w in step3},
        "truncated_at_vertices": p.max_path_vertices,
        "experiment": {
            "days_of_lab_work": 7,
            "pmol_per_oligonucleotide": float(p.pmol),
            "copies_per_oligonucleotide": float(p.pmol) * 1e-12 * AVOGADRO,
            "ligation": "50 pmol of each oligonucleotide, 5 units T4 DNA ligase, "
                        "4 hours at room temperature in 100 ul",
            "pcr": "35 cycles of 94 C 15 s / 30 C 60 s; graduated PCR 25 cycles at 40 C",
            "gel": "3% or 5% agarose in TBE with ethidium bromide",
            "beads": "biotinylated splints on streptavidin paramagnetic particles",
            "scaling": "the copies needed grow exponentially with the number of "
                       "vertices; Hartmanis (1995) put a 200-city instance at more "
                       "than the weight of the Earth in DNA",
        },
    }
