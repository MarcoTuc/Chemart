"""Aevol (Knibbe, Parsons & Beslon), book 18.1.2. Catalog id: aevol.

A population of digital organisms, each a circular double-stranded binary
chromosome, evolving on a toroidal grid. The whole chemistry is the decoding
of one genome into a metabolism:

- **transcription**: a promoter is a 22-base window within ``PROM_MAX_DIFF``
  mismatches of the consensus ``0101011001110010010110`` (on either strand);
  the number of mismatches d fixes the transcription level 1 - d/5. The RNA
  runs from the base after the promoter to the first terminator, an 11-base
  stem-loop (4 complementary pairs around a 3-base loop).
- **translation**: inside an RNA, a ribosome binding site is the
  Shine-Dalgarno motif ``011011``, a 4-base spacer and the START codon
  ``000``; the gene is read in 3-base codons from there to the first STOP
  codon ``001``.
- **folding**: the codons are an alphabet of six "amino acids", two per
  parameter: M0/M1, W0/W1, H0/H1 (START also codes for H0). Each kind is read
  as a Gray-coded binary number, normalised to [0, 1] by its own number of
  codons, and rescaled: mean m in [0, 1], half-width w in [0, w_max], height
  h in [-1, 1]. The protein is the triangular fuzzy set (m, w, h): a narrow
  tall triangle is a specialised enzyme, a wide flat one a polyvalent protein
  of low efficiency.
- **phenotype**: the triangles of all functional proteins, each scaled by its
  transcription level, are summed (activators and inhibitors separately,
  each clipped to +-1, then added and clipped below at 0) on a grid of
  ``env_sampling`` points of [0, 1].
- **fitness**: the metabolic error is the area between phenotype and the
  environmental target (a sum of gaussians); fitness = exp(-k * error).

Each generation every grid cell is refilled by an offspring of a parent drawn
from its 3x3 neighbourhood with probability proportional to fitness. The
offspring's genome undergoes local mutations (switch, small insertion, small
deletion) and chromosomal rearrangements (duplication, deletion,
translocation, inversion), each drawn Binomial(genome length, rate).

The model is a port of the aevol source (BASE_2 "standard" flavour), see the
catalog entry. The run is returned as the observed network of a population
run:

- ``P + V -> P + O``  the organism P of a neighbouring cell replicates and its
  offspring O takes the place of the cell's previous occupant V,
- ``P -> O``          the cell's own occupant replicates in place,
- ``G -> G + Prot...`` expression: the proteome a genotype folds into.

Species are genotypes (``G<length>-<8 hex of sha1>``, structure = the genome)
and proteins (``P<codons>-<8 hex of sha1>``, structure = the coding bits and
the folded triple m, w, h).
"""

from __future__ import annotations

import hashlib
import math
from collections import Counter

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from chemart.network import Network, Reaction, Species

# ---------------------------------------------------------------------------
# constants: aevol src/libaevol/macros.h and Promoter.h (BASE_2)

PROM_SEQ = "0101011001110010010110"
PROM_SIZE = 22
PROM_MAX_DIFF = 4
TERM_STEM_SIZE = 4
TERM_LOOP_SIZE = 3
TERM_SIZE = 2 * TERM_STEM_SIZE + TERM_LOOP_SIZE          # 11
SHINE_DAL_SEQ = "011011"
SHINE_DAL_SIZE = 6
SHINE_START_SPACER = 4
CODON_SIZE = 3
RBS_SIZE = SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE      # 13
DO_TRANSLATION_LOOP = SHINE_DAL_SIZE + SHINE_START_SPACER + 3 * CODON_SIZE   # 19
MAX_CODONS = 64 * 3                                       # translation cap

CODON_START, CODON_STOP = 0b000, 0b001
CODON_W0, CODON_W1 = 0b010, 0b011
CODON_M0, CODON_M1 = 0b100, 0b101
CODON_H0, CODON_H1 = 0b110, 0b111

X_MIN, X_MAX = 0.0, 1.0
Y_MIN, Y_MAX = 0.0, 1.0
H_MIN, H_MAX = -1.0, 1.0
W_MIN = 0.0

LEADING, LAGGING = 0, 1

#: aevol examples/basic/param.in: the classic three-gaussian target.
DEFAULT_TARGET = [[1.2, 0.52, 0.12], [-1.4, 0.5, 0.07], [0.3, 0.8, 0.03]]

_PROM = np.array([int(c) for c in PROM_SEQ], dtype=np.uint8)
_PROM_LAG = 1 - _PROM                                     # PROM_SEQ_LAG
#: SHINE_DAL_SEQ_LEAD "0110111111000": the 6 SD bases, 4 free spacer bases, START.
_RBS_IDX = np.array([0, 1, 2, 3, 4, 5, 10, 11, 12])
_RBS_LEAD = np.array([int(c) for c in SHINE_DAL_SEQ + "0000" + "000"], dtype=np.uint8)
_RBS_LAG = 1 - _RBS_LEAD


def bits(text: str) -> np.ndarray:
    bad = set(text) - {"0", "1"}
    if bad:
        raise ValueError(f"a genome is a string of 0s and 1s, got {sorted(bad)}")
    return np.frombuffer(text.encode(), dtype=np.uint8) - ord("0")


def text(g: np.ndarray) -> str:
    return "".join("1" if b else "0" for b in g)


def genotype_id(g: np.ndarray) -> str:
    """Identity of a genome: its length and a hash of its bases."""
    return f"G{g.size}-{hashlib.sha1(g.tobytes()).hexdigest()[:8]}"


def _windows(g: np.ndarray, size: int, backwards: bool = False) -> np.ndarray:
    """Circular windows: row p is g[p:p+size] (or g[p-size+1:p+1] backwards)."""
    if backwards:
        padded = np.concatenate([g[g.size - size + 1:], g])
    else:
        padded = np.concatenate([g, g[:size - 1]])
    return sliding_window_view(padded, size)


# ---------------------------------------------------------------------------
# signals on the chromosome (Promoter.cpp, TranscriptionTerminationSequence.cpp,
# TranslationInitiationSequence.cpp)

def promoter_distances(g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mismatches to the promoter consensus at every position, per strand."""
    lead = (_windows(g, PROM_SIZE) != _PROM).sum(axis=1)
    lag = (_windows(g, PROM_SIZE, backwards=True) != _PROM_LAG[::-1]).sum(axis=1)
    return lead.astype(np.int64), lag.astype(np.int64)


def terminators(g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Stem-loop terminators: 4 complementary pairs around a 3-base loop.

    A lagging terminator at q is a leading terminator at q - (TERM_SIZE - 1),
    since both test the same four pairs of positions.
    """
    w = _windows(g, TERM_SIZE)
    lead = np.ones(g.size, dtype=bool)
    for t in range(TERM_STEM_SIZE):
        lead &= w[:, t] != w[:, TERM_SIZE - 1 - t]
    return lead, np.roll(lead, TERM_SIZE - 1)


def rbs_sites(g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Shine-Dalgarno + 4-base spacer + START codon, per strand."""
    lead = (_windows(g, 13)[:, _RBS_IDX] == _RBS_LEAD[_RBS_IDX]).all(axis=1)
    back = _windows(g, 13, backwards=True)[:, 12 - _RBS_IDX]
    lag = (back == _RBS_LAG[_RBS_IDX]).all(axis=1)
    return lead, lag


class Rna:
    """One transcript: promoter position, strand, transcription level, length."""

    __slots__ = ("pos", "strand", "basal", "length", "start")

    def __init__(self, pos: int, strand: int, basal: float, length: int, start: int):
        self.pos, self.strand, self.basal, self.length, self.start = pos, strand, basal, length, start

    def is_coding(self, proteins) -> bool:
        return any(p.rna is self for p in proteins)


class Protein:
    """One gene product: its codons and the triangle (m, w, h) it folds into."""

    __slots__ = ("first_aa", "strand", "codons", "e", "m", "w", "h", "functional", "rna")

    def __init__(self, first_aa: int, strand: int, codons: list[int], e: float, rna: Rna):
        self.first_aa, self.strand, self.codons, self.e, self.rna = first_aa, strand, codons, e, rna
        self.m = self.w = self.h = 0.0
        self.functional = False

    @property
    def code(self) -> str:
        return "".join(f"{c:03b}" for c in self.codons)

    @property
    def id(self) -> str:
        return f"P{len(self.codons)}-{hashlib.sha1(self.code.encode()).hexdigest()[:8]}"

    def structure(self) -> str:
        return (f"code={self.code} codons={len(self.codons)} m={self.m:.6g} "
                f"w={self.w:.6g} h={self.h:.6g} functional={int(self.functional)}")


def transcribe(g: np.ndarray) -> list[Rna]:
    """Every promoter with the transcript it produces (Individual::do_trancription)."""
    L = g.size
    if L < PROM_SIZE:
        return []
    dist = promoter_distances(g)
    term = terminators(g)
    out: list[Rna] = []
    for strand in (LEADING, LAGGING):
        positions = np.flatnonzero(dist[strand] <= PROM_MAX_DIFF)
        stops = np.flatnonzero(term[strand])
        if positions.size == 0 or stops.size == 0:
            continue
        if strand == LEADING:                             # the first terminator after the promoter
            starts = (positions + PROM_SIZE) % L
            i = np.searchsorted(stops, starts)
            gaps = np.where(i == stops.size, stops[0] + L - starts,
                            stops[np.minimum(i, stops.size - 1)] - starts)
        else:                                             # ... scanning backwards
            starts = (positions - PROM_SIZE) % L
            i = np.searchsorted(stops, starts, side="right") - 1
            gaps = np.where(i < 0, starts + L - stops[-1], starts - stops[np.maximum(i, 0)])
        levels = 1.0 - dist[strand][positions] / (PROM_MAX_DIFF + 1)
        for pos, start, gap, basal in zip(positions.tolist(), starts.tolist(),
                                          gaps.tolist(), levels.tolist()):
            out.append(Rna(pos, strand, basal, gap + TERM_SIZE, start))
    return out


def _codon(g: np.ndarray, pos: int, strand: int) -> int:
    L = g.size
    if strand == LEADING:
        b = (g[pos % L], g[(pos + 1) % L], g[(pos + 2) % L])
    else:
        b = (1 - g[pos % L], 1 - g[(pos - 1) % L], 1 - g[(pos - 2) % L])
    return int(b[0]) * 4 + int(b[1]) * 2 + int(b[2])


def _read_gene(g: np.ndarray, first_aa: int, strand: int, budget: int) -> list[int] | None:
    """Codons from the first amino acid to the STOP codon, None if there is none."""
    L = g.size
    codons: list[int] = []
    cur = first_aa
    while budget >= CODON_SIZE:
        value = _codon(g, cur, strand)
        if value == CODON_STOP:
            return codons or None
        codons.append(value)
        cur = (cur + CODON_SIZE) % L if strand == LEADING else (cur - CODON_SIZE) % L
        budget -= CODON_SIZE
    return None


def fold(codons: list[int], w_max: float) -> tuple[float, float, float, bool]:
    """Codons -> (mean, half-width, height, functional) (Individual::translate_protein).

    Each of M, W and H is a Gray-coded binary number read from its own codons,
    normalised by their number and rescaled to its range.
    """
    values = {"m": 0, "w": 0, "h": 0}
    counts = {"m": 0, "w": 0, "h": 0}
    gray = {"m": 0, "w": 0, "h": 0}
    kind = {CODON_M0: ("m", 0), CODON_M1: ("m", 1), CODON_W0: ("w", 0), CODON_W1: ("w", 1),
            CODON_H0: ("h", 0), CODON_START: ("h", 0), CODON_H1: ("h", 1)}
    for c in codons[:MAX_CODONS]:
        if c not in kind:
            continue
        key, bit = kind[c]
        counts[key] += 1
        gray[key] ^= bit
        values[key] = values[key] * 2 + gray[key]
    m = values["m"] / (2 ** counts["m"] - 1) if counts["m"] else 0.5
    w = values["w"] / (2 ** counts["w"] - 1) if counts["w"] else 0.0
    h = values["h"] / (2 ** counts["h"] - 1) if counts["h"] else 0.5
    m = (X_MAX - X_MIN) * m + X_MIN
    w = (w_max - W_MIN) * w + W_MIN
    h = (H_MAX - H_MIN) * h + H_MIN
    functional = bool(counts["m"] and counts["w"] and counts["h"] and w != 0.0 and h != 0.0)
    return m, w, h, functional


def translate(g: np.ndarray, rnas: list[Rna], w_max: float) -> list[Protein]:
    """Proteins of every transcript; transcripts sharing a gene pool their levels."""
    L = g.size
    sites = [np.flatnonzero(s) for s in rbs_sites(g)]
    found: dict[tuple[int, int], Protein] = {}
    order: list[Protein] = []
    for rna in rnas:
        span = rna.length - DO_TRANSLATION_LOOP
        if rna.length < 21 or span <= 0:
            continue
        pool = sites[rna.strand]
        if pool.size == 0:
            continue
        # the transcript covers span positions from its first transcribed base
        if rna.strand == LEADING:
            low, high = rna.start, rna.start + span              # [low, high)
        else:
            low, high = rna.start - span + 1, rna.start + 1
        if 0 <= low and high <= L:
            chosen = pool[(pool >= low) & (pool < high)]
        else:
            chosen = pool[(pool >= low % L) | (pool < high % L)]
        offsets = (chosen - rna.start) % L if rna.strand == LEADING else (rna.start - chosen) % L
        keep = np.argsort(offsets, kind="stable")
        for offset, c in zip(offsets[keep].tolist(), chosen[keep].tolist()):
            first_aa = (c + RBS_SIZE) % L if rna.strand == LEADING else (c - RBS_SIZE) % L
            codons = _read_gene(g, first_aa, rna.strand, rna.length - offset - RBS_SIZE)
            if codons is None:
                continue
            key = (rna.strand, first_aa)
            if key in found:                     # same gene read from another promoter
                found[key].e += rna.basal
                continue
            prot = Protein(first_aa, rna.strand, codons, rna.basal, rna)
            prot.m, prot.w, prot.h, prot.functional = fold(codons, w_max)
            found[key] = prot
            order.append(prot)
    return order


# ---------------------------------------------------------------------------
# fuzzy sets (phenotype/fuzzy/Discrete_Double_Fuzzy.cpp)

def add_triangle(points: np.ndarray, mean: float, half_width: float, height: float) -> None:
    """Add the triangle (mean, half-width, height) to a sampled fuzzy set."""
    n = points.size - 1
    if abs(half_width) < 1e-15 or abs(height) < 1e-15:
        return
    x0, x1, x2 = mean - half_width, mean, mean + half_width
    a0 = min(max(math.ceil(x0 * n), 0), n)
    a1 = min(max(math.ceil(x1 * n), 0), n)
    i = np.arange(a0, a1)
    points[i] += ((i / n) - x0) / (x1 - x0) * height
    b_end = math.ceil(x2 * n)
    if b_end > n:
        points[n - 1] += height * ((x2 - 1.0) / (x2 - x1))
    j = np.arange(a1, min(b_end, n))
    points[j] += height * ((x2 - (j / n)) / (x2 - x1))


def geometric_area(points: np.ndarray) -> float:
    """Integral of |f| by the trapezium rule on the sampling grid."""
    n = points.size - 1
    return float(np.abs((points[:-1] + points[1:]) / (2 * n)).sum())


def phenotype(proteins: list[Protein], sampling: int) -> np.ndarray:
    """Sum of the protein triangles (Individual::compute_absolute_phenotype)."""
    activ = np.zeros(sampling)
    inhib = np.zeros(sampling)
    for prot in proteins:
        if not prot.functional or abs(prot.w) < 1e-15 or abs(prot.h) < 1e-15:
            continue
        add_triangle(activ if prot.h > 0 else inhib, prot.m, prot.w, prot.h * prot.e)
    np.minimum(activ, Y_MAX, out=activ)
    np.maximum(inhib, -Y_MAX, out=inhib)
    return np.maximum(activ + inhib, Y_MIN)


def target_function(gaussians, sampling: int) -> np.ndarray:
    """The environmental target, sampled and clipped to [0, 1] (PhenotypicTarget::build)."""
    x = np.arange(sampling) / (sampling - 1)
    y = np.zeros(sampling)
    for height, mean, width in gaussians:
        y += height * np.exp(-(x - mean) ** 2 / (2 * width * width))
    return np.clip(y, Y_MIN, Y_MAX)


def check_target(value) -> list[list[float]]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"target must be a non-empty list of [height, mean, width] gaussians, got {value!r}")
    out = []
    for item in value:
        if (not isinstance(item, (list, tuple)) or len(item) != 3
                or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in item)):
            raise ValueError(f"every target gaussian must be [height, mean, width], got {item!r}")
        if item[2] <= 0:
            raise ValueError(f"a target gaussian needs a positive width, got {item!r}")
        out.append([float(v) for v in item])
    return out


# ---------------------------------------------------------------------------
# mutations (biochemistry/Dna.cpp and mutations/mutators/DnaMutator.cpp)

def do_switch(g: np.ndarray, pos: int) -> np.ndarray:
    g[pos] ^= 1
    return g


def do_small_insertion(g: np.ndarray, pos: int, seq: np.ndarray) -> np.ndarray:
    return np.concatenate([g[:pos], seq, g[pos:]])


def do_small_deletion(g: np.ndarray, pos: int, nb: int) -> np.ndarray:
    if pos + nb <= g.size:
        return np.concatenate([g[:pos], g[pos + nb:]])
    return g[nb - g.size + pos:pos]                       # the deletion spans the origin


def _segment(g: np.ndarray, pos_1: int, pos_2: int) -> np.ndarray:
    if pos_1 < pos_2:
        return g[pos_1:pos_2]
    return np.concatenate([g[pos_1:], g[:pos_2]])


def do_duplication(g: np.ndarray, pos_1: int, pos_2: int, pos_3: int) -> np.ndarray:
    segment = _segment(g, pos_1, pos_2)
    return np.concatenate([g[:pos_3], segment, g[pos_3:]])


def do_deletion(g: np.ndarray, pos_1: int, pos_2: int) -> np.ndarray:
    if pos_1 < pos_2:
        return np.concatenate([g[:pos_1], g[pos_2:]])
    return g[pos_2:pos_1]


def do_inversion(g: np.ndarray, pos_1: int, pos_2: int) -> np.ndarray:
    return np.concatenate([g[:pos_1], 1 - g[pos_1:pos_2][::-1], g[pos_2:]])


def do_translocation(g: np.ndarray, pos_1: int, pos_2: int, pos_3: int, pos_4: int,
                     invert: bool) -> np.ndarray:
    """ABCDE -> ADCBE, ADB'C'E or AC'D'BE, after Dna::do_translocation."""
    pos_min = min(pos_1, pos_2, pos_3, pos_4)
    if not invert:
        if pos_min == pos_1:
            b, c, d, e = pos_1, pos_3, pos_2, pos_4
        elif pos_min == pos_2:
            b, c, d, e = pos_2, pos_4, pos_1, pos_3
        elif pos_min == pos_3:
            b, c, d, e = pos_3, pos_2, pos_4, pos_1
        else:
            b, c, d, e = pos_4, pos_1, pos_3, pos_2
        return np.concatenate([g[:b], g[d:e], g[c:d], g[b:c], g[e:]])
    if pos_min in (pos_1, pos_2):
        b, c, d, e = (pos_1, pos_3, pos_2, pos_4) if pos_min == pos_1 else (pos_2, pos_4, pos_1, pos_3)
        return np.concatenate([g[:b], g[d:e], 1 - g[b:c][::-1], 1 - g[c:d][::-1], g[e:]])
    b, c, d, e = (pos_3, pos_2, pos_4, pos_1) if pos_min == pos_3 else (pos_4, pos_1, pos_3, pos_2)
    return np.concatenate([g[:b], 1 - g[c:d][::-1], 1 - g[d:e][::-1], g[b:c], g[e:]])


MUTATION_TYPES = ("duplication", "deletion", "translocation", "inversion",
                  "switch", "small_insertion", "small_deletion")


class Mutator:
    """Draws the mutations of one replication, then applies them in order.

    Counts are Binomial(genome length, rate) per type; rearrangements are drawn
    (and their positions chosen) first, in a random order, then local mutations,
    each with the genome size left by the previous ones.
    """

    def __init__(self, rng, rates: dict, max_indel_size: int, limits: tuple[int, int]):
        self.rng = rng
        self.rates = rates
        self.max_indel_size = max_indel_size
        self.min_length, self.max_length = limits

    def draw(self, length: int) -> list[tuple]:
        events: list[tuple] = []
        size = length
        for group in (("duplication", "deletion", "translocation", "inversion"),
                      ("switch", "small_insertion", "small_deletion")):
            urn = [int(self.rng.binomial(size, self.rates[k])) for k in group]
            total = sum(urn)
            while total > 0:
                pick = int(self.rng.integers(total))
                total -= 1
                acc = 0
                for i, kind in enumerate(group):
                    acc += urn[i]
                    if pick < acc:
                        urn[i] -= 1
                        break
                event = self._event(kind, size)
                if event is not None:
                    events.append(event)
                    size = self._size_after(event, size)
        return events

    def _pair(self, size: int) -> tuple[int, int]:
        pos_1 = int(self.rng.integers(size))
        pos_2 = pos_1
        while pos_2 == pos_1:
            pos_2 = int(self.rng.integers(size))
        return pos_1, pos_2

    def _event(self, kind: str, size: int):
        rng = self.rng
        if kind in ("duplication", "deletion", "translocation", "inversion") and size == 1:
            return None
        if kind == "duplication":
            pos_1, pos_2 = self._pair(size)
            pos_3 = int(rng.integers(size))
            seglen = pos_2 - pos_1 if pos_1 < pos_2 else size - pos_1 + pos_2
            if size + seglen > self.max_length:
                return None
            return ("duplication", pos_1, pos_2, pos_3)
        if kind == "deletion":
            pos_1, pos_2 = self._pair(size)
            after = size - (pos_2 - pos_1) if pos_1 < pos_2 else pos_1 - pos_2
            if after < self.min_length:
                return None
            return ("deletion", pos_1, pos_2)
        if kind == "translocation":
            pos_1, pos_2 = sorted(self._pair(size))
            seglen = pos_2 - pos_1
            pos_3 = pos_1 + int(rng.integers(seglen))
            pos_4 = int(rng.integers(size - seglen))
            if pos_4 >= pos_1:
                pos_4 += seglen
            return ("translocation", pos_1, pos_2, pos_3, pos_4, bool(rng.integers(2) == 0))
        if kind == "inversion":
            pos_1, pos_2 = sorted(self._pair(size))
            return ("inversion", pos_1, pos_2)
        if kind == "switch":
            return ("switch", int(rng.integers(size)))
        nb = 1 if self.max_indel_size == 1 else 1 + int(rng.integers(self.max_indel_size))
        if kind == "small_insertion":
            pos = int(rng.integers(size))
            if size + nb > self.max_length:
                return None
            return ("small_insertion", pos, rng.integers(0, 2, nb, dtype=np.uint8))
        pos = int(rng.integers(size))
        if size - nb < self.min_length:
            return None
        return ("small_deletion", pos, nb)

    @staticmethod
    def _size_after(event, size: int) -> int:
        kind = event[0]
        if kind == "duplication":
            _, pos_1, pos_2, _ = event
            return size + (pos_2 - pos_1 if pos_1 < pos_2 else size - pos_1 + pos_2)
        if kind == "deletion":
            _, pos_1, pos_2 = event
            return size - (pos_2 - pos_1) if pos_1 < pos_2 else pos_1 - pos_2
        if kind == "small_insertion":
            return size + event[2].size
        if kind == "small_deletion":
            return size - event[2]
        return size


def apply_mutations(g: np.ndarray, events) -> np.ndarray:
    for event in events:
        kind, args = event[0], event[1:]
        if kind == "switch":
            g = do_switch(g.copy(), *args)
        elif kind == "small_insertion":
            g = do_small_insertion(g, *args)
        elif kind == "small_deletion":
            g = do_small_deletion(g, *args)
        elif kind == "duplication":
            g = do_duplication(g, *args)
        elif kind == "deletion":
            g = do_deletion(g, *args)
        elif kind == "inversion":
            g = do_inversion(g, *args)
        else:
            g = do_translocation(g, *args)
    return g


# ---------------------------------------------------------------------------
# the individual and the population

class Individual:
    """A genome with everything the chemistry reads off it."""

    __slots__ = ("genome", "gid", "rnas", "proteins", "phenotype", "metabolic_error", "fitness")

    def __init__(self, genome: np.ndarray, w_max: float, sampling: int,
                 target: np.ndarray, selection_pressure: float):
        self.genome = genome
        self.gid = genotype_id(genome)
        self.rnas = transcribe(genome)
        self.proteins = translate(genome, self.rnas, w_max)
        self.phenotype = phenotype(self.proteins, sampling)
        self.metabolic_error = geometric_area(self.phenotype - target)
        self.fitness = (0.0 if geometric_area(self.phenotype) == 0.0
                        else math.exp(-selection_pressure * self.metabolic_error))

    @property
    def functional(self) -> list[Protein]:
        return [p for p in self.proteins if p.functional]


class World:
    """The grid of organisms, its target and its evolutionary loop."""

    def __init__(self, rng, *, width=4, height=4, w_max=0.033333333, selection_pressure=1000.0,
                 env_sampling=300, target=None, point_mutation_rate=1e-5,
                 small_insertion_rate=1e-5, small_deletion_rate=1e-5, duplication_rate=1e-5,
                 deletion_rate=1e-5, translocation_rate=1e-5, inversion_rate=1e-5,
                 max_indel_size=6, min_genome_length=1, max_genome_length=10_000_000,
                 selection_patch_size=3):
        self.rng = rng
        self.width, self.height = width, height
        self.size = width * height
        self.w_max = w_max
        self.selection_pressure = selection_pressure
        self.sampling = env_sampling
        self.gaussians = DEFAULT_TARGET if target is None else target
        self.target = target_function(self.gaussians, env_sampling)
        self.target_area = geometric_area(self.target)
        self.patch = selection_patch_size
        self.mutator = Mutator(
            rng,
            {"switch": point_mutation_rate, "small_insertion": small_insertion_rate,
             "small_deletion": small_deletion_rate, "duplication": duplication_rate,
             "deletion": deletion_rate, "translocation": translocation_rate,
             "inversion": inversion_rate},
            max_indel_size, (min_genome_length, max_genome_length))
        self.known: dict[str, Individual] = {}
        self.cells: list[str] = []
        self.parent: dict[str, str | None] = {}
        self.born: Counter = Counter()
        self.events: dict[tuple, list] = {}
        self.mutations: Counter = Counter()
        self.generation = 0

    # --- evaluation -------------------------------------------------------
    def evaluate(self, genome: np.ndarray) -> Individual:
        gid = genotype_id(genome)
        indiv = self.known.get(gid)
        if indiv is None:
            indiv = Individual(genome, self.w_max, self.sampling, self.target, self.selection_pressure)
            self.known[gid] = indiv
        return indiv

    def random_individual(self, length: int, attempts: int = 100000) -> Individual:
        """A random genome that does better than a flat phenotype (Individual::make_random)."""
        for _ in range(attempts):
            indiv = self.evaluate(self.rng.integers(0, 2, length, dtype=np.uint8))
            if round((indiv.metabolic_error - self.target_area) * 1e6) / 1e6 < 0.0:
                return indiv
        raise ValueError(
            f"no random genome of {length} bases beat a flat phenotype in {attempts} draws; "
            "use a longer genome_length or an easier target")

    def seed(self, indiv: Individual) -> None:
        """Fill every cell with a clone (ExperimentCreator::create_random_clonal_population)."""
        self.cells = [indiv.gid] * self.size
        self.parent.setdefault(indiv.gid, None)
        self.born[indiv.gid] += self.size

    # --- selection --------------------------------------------------------
    def neighbourhood(self, cell: int) -> list[int]:
        """The selection patch around a cell of the torus (Grid::apply_offset)."""
        x, y = cell // self.height, cell % self.height
        half = self.patch // 2
        return [((x + dx) % self.width) * self.height + ((y + dy) % self.height)
                for dx in range(-half, half + 1) for dy in range(-half, half + 1)]

    def select(self, cell: int) -> int:
        patch = self.neighbourhood(cell)
        fit = np.array([self.known[self.cells[c]].fitness for c in patch])
        total = fit.sum()
        probs = np.full(fit.size, 1.0 / fit.size) if total == 0.0 else fit / total
        return patch[int(self.rng.choice(fit.size, p=probs))]

    # --- one generation ---------------------------------------------------
    def _event(self, reactants: list[str], products: list[str]) -> None:
        lhs, rhs = Counter(reactants), Counter(products)
        key = (frozenset(lhs.items()), frozenset(rhs.items()))
        entry = self.events.get(key)
        if entry is None:
            self.events[key] = entry = [dict(lhs), dict(rhs), 0]
        entry[2] += 1

    def step(self) -> None:
        self.generation += 1
        offspring: list[str] = []
        for cell in range(self.size):
            source = self.select(cell)
            parent = self.known[self.cells[source]]
            events = self.mutator.draw(parent.genome.size)
            if events:
                child = self.evaluate(apply_mutations(parent.genome, events))
                for event in events:
                    self.mutations[event[0]] += 1
            else:
                child = parent
            self.parent.setdefault(child.gid, parent.gid)
            self.born[child.gid] += 1
            offspring.append(child.gid)
            if source == cell:
                self._event([parent.gid], [child.gid])
            else:
                self._event([parent.gid, self.cells[cell]], [parent.gid, child.gid])
        self.cells = offspring

    def sample(self) -> dict:
        indivs = [self.known[g] for g in self.cells]
        best = max(indivs, key=lambda i: i.fitness)
        return {
            "generation": self.generation,
            "best_fitness": best.fitness,
            "best_metabolic_error": best.metabolic_error,
            "mean_fitness": float(np.mean([i.fitness for i in indivs])),
            "mean_metabolic_error": float(np.mean([i.metabolic_error for i in indivs])),
            "mean_genome_length": float(np.mean([i.genome.size for i in indivs])),
            "mean_functional_proteins": float(np.mean([len(i.functional) for i in indivs])),
            "genotypes": len(set(self.cells)),
        }

    def run(self, generations: int) -> list[dict]:
        history = [self.sample()]
        for _ in range(generations):
            self.step()
            history.append(self.sample())
        return history

    def best(self) -> Individual:
        return max((self.known[g] for g in self.cells), key=lambda i: i.fitness)


# ---------------------------------------------------------------------------

def generate(p, rng):
    if p.selection_patch_size % 2 == 0:
        raise ValueError(f"selection_patch_size must be odd, got {p.selection_patch_size}")
    if p.min_genome_length > p.max_genome_length:
        raise ValueError(
            f"min_genome_length {p.min_genome_length} exceeds max_genome_length {p.max_genome_length}")
    if not p.min_genome_length <= p.genome_length <= p.max_genome_length:
        raise ValueError(
            f"genome_length {p.genome_length} is outside "
            f"[{p.min_genome_length}, {p.max_genome_length}]")
    gaussians = check_target(p.target)

    world = World(
        rng, width=p.grid_width, height=p.grid_height, w_max=p.w_max,
        selection_pressure=p.selection_pressure, env_sampling=p.env_sampling, target=gaussians,
        point_mutation_rate=p.point_mutation_rate, small_insertion_rate=p.small_insertion_rate,
        small_deletion_rate=p.small_deletion_rate, duplication_rate=p.duplication_rate,
        deletion_rate=p.deletion_rate, translocation_rate=p.translocation_rate,
        inversion_rate=p.inversion_rate, max_indel_size=p.max_indel_size,
        min_genome_length=p.min_genome_length, max_genome_length=p.max_genome_length,
        selection_patch_size=p.selection_patch_size)
    ancestor = world.random_individual(p.genome_length)
    world.seed(ancestor)
    history = world.run(p.generations)

    # species: the genotypes that lived, and the proteins they fold into
    final = Counter(world.cells)
    lineage = [g for g in world.known if world.born[g]]
    species, reactions, proteome = [], [], {}
    for gid in lineage:
        indiv = world.known[gid]
        species.append(Species(gid, structure=text(indiv.genome)))
        proteome[gid] = Counter(prot.id for prot in indiv.proteins)
    proteins: dict[str, Protein] = {}
    for gid in lineage:
        for prot in world.known[gid].proteins:
            proteins.setdefault(prot.id, prot)
    species += [Species(pid, structure=prot.structure()) for pid, prot in proteins.items()]

    for gid in lineage:                                   # expression: genome -> genome + proteome
        if proteome[gid]:
            products = dict(proteome[gid])
            products[gid] = products.get(gid, 0) + 1
            reactions.append(Reaction({gid: 1}, products, count=int(world.born[gid])))
    for lhs, rhs, n in world.events.values():             # replication and replacement
        reactions.append(Reaction(lhs, rhs, count=n))

    best = world.best()
    return Network(
        species=species,
        reactions=reactions,
        status="observed",
        initial_state={ancestor.gid: float(world.size)},
        extras={
            "space": {
                "lattice": "torus", "width": p.grid_width, "height": p.grid_height,
                "neighbourhood": f"{p.selection_patch_size}x{p.selection_patch_size} selection patch",
                "cell_index": "x * height + y",
                "final_grid": list(world.cells),
            },
            "final_state": {g: float(n) for g, n in sorted(final.items())},
            "genotype_parent": {g: world.parent[g] for g in lineage},
            "proteome": {g: dict(sorted(proteome[g].items())) for g in lineage if proteome[g]},
            "proteins": {pid: {"m": prot.m, "w": prot.w, "h": prot.h,
                               "codons": len(prot.codons), "functional": prot.functional}
                         for pid, prot in proteins.items()},
            "events": ("P + V -> P + O: the organism P of a neighbouring cell replicates and its "
                       "offspring O replaces the cell's previous occupant V; P -> O: the cell's own "
                       "occupant replicates in place; G -> G + proteins: expression of a genotype"),
            "analysis": {
                "generations": p.generations,
                "population": world.size,
                "genotypes_seen": len(lineage),
                "mutations": {k: int(world.mutations[k]) for k in MUTATION_TYPES},
                "phenotypic_target": {
                    "gaussians": gaussians,
                    "sampling": p.env_sampling,
                    "area": world.target_area,
                    "points": [float(v) for v in world.target],
                },
                "best": {
                    "genotype": best.gid,
                    "fitness": best.fitness,
                    "metabolic_error": best.metabolic_error,
                    "genome_length": int(best.genome.size),
                    "rnas": len(best.rnas),
                    "proteins": len(best.proteins),
                    "functional_proteins": len(best.functional),
                    "phenotype": [float(v) for v in best.phenotype],
                },
                "ancestor": {
                    "genotype": ancestor.gid,
                    "fitness": ancestor.fitness,
                    "metabolic_error": ancestor.metabolic_error,
                    "functional_proteins": len(ancestor.functional),
                },
                "per_generation": history,
            },
        },
    )
