"""Algorithmic Chemistry GP: a register-machine program read as a reaction multiset (book 16.6).

Catalog id: acgp. Each instruction r3 <- r1 op r2 is the reaction r1 + r2 -> r3
(eq. 16.9). The program is a random multiset of n_instructions instructions;
its operators are kept in extras.program, aligned with the reactions.
"""

from chemart.helpers.explicit import network, term

OPERATORS = ("+", "-", "*", "/")


def _registers(name, value, n):
    if not isinstance(value, list) or not all(isinstance(r, int) and not isinstance(r, bool) and 0 <= r < n for r in value):
        raise ValueError(f"{name} must be a list of register indices in 0..{n - 1}, got {value!r}")
    return list(dict.fromkeys(value))


def generate(p, rng):
    n = p.n_registers
    inputs = _registers("input_registers", p.input_registers, n)
    outputs = _registers("output_registers", p.output_registers, n)
    if set(inputs) & set(outputs):
        raise ValueError(f"input registers are read-only, so they cannot be outputs: {sorted(set(inputs) & set(outputs))}")
    writable = [r for r in range(n) if r not in inputs]
    if not writable:
        raise ValueError("every register is an input; at least one register must be writable")

    R = [f"r{i}" for i in range(n)]
    reactions, program = [], []
    for _ in range(p.n_instructions):
        a, b = (int(v) for v in rng.integers(0, n, 2))
        dst = writable[int(rng.integers(0, len(writable)))]
        op = OPERATORS[int(rng.integers(0, len(OPERATORS)))]
        lhs = term(2, R[a]) if a == b else f"{R[a]} + {R[b]}"
        reactions.append((f"{lhs} -> {R[dst]}", None))
        program.append(f"{R[dst]} = {R[a]} {op} {R[b]}")
    return network(
        reactions,
        species=R,
        extras={
            "program": program,
            "input_registers": [R[i] for i in inputs],
            "output_registers": [R[i] for i in outputs],
        },
    )
