#!/usr/bin/env python3
"""Survey the Chemart catalog: filter by capability, family, kind, fidelity.

    uv run python .claude/skills/chemart/scripts/survey.py --provides rate-constants
    uv run python .claude/skills/chemart/scripts/survey.py --provides energies mass-conservation
    uv run python .claude/skills/chemart/scripts/survey.py --family origin-of-life --verbose
    uv run python .claude/skills/chemart/scripts/survey.py --constructive --generate --limit 10

--provides takes several tags and requires all of them, which is the question
the plain CLI cannot answer.

--generate actually builds each network with default parameters and reports
real sizes and computed capabilities instead of catalog claims. Slower
(roughly a second per chemistry) but it is the honest version: the catalog's
`provides` is a claim about the chemistry, while a generated network tells you
what the defaults really hand you.
"""

from __future__ import annotations

import argparse
import sys

from chemart.catalog import PROVIDES, active, load

FIDELITIES = ("book", "book+decisions", "reconstructed")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--provides", nargs="+", metavar="TAG",
                    help=f"require all of these; one of {', '.join(sorted(PROVIDES))}")
    ap.add_argument("--family")
    ap.add_argument("--kind", choices=("generator", "formalism", "framework", "analysis", "wet"))
    ap.add_argument("--fidelity", choices=FIDELITIES)
    ap.add_argument("--constructive", action="store_true")
    ap.add_argument("--fixed", action="store_true", help="only non-constructive entries")
    ap.add_argument("--generate", action="store_true",
                    help="build each network and report real sizes")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--verbose", action="store_true", help="show the one-line summary too")
    args = ap.parse_args(argv)

    if args.provides:
        unknown = sorted(set(args.provides) - PROVIDES)
        if unknown:
            print(f"error: unknown capability {unknown}. Valid: {', '.join(sorted(PROVIDES))}",
                  file=sys.stderr)
            return 2

    rows = []
    for c in active(load()):
        if args.provides and not set(args.provides) <= set(c.provides):
            continue
        if args.family and c.family != args.family:
            continue
        if args.kind and c.kind != args.kind:
            continue
        if args.fidelity and c.fidelity != args.fidelity:
            continue
        if args.constructive and not c.constructive:
            continue
        if args.fixed and c.constructive:
            continue
        rows.append(c)

    if args.limit:
        rows = rows[: args.limit]

    if not rows:
        print("no chemistries match; relax a filter (try --provides with fewer tags)")
        return 0

    if not args.generate:
        width = max(len(c.id) for c in rows)
        print(f"{len(rows)} chemistries\n")
        print(f"{'id':<{width}}  {'family':<22} {'kind':<10} {'fidelity':<15} constructive")
        print("-" * (width + 62))
        for c in rows:
            print(f"{c.id:<{width}}  {c.family:<22} {c.kind:<10} "
                  f"{str(c.fidelity):<15} {'yes' if c.constructive else 'no'}")
            if args.verbose:
                print(f"{'':<{width}}  {', '.join(c.provides)}")
                if c.phenomena:
                    print(f"{'':<{width}}  phenomena: {'; '.join(c.phenomena[:2])}")
        return 0

    import chemart

    width = max(len(c.id) for c in rows)
    print(f"{len(rows)} chemistries, generating with defaults (seed={args.seed})\n")
    print(f"{'id':<{width}}  {'species':>8} {'reactions':>10}  {'status':<10} computed provides")
    print("-" * (width + 70))
    for c in rows:
        try:
            net = chemart.generate_network(c.id, seed=args.seed)
        except Exception as err:                                  # noqa: BLE001
            print(f"{c.id:<{width}}  {'—':>8} {'—':>10}  {'error':<10} "
                  f"{type(err).__name__}: {err}")
            continue
        print(f"{c.id:<{width}}  {len(net.species):>8} {len(net.reactions):>10}  "
              f"{net.status:<10} {', '.join(net.provides)}")
        if args.verbose:
            claimed = set(c.provides) - set(net.provides)
            if claimed:
                print(f"{'':<{width}}  claimed but not in this network: {', '.join(sorted(claimed))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
