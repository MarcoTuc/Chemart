"""Command line interface over `chemart.api`.

    chemart list
    chemart describe matrix-chemistry
    chemart generate matrix-chemistry -p N=4 --seed 0 --format json

Parameter values are parsed as JSON when possible (``-p N=4``,
``-p food_set='["a","b"]'``) and taken as plain strings otherwise.
"""

from __future__ import annotations

import argparse
import json
import sys

from chemart import api


def _value(text: str):
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="chemart", description="A mart of artificial chemistries.")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="list every chemistry (JSON)")
    describe = sub.add_parser("describe", help="describe one chemistry (JSON)")
    describe.add_argument("chemistry")
    generate = sub.add_parser("generate", help="generate a reaction network")
    generate.add_argument("chemistry")
    generate.add_argument("-p", "--param", action="append", default=[], metavar="NAME=VALUE")
    generate.add_argument("--seed", type=int)
    generate.add_argument("--format", choices=("summary", "text", "json"), default="summary")
    args = parser.parse_args(argv)

    try:
        if args.command == "list":
            print(json.dumps(api.list_chemistries(), indent=2, ensure_ascii=False))
        elif args.command == "describe":
            print(json.dumps(api.describe_chemistry(args.chemistry), indent=2, ensure_ascii=False))
        else:
            params = {}
            for item in args.param:
                name, sep, raw = item.partition("=")
                if not sep:
                    parser.error(f"--param expects NAME=VALUE, got {item!r}")
                params[name] = _value(raw)
            net = api.generate_network(args.chemistry, args.seed, **params)
            if args.format == "json":
                print(json.dumps(net.to_dict(), indent=2, ensure_ascii=False))
            elif args.format == "text":
                print(net.to_text())
            else:
                print(net.summary())
    except (ValueError, NotImplementedError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
