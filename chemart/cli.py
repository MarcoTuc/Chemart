"""Command line interface over `chemart.api` and `chemart.hub`.

    chemart list                 # the chemistry catalog (--all adds the archive)
    chemart describe matrix-chemistry
    chemart generate matrix-chemistry -p N=4 --seed 0 --format json
    chemart simulate brusselator --t-end 40                      # rate equations
    chemart simulate brusselator --method ssa --volume 100 --seed 1 --format csv
    chemart simulate kauffman-autocatalytic-sets --x0 1 --rates '{"dist": "lognormal", "mean": 0, "sigma": 1}'

Chemart Hub:

    chemart login                       paste a token from <hub>/settings/tokens
    chemart whoami | logout
    chemart search "autocatalysis" --provides rate-constants
    chemart describe alice/my-chem
    chemart generate alice/my-chem --trust-remote-code --revision 3f2a9c1
    chemart download bob/ecoli-core     prints the local folder
    chemart new my-chem                 a generator skeleton to edit
    chemart check my-chem               the checks `push` runs, without pushing
    chemart push my-chem [alice/my-chem]
    chemart push-network net.json alice/my-network

Parameter values are parsed as JSON when possible (``-p N=4``,
``-p food_set='["a","b"]'``) and taken as plain strings otherwise.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _value(text: str):
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def _print_json(data) -> None:
    print(json.dumps(data, indent=2, ensure_ascii=False))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="chemart", description="A mart of artificial chemistries.")
    sub = parser.add_subparsers(dest="command", required=True)
    listing = sub.add_parser("list", help="list every built-in chemistry (JSON)")
    listing.add_argument("--all", action="store_true", help="also list archived entries")

    describe = sub.add_parser("describe", help="describe one chemistry (JSON)")
    describe.add_argument("chemistry", help="catalog id, or namespace/name on the hub")
    describe.add_argument("--revision")

    generate = sub.add_parser("generate", help="generate a reaction network")
    generate.add_argument("chemistry", help="catalog id, or namespace/name on the hub")
    generate.add_argument("-p", "--param", action="append", default=[], metavar="NAME=VALUE")
    generate.add_argument("--seed", type=int)
    generate.add_argument("--format", choices=("summary", "text", "json"), default="summary")
    generate.add_argument("--revision", help="hub only: 'main' or a commit id")
    generate.add_argument("--trust-remote-code", action="store_true",
                          help="hub only: allow the repo's own Python code to run on this machine")

    simulate = sub.add_parser("simulate", help="simulate a chemistry's network (ODE or SSA)")
    simulate.add_argument("chemistry", help="catalog id, or namespace/name on the hub")
    simulate.add_argument("-p", "--param", action="append", default=[], metavar="NAME=VALUE",
                          help="chemistry parameter")
    simulate.add_argument("--seed", type=int, help="seeds the chemistry, the rate draws and the SSA")
    simulate.add_argument("--method", choices=("ode", "ssa"), default="ode")
    simulate.add_argument("--t-end", type=float, default=40.0)
    simulate.add_argument("--points", type=int, default=200)
    simulate.add_argument("--volume", type=float, default=1.0, help="ssa only: counts = amount x volume")
    simulate.add_argument("--rates", metavar="SPEC",
                          help="a number, a JSON distribution/table, or a .json/.csv file")
    simulate.add_argument("--x0", metavar="SPEC", help="initial state, as for --rates (by species)")
    simulate.add_argument("--fill-only", action="store_true", help="only rate reactions without a rate")
    simulate.add_argument("--species", nargs="*", help="columns to print (default: all)")
    simulate.add_argument("--format", choices=("table", "csv", "json"), default="table")
    simulate.add_argument("--revision", help="hub only: 'main' or a commit id")
    simulate.add_argument("--trust-remote-code", action="store_true",
                          help="hub only: allow the repo's own Python code to run on this machine")

    login = sub.add_parser("login", help="save an API token for the hub")
    login.add_argument("--token", help="default: prompt")
    sub.add_parser("logout", help="forget the saved token")
    sub.add_parser("whoami", help="show the logged-in hub user")

    search = sub.add_parser("search", help="search the hub")
    search.add_argument("query", nargs="?")
    search.add_argument("--type", dest="repo_type", choices=("generator", "network"))
    search.add_argument("--family")
    search.add_argument("--provides", action="append", default=[])
    search.add_argument("--tag")
    search.add_argument("--author")
    search.add_argument("--sort", default="trending")
    search.add_argument("--limit", type=int, default=30)
    search.add_argument("--json", action="store_true")

    download = sub.add_parser("download", help="download a hub repo; prints the local folder")
    download.add_argument("repo")
    download.add_argument("--revision")

    new = sub.add_parser("new", help="write a generator skeleton to edit and push")
    new.add_argument("folder")
    new.add_argument("--id", help="chemistry id (default: the folder name)")
    new.add_argument("--name")

    check = sub.add_parser("check", help="run the checks `push` runs, without pushing")
    check.add_argument("folder")

    push = sub.add_parser("push", help="share a generator folder on the hub")
    push.add_argument("folder")
    push.add_argument("repo", nargs="?", help="namespace/name (default: <you>/<entry id>)")
    push.add_argument("-m", "--message")

    push_net = sub.add_parser("push-network", help="share a network JSON file on the hub")
    push_net.add_argument("file")
    push_net.add_argument("repo")
    push_net.add_argument("-m", "--message")
    push_net.add_argument("--title")
    push_net.add_argument("--description")
    push_net.add_argument("--license")
    push_net.add_argument("--tag", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)

    from chemart import api
    from chemart.hub._http import HubError

    try:
        return _run(parser, args, api) or 0
    except HubError as err:
        print(f"error: {err}", file=sys.stderr)
        return 3
    except (ValueError, NotImplementedError, ImportError, FileExistsError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


def _run(parser, args, api) -> int | None:
    if args.command == "list":
        _print_json(api.list_chemistries(args.all))
    elif args.command == "describe":
        _print_json(api.describe_chemistry(args.chemistry, args.revision))
    elif args.command in ("generate", "simulate"):
        params = {}
        for item in args.param:
            name, sep, raw = item.partition("=")
            if not sep:
                parser.error(f"--param expects NAME=VALUE, got {item!r}")
            params[name] = _value(raw)
        net = api.generate_network(args.chemistry, args.seed, revision=args.revision,
                                   trust_remote_code=args.trust_remote_code, **params)
        if args.command == "simulate":
            return _simulate(net, args)
        if args.format == "json":
            _print_json(net.to_dict())
        elif args.format == "text":
            print(net.to_text())
        else:
            print(net.summary())
    else:
        return _hub_command(args)
    return None


def _spec(text: str | None):
    """A --rates/--x0 value: a file path, JSON, or a plain string."""
    if text is None:
        return None
    if text.endswith((".json", ".csv")):
        return Path(text)
    return _value(text)


def _simulate(net, args) -> int:
    from chemart import simulate

    options = dict(rates=_spec(args.rates), x0=_spec(args.x0), points=args.points, seed=args.seed,
                   fill_only=args.fill_only)
    try:
        if args.method == "ode":
            traj = simulate.ode(net, args.t_end, **options)
        else:
            traj = simulate.ssa(net, args.t_end, volume=args.volume, **options)
    except simulate.NotSimulable as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    if args.format == "json":
        _print_json(traj.to_dict())
        return 0
    names, t, X = traj.array(args.species)
    if args.format == "csv":
        print("t," + ",".join(names))
        for k in range(len(t)):
            print(f"{t[k]}," + ",".join(repr(float(v)) for v in X[k]))
        return 0
    print(traj.network.summary())
    if traj.settings.get("stopped"):
        print(f"stopped early: {traj.settings['stopped']}")
    width = max([10, *map(len, names)])
    print(f"\n{'t':>10}  " + "  ".join(f"{n:>{width}}" for n in names))
    for k in range(0, len(t), max(1, len(t) // 10)):
        print(f"{t[k]:>10.4g}  " + "  ".join(f"{v:>{width}.4g}" for v in X[k]))
    return 0


def _hub_command(args) -> int | None:
    from chemart import hub
    from chemart.hub import _config

    if args.command == "login":
        name = hub.login(args.token)
        print(f"logged in to {_config.hub_url()} as {name}")
    elif args.command == "logout":
        hub.logout()
        print(f"logged out of {_config.hub_url()}")
    elif args.command == "whoami":
        me = hub.whoami()
        orgs = f" (orgs: {', '.join(me['orgs'])})" if me["orgs"] else ""
        print(f"{me['name']}{orgs} on {_config.hub_url()}")
    elif args.command == "search":
        repos = hub.search(args.query, repo_type=args.repo_type, family=args.family,
                           provides=args.provides, tag=args.tag, author=args.author,
                           sort=args.sort, limit=args.limit)
        if args.json:
            _print_json(repos)
        for r in [] if args.json else repos:
            code = "  [code]" if r["has_code"] else ""
            print(f"{r['id']:40s} {r['repo_type']:9s} ♥{r['likes']:<4d} ↓{r['downloads']:<6d} {r['title']}{code}")
    elif args.command == "download":
        print(hub.snapshot_download(args.repo, args.revision))
    elif args.command == "new":
        for path in hub.new(args.folder, args.id, args.name):
            print(f"wrote {path}")
        print(f"next: edit them, then `chemart check {args.folder}` and `chemart push {args.folder}`")
    elif args.command == "check":
        problems = hub.check(args.folder)
        for p in problems:
            print(f"  - {p}")
        print("ok: ready to push" if not problems else f"{len(problems)} problem(s)")
        return 1 if problems else 0
    elif args.command == "push":
        result = hub.push_generator(args.folder, args.repo, message=args.message)
        _report_push(result)
    elif args.command == "push-network":
        from chemart.network import Network

        net = Network.from_dict(json.loads(Path(args.file).read_text()))
        result = hub.push_network(net, args.repo, message=args.message, title=args.title,
                                  description=args.description, license=args.license,
                                  tags=args.tag or None)
        _report_push(result)
    return None


def _report_push(result: dict) -> None:
    if "commit" not in result:                    # a static hub: a pull request, or the steps for one
        from chemart.hub import _pr

        print(_pr.describe(result))
        return
    if result.get("unchanged"):
        print(f"no changes; {result['url']} is at {result['commit'][:12]}")
    else:
        print(f"pushed {result['commit'][:12]} -> {result['url']}")


if __name__ == "__main__":
    raise SystemExit(main())
