"""`chemart-hub`: run and administer a Chemart Hub.

    chemart-hub init                         create the data folder and database
    chemart-hub serve [--host --port]        run the site and API
    chemart-hub pit [--port]                 the simulation pit, on this machine only
    chemart-hub create-user NAME [--admin]   prompts for a password
    chemart-hub set-admin NAME [--off]       grant (or revoke) superadmin
    chemart-hub set-password NAME            prompts for the new password
    chemart-hub create-org NAME --owner USER
    chemart-hub add-member ORG USER
    chemart-hub token USER [--name --scope]  print a new API token
    chemart-hub seed [--only a,b]            publish the built-in catalog as chemart/*
    chemart-hub gc                           delete stored files no commit uses (orphans
                                             of interrupted uploads); the archive is untouched

Every command takes --data-dir (default $CHEMART_HUB_DATA or ./hub-data).
"""

from __future__ import annotations

import argparse
import getpass
import os
import sys

from chemart_hub import db, service
from chemart_hub.config import Settings
from chemart_hub.errors import ApiError


def _settings(args) -> Settings:
    return Settings.from_env(data_dir=args.data_dir, public_url=getattr(args, "public_url", None))


def _conn(settings: Settings):
    db.migrate(settings.db_path)
    return db.connect(settings.db_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="chemart-hub", description="Run and administer a Chemart Hub.")
    parser.add_argument("--data-dir", help="data folder (default $CHEMART_HUB_DATA or ./hub-data)")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("init", help="create the data folder and database")

    serve = sub.add_parser("serve", help="run the site and API")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8000)
    serve.add_argument("--public-url", help="URL the hub is reached at (default http://HOST:PORT)")

    pit = sub.add_parser("pit", help="run the simulation pit on this machine (localhost only)")
    pit.add_argument("--port", type=int, default=8765)
    pit.add_argument("--no-browser", action="store_true", help="do not open a browser tab")

    user = sub.add_parser("create-user", help="create a user (prompts for the password)")
    user.add_argument("name")
    user.add_argument("--email")
    user.add_argument("--admin", action="store_true")
    user.add_argument("--password", help="non-interactive; prefer the prompt")

    set_admin = sub.add_parser("set-admin", help="make a user a superadmin (full power on the site)")
    set_admin.add_argument("name")
    set_admin.add_argument("--off", action="store_true", help="revoke superadmin instead")

    set_password = sub.add_parser("set-password", help="set a user's password (prompts)")
    set_password.add_argument("name")
    set_password.add_argument("--password", help="non-interactive; prefer the prompt")

    org = sub.add_parser("create-org", help="create an organisation")
    org.add_argument("name")
    org.add_argument("--owner", required=True)

    member = sub.add_parser("add-member", help="add a user to an organisation")
    member.add_argument("org")
    member.add_argument("user")

    token = sub.add_parser("token", help="print a new API token for a user")
    token.add_argument("user")
    token.add_argument("--name", default="cli")
    token.add_argument("--scope", choices=("read", "write"), default="write")

    seed = sub.add_parser("seed", help="publish the built-in catalog as the official chemart/ repos")
    seed.add_argument("--url", help="hub URL (default $CHEMART_HUB_URL)")
    seed.add_argument("--token", help="write token of a chemart org member (default: saved login)")
    seed.add_argument("--only", help="comma-separated catalog ids")

    sub.add_parser("gc", help="delete orphaned stored files and compact the database")

    build = sub.add_parser("build-static", help="build the hub as a static site from a registry repository")
    build.add_argument("registry", help="a clone of the registry (a git repository with repos/<ns>/<name>/)")
    build.add_argument("out", help="output folder (replaced)")
    build.add_argument("--site-url", required=True, help="where the site will live, e.g. https://you.github.io/chemart-hub")
    build.add_argument("--registry-repo", required=True, help="the registry on GitHub, as owner/repo")
    build.add_argument("--branch", default="main", help="the registry's default branch")

    init_reg = sub.add_parser("init-registry", help="start a registry repository for a hub on GitHub Pages")
    init_reg.add_argument("folder")
    init_reg.add_argument("--owner", required=True, help="your GitHub login: the first maintainer")
    init_reg.add_argument("--site-url", required=True, help="where the site will live, e.g. https://you.github.io/chemart-hub")
    init_reg.add_argument("--chemart-repo", default="MarcoTuc/Chemart", help="where the workflows get Chemart from")
    init_reg.add_argument("--chemart-ref", default="chemart-library")
    init_reg.add_argument("--branch", default="main", help="the registry's default branch")
    init_reg.add_argument("--no-official", action="store_true", help="do not stock the official shelf")

    sync = sub.add_parser("sync-official", help="write repos/chemart/ in a registry from the installed catalog")
    sync.add_argument("registry")

    check = sub.add_parser("validate-pr", help="check a registry pull request (runs the contributed code)")
    check.add_argument("registry")
    check.add_argument("--base", required=True, help="the base commit or ref")
    check.add_argument("--head", default="HEAD", help="the pull request's head")
    check.add_argument("--author", required=True, help="the pull request author's GitHub login")
    check.add_argument("--no-contract", action="store_true", help="do not import and run generator code")

    args = parser.parse_args(argv)
    if args.command in ("build-static", "init-registry", "sync-official", "validate-pr"):
        return _registry_command(args)
    settings = _settings(args)
    try:
        if args.command == "init":
            _conn(settings).close()
            print(f"initialised {settings.data_dir}")
        elif args.command == "pit":
            import threading
            import webbrowser

            import uvicorn

            from chemart_hub.pit import create_pit_app

            url = f"http://127.0.0.1:{args.port}/"
            print(f"the pit is at {url} (Ctrl+C to stop)")
            if not args.no_browser:
                threading.Timer(1.0, webbrowser.open, [url]).start()
            uvicorn.run(create_pit_app(), host="127.0.0.1", port=args.port, log_level="warning")
        elif args.command == "serve":
            import uvicorn

            from chemart_hub.app import create_app

            if not getattr(args, "public_url", None) and "CHEMART_HUB_PUBLIC_URL" not in os.environ:
                settings.public_url = f"http://{args.host}:{args.port}"
            uvicorn.run(create_app(settings), host=args.host, port=args.port)
        elif args.command == "create-user":
            password = args.password or getpass.getpass(f"password for {args.name}: ")
            with _conn(settings) as c:
                row = service.create_user(c, args.name, password, args.email, admin=args.admin)
            print(f"created user {row['name']}{' (admin)' if args.admin else ''}")
        elif args.command == "set-admin":
            with _conn(settings) as c:
                row = service.account(c, args.name)
                if row is None or row["is_org"]:
                    raise ApiError(404, f"user {args.name!r} not found")
                if args.off and c.execute(
                        "SELECT COUNT(*) FROM accounts WHERE is_admin = 1 AND is_org = 0 AND id != ?",
                        (row["id"],)).fetchone()[0] == 0:
                    raise ApiError(409, "this is the last superadmin")
                c.execute("UPDATE accounts SET is_admin = ?, disabled = 0 WHERE id = ?", (int(not args.off), row["id"]))
            print(f"{row['name']} is {'no longer ' if args.off else ''}a superadmin")
        elif args.command == "set-password":
            from chemart_hub import auth

            password = args.password or getpass.getpass(f"new password for {args.name}: ")
            auth.check_password_strength(password)
            with _conn(settings) as c:
                row = service.account(c, args.name)
                if row is None or row["is_org"]:
                    raise ApiError(404, f"user {args.name!r} not found")
                c.execute("UPDATE accounts SET password_hash = ? WHERE id = ?", (auth.hash_password(password), row["id"]))
                c.execute("DELETE FROM sessions WHERE user_id = ?", (row["id"],))
            print(f"password set for {row['name']}")
        elif args.command == "create-org":
            with _conn(settings) as c:
                row = service.create_org(c, args.name, args.owner, allow_reserved=True)
            print(f"created organisation {row['name']} owned by {args.owner}")
        elif args.command == "add-member":
            with _conn(settings) as c:
                service.add_member(c, args.org, args.user)
            print(f"added {args.user} to {args.org}")
        elif args.command == "token":
            from chemart_hub import auth

            with _conn(settings) as c:
                row = service.account(c, args.user)
                if row is None or row["is_org"]:
                    raise ApiError(404, f"user {args.user!r} not found")
                print(auth.create_token(c, row["id"], args.name, args.scope))
        elif args.command == "gc":
            from chemart_hub import curation
            from chemart_hub.storage import BlobStore

            with _conn(settings) as c:
                result = curation.collect_garbage(c, BlobStore(settings.blobs_dir))
            print(f"removed {result['files_removed']} orphaned files, freed {result['bytes_freed']} bytes")
        elif args.command == "seed":
            from chemart_hub.seed import seed as run_seed

            only = [s for s in (args.only or "").split(",") if s] or None
            return run_seed(url=args.url, token=args.token, only=only)
    except ApiError as err:
        print(f"error: {err.message}", file=sys.stderr)
        for p in err.problems:
            print(f"  - {p}", file=sys.stderr)
        return 2
    return 0


def _registry_command(args) -> int:
    """The commands that work on a registry repository; they need no hub data folder."""
    from chemart_hub import static

    if args.command == "build-static":
        report = static.build(args.registry, args.out, site_url=args.site_url,
                              registry_repo=args.registry_repo, branch=args.branch)
        print(f"built {len(report.repos)} repos, {report.commits} commits, {report.pages} pages, "
              f"{report.files} raw files into {args.out}")
        for p in report.problems:
            print(f"  - {p}", file=sys.stderr)
        return 0
    if args.command == "init-registry":
        written = static.init_registry(args.folder, owner=args.owner, site_url=args.site_url,
                                       chemart_repo=args.chemart_repo, chemart_ref=args.chemart_ref,
                                       branch=args.branch, official=not args.no_official)
        print(f"wrote {len(written)} files" + ("" if args.no_official else " and the official shelf")
              + f" into {args.folder}")
        print("next: commit it, push it to GitHub, and in Settings > Pages set Source to 'GitHub Actions'")
        return 0
    if args.command == "sync-official":
        changed = static.sync_official(args.registry)
        print(f"{len(changed)} official repo(s) changed" + (f": {', '.join(changed)}" if changed else ""))
        return 0
    problems = static.validate_pr(args.registry, args.base, args.head, args.author,
                                  run_contract=not args.no_contract)
    for p in problems:
        print(f"  - {p}")
    print("ok: the pull request can be merged" if not problems else f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
