"""Publish Chemart's built-in catalog as the official ``chemart/<id>`` repos.

Each repo is *builtin-backed*: its chemart.yaml is the catalog entry plus
``hub.builtin: <id>``, and it carries no code, because the generator ships
with the chemart library. Loading ``chemart/brusselator`` therefore runs the
installed ``brusselator`` generator and needs no trust_remote_code.

Seeding goes through the public client, like any other push, so it exercises
the real upload path. Re-running it only commits repos whose files changed.
"""

from __future__ import annotations

import json
import os
import sys

import yaml


class _NoAliases(yaml.SafeDumper):
    def ignore_aliases(self, data):  # the hub refuses YAML anchors and aliases
        return True


def official_files(chemistry_id: str) -> dict[str, bytes]:
    from chemart import api, catalog

    entry = next(c for c in catalog.load() if c.id == chemistry_id)
    raw = yaml.safe_load((catalog.CATALOG_DIR / entry.source_file).read_text())["chemistries"][0]
    hub = {"repo_type": "generator", "builtin": entry.id, "tags": ["banzhaf-yamamoto"]}
    chemart_yaml = yaml.dump({"hub": hub, "chemistries": [raw]}, Dumper=_NoAliases,
                             sort_keys=False, allow_unicode=True, width=100)

    preview = api.generate_network(entry.id, seed=1)
    # The page already shows the entry's intuition and (S, R, A); the card
    # says where the chemistry comes from and how to load it.
    lines = [f"# {entry.name}", ""]
    lines += [
        "One of the artificial chemistries of Banzhaf & Yamamoto, "
        f"*Artificial Chemistries* (MIT Press, 2015), §{entry.book}"
        + (f"; first proposed by {entry.origin}." if entry.origin else "."),
        "",
        "Its generator ships with the `chemart` library, so loading it runs no code "
        "from the hub and needs no `trust_remote_code`.",
        "",
        "```python",
        "import chemart",
        "",
        f'net = chemart.generate_network("chemart/{entry.id}", seed=0)'
        f'   # same as "{entry.id}"',
        "print(net.summary())",
        "```",
        "",
    ]
    return {
        "chemart.yaml": chemart_yaml.encode(),
        "README.md": "\n".join(lines).encode(),
        "preview.json": json.dumps(preview.to_dict(), ensure_ascii=False).encode(),
    }


def seed(url: str | None = None, token: str | None = None, only: list[str] | None = None,
         namespace: str = "chemart") -> int:
    if url:
        os.environ["CHEMART_HUB_URL"] = url
    if token:
        os.environ["CHEMART_HUB_TOKEN"] = token
    from chemart import __version__, catalog, hub

    me = hub.whoami()
    if namespace not in me["orgs"] and me["name"] != namespace and not me.get("is_admin"):
        print(f"error: {me['name']} is not a member of the {namespace!r} organisation "
              f"(chemart-hub create-org {namespace} --owner {me['name']})", file=sys.stderr)
        return 2
    ids = [c.id for c in catalog.load() if c.implemented]
    if only:
        unknown = sorted(set(only) - set(ids))
        if unknown:
            print(f"error: not built-in chemistries: {', '.join(unknown)}", file=sys.stderr)
            return 2
        ids = [i for i in ids if i in only]
    changed = 0
    for i, cid in enumerate(ids, 1):
        result = hub.upload_files(f"{namespace}/{cid}", official_files(cid), repo_type="generator",
                                  message=f"chemart {__version__}: {cid}")
        state = "unchanged" if result.get("unchanged") else f"-> {result['commit'][:12]}"
        changed += not result.get("unchanged")
        print(f"[{i:3d}/{len(ids)}] {namespace}/{cid:36s} {state}")
    print(f"{changed} of {len(ids)} official repos updated")
    return 0
