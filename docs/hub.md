# Chemart Hub

The catalog holds the 98 chemistries of the book. The **Chemart Hub** is where
everything else goes: chemistries people invent or reconstruct, and reaction
networks they want to hand to someone else. It is a website you can browse,
and `chemart` talks to it directly, so a shared chemistry loads with the same
call as a built-in one:

```python
import chemart

chemart.generate_network("brusselator")                    # built in
chemart.generate_network("chemart/brusselator")            # the same, from the hub's official shelf
chemart.load_network("ada/brusselator-b35")                # a network someone shared
chemart.generate_network("ada/hypercycle-lite",            # a chemistry someone shared
                         revision="f2634dc6b1d4", trust_remote_code=True)
```

## Two kinds of repo

Everything on the hub lives in a **repo** named `namespace/name`, where the
namespace is a user or an organisation.

| | a **chemistry** (generator repo) | a **network** (network repo) |
|---|---|---|
| holds | a catalog entry and the Python that builds its network | one reaction network, as plain JSON |
| you load it with | `chemart.generate_network(id, ...)`, choosing your own parameters and seed | `chemart.load_network(id)` |
| runs code on your machine | yes: needs `trust_remote_code=True` | never |
| files | `chemart.yaml`, `generator.py` (+ helper `.py` files), `preview.json`, `README.md` | `network.json`, `README.md`, optional `chemart.yaml` |

The official `chemart/` organisation holds the 98 book chemistries. Those repos
are **built-in backed**: they contain the catalog entry but no code, because
the generator ships with the library. Loading `chemart/tierra` runs your
installed `tierra` generator and needs no trust flag.

## Loading

```python
chemart.describe_chemistry("ada/hypercycle-lite")   # entry, parameter schema, hub info; runs nothing
chemart.hub.search("autocatalysis", provides=["rate-constants"])
chemart.hub.repo_info("ada/hypercycle-lite")        # head commit, likes, downloads, facets
chemart.hub.snapshot_download("ada/brusselator-b35")  # every file, as a local folder
```

Every repo is versioned. Each push makes an immutable **commit**, and `main`
points at the latest one. `revision=` accepts `"main"`, a full commit id, or a
unique prefix of at least 7 characters; `"ada/x@3f2a9c1"` is shorthand for
`revision="3f2a9c1"`.

A network generated from the hub records exactly where it came from:

```python
net = chemart.generate_network("ada/hypercycle-lite", trust_remote_code=True, k=2.0)
net.chemistry   # 'ada/hypercycle-lite@f2634dc6b1d4…' (the full commit id)
chemart.generate_network(net.chemistry, seed=net.seed, trust_remote_code=True, **net.params) == net   # True
```

## Trusting code

A chemistry's `generator.py` is ordinary Python. It runs on **your** machine
with **your** permissions, so Chemart will not import it unless you say so:

```python
chemart.generate_network("ada/hypercycle-lite")
# ValueError: ada/hypercycle-lite runs its own Python code on your machine (generator.py
# at commit f2634dc6b1d4). Read it first: http://…/ada/hypercycle-lite/blob/…/generator.py
# If you trust it, pass trust_remote_code=True, and pin the version you read with
# revision='f2634dc6b1d4' so later pushes cannot change what runs.
```

This is the same contract as `transformers`. Keep three things in mind:

- **Read, then pin.** Open the code on the repo's *Files* tab, then pass the
  commit you read as `revision`. Without a pin, the next push to `main` runs
  too, and Chemart warns you whenever that would be new code.
- **Nothing runs until then.** `describe_chemistry`, `search`, `repo_info` and
  `snapshot_download` never import anything. The hub itself never runs uploaded
  code either: it validates uploads by parsing them, never by executing them.
- **Tool calls never trust.** `chemart.call_tool` refuses to pass
  `trust_remote_code`, so an LLM driving Chemart through its tools cannot be
  talked into running a stranger's code.

A repo can declare extra modules it imports (`hub.requires`). Chemart checks
that they are installed and tells you what is missing. It never installs
anything itself.

## Sharing

Create an account on the hub, make a token under **Settings → Tokens**, then:

```bash
uv run chemart login          # paste the token; stored per hub in ~/.cache/chemart/tokens.json
uv run chemart whoami
```

### A network

```python
net = chemart.generate_network("brusselator", seed=1, b=3.5)
net.push_to_hub("you/brusselator-b35",
                title="Brusselator past the Hopf point",
                description="b = 3.5 > 1 + a^2: a limit cycle.",
                tags=["oscillator", "limit-cycle"], license="CC-BY-4.0")
```

or, from a file: `uv run chemart push-network net.json you/my-network --title "..."`.

### A chemistry

```bash
uv run chemart new my-chem          # a working skeleton: chemart.yaml, generator.py, README.md
# edit my-chem/chemart.yaml (the catalog entry) and my-chem/generator.py
uv run chemart check my-chem        # the checks the hub will apply, run locally
uv run chemart push my-chem         # -> you/my-chem
```

`chemart.yaml` is a catalog file, with the same schema as the 98 built-ins
(`catalog/SCHEMA.md` in the repository), plus a `hub:` block:

```yaml
hub:
  repo_type: generator
  license: MIT
  tags: [autocatalysis]
  requires: [networkx]            # optional: modules your code imports
  requires_chemart: ">=0.1"       # optional
chemistries:
- id: my-chem                     # must equal the repo name
  name: My chemistry
  intuition: >
    What the molecules stand for, what a reaction does, and why it is interesting.
  family: core
  kind: generator
  constructive: false
  S: {repr: "...", definition: explicit}
  R: {definition: explicit, arity: [1, 2], scheme: "..."}
  A: {reactor: [ode, ssa], dilution: none}
  params:
    - {name: k, type: float, default: 1.0, min: 0, role: kinetic, meaning: "..."}
  provides: [topology, stoichiometry, rate-constants]
  fidelity: original              # or book / book+decisions / reconstructed
```

`generator.py` defines `generate(p, rng)` exactly like a built-in chemistry
module, and may use `chemart.expand`, `chemart.soup` and `chemart.helpers`.
Other `.py` files in the folder are imported relatively (`from .rules import x`).

**What a push checks.** `chemart push` runs [the contract](contributing.md)
on your machine: defaults finish in under five seconds, the network
round-trips through JSON, the same seed gives the same network, and the entry
claims everything the network contains in `provides`. It then writes
`preview.json` (the network at default parameters) and uploads. The hub
re-checks the entry and the preview by parsing them; it does not run your code.
Relative to the built-in catalog, hub entries may omit `book` and may use
`fidelity: original`; `seed`, `revision` and `trust_remote_code` cannot be
parameter names.

## Cache and offline use

Downloads are cached under `~/.cache/chemart/hub/<hub>/<namespace>--<name>/`,
one folder per commit. Every file is checked against its sha256 before it is
stored. A cached full commit id never touches the network. If the hub cannot
be reached, `main` falls back to the last commit you resolved, with a warning.

| variable | meaning |
|---|---|
| `CHEMART_HUB_URL` | the hub to talk to (default `http://127.0.0.1:8000`, a local hub) |
| `CHEMART_HOME` | cache and credentials (default `~/.cache/chemart`) |
| `CHEMART_HUB_TOKEN` | a token, overriding the one saved by `chemart login` |
| `CHEMART_HUB_OFFLINE=1` | never touch the network; use the cache only |

## Running a hub

The server lives in this repository as the `chemart-hub` package (`hub/`). It
is a FastAPI app with a SQLite database and a content-addressed file store.
Everything sits in one data folder.

```bash
uv sync --all-packages
export CHEMART_HUB_DATA=hub-data
uv run --package chemart-hub chemart-hub init
uv run --package chemart-hub chemart-hub create-user you --admin     # prompts for a password
uv run --package chemart-hub chemart-hub create-org chemart --owner you
uv run --package chemart-hub chemart-hub serve                       # http://127.0.0.1:8000

# in another terminal: stock the official shelf with the 98 built-ins
uv run --package chemart-hub chemart-hub token you > /tmp/token
uv run --package chemart-hub chemart-hub seed --token "$(cat /tmp/token)"
```

Set `CHEMART_HUB_PUBLIC_URL` to the address people reach it at, and
`CHEMART_HUB_SECURE=1` once it is behind TLS. The JSON API is documented at
`/api/docs`.

### Curating a hub by hand

A **superadmin** can change anything on the site from the browser. Make one
from the command line (this also works as a way back in if you lose access):

```bash
uv run --package chemart-hub chemart-hub set-admin you        # --off to revoke
uv run --package chemart-hub chemart-hub set-password you     # prompts
```

Superadmins get an **Admin** link in the header, which opens the panel at `/admin`:

| tab | what you can do |
|---|---|
| Repos | **feature** a repo (it appears on the home page's featured shelf), **hide** it (gone from listings and search, still loadable by id), open its editor |
| Users & orgs | grant or revoke superadmin, suspend an account (no log-in, tokens stop working, repos stay), reset a password, edit a profile or an organisation's members, delete an account with its repos, create organisations |
| Site text | the home page headline and subtitle, shelf headings, the Share page introduction, the footer, a site-wide **announcement banner** and the **About** page (markdown, sanitised) |
| Archive | everything that was deleted: **restore** it, or **purge** it for good (see below) |
| Log | every action a superadmin took on someone else's things |

Every repo also has an **✎ Edit** tab, open to its owner, its organisation's
members and superadmins. There you can change the card (title, description,
tags, license), edit, add, rename or delete any file in the browser, restore an
older version from the *Commits* tab, and delete the repo. Superadmins also get
feature, hide and **transfer** there. Browser edits are ordinary commits: they
are checked by the same rules as a push (an edit that would break the repo is
refused, with the reasons), they appear in the history, and they can be
undone with *Restore*. The hub cannot re-run a chemistry's code, so an edit to
a `.py` file does not refresh its `preview.json`; push with `chemart push` for
that.

### Deleting, the archive, and purging

Deleting never destroys anything by itself. When anyone deletes a repo (from
its Edit tab, the API or `chemart.hub.delete_repo`), or a superadmin deletes an
account, it moves to the **archive**:

- it disappears from the site, search and the API, and can no longer be
  loaded, even by a pinned commit id;
- every row and file is kept, and **Restore** in *Admin → Archive* brings it
  back exactly as it was (same versions, likes and downloads);
- the owner can reuse a deleted repo's name straight away. A deleted
  account's name stays reserved, so nobody can take over its namespace.

Only **Purge**, in the archive, removes things for good, after you type the
name to confirm (or `purge everything` for the whole archive). A purge deletes
the database rows. It also deletes from disk every stored file that no
remaining repo uses. Files shared with other repos stay, because they still
need them. Finally it compacts the database, so the purged rows leave no trace
in its free pages. The admin log keeps one line saying what was purged, and
when. Copies people already downloaded to their own machines are beyond the
hub's reach.

**Clean up orphaned files** (or `chemart-hub gc`) removes only files that no
commit uses at all, left behind by interrupted uploads, once they are an hour
old. It never touches the archive.
