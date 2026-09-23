# Chemart Hub

The catalog holds the chemistries of the book. The **Chemart Hub** is where
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

| | a **chemistry** (code repo) | a **network** (network repo) |
|---|---|---|
| holds | a catalog entry and the Python that builds its network | one reaction network, as plain JSON |
| you load it with | `chemart.generate_network(id, ...)`, choosing your own parameters and seed | `chemart.load_network(id)` |
| runs code on your machine | yes: needs `trust_remote_code=True` | never |
| files | `chemart.yaml`, `generator.py` (+ helper `.py` files), `preview.json`, `README.md` | `network.json`, `README.md`, optional `chemart.yaml` |

The official `chemart/` organisation holds the book chemistries. Those repos
are **built-in backed**: they contain the catalog entry but no code, because
the generator ships with the library. Loading `chemart/gamma` runs your
installed `gamma` generator and needs no trust flag.

The public hub is <https://marcotuc.github.io/chemart-hub/>, built from the
[MarcoTuc/chemart-hub](https://github.com/MarcoTuc/chemart-hub) registry; it
is `chemart`'s default. See [A hub on GitHub Pages](#a-hub-on-github-pages).

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

On the public hub, sharing opens a **pull request** on its registry. Sign in
to GitHub once with the [GitHub CLI](https://cli.github.com/) (`gh auth login`)
and use the commands below as they are: your repos go under your GitHub login,
and they appear on the site once the pull request is merged. Without the
GitHub CLI, `chemart push` writes the folder and prints the steps to open the
pull request by hand. `chemart login` and `whoami` are for a hub you run
yourself.

On a hub you run yourself, create an account, make a token under
**Settings → Tokens**, then:

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

`chemart.yaml` is a catalog file, with the same schema as the built-ins
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
| `CHEMART_HUB_URL` | the hub to talk to (default `https://marcotuc.github.io/chemart-hub`, the public hub; `http://127.0.0.1:8000` for one you run with `chemart-hub serve`) |
| `CHEMART_HOME` | cache and credentials (default `~/.cache/chemart`) |
| `CHEMART_HUB_TOKEN` | a token, overriding the one saved by `chemart login` |
| `CHEMART_HUB_OFFLINE=1` | never touch the network; use the cache only |

## The simulation pit

The pit is a local app for trying chemistries out: pick one, run it, and watch
what happens. It ships in the `chemart-hub` package but has nothing to do with
a hub's repos. It has no database and no accounts, and it runs only on your
machine.

```bash
uv sync --all-packages
uv run chemart-hub pit                  # opens http://127.0.0.1:8765
uv run chemart-hub pit --port 9000 --no-browser
```

The page follows the chemistry's type and faces:

- **Given and generator chemistries.**
    - A form for the chemistry's parameters.
    - Editors for rates and initial amounts: the chemistry's own, a constant,
      a distribution, a JSON table or a file.
    - A choice of rate equations or a stochastic path. The species are
      plotted as the result streams in.
    - **Measure the network**, for its cheap and moderate measures.
    - A sweep of one argument over values and seeds, with a measure plotted
      against it.
- **Chemistries with an evolve face.** The process runs live. The page plots
  the most abundant species, the chemistry's observables and the measures
  you choose to track, against its clock, and a network window sets how many
  frames a tracked network measure sees. Tick **run until I stop** and it
  keeps going until you press **Stop**; either way what ran is kept, and the
  download link appears as soon as the run starts. The charts hold the last
  3000 points and the download the last 5000 frames.

Each finished run can be downloaded as a
[Trajectory](reference/record.md#the-trajectory-record) JSON file; the pit
keeps the last five. The same work can be scripted with `chemart simulate`,
`chemart evolve` and `chemart measure` (see [API and CLI](reference/api.md#command-line)).

The pit binds to 127.0.0.1 only. It refuses requests addressed to any other
host, and runs nothing for a request without its own header, so another
website open in your browser cannot start runs on your machine. The static
build of a hub never includes it.

## A hub on GitHub Pages

A hub does not need a server. It can be a GitHub repository, the
**registry**, from which a workflow builds a static site on GitHub Pages:
the same pages as a live hub, served as plain files, free. Contributions are
pull requests; when one is merged, the site rebuilds.

| | a live hub (`chemart-hub serve`) | a hub on GitHub Pages |
|---|---|---|
| where repos live | the hub's database and file store | `repos/<namespace>/<name>/` folders in the registry |
| sharing | `chemart push` uploads | `chemart push` opens a pull request |
| accounts | hub accounts and tokens | GitHub accounts; organisations in `namespaces.yaml` |
| checks | on upload | on every pull request (`validate.yml`) |
| revisions | one per push | one per merged commit that changed the folder |
| likes, downloads, browser editing, admin panel | yes | no |

Loading is the same on both. Point the client at the site and nothing else
changes; it recognises a static hub by its `hub.json`:

```bash
export CHEMART_HUB_URL=https://you.github.io/chemart-hub
```

**Sharing** runs `chemart push` / `push_to_hub` as usual. With the
[GitHub CLI](https://cli.github.com/) signed in (`gh auth login`), chemart
forks the registry, commits your folder on a new branch and opens the pull
request. Without it, it writes the folder into `./chemart-pull-request/` and
prints the steps. Your namespace is your GitHub login.

**Starting a registry:**

```bash
uv sync --all-packages
uv run --package chemart-hub chemart-hub init-registry chemart-hub \
    --owner <your GitHub login> --site-url https://<you>.github.io/chemart-hub
cd chemart-hub && git init -b main && git add -A && git commit -m "Start the hub"
# create the GitHub repository, push, then Settings > Pages > Source: GitHub Actions
```

That writes a README, a contributing guide, `namespaces.yaml` (you are the
first maintainer), `site.yaml`, the official shelf under `repos/chemart/`, and
three workflows:

- **validate**: on every pull request, `chemart-hub validate-pr` checks that the
  changed repos are valid, that chemistries pass `chemart check` (running the
  contributed code on a throwaway runner), and that the author changed only
  their own namespace or an organisation that lists them.
- **publish**: on every push to the default branch, `chemart-hub build-static`
  builds the site and deploys it to Pages. It loads the registry's history into
  a throwaway hub and saves every public page of the real app, so each commit
  that changed a repo becomes a revision with a stable id, and a pinned
  `revision=` keeps loading the same files.
- **sync-official**: weekly, `chemart-hub sync-official` regenerates
  `repos/chemart/` from the chemart library (archived entries excluded) and
  opens a pull request if anything changed.

To try a build locally:

```bash
uv run --package chemart-hub chemart-hub build-static ./chemart-hub ./site \
    --site-url http://127.0.0.1:8765/chemart-hub --registry-repo you/chemart-hub
mkdir -p serve && ln -s "$PWD/site" serve/chemart-hub
python -m http.server 8765 --directory serve       # http://127.0.0.1:8765/chemart-hub/
```

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

# in another terminal: stock the official shelf with the built-ins
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
