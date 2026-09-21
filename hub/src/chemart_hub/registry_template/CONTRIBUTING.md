# Contributing to the Chemart Hub

Everything on this hub arrives as a pull request. Your repos live under
`repos/<your GitHub name>/<name>/`, and a check makes sure they are valid and
that you only change your own namespace. Once a maintainer merges your pull
request, the site rebuilds and your repo is on the shelf.

## The easy way: let chemart open the pull request

Install the [GitHub CLI](https://cli.github.com/) and sign in once with
`gh auth login`. Then point chemart at this hub:

```bash
export CHEMART_HUB_URL={{SITE_URL}}
```

**A chemistry** (a catalog entry plus the Python that builds its network):

```bash
chemart new my-chem           # a working skeleton in ./my-chem
# edit my-chem/chemart.yaml and my-chem/generator.py
chemart check my-chem         # the same checks the registry runs
chemart push my-chem          # forks the registry and opens a pull request
```

**A network** (plain data, no code):

```python
import chemart

net = chemart.generate_network("chemart/brusselator", seed=1)
net.push_to_hub("<your GitHub name>/brusselator-snapshot",
                title="A Brusselator snapshot", tags=["oscillator"])
```

Push again to update: each merged change becomes a new revision, and older
revisions stay loadable by their commit id.

## By hand

Without the GitHub CLI, `chemart push` writes the folder into
`./chemart-pull-request/` and prints the steps. They are the usual ones: fork
this repository, add your folder as `repos/<your GitHub name>/<name>/` on a new
branch, and open a pull request.

A chemistry folder holds `chemart.yaml` (a `hub:` block and one catalog entry),
`generator.py` (a top-level `generate(p, rng)`), any helper `.py` files,
`preview.json` (the network at default parameters: `chemart push` makes it for
you) and a `README.md`. A network folder holds `network.json`, an optional
`chemart.yaml` with a `hub:` block, and a `README.md`.

## What the check does

- The folder must form a valid repo: the same format rules as a live Chemart
  Hub (file types and sizes, `chemart.yaml`, `preview.json` made with the
  default parameters, no code in network repos).
- A chemistry must pass `chemart check`: the entry validates, the defaults run
  in seconds, the network is plain JSON, and the same seed gives the same
  network. This runs your code on a throwaway GitHub runner.
- You may change only `repos/<your GitHub name>/`, or the folder of an
  organisation that lists you in `namespaces.yaml`. Everything else, including
  `repos/chemart/`, the workflows, `namespaces.yaml` and `site.yaml`, is for
  maintainers: the members of the `chemart` organisation.

## Organisations

To share under an organisation's name, a maintainer adds it to
`namespaces.yaml` with its members:

```yaml
orgs:
  my-lab:
    fullname: My Lab
    members: [alice, bob]
```

## Trust

The site never runs anyone's code, and neither does loading a network.
Running a shared chemistry runs its `generator.py` on your machine, so
`chemart.generate_network` asks for `trust_remote_code=True`: read the code
first, and pin `revision=` to the commit you read.
