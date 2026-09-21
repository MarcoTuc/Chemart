# Chemart Hub registry

This repository is a [Chemart Hub](https://github.com/{{CHEMART_REPO}}): a shelf of
artificial chemistries and reaction networks that anyone can load in one line
of Python. The website at **{{SITE_URL}}** is built from it, and rebuilds
every time a pull request is merged.

```python
import chemart

net = chemart.generate_network("chemart/brusselator", seed=0)   # a book chemistry
net = chemart.load_network("<someone>/<network>")                # a shared network
```

Point `chemart` at this hub with `CHEMART_HUB_URL={{SITE_URL}}`.

## What is in here

| path | what it holds |
|---|---|
| `repos/<namespace>/<name>/` | one hub repo per folder: a chemistry (`chemart.yaml`, `generator.py`, `preview.json`, `README.md`) or a network (`network.json`, `chemart.yaml`, `README.md`) |
| `repos/chemart/` | the official shelf: the chemistries of Banzhaf & Yamamoto's *Artificial Chemistries*, kept in step with the `chemart` library by a scheduled workflow |
| `namespaces.yaml` | organisations and their members, and optional profile texts |
| `site.yaml` | featured and hidden repos, and the site's texts |

Every commit that changes a repo's folder becomes a revision of that repo, so
`revision="<commit>"` pins exactly what you read.

## Sharing yours

See [CONTRIBUTING.md](CONTRIBUTING.md). In short:

```bash
chemart new my-chem           # a working skeleton
chemart check my-chem         # the checks this registry runs
chemart push my-chem          # opens a pull request (needs `gh auth login` once)
```

## How the site is built

`.github/workflows/publish.yml` runs `chemart-hub build-static` on every push to
the default branch and deploys the result to GitHub Pages. It renders the
Chemart Hub's own pages, so the site looks and works like a live hub, without
a server. `.github/workflows/validate.yml` checks every pull request, and
`.github/workflows/sync-official.yml` refreshes `repos/chemart/` weekly.
