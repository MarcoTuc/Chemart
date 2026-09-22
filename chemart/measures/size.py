"""A. Size and composition: reference values for normalising, not axes."""

from __future__ import annotations

from collections import Counter

from chemart.measures import register


@register("n_species", "A")
def n_species(ctx) -> int:
    """Number of species."""
    return len(ctx.net.species)


@register("n_reactions", "A")
def n_reactions(ctx) -> int:
    """Number of reactions."""
    return len(ctx.net.reactions)


@register("reaction_density", "A")
def reaction_density(ctx) -> float:
    """Reactions per species, r / n."""
    return len(ctx.net.reactions) / len(ctx.net.species)


@register("arity", "A", needs="S")
def arity(ctx) -> dict[str, float]:
    """Share of reactions by molecularity, "reactants->products" (e.g. "2->1": 0.4),
    counting molecules with their multiplicity."""
    counts = Counter(f"{sum(r.reactants.values())}->{sum(r.products.values())}" for r in ctx.net.reactions)
    total = sum(counts.values())
    return {k: n / total for k, n in sorted(counts.items())}


@register("mean_reactants", "A", needs="S")
def mean_reactants(ctx) -> float:
    """Mean number of reactant molecules per reaction."""
    return sum(sum(r.reactants.values()) for r in ctx.net.reactions) / len(ctx.net.reactions)


@register("mean_products", "A", needs="S")
def mean_products(ctx) -> float:
    """Mean number of product molecules per reaction."""
    return sum(sum(r.products.values()) for r in ctx.net.reactions) / len(ctx.net.reactions)


@register("reversible_fraction", "A", needs="S")
def reversible_fraction(ctx) -> float:
    """Share of reactions whose exact reverse is also in the network."""
    keys = {(frozenset(r.reactants.items()), frozenset(r.products.items())) for r in ctx.net.reactions}
    back = sum(1 for r in ctx.net.reactions
               if (frozenset(r.products.items()), frozenset(r.reactants.items())) in keys)
    return back / len(ctx.net.reactions)


@register("catalysed_fraction", "A")
def catalysed_fraction(ctx) -> float:
    """Share of reactions with a catalyst: a species on both sides."""
    return sum(1 for r in ctx.net.reactions if r.catalysts) / len(ctx.net.reactions)
