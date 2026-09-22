"""Prime number chemistry: book 2.5.2, appendix NumberChem.py, Banzhaf et al. 1996 [72]."""

from collections import Counter
from statistics import mean

import pytest

from chemart import evolve, generate_network
from chemart.chemistries.prime_number_chemistry import is_prime

ID = "prime-number-chemistry"


def value(species_id):
    return int(species_id[1:])


def consumed(reaction):
    """Species whose count drops in the reaction."""
    return {s for s, n in reaction.reactants.items() if reaction.products.get(s, 0) < n}


def assert_division_rule(reaction):
    """Every reaction is s1 + s2 -> s1 + s2/s1 with s1 < s2 and s1 | s2 (book eq. 2.43)."""
    lhs = sorted(value(s) for s, n in reaction.reactants.items() for _ in range(n))
    rhs = sorted(value(s) for s, n in reaction.products.items() for _ in range(n))
    s1, s2 = lhs
    assert s1 < s2 and s2 % s1 == 0
    assert rhs == sorted([s1, s2 // s1])


def division_fixed_point(seed):
    """Independent closure: add quotients b // a until nothing new appears."""
    known = set(seed)
    while True:
        new = {b // a for a in known for b in known if a < b and b % a == 0} - known
        if not new:
            return known
        known |= new


def test_closure_of_12_2_3():
    net = generate_network(ID, numbers=[12, 2, 3])
    assert net.status == "complete"
    assert sorted(value(s.id) for s in net.species) == [2, 3, 4, 6, 12]
    lines = set(net.to_text().splitlines())
    assert {"n2 + n12 -> n2 + n6", "n3 + n12 -> n3 + n4", "n2 + n4 -> 2 n2"} <= lines
    for r in net.reactions:
        assert_division_rule(r)


def test_closure_contains_only_numbers_reachable_by_division():
    net = generate_network(ID, seed=3)
    seed = [value(s) for s in net.extras["seed"]]
    numbers = {value(s.id) for s in net.species}
    assert numbers == division_fixed_point(seed)
    # every product divides some seed number: the closure is finite
    assert all(any(m % n == 0 for m in seed) for n in numbers)
    assert net.status == "complete"


def test_primes_are_exactly_the_non_reactive_numbers():
    """Book 2.5.2: 'the reaction rule precisely determines primes as nonreactive'."""
    net = generate_network(ID, seed=5)
    numbers = {value(s.id) for s in net.species}
    eaten = {value(s) for r in net.reactions for s in consumed(r)}
    assert not any(is_prime(n) for n in eaten)
    divisible = {b for b in numbers if any(a < b and b % a == 0 for a in numbers)}
    assert eaten == divisible
    assert set(net.extras["primes"]) == {s.id for s in net.species if is_prime(value(s.id))}


def test_closure_truncates_on_budget():
    net = generate_network(ID, numbers=[2 ** 10, 2], max_species=4)
    assert net.status == "truncated" and len(net.species) <= 4


def test_soup_follows_numberchem():
    traj = evolve(ID, seed=11)
    net = traj.network
    assert net.status == "observed"
    start = net.initial_state
    assert sum(start.values()) == 100
    assert all(2 <= value(s) <= 1000 for s in start)
    # divisor is a catalyst and both molecules are reinjected: population stays M
    assert sum(net.extras["final_state"].values()) == 100
    assert sum(r.count for r in net.reactions) == net.extras["analysis"]["effective_collisions"] > 0
    for r in net.reactions:
        assert_division_rule(r)
        assert not any(is_prime(value(s)) for s in consumed(r))
    fraction = traj.series("prime_fraction")
    assert len(fraction) == 10000 // 100 + 1
    assert traj.times()[1] == 100.0
    # primes are never consumed, so the prime count never decreases
    assert all(b >= a for a, b in zip(fraction, fraction[1:]))


def test_prime_fraction_grows_in_the_book_default_run():
    """Book 2.5.2 / figs. 2.8-2.9: M = 100 from [2, 1000], 10000 iterations; primes emerge and stay."""
    runs = [evolve(ID, seed=s).series("prime_fraction") for s in range(10)]
    assert all(f[-1] > f[0] for f in runs)
    assert mean(f[0] for f in runs) < 0.3
    assert mean(f[-1] for f in runs) > 0.95
    # grows through the run, not only at the end
    assert all(f[10] <= f[50] <= f[-1] for f in runs)


def test_constructive_only_with_a_wide_interval():
    """Book 2.5.2: each of 2..101 gives no new numbers; 100 draws from [2, 1000] do."""
    closed = evolve(ID, seed=2, numbers=list(range(2, 102)), iterations=2000)
    assert closed.network.extras["analysis"]["new_numbers"] == []
    assert closed.series("prime_fraction")[-1] > closed.series("prime_fraction")[0]
    open_ = evolve(ID, seed=2).network
    assert len(open_.extras["analysis"]["new_numbers"]) > 0


def test_no_prime_factorisation():
    """Book 2.5.2: divisors are catalysts, so no copies are made. Each molecule only ever
    becomes a divisor of itself, so a prime p cannot outnumber the initial multiples of p."""
    net = evolve(ID, seed=4).network
    start = Counter({value(s): c for s, c in net.initial_state.items()})
    for s, copies in net.extras["final_state"].items():
        p = value(s)
        if is_prime(p):
            assert copies <= sum(c for n, c in start.items() if n % p == 0)


def test_soup_size_matters():
    """Paper [72] fig. 6 (reactor algorithm II, maxn = 10000): small soups run into dead
    ends with non-primes left, large soups end almost all-prime. Shortened to 200 generations."""
    def final(M, s):
        return evolve(ID, seed=s, M=M, maxn=10000, iterations=200 * M).series("prime_fraction")[-1]

    small = [final(15, s) for s in range(6)]
    large = [final(150, s) for s in range(6)]
    assert mean(small) < 0.5
    assert mean(large) > mean(small) + 0.3
    assert sum(f == 1.0 for f in small) == 0


@pytest.mark.parametrize(
    "run, given, message",
    [
        (generate_network, dict(minn=500, maxn=100), "maxn must be >= minn"),
        (evolve, dict(minn=500, maxn=100), "maxn must be >= minn"),
        (generate_network, dict(numbers=[4, 1]), "integers >= 2"),
        (evolve, dict(numbers=[7]), "at least 2 molecules"),
    ],
)
def test_rejects_inconsistent_parameters(run, given, message):
    with pytest.raises(ValueError, match=message):
        run(ID, **given)
