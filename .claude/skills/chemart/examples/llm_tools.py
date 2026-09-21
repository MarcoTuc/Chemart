#!/usr/bin/env python3
"""Drive Chemart from an LLM tool loop.

    uv run python .claude/skills/chemart/examples/llm_tools.py

Chemart exposes its whole interface as three JSON-Schema tools, so a model can
discover and generate chemistries without any bespoke glue. `call_tool` returns
JSON-ready data, and `generate_network` comes back as the plain record dict.
"""

import json

import chemart

# ------------------------------------------------------------ the tool specs
tools = chemart.tool_definitions()
print(f"{len(tools)} tools: {[t['name'] for t in tools]}\n")

for t in tools:
    print(f"--- {t['name']}")
    print(f"    {t['description'][:96]}…")
    print(f"    input: {json.dumps(t['input_schema'])[:120]}…")

# The chemistry argument is an enum of implemented ids, so a model cannot ask
# for one that does not exist.
gen = next(t for t in tools if t["name"] == "generate_network")
enum = gen["input_schema"]["properties"]["chemistry"].get("enum", [])
print(f"\ngenerate_network accepts {len(enum)} chemistry ids as an enum")

# ------------------------------------------------------------- executing calls
listing = chemart.call_tool("list_chemistries")
print(f"\nlist_chemistries -> {len(listing)} rows, first: {listing[0]['id']}")

described = chemart.call_tool("describe_chemistry", {"chemistry": "gamma"})
print(f"describe_chemistry('gamma') -> keys: {sorted(described)[:7]}…")

record = chemart.call_tool("generate_network", {
    "chemistry": "brusselator",
    "seed": 1,
    "params": {"b": 3.0},
})
print(f"generate_network -> {len(record['species'])} species, "
      f"{len(record['reactions'])} reactions, status={record['status']}")

# Everything is plain JSON, so it can go straight back into a model context.
print(f"\nserialises cleanly: {len(json.dumps(record))} bytes")

# The tool path and the direct path agree exactly — useful when mixing both.
direct = chemart.generate_network("brusselator", seed=1, b=3.0).to_dict()
assert record == direct
print("call_tool result == direct API result")

# ------------------------------------------------------- a minimal agent loop
def handle(name, arguments=None):
    """What a tool-use handler looks like. Errors are actionable by design."""
    try:
        return {"ok": True, "result": chemart.call_tool(name, arguments or {})}
    except (ValueError, NotImplementedError) as err:
        # e.g. unknown chemistry -> suggests close ids; bad param -> names it
        return {"ok": False, "error": str(err)}

print("\nerror handling:")
print(" ", handle("describe_chemistry", {"chemistry": "brusselater"})["error"])
print(" ", handle("generate_network", {"chemistry": "brusselator",
                                       "params": {"b": "not a number"}})["error"])

# A practical hint for agent loops: keep records small before returning them to
# a model. Several chemistries carry genomes or sequences in Species.structure,
# and a default network can run to megabytes of JSON.
big = chemart.generate_network("automata-reaction", seed=1).to_dict()
size = len(json.dumps(big))
print(f"\nautomata-reaction default record is {size/1e6:.1f} MB of JSON — summarise before "
      f"returning records like this to a model")
print("  e.g. drop structures:",
      len(json.dumps({**big, "species": [{"id": s["id"]} for s in big["species"]]}))/1e6,
      "MB")
