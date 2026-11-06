"""Protein preparation: coarse-grained Martini proteins for the membrane builder.

You either upload a protein that is already coarse-grained (structure + `.itp` topology), or
start from an atomistic structure, which is cleaned, optionally trimmed and mutated, and then
converted with martinize2. Every protein gets a name you choose; its files are handed to the
membrane builder as `<name>.pdb` and `<name>.itp`.
"""

