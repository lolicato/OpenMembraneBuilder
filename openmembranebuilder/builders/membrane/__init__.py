"""The membrane builder: from composition and box to ready-to-simulate systems.

It collects the lipid composition (by hand or from a template), the box, the water and ions,
and the placement of the proteins, shows a preview, and then builds one or more systems with
COBY. Each system gets its own folder with `system.gro`, `system.pdb` and `topol.top`; the
force field, simulation settings and run scripts are copied next to them.
"""

