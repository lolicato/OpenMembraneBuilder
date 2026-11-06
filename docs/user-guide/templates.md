# Membrane templates

A membrane template is a small text file with one or more lipid compositions. Use it to keep
compositions you need often, to share them, or to build several compositions at once.

Templates can be

- **predefined** – shipped with OpenMembraneBuilder and offered under *Use predefined template*;
- **uploaded** – your own file under *Upload template* (`.csv` or `.ob`, the same format).

## Format

```text title="plasma_membranes.ob"
# PlasmaMembrane

;resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
POPC,0.70,0.30,0.61,0.61
DOPC,0.30,0.70,

# CholesterolRich

;resname,upper_leaflet,lower_leaflet,apl_upper,apl_lower
POPC,0.50,0.40,0.61,0.61
CHOL,0.50,0.60,0.40,0.40
```

- **`# Name`** starts a composition (a *block*). Spaces in the name become `_`.
  A file without any `# Name` line is a single composition.
- **One lipid per line**: `resname, upper, lower[, apl_upper, apl_lower]`.
- **`;`** starts a comment; everything after it on the line is ignored. Empty lines are
  ignored too.

| Column | Meaning |
|---|---|
| `resname` | lipid residue name; must exist in the chosen force field |
| `upper_leaflet` | fraction (or number) of this lipid in the upper leaflet |
| `lower_leaflet` | fraction (or number) in the lower leaflet |
| `apl_upper` | area per lipid in the upper leaflet, nm² – optional, default **0.6** |
| `apl_lower` | area per lipid in the lower leaflet, nm² – optional, default **0.6** |

A line has either **3** columns (no areas) or **5** (areas may be left empty: `DOPC,0.30,0.70,,`).
Empty leaflet values count as 0.

## Fractions or numbers

The first lipid line decides how the whole template is read:

=== "Fractions"

    ```text
    POPC,0.70,0.30
    DOPC,0.30,0.70
    ```

    Values between 0 and 1; each leaflet should add up to 1.

=== "Numbers of lipids"

    ```text
    POPC,70,30
    DOPC,30,70
    ```

    Whole numbers: the exact number of lipids per leaflet. *Give absolute lipid values* is
    switched on automatically.

!!! warning "Do not mix"
    Use either fractions or numbers in one template.

!!! note "Whole numbers on the first line"
    A first line with only whole numbers – also `1.0` and `0` – is read as *numbers of lipids*.
    For a fraction template that starts with such a line (e.g. a pure `POPC,1.0,1.0`), untick
    **Give absolute lipid values** after loading it.

## Area per lipid

The area per lipid sets how much membrane area each lipid takes, and therefore how densely the
leaflet is packed (and, with fractions, how many lipids fit into the box).

| Smaller area | Larger area |
|---|---|
| tighter packing, more lipids | looser packing, fewer lipids |

Typical Martini values: PC lipids ≈ 0.6–0.65 nm², cholesterol ≈ 0.4 nm².

## Using a template

=== "In the app"

    Choose the template under *Use predefined template* or *Upload template*. Its first
    composition fills the lipid table, where you can still change values. The template file
    is copied to `user_inputs/` of the build.

=== "Several compositions at once"

    To build **every** composition of a template in one go, list them in `config.json` and
    use the command line. `entries` then holds one block per composition:

    ```json
    "entries": {
        "PlasmaMembrane":  [["POPC", 0.7, 0.3, 0.61, 0.61], ["DOPC", 0.3, 0.7, 0.6, 0.6]],
        "CholesterolRich": [["POPC", 0.5, 0.4, 0.61, 0.61], ["CHOL", 0.5, 0.6, 0.4, 0.4]]
    }
    ```

    ```bash
    omb cli config.json
    ```

    Each composition gets its own systems, named after the block:

    ```text
    systems/PlasmaMembrane_R0001/   systems/CholesterolRich_R0001/
    systems/PlasmaMembrane_R0002/   systems/CholesterolRich_R0002/
    ```

## Adding a predefined template

Put the file into the package's `resources/membrane_templates/` folder (in a source checkout:
`openmembranebuilder/resources/membrane_templates/`). Its file name (without `.csv`/`.ob`)
appears in the list.
