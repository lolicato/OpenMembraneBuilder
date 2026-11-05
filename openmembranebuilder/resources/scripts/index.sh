#!/usr/bin/env bash

# Writes the make_ndx input for the index groups used by the mdp files
# (MEMB, SOLV, SOLU and, with proteins, PROT), from the build's config.json.
# Works with the wizard's config.json ({"martinize": …, "membrane": {…}}) and with a
# flat membrane config; "entries" can be one lipid list or one list per template.
make_index_input() {
    local json_file="$1"
    local out_file="${2:-inp.inp}"

    local cfg='(.membrane // .)'

    # lipids with a non-zero amount in at least one leaflet / template
    local lipids
    lipids=$(jq -r "$cfg"' | .entries
        | (if type == "object" then [.[][]] else . end)
        | map(select(((.[1] | tonumber? // 0) > 0) or ((.[2] | tonumber? // 0) > 0)) | .[0])
        | unique | .[]' "$json_file" \
        | sed 's/^/r /' \
        | paste -sd'|' - \
        | sed 's/|/ | /g')

    local solvation
    solvation=$(jq -r "$cfg"' | .solvation // ""' "$json_file")

    local solvents
    solvents=$(echo "$solvation" \
        | grep -oE 'solv:[^ ]+' \
        | awk -F: '{print "r "$2}' \
        | paste -sd'|' - \
        | sed 's/|/ | /g')

    local pos neg
    pos=$(echo "$solvation" | grep -oE 'pos:[^ ]+' | cut -d: -f2 || true)
    neg=$(echo "$solvation" | grep -oE 'neg:[^ ]+' | cut -d: -f2 || true)

    local solv_sel="$solvents"
    [ -n "$pos" ] && solv_sel="$solv_sel | r $pos"
    [ -n "$neg" ] && solv_sel="$solv_sel | r $neg"

    # proteins: anything inserted on the protein page
    local n_prot
    n_prot=$(jq -r "$cfg"' | (.protein_files // {}) | length' "$json_file")

    {
        echo "del 1-100"

        # Group 1: membrane
        echo "$lipids"
        echo "name 1 MEMB"

        # Group 2: solvent + ions
        echo "$solv_sel"
        echo "name 2 SOLV"

        if [ "${n_prot:-0}" -gt 0 ]; then
            # Group 3: protein
            echo "0 &! 1 &! 2"
            echo "name 3 PROT"

            # Group 4: protein + membrane
            echo "1 | 3"
            echo "name 4 SOLU"
        else
            # Group 3: membrane only
            echo "1"
            echo "name 3 SOLU"
        fi

        echo "q"
    } > "$out_file"
}
