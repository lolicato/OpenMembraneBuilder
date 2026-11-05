#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$RUN_DIR/.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/scripts"
MDP_DIR="$ROOT_DIR/mdp"
CONFIG_FILE="$ROOT_DIR/config.json"

source "$SCRIPT_DIR/index.sh"

STEPS=("em" "eq1" "eq2" "eq3" "md")

for system_dir in ./*; do
    [ -d "$system_dir" ] || continue

    echo "========================================"
    echo "Running system: $system_dir"
    echo "========================================"

    cd "$system_dir"

    input_gro="system.gro"

    for step in "${STEPS[@]}"; do
        mdp_file="$MDP_DIR/${step}.mdp"
        tpr_file="${step}.tpr"
        output_gro="${step}.gro"

        echo ""
        echo ">>> grompp: $step"

        if [ "$step" = "em" ]; then
            gmx grompp \
                -f "$mdp_file" \
                -c "$input_gro" \
                -p topol.top \
                -o "$tpr_file" \
                -maxwarn 1

        elif [ "$step" = "md" ]; then
            gmx grompp \
                -f "$mdp_file" \
                -c "$input_gro" \
                -p topol.top \
                -n index.ndx \
                -o "$tpr_file" \
                -maxwarn 1

        else
            gmx grompp \
                -f "$mdp_file" \
                -c "$input_gro" \
                -p topol.top \
                -n index.ndx \
                -o "$tpr_file" \
                -r "$input_gro" \
                -maxwarn 1
        fi

        echo ">>> mdrun: $step"

        gmx mdrun \
            -deffnm "$step" \
            -v

        input_gro="$output_gro"

        if [ "$step" = "em" ]; then
            echo ">>> Creating index.ndx"

            make_index_input "$CONFIG_FILE"

            gmx make_ndx \
                -f system.gro \
                -o index.ndx \
                < inp.inp

            rm -f inp.inp
        fi
    done

    cd "$RUN_DIR"
done

echo ""
echo "All simulations completed."
