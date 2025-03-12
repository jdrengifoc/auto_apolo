#!/bin/bash

# Ensure correct usage
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <param_dict> <use_all_combinations>"
    exit 1
fi

declare -A ff_PARAMS
eval "declare -A PARAMS="${1#*=}

ff_USE_ALL_COMBINATIONS="$2"
ff_COMBINATIONS=()

if [ "$ff_USE_ALL_COMBINATIONS" == "TRUE" ]; then
    # Extract param names
    PARAM_KEYS=("${!ff_PARAMS[@]}")
    NUM_PARAMS=${#PARAM_KEYS[@]}

    # Recursive function to generate Cartesian product
    cartesian_product() {
        local depth="$1"
        local current="$2"

        if [ "$depth" -eq "$NUM_PARAMS" ]; then
            ff_COMBINATIONS+=("$current")
            return
        fi

        local param="${PARAM_KEYS[depth]}"
        for value in ${ff_PARAMS[$param]}; do
            cartesian_product "$((depth + 1))" "$current,$value"
        done
    }

    cartesian_product 0 ""
else
    # Use ntuples as they are
    local IFS=$'\n'
    for tuple in ${ff_PARAMS[@]}; do
        ff_COMBINATIONS+=("$tuple")
    done
fi

# Print results
for combo in "${ff_COMBINATIONS[@]}"; do
    echo "$combo"
done
