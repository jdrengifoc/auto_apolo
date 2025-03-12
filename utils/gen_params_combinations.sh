#!/bin/bash

# Verify that the project directory is passed as an argument
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 <ruta_al_proyecto>"
    exit 1
fi

PROYECTO_DIR="$1"

# Verify that the folder exists
if [ ! -d "$PROYECTO_DIR" ]; then
    echo "ERROR: La carpeta '$PROYECTO_DIR' no existe."
    exit 1
fi

# Find the first `.config` and `.R` files in the project folder.
CONFIG_FILE=$(find "$PROYECTO_DIR" -maxdepth 1 -type f -name "*.config" | head -n 1)

# Verify that the config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: El archivo de configuración '$CONFIG_FILE' no existe."
    exit 1
fi

# Load configuration
source "$CONFIG_FILE"

COMBINATIONS=()
PARAM_KEYS=("${!aa_PARAMS[@]}")
NUM_PARAMS=${#PARAM_KEYS[@]}

if [ "$aa_USE_ALL_COMBINATIONS" = true ]; then
    # Generate full Cartesian product
    cartesian_product() {
        local depth="$1"
        local current="$2"

        if [ "$depth" -eq "$NUM_PARAMS" ]; then
            COMBINATIONS+=("${current:1}")  # Remove leading comma
            return
        fi

        local param="${PARAM_KEYS[depth]}"
        for value in ${aa_PARAMS[$param]}; do
            cartesian_product "$((depth + 1))" "$current,$value"
        done
    }

    cartesian_product 0 ""
else
    # Generate ntuple combinations (one per index)
    local IFS=$'\n'
    param_values=()

    for param in "${PARAM_KEYS[@]}"; do
        param_values+=("${aa_PARAMS[$param]}")
    done

    num_values=$(wc -w <<< "${aa_PARAMS[${PARAM_KEYS[0]}]}")

    for ((i = 0; i < num_values; i++)); do
        tuple=""
        for param in "${PARAM_KEYS[@]}"; do
            values=(${aa_PARAMS[$param]})
            tuple+="${values[i]},"
        done
        COMBINATIONS+=("${tuple%,}") # Remove trailing comma
    done
fi

# Print results
for combo in "${COMBINATIONS[@]}"; do
    echo "$combo"
done