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

# Find and delete all directories inside the given project directory
find "$PROYECTO_DIR" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +

echo "✅ All subdirectories inside '$PROYECTO_DIR' have been deleted."