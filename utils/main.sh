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

bash "$(dirname "$0")/clean_job_folder.sh" "$PROYECTO_DIR"
bash "$(dirname "$0")/gen_jobs.sh" "$PROYECTO_DIR"
bash "$(dirname "$0")/run_jobs.sh" "$PROYECTO_DIR"