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

# Find the first `.config` file in the project folder.
CONFIG_FILE=$(find "$PROYECTO_DIR" -maxdepth 1 -type f -name "*.config" | head -n 1)

# Verify that the config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: No se encontró un archivo de configuración en '$PROYECTO_DIR'."
    exit 1
fi

# Load configuration
source "$CONFIG_FILE"

# Verify that the job folder exists
if [ -z "$aa_JOBS_folder" ]; then
    echo "ERROR: La variable 'aa_JOBS_folder' no está definida en el archivo de configuración."
    exit 1
fi

JOBS_DIR="$PROYECTO_DIR/$aa_JOBS_folder"

if [ ! -d "$JOBS_DIR" ]; then
    echo "ERROR: La carpeta de trabajos '$JOBS_DIR' no existe."
    exit 1
fi

# Run each .sh file in the job folder using sbatch
for job_script in "$JOBS_DIR"/*.sh; do
    if [ -f "$job_script" ]; then
        echo "Ejecutando: sbatch $job_script"
        sbatch "$job_script"
    fi
done

echo "✅ Todos los scripts en '$JOBS_DIR' han sido enviados a Slurm."
