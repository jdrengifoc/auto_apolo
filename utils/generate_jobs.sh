# Load the config
source "$CONFIG_FILE"

# Convert arrays from config to usable arrays in the script
eval "PARAM1_LIST=(${aa_param1[@]})"
eval "PARAM2_LIST=(${aa_param2[@]})"

echo "Param1 list: ${PARAM1_LIST[@]}"
echo "Param2 list: ${PARAM2_LIST[@]}"

# Check if arrays are empty
if [ ${#PARAM1_LIST[@]} -eq 0 ] || [ ${#PARAM2_LIST[@]} -eq 0 ]; then
  echo "Error: PARAM1_LIST or PARAM2_LIST is empty"
  exit 1
fi

# Generate jobs based on combinations
if [ "$aa_USE_ALL_COMBINATIONS" = true ]; then
  for p1 in "${PARAM1_LIST[@]}"; do
    for p2 in "${PARAM2_LIST[@]}"; do
      echo "Creating job for param1=$p1, param2=$p2"
      # Your job creation logic here
    done
  done
else
  for ((i=0; i<${#PARAM1_LIST[@]}; i++)); do
    echo "Creating paired job for param1=${PARAM1_LIST[i]}, param2=${PARAM2_LIST[i]}"
    # Your job creation logic here
  done
fi