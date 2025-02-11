#!/bin/sh

trap "echo Exited!; exit;" SIGINT SIGTERM

MAX_RUNTIME=10
MAX_SIZE="-20k"
FILENAME="*"
TESTDATA_DIR="../Testdata/Thingi10K/"
SD_BINARY="./build/sd"
# declare -a ALGORITHMS=("tet" "boundary" "boundary-lp")
declare -a ALGORITHMS=("tet" "boundary")

# Collect the files
FILES=$(find $TESTDATA_DIR -size $MAX_SIZE -name "${FILENAME}.*" \( -iname "*.vtk" -or -iname "*.stl" \))

ROW="File"
for ALGORITHM in "${ALGORITHMS[@]}"
do
  ROW="$ROW,$ALGORITHM"
done
echo $ROW > results.csv

# Run the benchmark
for FILE in $FILES
do
  echo "Testing on $FILE"
  ROW="$FILE"
  for ALGORITHM in "${ALGORITHMS[@]}"
  do
    OUTPUT=$(timeout -f $MAX_RUNTIME $SD_BINARY $FILE -a $ALGORITHM -b | grep "components" || echo "Timeout")
    echo "$ALGORITHM: $OUTPUT"
    ROW="$ROW,$OUTPUT"
  done
  echo $ROW >> results.csv
  echo ""
done
