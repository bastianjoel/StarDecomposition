#!/bin/sh

trap "echo Exited!; exit;" SIGINT SIGTERM

MAX_RUNTIME=10
MAX_SIZE="-2000k"
FILENAME="rest*"
TESTDATA_DIR="../Testdata/"
SD_BINARY="./build/sd"
declare -a ALGORITHMS=("tet" "boundary" "boundary-lp")

# Collect the files
FILES=$(find $TESTDATA_DIR -size $MAX_SIZE -name "${FILENAME}.*" \( -iname "*.vtk" -or -iname "*.stl" \))

# Run the benchmark
for FILE in $FILES
do
  for ALGORITHM in "${ALGORITHMS[@]}"
  do
    echo "Running $FILE with $ALGORITHM"
    timeout -f $MAX_RUNTIME $SD_BINARY $FILE -a $ALGORITHM || echo "Timeout on $FILE with $ALGORITHM"
  done
done
