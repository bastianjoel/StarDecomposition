#!/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

THREADS=$(nproc)
MAX_RUNTIME=30m
MAX_SIZE="-10M"
FILENAME=$2
TESTDATA_DIR=$1
SD_BINARY="./build/sd"
NUM_REPEATS=10
declare -a ALGORITHMS=("boundary" "boundary-lp" "tet")
# declare -a ALGORITHMS=("tet")

# Collect the files
FILES=$(find $TESTDATA_DIR -size $MAX_SIZE -name "${FILENAME}.*" \( -iname "*.vtk" -or -iname "*.off" -or -iname "*.stl" \) -printf '%s\t%p\n' | sort -n | cut -f2-)

rm -rf .output || true
mkdir .output

IT=0
# Run the benchmark
for FILE in $FILES
do
  echo "Testing on $FILE"
  ROW="$FILE"
  for i in $( eval echo {1..$NUM_REPEATS} )
  do
    for ALGORITHM in "${ALGORITHMS[@]}"
    do
      while [ $(jobs -p | wc -l) -ge ${THREADS} ] ; do sleep 1 ; done
      (
        timeout --foreground "$MAX_RUNTIME" "$SD_BINARY" "$FILE" -a "$ALGORITHM" -b -f -s "$i"

        status=$?

        if [ $status -eq 124 ]; then
          echo "$FILE,$ALGORITHM,TIMEOUT,0,,,,"
        elif [ $status -ne 0 ]; then
          echo "$FILE,$ALGORITHM,ERROR,0,,,,"
        fi
      ) > ".output/${ALGORITHM}_$(date +%s%N).txt" &
    done
  done
done

wait
echo "filename,algorithm,time,result_feasible,result_components,seed,result_cell_counts,result_boundary_face_counts" > results.csv
cat .output/*.txt >> results.csv
rm -rf .output
