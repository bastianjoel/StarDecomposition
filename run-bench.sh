#!/bin/bash

trap "echo Exited!; exit;" SIGINT SIGTERM

THREADS=7
MAX_RUNTIME=30m
MAX_SIZE="-5000k"
FILENAME="rest"
TESTDATA_DIR="../Testdata/Locally-Injective-Mapping/3D_Parameterization/"
SD_BINARY="./build/sd"
# declare -a ALGORITHMS=("tet" "boundary" "boundary-lp")
declare -a ALGORITHMS=("tet")

# Collect the files
FILES=$(find $TESTDATA_DIR -size $MAX_SIZE -name "${FILENAME}.*" \( -iname "*.vtk" -or -iname "*.stl" \))

rm -rf .output || true
mkdir .output

IT=0
# Run the benchmark
for FILE in $FILES
do
  echo "Testing on $FILE"
  ROW="$FILE"
  for ALGORITHM in "${ALGORITHMS[@]}"
  do
    while [ $(jobs | wc -l) -ge ${THREADS} ] ; do sleep 1 ; done
    (timeout -f $MAX_RUNTIME $SD_BINARY $FILE -a $ALGORITHM -b || echo "$FILE,$ALGORITHM,TIMEOUT,0,,") > .output/${ALGORITHM}_$(date +%s%N).txt &
  done
done

wait
echo "filename,algorithm,time,result_components,result_cell_counts,result_boundary_face_counts" > results.csv
cat .output/*.txt >> results.csv
rm -rf .output
