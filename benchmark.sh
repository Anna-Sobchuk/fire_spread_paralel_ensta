#!/bin/bash
# Usage: ./benchmark.sh [grid_size] [wind_x] [wind_y]

GRID_SIZE=${1:-256}      # Default 256x256 grid
WIND_X=${2:-5}           # Default wind X=5
WIND_Y=${3:-0}           # Default wind Y=0
REPEAT=${4:-5}
OUTPUT="results.csv"     # Output file

echo "Threads,TotalTime,UpdateTime" > $OUTPUT

for threads in 1 2 4 8; do
    total_time=0
    update_time=0
    echo "Benchmarking with $threads threads..."
    
    for i in $(seq $REPEAT); do
        # Run simulation with timing output
        data=$(OMP_NUM_THREADS=$threads ./simulation.exe \
            -n $GRID_SIZE \
            --wind=$WIND_X,$WIND_Y \
            --start=$((GRID_SIZE/2)),$((GRID_SIZE/2)) 2>&1 | \
            grep "TimeData: ")
        
        step_total=$(echo $data | awk '{print $2}')
        step_update=$(echo $data | awk '{print $3}')

        total_time=$(echo "$total_time + $step_total" | bc)
        update_time=$(echo "$update_time + $step_update" | bc)
    done

    avg_total=$(echo "scale=4; $total_time / $REPEAT" | bc)
    avg_update=$(echo "scale=4; $update_time / $REPEAT" | bc)
    echo "$threads,$avg_total,$avg_update" >> $OUTPUT
done

echo "Results saved to $OUTPUT"