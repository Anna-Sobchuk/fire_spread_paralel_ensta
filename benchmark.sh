#!/bin/bash
# Usage: ./benchmark.sh [grid_size] [wind_x] [wind_y] [repeat]

GRID_SIZE=${1:-256}      # Grille par défaut 256x256
WIND_X=${2:-5}           # Vent X=5 par défaut
WIND_Y=${3:-0}           # Vent Y=0 par défaut
REPEAT=${4:-5}           # 5 répétitions par défaut
OUTPUT="results.csv"     # Fichier de sortie

echo "Processes,TotalTime(s),UpdateTime(s)" > $OUTPUT

for processes in 1 2 4 8; do
    total_time=0
    update_time=0
    echo "Benchmark avec $processes processus MPI..."
    
    for i in $(seq $REPEAT); do
        # Exécution avec MPI et récupération des temps
        data=$(mpirun -np $processes ./simulation \
            -n $GRID_SIZE \
            --wind=$WIND_X,$WIND_Y \
            --start=$((GRID_SIZE/2)),$((GRID_SIZE/2)) 2>&1 | \
            grep "PERF_DATA:" | awk -F':' '{print $2}')
        
        # Extraction des temps
        step_total=$(echo $data | awk '{print $1}')
        step_update=$(echo $data | awk '{print $2}')

        total_time=$(echo "$total_time + $step_total" | bc)
        update_time=$(echo "$update_time + $step_update" | bc)
    done

    # Calcul des moyennes
    avg_total=$(echo "scale=4; $total_time / $REPEAT" | bc)
    avg_update=$(echo "scale=4; $update_time / $REPEAT" | bc)
    
    # Enregistrement dans le CSV
    echo "$processes,$avg_total,$avg_update" >> $OUTPUT
done

echo "Résultats enregistrés dans $OUTPUT"