#!/bin/bash
#MSUB -r hop           # Nom du job                
#MSUB -n 1                  # Reservation de 4 processus au total
#MSUB -N 1                  # Les processus sont répartis sur 2 noeuds
#MSUB -T 86400                # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout_hop          # Sortie standard
#MSUB -e stderr_hop          # Sortie d'erreur       
#MSUB -p gen2191      # Allocation
#MSUB -q mono                # sélection de la queue GPU (groupe genci ou challenge uniquement)

set -x
export DATE=`date +%F_%Hh%M`
./hop -in /scratch/cont003/pocvirk/localgroup_64proc/output_00015/part_00015.out -p 1. > run_hop$DATE.log
