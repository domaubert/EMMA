#!/bin/bash
#MSUB -r hop           # Nom du job                
#MSUB -n 1                  # Reservation de 4 processus au total
#MSUB -N 1                  # Reservation de 4 processus au total
#MSUB -T 86400                # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout_hop          # Sortie standard
#MSUB -e stderr_hop          # Sortie d'erreur       
#MSUB -p gen6667      # Allocation
#MSUB -q gpu                # sélection de la queue GPU (groupe genci ou challenge uniquement)

set -x
export DATE=`date +%F_%Hh%M`
yorick -batch split_1024.i > run_test$DATE.log
