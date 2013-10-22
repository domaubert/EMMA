#!/bin/bash
#MSUB -r lidau50           # Nom du job                
#MSUB -n 64                  # Reservation de 4 processus au total
#MSUB -N 32                  # Les processus sont répartis sur 2 noeuds
#MSUB -T 86400                # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout_radramses          # Sortie standard
#MSUB -e stderr_radramses          # Sortie d'erreur       
#MSUB -p gen2191      # Allocation
#MSUB -q gpu                # sélection de la queue GPU (groupe genci ou challenge uniquement)

set -x
export DATE=`date +%F_%Hh%M`
mpirun ./ramses_1024_50 ./gal_radramses.nml > run$DATE.log
