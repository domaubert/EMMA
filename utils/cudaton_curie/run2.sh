#!/bin/bash
#MSUB -r atonclues           # Nom du job                
#MSUB -n 8                   # Reservation de 4 processus au total
#MSUB -N 8                   # Les processus sont répartis sur 2 noeuds
#MSUB -T 15000               # Limite de temps elapsed du job ici 600s      
#MSUB -o stdout_aton         # Sortie standard
#MSUB -e stderr_aton         # Sortie d'erreur       
#MSUB -q hybrid              # queue hybride
#MSUB -A gch0003             # ID projet

set -x
cd ${BRIDGE_MSUB_PWD}
ccc_mprun ./a.out >test.log
