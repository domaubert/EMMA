Pour lancer

mpirun -np 8 ./fof

le fof va chercher automatiquement le fichier fof.in qui contient les
parametres d'execution

//-- CONTENU de FOF.in

./halo13      #nom de base des sorties
EMM           # ne pas changer, dit a fof d'utiliser le format EMMA
./data_orage/00013/part/part.00013 #les fichiers a analyser
dummyname
dummyname
0 !group size should be set to zero
0.2 !percolation parameter
10 !min mass
100000000 !mass max
.false. !star ?
.false. ! metal ?
.true. !output apres lecture
.true.  !detection des structures
.false. !lire depuis cube
.true. !timings
