#moduli da caricare
module load profile/advanced
module load autoload openmpi/1.6.3--gnu--4.5.2
module load lapack/3.3.1--gnu--4.5.2
module load blas/2007--gnu--4.5.2

#pulizia dell'esistente
make clean

#compilazione
make

