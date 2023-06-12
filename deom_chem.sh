docker run -it --workdir="/home/codes/$1" --mount type=bind,source=$PWD/$1,target=/home/codes/$1 dhem/deom_mpi:chem
# docker run -it --workdir="/home/codes/$1" --mount type=bind,source=$PWD/$1,target=/home/codes/$1 dhem/deom_mpi:conda
# docker run -it --platform linux/arm64 --workdir="/home/codes/$1" --mount type=bind,source=$PWD/$1,target=/home/codes/$1 dhem/deom_mpi:arm
