#!/bin/bash


# ----------------------------------------------------------------------------

sudo apt-get update
sudo apt-get upgrade -y

sudo apt-get install -y gcc g++  cmake make build-essential unzip m4 libxml2-dev python3 python2.7 python-numpy python-scipy python3-numpy python3-scipy python-pmw freeglut3-dev libpng-dev libfreetype6-dev libglew-dev gabedit qtiplot r-base grace csh python-matplotlib   libblas-dev libatlas-base-dev libreadline7 libreadline-dev libgfortran3 openssh-client openssh-server dssp gnuplot ffmpeg python-gtk2 python-pyqt5

sudo apt-get install csh gfortran flex bison    # For AmberTools19
sudo apt-get install python3-setuptools # for propka31
sudo apt-get install -y gcc-5 g++-5 nvidia-cuda-toolkit 

# ----------------------------------------------------------------------------

n=2
download=1

mkdir local
cd local

# FFTW - Fastest Fourier Transform in the West ; http://www.fftw.org/
v='3.3.8'
if [ $download -eq 1 ]; then wget http://www.fftw.org/fftw-"$v".tar.gz ; fi
tar -xzvf fftw-"$v".tar.gz 
cd fftw-"$v"/
./configure --with-pic --enable-float --enable-threads --enable-sse2
make -j "$n"
sudo make install -j "$n"
make distclean
./configure --with-pic --enable-threads --enable-sse2
make -j "$n"
sudo make install -j "$n"
cd ../



# OpenMPI
v='4.0'
vv='4.0.0'
if [ $download -eq 1 ]; then wget  https://www.open-mpi.org/software/ompi/v"$v"/downloads/openmpi-"$vv".tar.gz ; fi
tar -xzvf openmpi-"$vv".tar.gz
cd openmpi-"$vv"
./configure --enable-static --disable-shared
make -j "$n"
sudo make install -j "$n"
cd ../



# GROMACS --------------------------------------------------------------------
v='2019.4'
if [ $download -eq 1 ]; then wget  ftp://ftp.gromacs.org/pub/gromacs/gromacs-"$v".tar.gz ; fi
tar -xzvf gromacs-"$v".tar.gz 
cd gromacs-"$v"/

# "default" build
mkdir build
cd build
cmake .. -DSHARED_LIBS_DEFAULT=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=YES -DGMX_BUILD_OWN_FFTW=OFF -DGMX_DEFAULT_SUFFIX=OFF -DGMX_MPI=OFF -DGMX_GPU=OFF -DGMX_DOUBLE=OFF -DGMX_BUILD_MDRUN_ONLY=OFF -DCMAKE_INSTALL_PREFIX=/home/"$USER"/local/gromacs-"$v"/
make -j "$n"
make install -j "$n"
sudo ln -s /home/"$USER"/local/gromacs-"$v"/bin/gmx /usr/local/bin/gmx-"$v"
cd ..

# double precision
mkdir build_double
cd build_double
cmake .. -DSHARED_LIBS_DEFAULT=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=YES -DGMX_BUILD_OWN_FFTW=OFF -DGMX_DEFAULT_SUFFIX=ON -DGMX_MPI=OFF -DGMX_GPU=OFF -DGMX_DOUBLE=ON -DGMX_BUILD_MDRUN_ONLY=OFF -DCMAKE_INSTALL_PREFIX=/home/"$USER"/local/gromacs-"$v"/
make -j "$n"
make install -j "$n"
sudo ln -s /home/"$USER"/local/gromacs-"$v"/bin/gmx_d /usr/local/bin/gmx_d-"$v"
cd ..

# gpus
mkdir build_gpu
cd build_gpu
CC=gcc-5 CXX=g++-5 cmake .. -DSHARED_LIBS_DEFAULT=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=YES -DGMX_BUILD_OWN_FFTW=OFF -DGMX_DEFAULT_SUFFIX=ON -DGMX_MPI=ON -DGMX_GPU=ON -DGMX_DOUBLE=OFF -DGMX_BUILD_MDRUN_ONLY=ON -DCMAKE_INSTALL_PREFIX=/home/"$USER"/local/gromacs-"$v"/
make -j "$n"
make install -j "$n"
sudo ln -s /home/"$USER"/local/gromacs-"$v"/bin/mdrun_mpi /usr/local/bin/mdrun_mpi-"$v"
cd ..


sudo ln -s /usr/local/bin/gmx-"$v" /usr/local/bin/gmx
sudo ln -s /usr/local/bin/gmx_d-"$v" /usr/local/bin/gmx_d


cd ..



# ----------------------------------------------------------------------------


# DSSP
echo '' >> ~/.bashrc
echo 'export DSSP=/usr/bin/dssp' >> ~/.bashrc
source ~/.bashrc


# MGL Tools - Molecular Graphics Laboratory Tools ; http://mgltools.scripps.edu/
v='x86_64Linux2_1.5.6'
if [ $download -eq 1 ]; then wget  http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_"$v".tar.gz ; fi
tar -xzvf mgltools_"$v".tar.gz
cd mgltools_"$v"/
./install.sh
dir=$(pwd)
echo 'source '"$dir"'/initMGLtools.sh' >> ~/.bashrc
source initMGLtools.sh
cd ../






# AutoDock Vina
v='1_1_2_linux_x86'
if [ $download -eq 1 ]; then wget  http://vina.scripps.edu/download/autodock_vina_"$v".tgz ; fi
tar -xzvf autodock_vina_"$v".tgz 
sudo ln -s /home/"$USER"/local/autodock_vina_"$v"/bin/vina /usr/local/bin/







# PyMOL
if [ $download -eq 1 ]; then wget  https://github.com/msgpack/msgpack-c/archive/master.zip ; fi
mv master.zip msgpack-c-master.zip
unzip msgpack-c-master.zip 
cd msgpack-c-master/
cmake .
make -j "$n"
sudo make install
cd ../

v='2.1.0'
if [ $download -eq 1 ]; then wget  https://downloads.sourceforge.net/project/pymol/pymol/1.8/pymol-v"$v".tar.bz2 ; fi
tar -xjvf pymol-v"$v".tar.bz2
mv pymol/ pymol-v"$v"/
cd pymol-v"$v"/
python2.7 setup.py build install --home=/home/"$USER"/local/pymol-v"$v"/ --install-lib=/home/"$USER"/local/pymol-v"$v"/modules/ --install-scripts=/home/"$USER"/local/pymol-v"$v"/

echo '' >> ~/.bashrc
echo 'export PYMOL_PATH=/home/'"$USER"'/local/pymol-v'"$v"'/' >> ~/.bashrc
echo 'export PYTHONPATH=$PYTHONPATH:$PYMOL_PATH/modules:' >> ~/.bashrc
source ~/.bashrc
sudo ln -s /home/"$USER"/local/pymol-v"$v"/pymol /usr/local/bin/
cd ../









# AmberTools; http://ambermd.org/
vat='19'  # AmberTools version
va='18'   # Amber version
tar -xjvf AmberTools"$vat".tar.bz2 
cd amber"$va"
dir=$(pwd)
path_python=$(whereis python3.5 | awk '{print $2}')
export AMBERHOME="$dir"
./update_amber --update
./configure -noX11 --with-python "$path_python" gnu
source amber.sh
echo 'source '"$dir"'/amber.sh' >> ~/.bashrc
make install -j "$n"
sudo ln -s /home/"$USER"/local/amber"$va"/bin/* /usr/local/bin/
cd ..










# Coot
./ccp4-7.0-setup-linux64 
echo 'alias coot=/home/'"$USERS"'/local/ccp4/ccp4-7.0/bin/coot' >> ~/.bashrc
source ~/.bashrc


# Phenix
v='1.14-3260'
tar -xzvf phenix-installer-"$v"-intel-linux-2.6-x86_64-centos6.tar.gz
cd phenix-installer-"$v"-intel-linux-2.6-x86_64-centos6/
./install --prefix=/home/"$USER"/local/ --nproc="$n" --openmp --python_static
echo 'source /home/'"$USER"'/local/phenix-'"$v"'/phenix_env.sh' >> ~/.bashrc
source ~/.bashrc
cd ../
mv phenix-installer-"$v"-intel-linux-2.6-x86_64-centos6/ phenix-"$v"/

# XDS
tar -xzvf XDS-INTEL64_Linux_x86_64.tar.gz 
sudo ln -s /home/"$USER"/local/XDS-INTEL64_Linux_x86_64/* /usr/local/bin/

# autoPROC
mkdir autoproc
mv ../GPhL_autoPROC_snapshot_20180515_install.sh .
mv ../GPhL_autoPROC_snapshot_20180515.linux64.tar .
mv ../GPhL_autoPROC_snapshot_20180515_licence .licence
chmod +x GPhL_autoPROC_snapshot_20180515_install.sh 
echo 'source /home/'"$USER"'/local/autoproc/setup.sh' > ~/.bashrc
source ~/.bashrc
sudo ln -s /home/"$USER"/local/XDS-INTEL64_Linux_x86_64/* /usr/local/bin/
# 	process -checkdeps

# XCHEM
unzip XChemExplorer-master.zip 
cd XChemExplorer-master/
dir=$(pwd)
cat XChemExplorer_dmd.sh | grep -v 'module unload ccp4' | grep -v 'source /dls/science/groups/i04-1/software/pandda-update/ccp4/ccp4-7.0/bin/ccp4.setup-sh' 

echo '#!/bin/bash'                                           > XChemExplorer_dmd.sh
echo ''                                                     >> XChemExplorer_dmd.sh
echo 'export XChemExplorer_DIR='"$dir"                      >> XChemExplorer_dmd.sh
echo 'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh' >> XChemExplorer_dmd.sh
echo ''                                                     >> XChemExplorer_dmd.sh
echo 'ccp4-python $XChemExplorer_DIR/XChemExplorer.py'      >> XChemExplorer_dmd.sh

sudo ln -s /home/"$USER"/local/XChemExplorer-master/XChemExplorer  /usr/local/bin/













# MOPAC - Molecular Orbital PACkage ; http://openmopac.net/
v='2016'
mkdir mopac
cd mopac/
unzip ../MOPAC"$v"_for_Linux_64_bit.zip
chmod +x MOPAC"$v".exe
dir=$(pwd)
echo '' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":'"$dir" >> ~/.bashrc
echo 'export MOPAC_LICENSE='"$dir" >> ~/.bashrc
source ~/.bashrc 
sudo ln -s /home/"$USER"/local/mopac/MOPAC"$v".exe /usr/local/bin/mopac
./MOPAC"$v".exe 11902245a29278647
cd ../











# Propka ; http://propka.ki.ku.dk/
sudo apt-get install apbs pdb2pqr
sudo rm -f /usr/bin/propka # this is created by apbs when installed via sudo apt-get
v='3.1'
vv='31'
if [ $download -eq 1 ]; then wget https://github.com/jensengroup/propka-"$v"/archive/master.zip ; mv master.zip propka-"$v"-master.zip ; fi
unzip propka-"$v"-master.zip
cd propka-"$v"-master/
python3 setup.py install --user --install-scripts /home/"$USER"/local/propka-"$v"-master/
sudo ln -s /home/"$USER"/local/propka-"$v"-master/propka"$vv" /usr/local/bin/propka"$v"

cd ..


# APBS - Adaptive Poisson-Boltzmann Solver ; http://www.poissonboltzmann.org/apbs; http://sourceforge.net/projects/apbs/
v='1.5-linux64'
if [ $download -eq 1 ]; then wget https://downloads.sourceforge.net/project/apbs/apbs/apbs-1.5/APBS-"$v".tar.gz ; fi
sudo ln -s /lib/x86_64-linux-gnu/libreadline.so.7 /lib/x86_64-linux-gnu/libreadline.so.6
tar -xzvf APBS-"$v".tar.gz
cd APBS-"$v"/
dir=$(pwd)
echo 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":'"$dir"'/lib' >> ~/.bashrc
source ~/.bashrc 
sudo ln -s /home/"$USER"/local/APBS-"$v"/bin/apbs /usr/local/bin/
cd ..

#PDB2PQR ; https://sourceforge.net/projects/pdb2pqr/
v='2.1.0'
tar -xzvf pdb2pqr-linux-bin64-"$v".tar.gz
cd pdb2pqr-linux-bin64-"$v"/
sudo ln -s /home/"$USER"/local/pdb2pqr-linux-bin64-"$v"/pdb2pqr /usr/local/bin/
cd ../









# GAMESS
tar -xavf gamess-current.tar.gz
cd gamess/
echo -e "\n""linux64""\n""/home/""$USER""/local/gamess""\n""/home/""$USER""/local/gamess""\n""00""\n""gfortran""\n""4.8""\n""\n""atlas""\n"/usr/lib/atlas-base/"\n""\n""\n""sockets""\n""no""\n" | ./config

cd ddi
./compddi
mv ddikick.x ../
cd ../

sudo sh -c "echo 'kernel.shmmax = 3064372224' >> /etc/sysctl.conf"
sudo sh -c "echo 'kernel.shmall = 748137' >> /etc/sysctl.conf"
sudo sysctl -p

./compall

./lked gamess 00

mv rungms original_rungms
pc_name=$(hostname)
max_cores=$(grep -c ^processor /proc/cpuinfo)
sed 's/set SCR=\/scr\/$USER/set SCR=`pwd`/g' original_rungms | sed 's/set USERSCR=\/u1\/$USER\/scr/set USERSCR=`pwd`/g' | sed 's/set GMSPATH=\/u1\/mike\/gamess/set GMSPATH=\/home\/$USER\/programs\/gamess\//g' | sed 's/switch (`hostname`)/switch (`hostname`)\n         case '"$pc_name"':\n            set NNODES=1\n            set NCPUS=$NCPUS\n            set HOSTLIST=(`hostname`:cpus=$NCPUS)\n            breaksw/g' > rungms
chmod +x rungms









