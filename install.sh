sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics pip 
sudo dist-upgrade
pip install --upgrade scipy
echo "sorry for the sudo but gmsh2.12 needs to be installed"
wget http://gmsh.info/bin/Linux/gmsh-2.12.0-Linux64.tgz
tar zxvf gmsh-2.12.0-Linux64.tgz
sudo cp gmsh-2.12.0-Linux/bin/gmsh /usr/bin/
sudo cp -r gmsh-2.12.0-Linux/share/doc/gmsh /usr/share/doc/
sudo cp gmsh-2.12.0-Linux/share/man/man1/gmsh.1 /usr/share/man/man1/
rm gmsh-2.12.0-Linux64.tgz

