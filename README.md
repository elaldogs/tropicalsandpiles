# Tropical Sandpiles
Simulations of tropical sandpiles and tropical curves.
 - parallelsandpile outputs a file called grid.dat which contains a representation of the final state of the sandpile. It requires MPI to be compiled and run it.
 - visualizegrid reads grid.dat and displays the final state of the sandpile.

These programs are part of work with Nikita Kalinin, Ernesto Lupercio, Mikhail Shkolnikov, Vera Kalinina, and Yulieth Prieto.


# Manual to tropical (linearized) sandpile model:
- To run tropical sandpiles compile:
g++ -std=c++11 -O3 linearsandpile.cpp -o linearsandpile
- and run (with default parameters):
./linearsandpile
- this will produce a bunch of files in the folder tsandpile/
- file power1000_900_2.txt contains the avalanche sizes, with the sign '-' if the corresponding avalanche touched the boundary
- To visualize the corresponding tropical curve type:
python vizualiselinearsand.py
- this will show you a picture of the tropical curve, with blue points indicating the positions of initial unstable points.


# Manual to usual sandpiles
- To run a sandpile on the square with side=100 and 10^6 points, with initial background 3, type
python powerlaw_sandpile.py

- this produces a file in the folder classical_sandpile/


