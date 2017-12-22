In this repository we have placed all the code and experimental data associated to the paper Self-organized criticality, pattern emergence, and Tropical Geometry by N. Kalinin, A Guzmán-Sáenz, Y. Prieto, M. Shkolnokov, V. Kalinina and E. Lupercio.

In the paper we propose a new Mathematical tool to study self-organized criticality: Tropical Geometry. As evidence for this claim, we show that the statistics of avalanches in the continuous Tropical sandpile model obeys a power law. We also explain how Tropical Geometry and this model in particular are useful in the study of proportional growth phenomena. These programs are part of work with Nikita Kalinin, Ernesto Lupercio, Mikhail Shkolnikov, Vera Kalinina, and Yulieth Prieto.

# Tropical Sandpiles
Simulations of tropical sandpiles and tropical curves.
 - parallelsandpile outputs a file called grid.dat which contains a representation of the final state of the sandpile. 
 It requires MPI to be compiled and to run it.
 - visualizegrid reads grid.dat and displays the final state of the sandpile.


# Manual to tropical (linearized) sandpile model:
- To run tropical sandpiles compile:
g++ -std=c++11 -O3 linearsandpile.cpp -o linearsandpile
- and run (with default parameters):
./linearsandpile
- this will produce a bunch of files in the folder tsandpile/
- files power1000_900_2.txt contains the avalanche sizes, with the sign '-' if the corresponding avalanche touched the boundary
- To visualize the corresponding tropical curve type:
python vizualiselinearsand.py
- this will show you a picture of the tropical curve, with blue points indicating the positions of initial unstable points.


To see power law:
in the tsandpile/power1000_100000_n.txt files the data for avalanches is stored (n is the seed = 82,83,84,85,86,88,89,90)
which gives 0.9 as the critical exponent
size = 1000, 100000 stands for the total number of experiments.
In the file we have a bit more than 10000 of them performed (for each seed).

in R do the following:

f82<-c(…copypaste from the corresponding file)
... the same for other seeds
f<-c(f82,f83,f84,f85,f86,f88,f89,f90)

g<-subset(f,f>0)
g<-subset(g,g<0.99)

r<-hist(g, breaks=100);
plot(r$breaks[-1], r$counts, log='xy', type='p',xlab='log(s)',ylab='log(frequency)')
abline(a = 2.2, b = -0.9, col = "green")
mtext('y=-0.9x+2.2',line=-5)



# Manual to usual sandpiles
- To run a sandpile on the square with side=100 and 10^6 points, with initial background 3, type
python powerlaw_sandpile.py

- this produces a file in the folder classical_sandpile/
- file readme shows how in R you can obtain that frequency = c(size of avalanche)^(-1.15).


