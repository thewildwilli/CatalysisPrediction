# Catalysis Prediction

Catalysis is one of the main components of RAF (Reflexively
Autocatalytic and Food-generated) sets<sup>1</sup> to model pre-biotic self-sustaining
chemical systems, which are thought to be precursors of life. Yet, predicting if
a molecule will catalyse a reaction is an unsolved problem.

This project aims at predicting a particular case of catalysis: that of double docking.
If molecule A can bind to both molecules B and C, then A may bring B and C to a
favourable configuration and catalyse a reaction between them.

To this end, a fast molecular docking algorithm was developed for small molecules.
Molecular docking is a well-known problem, but only for cases when at least one of the
molecules is very large - e.g. a protein. This novel algorithm focuses on the particular
and largely unexplored case of small molecule to small molecule docking.

This project was initially developed as part of my Master's dissertation, available [here](https://1drv.ms/b/s!Ala0ZyC71s40gqpPv23y10dBuhyuxw).

## Fast Molecular Docking

At the present stage, this project provides a fast molecular docking algorithm. These videos show the docking program working:

### DNA Hexamer example
[<img src="/doc/DnaExample.png" title="DNA hexamer example" alt="DNA hexamer example video" width="40%"/>](https://www.youtube.com/watch?v=TiqyNjZQtF0)

### Protein-ligand example
Although this project does not target large molecules specificallty, experiments were carried out for comparability.
[<img src="/doc/ProteinExample.png" title="Protein-ligand example" alt="Protein-ligand example video" width="40%"/>](https://www.youtube.com/watch?v=t_llyh4zadE)


## Running the program

You can download a binary release. You will need Java 1.8 and Scala 2.11 installed.

## 

<hr/>
<sup>1</sup> M. Steel and W. Hordijk, “Detecting autocatalytic, self-sustaining sets in chemical reaction systems,” Journal of theoretical biology,vol. 277, no. 4, pp. 451-461, 2004

<sup>2</sup> Ernesto Ocampo, "A Fast Double Docking Algorithm for Catalysis Prediction", 2016. Dissertation submitted for the degree of Master of Science in Computer Science, University of Oxford, under the supervision of Jotun Hein and Peter Jeavons.
