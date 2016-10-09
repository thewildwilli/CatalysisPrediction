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

You can [download a binary release] (releases/). You will need Java 1.8 and Scala 2.11.7 installed. Extract the zip and run:

`java -jar CatalysisPrediction.jar [args]`

You can see the list of program arguments [here](wiki/Program-arguments). 

The binaries ship with some test data (`test` directory), so you can get started right away. For example, try:

```
java -jar CatalysisPrediction.jar -dir test/data/ -a 3HTB/3HTB_protein.pdb -b 3HTB/3HTB_ligand.pdb -out 3HTB/3htb_docked.mol2 -docker forcevector --ignoreAhydrogens -threshold 1.0e-5 -surface 1.4 -permeability 0.90 -balance 1,0,1,0
```

You may also want [to edit your `viewinit.txt`](wiki/viewinit.txt) file to adjust how molecules are visualised.

## Compiling

Intellij Idea with Scala plug-in was used for development and is the recommended choice. Clone the repository or download sources. You may need to check that the Scala library is correctly linked to the project.

## Libraries

- [Jmol](http://jmol.sourceforge.net/) is used for molecule 3D visualisation. 
- [Breeze](http://www.scalanlp.org/) is used for linear algebra operations.
- ThreadCSO<sup>3</sup> is used for concurrent programming.

Additionally, [SimRNA](http://genesilico.pl/software/stand-alone/simrna), [RNA Composer](http://rnacomposer.cs.put.poznan.pl/), [Make-NA](http://structure.usc.edu/make-na/server.html), and [OpenBabel](http://openbabel.org/) were used to generate some of the test data, but are not dependencies of the code.

<hr/>
<sup>1</sup> M. Steel and W. Hordijk, “Detecting autocatalytic, self-sustaining sets in chemical reaction systems,” Journal of theoretical biology,vol. 277, no. 4, pp. 451-461, 2004

<sup>2</sup> Ernesto Ocampo, "A Fast Double Docking Algorithm for Catalysis Prediction", 2016. Dissertation submitted for the degree of Master of Science in Computer Science, University of Oxford, under the supervision of Jotun Hein and Peter Jeavons.

<sup>3</sup> B. Sufrin, “Communicating Scala Objects.,” in CPA, 2008
