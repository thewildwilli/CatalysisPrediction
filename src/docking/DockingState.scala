package docking
import opt.State
import model.Molecule

// Created by Ernesto on 23/05/2016.
abstract class DockingState(val a: Molecule, val b: Molecule) extends State {

}
