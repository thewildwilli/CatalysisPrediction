package testprogs

import docking.dockscore._
import docking.docksearch._
import io.Mol2Writer
import io.threadcso._
import jmolint.JmolCmds._
import jmolint.{JmolCmds, JmolFrame, JmolMoleculeReader, JmolPanel}
import model.Molecule
import opt.Action
import prog.DockMain


object H20IceTest {


  def main(args: Array[String]): Unit = {
    DockMain.parseArgs(args)
    var pathA = DockMain.DockArgs.pathA
    val pathB = DockMain.DockArgs.pathB
    for (i <- 1 until 100){
      DockMain.doMainDock(pathA, pathB, DockMain.DockArgs.pathOut)

      // Only B is written to pathOut, combine them:
      val molA: Molecule = JmolMoleculeReader.read(DockMain.jmolPanel, 0)
      val molB: Molecule = JmolMoleculeReader.read(DockMain.jmolPanel, 1)
      molA.importM(molB)

      val pathOut = DockMain.DockArgs.pathOut + s"-$i"
      new Mol2Writer(pathOut).write(molA)

      pathA = pathOut
    }

  }



}