package prog

import docking.docksearch.forcevector.ForceVectorDocker
import jmolint.{JmolFrame, JmolMoleculeReader}
import model.Molecule
import prog.DockMain.{getFFParams, parseArgs}

/**
  *  This program scores a configuration and exits
  */
object Score {
  def main(args: Array[String]): Unit = {
    val dockArgs = parseArgs(args)
    if (dockArgs.dockerName != "forcevector")
      throw new Exception("Only forcevector currently supported")

    val visible = !dockArgs.exit
    val jmolPanel = new JmolFrame(500, 500, false, visible).getPanel
    jmolPanel.openFiles(List(dockArgs.fullPathA, dockArgs.fullPathB), dockArgs.pdbAddHydrogens)

    val molA: Molecule = JmolMoleculeReader.read(jmolPanel, 0)
    val molB: Molecule = JmolMoleculeReader.read(jmolPanel, 1)

    val ffparams = getFFParams(molA, molB, dockArgs)
    val score = new ForceVectorDocker(ffparams).score(molA, molB)

    println(s"Configuration score: $score")
    if (dockArgs.exit)
      sys.exit(0)
  }
}
