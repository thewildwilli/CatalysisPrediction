package testprogs

import io.Mol2Writer
import jmolint.{JmolFrame, JmolMoleculeReader}

object Mol2WriteTest {
  val frame = new JmolFrame(1, 1, false)
  frame.setVisible(false)

  val jmolPanel = frame.getPanel

  def main(args: Array[String]): Unit = {
    jmolPanel.openFiles(List(args(0)))
    val mol = JmolMoleculeReader.read(jmolPanel, 0)
    new Mol2Writer(args(1)).write(mol)
    sys.exit(0)
  }

}
