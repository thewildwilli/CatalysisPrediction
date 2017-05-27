package docking
import io.threadcso._
import opt.Action

/** Docking algorithms use a DockLog instance to
  * report their actions
  */
class DockLog(val log: ![Any], val enabled: Boolean = true) {
  private def send(contents: Any) =
    if (enabled) log ! contents

  def save = send("save")
  def reset = send("reset")
  def action(t: Action) = send(t)
  def other(thing: Any) = send(thing)  // used to report anything else
}

object DockLog {
  def dummy = new DockLog(null, false)
}
