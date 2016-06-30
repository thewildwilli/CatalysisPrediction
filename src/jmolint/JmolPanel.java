package jmolint;// Created by Ernesto on 15/06/2016.

import javax.swing.JPanel;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;

import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolSimpleViewer;
import org.jmol.adapter.smarter.SmarterJmolAdapter;

public class JmolPanel extends JPanel {
    private static final long serialVersionUID = -3661941083797644242L;
    private JmolSimpleViewer viewer;
    private JmolAdapter adapter;
    public JmolPanel() {
        adapter = new SmarterJmolAdapter();
        viewer = JmolSimpleViewer.allocateSimpleViewer(this, adapter);
    }

    public JmolSimpleViewer getViewer() {
        return viewer;
    }

    public void executeCmd(String... cmds){
        for (String cmd: cmds)
            viewer.evalString(cmd);
    }

    private final Dimension currentSize = new Dimension();
    private final Rectangle rectClip = new Rectangle();

    public void paint(Graphics g) {
        getSize(currentSize);
        g.getClipBounds(rectClip);
        viewer.renderScreenImage(g, currentSize.width, currentSize.height);
    }
}
