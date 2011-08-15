/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * JDialogXmippFilesList.java
 *
 * Created on Aug 2, 2011, 5:56:47 PM
 */
package browser.filebrowsers;

import browser.ICONS_MANAGER;
import browser.LABELS;
import ij.IJ;
import java.awt.BorderLayout;
import java.io.OutputStreamWriter;
import java.net.InetAddress;
import java.net.Socket;

/**
 *
 * @author Juanjo Vega
 */
public class JDialogXmippFilesList extends javax.swing.JDialog {//Frame {

    final static String EOT = "__END__";
    JPanelXmippBrowser panelXmippBrowser;
    int port;

    /** Creates new form JDialogXmippFilesList */
    public JDialogXmippFilesList(String directory, int port) {
        this(directory, port, false);
    }

    public JDialogXmippFilesList(String directory, int port, String expression) {
        this(directory, port, expression, false);
    }

    public JDialogXmippFilesList(String directory, int port, boolean singleSelection) {
        this(directory, port, "", singleSelection);
    }

    public JDialogXmippFilesList(String directory, int port, String expression, boolean singleSelection) {
        super();//

        this.port = port;

        setTitle(LABELS.TITLE_XMIPP_FILE_SELECTOR);

        initComponents();

        panelXmippBrowser = new JPanelXmippBrowser(directory, expression);
        panelXmippBrowser.setSingleSelection(singleSelection);

        add(panelXmippBrowser, BorderLayout.CENTER);

        pack();
        setLocationRelativeTo(null);
    }

    protected void button1Clicked() {
        if (sendSelectedFiles()) {
            dispose();
        }
    }

    protected boolean sendSelectedFiles() {
        return send(panelXmippBrowser.getSelectedValues());
    }

    protected boolean send(Object items[]) {
        try {
            Socket socket = new Socket(InetAddress.getByName("127.0.0.1"), port);

            // Get streams.
            OutputStreamWriter output = new OutputStreamWriter(socket.getOutputStream());
            output.flush();

            if (items != null) {
                for (int i = 0; i < items.length; i++) {
                    output.write(items[i].toString() + "\n");
                }
            }

            output.write(EOT);
            output.flush();

            // Closes connection.
            output.close();
            socket.close();

            return true;
        } catch (Exception ex) {
            IJ.error(ex.getMessage());
        }

        return false;
    }

    protected void button2Clicked() {
        cancel();
    }

    void cancel() {
        if (send(null)) {
            dispose();
        }
    }

    protected void goParent() {
        panelXmippBrowser.goParent();
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jToolBar = new javax.swing.JToolBar();
        jbParent = new javax.swing.JButton();
        jbRefresh = new javax.swing.JButton();
        jpButtons = new javax.swing.JPanel();
        jbOk = new javax.swing.JButton();
        jbCancel = new javax.swing.JButton();

        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });

        jToolBar.setFloatable(false);
        jToolBar.setRollover(true);

        jbParent.setIcon(ICONS_MANAGER.PARENT_DIRECTORY);
        jbParent.setText(LABELS.BUTTON_PARENT_DIRECTORY);
        jbParent.setFocusable(false);
        jbParent.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbParent.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbParent.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbParentActionPerformed(evt);
            }
        });
        jToolBar.add(jbParent);

        jbRefresh.setIcon(ICONS_MANAGER.REFRESH_DIRECTORY);
        jbRefresh.setText(LABELS.BUTTON_REFRESH_DIRECTORY);
        jbRefresh.setFocusable(false);
        jbRefresh.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jbRefresh.setVerticalTextPosition(javax.swing.SwingConstants.BOTTOM);
        jbRefresh.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbRefreshActionPerformed(evt);
            }
        });
        jToolBar.add(jbRefresh);

        getContentPane().add(jToolBar, java.awt.BorderLayout.NORTH);

        jbOk.setText(LABELS.BUTTON_OK);
        jbOk.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbOkActionPerformed(evt);
            }
        });
        jpButtons.add(jbOk);

        jbCancel.setText(LABELS.BUTTON_CANCEL);
        jbCancel.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jbCancelActionPerformed(evt);
            }
        });
        jpButtons.add(jbCancel);

        getContentPane().add(jpButtons, java.awt.BorderLayout.SOUTH);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jbRefreshActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbRefreshActionPerformed
        panelXmippBrowser.refreshCurrentDirectory();
}//GEN-LAST:event_jbRefreshActionPerformed

    private void jbOkActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbOkActionPerformed
        button1Clicked();
}//GEN-LAST:event_jbOkActionPerformed

    private void jbCancelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbCancelActionPerformed
        button2Clicked();
}//GEN-LAST:event_jbCancelActionPerformed

    private void jbParentActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jbParentActionPerformed
        goParent();
    }//GEN-LAST:event_jbParentActionPerformed

    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        cancel();
    }//GEN-LAST:event_formWindowClosing
    // Variables declaration - do not modify//GEN-BEGIN:variables
    javax.swing.JToolBar jToolBar;
    javax.swing.JButton jbCancel;
    javax.swing.JButton jbOk;
    javax.swing.JButton jbParent;
    javax.swing.JButton jbRefresh;
    private javax.swing.JPanel jpButtons;
    // End of variables declaration//GEN-END:variables
}
