/*
 * $Id$
 * $URL$
 * ---------------------------------------------------------------------
 * This file is part of SBVC, the systems biology visualizer and
 * converter. This tools is able to read a plethora of systems biology
 * file formats and convert them to an internal data structure.
 * These files can then be visualized, either using a simple graph
 * (KEGG-style) or using the SBGN-PD layout and rendering constraints.
 * Some currently supported IO formats are SBML (+qual, +layout), KGML,
 * BioPAX, SBGN, etc. Please visit the project homepage at
 * <http://www.cogsys.cs.uni-tuebingen.de/software/SBVC> to obtain the
 * latest version of SBVC.
 *
 * Copyright (C) 2012-2012 by the University of Tuebingen, Germany.
 *
 * SBVC is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation. A copy of the license
 * agreement is provided in the file named "LICENSE.txt" included with
 * this software distribution and also available online as
 * <http://www.gnu.org/licenses/lgpl-3.0-standalone.html>.
 * ---------------------------------------------------------------------
 */
package de.zbit.sbvc.gui;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.prefs.BackingStoreException;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JTabbedPane;
import javax.swing.JToolBar;
import javax.swing.KeyStroke;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.sbml.jsbml.SBMLDocument;

import de.zbit.AppConf;
import de.zbit.graph.RestrictedEditMode;
import de.zbit.graph.gui.TranslatorPanel;
import de.zbit.gui.BaseFrame;
import de.zbit.gui.GUIOptions;
import de.zbit.gui.GUITools;
import de.zbit.gui.JLabeledComponent;
import de.zbit.gui.JTabbedLogoPane;
import de.zbit.gui.actioncommand.ActionCommand;
import de.zbit.gui.prefs.FileSelector;
import de.zbit.gui.prefs.PreferencesPanel;
import de.zbit.io.filefilter.SBFileFilter;
import de.zbit.kegg.Translator;
import de.zbit.kegg.io.KEGG2jSBML;
import de.zbit.kegg.io.KEGGtranslatorIOOptions.Format;
import de.zbit.sbml.io.OpenedFile;
import de.zbit.sbvc.SBVC;
import de.zbit.sbvc.io.SBVCIOOptions;
import de.zbit.util.StringUtil;
import de.zbit.util.prefs.SBPreferences;
import de.zbit.util.prefs.SBProperties;
import de.zbit.util.progressbar.AbstractProgressBar;

/**
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class SBVCUI extends BaseFrame implements ActionListener, KeyListener, ItemListener {

  public static enum Action implements ActionCommand {
    // /**
    // * {@link Action} for LaTeX export.
    // */
    // TO_LATEX,
    /**
     * {@link Action} for new instances of {@link AbstractProgressBar}s. This
     * will display the progress in the {@link #statusBar}.
     */
    NEW_PROGRESSBAR,
    /**
     * {@link Action} for downloading KGMLs.
     */
    DOWNLOAD_PATHWAY,
    /**
     * This is coming from {@link RestrictedEditMode#OPEN_PATHWAY} and must be
     * renamed accordingly. The source is a kegg pathway id that should be
     * opened as new tab, when this action is fired. This is an invisible
     * action.
     */
    OPEN_PATHWAY,
    /**
     * Invisible {@link Action} that should be performed, whenever an
     * translation is done.
     */
    TRANSLATION_DONE;
    /**
     * Invisible {@link Action} that should be performed, whenever a file has
     * been dropped on this panel.
     */
    // FILE_DROPPED

    /*
     * (non-Javadoc)
     * 
     * @see de.zbit.gui.ActionCommand#getName()
     */
    public String getName() {
      switch (this) {
        // case TO_LATEX:
        // return "Export to LaTeX";
        case DOWNLOAD_PATHWAY:
          return "Download pathway";

        default:
          return StringUtil.firstLetterUpperCase(toString().toLowerCase().replace('_', ' '));
      }
    }

    /*
     * (non-Javadoc)
     * 
     * @see de.zbit.gui.ActionCommand#getToolTip()
     */
    public String getToolTip() {
      switch (this) {
        // case TO_LATEX:
        // return "Converts the currently opened model to a LaTeX report file.";
        case DOWNLOAD_PATHWAY:
          return "Downloads pathways from online servers.";
        default:
          return "";
      }
    }
  }

  /**
   * Generated serial version identifier.
   */
  private static final long serialVersionUID = 6631262606716052915L;

  //TODO: icon erstellen
//  static {
//    
//    String iconPaths[] = { "SBVCIcon_16.png", "SBVCIcon_32.png",
//        "SBVCIcon_48.png", "SBVCIcon_128.png", "SBVCIcon_256.png" };
//    for (String path : iconPaths) {
//      URL url = TranslatorUI.class.getResource("img/" + path);
//      if (url != null) {
//        UIManager.put(path.substring(0, path.lastIndexOf('.')), new ImageIcon(url));
//      }
//    }
//    try {
//      org.sbml.tolatex.gui.LaTeXExportDialog.initImages();
//    } catch (Throwable t) {
//      // Also allow SBVC to compile without
//      // SBML2LaTeX !
//    }
//  }

  /**
   * Default directory path's for saving and opening files. Only init them once.
   * Other classes should use these variables.
   */
  public static String openDir, saveDir;
  /**
   * This is where we place all the converted models.
   */
  private JTabbedPane tabbedPane;
  /**
   * preferences is holding all project specific preferences
   */
  private SBPreferences prefsIO;

  /**
 * 
 */
  public SBVCUI() {
    this(null);
  }

  public SBVCUI(AppConf appConf) {
    super(appConf);

    // init preferences
    initPreferences();
    File file = new File(prefsIO.get(SBVCIOOptions.INPUT));
    openDir = file.isDirectory() ? file.getAbsolutePath() : file.getParent();
    file = new File(prefsIO.get(SBVCIOOptions.OUTPUT));
    saveDir = file.isDirectory() ? file.getAbsolutePath() : file.getParent();

    //TODO: Icon erstellen
    // Depending on the current OS, we should add the following image
    // icons: 16x16, 32x32, 48x48, 128x128 (MAC), 256x256 (Vista).
//    int[] resolutions = new int[] { 16, 32, 48, 128, 256 };
//    List<Image> icons = new LinkedList<Image>();
//    for (int res : resolutions) {
//      Object icon = UIManager.get("SBVCIcon_" + res);
//      if ((icon != null) && (icon instanceof ImageIcon)) {
//        icons.add(((ImageIcon) icon).getImage());
//      }
//    }
//    setIconImages(icons);
  }

  /**
   * Init preferences, if not already done.
   */
  private void initPreferences() {
    if (prefsIO == null) {
      prefsIO = SBPreferences.getPreferencesFor(SBVCIOOptions.class);
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#createJToolBar()
   */
  protected JToolBar createJToolBar() {
    initPreferences();
    // final JPanel r = new JPanel(new VerticalLayout());
    final JToolBar r = new JToolBar("Translate new file", JToolBar.HORIZONTAL);

    JComponent jc = PreferencesPanel.createJComponentForOption(SBVCIOOptions.INPUT,
        prefsIO, this);
    // Allow a change of Focus (important!)
    if (jc instanceof FileSelector)
      ((FileSelector) jc).removeInputVerifier();
    r.add(jc);
    r.add(PreferencesPanel.createJComponentForOption(SBVCIOOptions.FORMAT, prefsIO, this));

    // Button and action
    JButton ok = new JButton("Translate now!", UIManager.getIcon("ICON_GEAR_16"));
    ok.setToolTipText(StringUtil
        .toHTMLToolTip("Starts the conversion of the input file to the selected output format and displays the result on this workbench."));
    ok.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        // Get selected file and format
        File inFile = getInputFile(r);
        String format = getOutputFileFormat(r);

        // Translate
        createNewTab(inFile, format);
      }
    });
    r.add(ok);

    GUITools.setOpaqueForAllElements(r, false);
    return r;
  }

  /**
   * Searches for any JComponent with "TranslatorOptions.FORMAT.getOptionName()"
   * on it and returns the selected format. Use it e.g. with
   * {@link #translateToolBar}.
   * 
   * @param r
   * @return String - format.
   */
  private String getOutputFileFormat(JComponent r) {
    String format = null;
    for (Component c : r.getComponents()) {
      if (c.getName() == null) {
        continue;
      } else if (c.getName().equals(SBVCIOOptions.FORMAT.getOptionName())
          && (JLabeledComponent.class.isAssignableFrom(c.getClass()))) {
        format = ((JLabeledComponent) c).getSelectedItem().toString();
        break;
      }
    }
    return format;
  }

  /**
   * Searches for any JComponent with "TranslatorOptions.INPUT.getOptionName()"
   * on it and returns the selected file. Use it e.g. with
   * {@link #translateToolBar}.
   * 
   * @param r
   * @return File - input file.
   */
  private File getInputFile(JComponent r) {
    File inFile = null;
    for (Component c : r.getComponents()) {
      if (c.getName() == null) {
        continue;
      } else if (c.getName().equals(SBVCIOOptions.INPUT.getOptionName())
          && (FileSelector.class.isAssignableFrom(c.getClass()))) {
        try {
          inFile = ((FileSelector) c).getSelectedFile();
        } catch (IOException e1) {
          GUITools.showErrorMessage(r, e1);
          e1.printStackTrace();
        }
      }
    }
    return inFile;
  }

  /**
   * A method to set the value of a currently displayed {@link FileSelector}
   * corresponding to the <code>SBVCIOOptions.INPUT</code> option.
   * 
   * @param r
   *          the JComponent on which the component for the mentioned optioned
   *          is placed.
   * @param file
   *          the file to set
   */
  private void setInputFile(JComponent r, File file) {
    for (Component c : r.getComponents()) {
      if (c.getName() == null) {
        continue;
      } else if (c.getName().equals(SBVCIOOptions.INPUT.getOptionName())
          && (FileSelector.class.isAssignableFrom(c.getClass()))) {
        ((FileSelector) c).setSelectedFile(file);
      }
    }
  }

  /**
   * Translate and create a new tab.
   * 
   * @param inFile
   * @param format
   */
  private void createNewTab(File inFile, String format) {
    // Check input
    if (!SBVCIOOptions.INPUT.getRange().isInRange(inFile)) {
      String message = "The given file is no valid input file.";
      if (inFile != null) {
        message = '\'' + inFile.getName() + "' is no valid input file.";
      }
      JOptionPane.showMessageDialog(this, message, System.getProperty("app.name"),
          JOptionPane.WARNING_MESSAGE);
    } else {
      Format f = null;
      try {
        f = Format.valueOf(format);
      } catch (Throwable exc) {
        exc.printStackTrace();
        JOptionPane.showMessageDialog(this, '\'' + format + "' is no valid output format.",
            System.getProperty("app.name"), JOptionPane.WARNING_MESSAGE);
      }
      if (f != null) {
        // Tanslate and add tab.
        try {
          openDir = inFile.getParent();
          addSBVCTab(SBVCPanelTools.createPanel(inFile, f, this));
        } catch (Exception e1) {
          GUITools.showErrorMessage(this, e1);
        }
      }
    }
  }

  /**
   * Adds a new {@link TranslatorPanel} to this {@link #tabbedPane} and changes
   * the selection to this new panel.
   * 
   * @param tp
   */
  public void addSBVCTab(TranslatorPanel<?> tp) {
    addSBVCTab(null, tp);
  }

  /**
   * Adds a new {@link TranslatorPanel} to this {@link #tabbedPane} and changes
   * the selection to this new panel.
   * 
   * @param tabName
   *          name for the tab
   * @param tp
   */
  public void addSBVCTab(String tabName, TranslatorPanel<?> tp) {
    try {
      if (tabName == null)
        tabName = tp.getTitle();
      tabbedPane.addTab(tabName, tp);
      tabbedPane.setSelectedIndex(tabbedPane.getTabCount() - 1);
    } catch (Exception e1) {
      GUITools.showErrorMessage(this, e1);
    }
  }

  /*
   * (non-Javadoc)
   * 
   * @see
   * java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
   */
  public void actionPerformed(ActionEvent e) {
    try {
      Action action = Action.valueOf(e.getActionCommand());
      switch (action) {
        case TRANSLATION_DONE:
          TranslatorPanel<?> source = (TranslatorPanel<?>) e.getSource();
          int index = tabbedPane.indexOfComponent(source);
          if (index >= 0) {// ELSE: User closed the tab before completion
            if (e.getID() != JOptionPane.OK_OPTION) {
              // If translation failed, remove the tab. The error
              // message has already been issued by the translator.
              tabbedPane.removeTabAt(index);
            } else {
              // Do not change title here. Initial title is mostly
              // better than this one ;-)
              // tabbedPane.setTitleAt(index, source.getTitle());
            }
          }
          getStatusBar().hideProgress();
          updateButtons();
          break;
        /*
         * Moved to BaseFrame. case FILE_DROPPED: String format =
         * getOutputFileFormat(toolBar); if ((format == null) ||
         * (format.length() < 1)) { break; } createNewTab(((File)
         * e.getSource()), format); break;
         */
        // case TO_LATEX:
        // writeLaTeXReport();
        // break;
        case DOWNLOAD_PATHWAY:
          TranslatePathwayDialog.showAndEvaluateDialog(tabbedPane, this, (Format) null);           
          break;
        case NEW_PROGRESSBAR:
          getStatusBar().showProgress((AbstractProgressBar) e.getSource());
          break;
        case OPEN_PATHWAY:
          try {
            addSBVCTab(e.getSource().toString(),
                SBVCPanelTools.createPanel(e.getSource().toString(), Format.GraphML, this));
          } catch (Exception e1) {
            GUITools.showErrorMessage(this, e1);
          }
          break;
        default:
          System.out.println(action);
          break;
      }
    } catch (Throwable exc) {
      GUITools.showErrorMessage(this, exc);
    }
  }

  // /**
  // * @param object
  // */
  // private void writeLaTeXReport() {
  // TranslatorPanel o = getCurrentlySelectedPanel();
  // if (o != null) {
  // o.writeLaTeXReport(null);
  // }
  // }

  /**
   * Closes the tab at the specified index.
   * 
   * @param index
   * @return true, if the tab has been closed.
   */
  private boolean closeTab(int index) {
    if ((index < 0) || (index >= tabbedPane.getTabCount())) {
      return false;
    }
    Component comp = tabbedPane.getComponentAt(index);
    String title = tabbedPane.getTitleAt(index);
    if ((title == null) || (title.length() < 1)) {
      title = "the currently selected document";
    }

    // Check if document already has been saved
    if ((comp instanceof TranslatorPanel<?>) && !((TranslatorPanel<?>) comp).isSaved()) {
      if ((JOptionPane.YES_OPTION != JOptionPane.showConfirmDialog(this, StringUtil.toHTMLToolTip(
          "Do you really want to close %s%s%s without saving?", KEGG2jSBML.quotStart, title,
          KEGG2jSBML.quotEnd), "Close selected document", JOptionPane.YES_NO_OPTION))) {
        return false;
      }
    }

    // Close the document.
    tabbedPane.removeTabAt(index);
    updateButtons();
    return true;
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#openFile(java.io.File[])
   */
  public File[] openFile(File... files) {
    boolean askOutputFormat = false;

    // Ask input file
    if ((files == null) || (files.length < 1)) {
      files = GUITools.openFileDialog(this, openDir, false, true, JFileChooser.FILES_ONLY,
          SBFileFilter.createAllFileFilter()); //TODO: change file filter

      askOutputFormat = true;
    }
    if ((files == null) || (files.length < 1)) {
      return files;
    } else {
      // Set value to box if it does not yet contain a valid value
      if (getInputFile(toolBar) == null) {
        setInputFile(toolBar, files[0]);
      }
    }

    // Ask output format
    String format = getOutputFileFormat(toolBar);
    if (askOutputFormat || (format == null) || (format.length() < 1)) {
      JLabeledComponent outputFormat = (JLabeledComponent) PreferencesPanel
          .createJComponentForOption(SBVCIOOptions.FORMAT, prefsIO, null);
      outputFormat.setTitle("Please select the output format");
      JOptionPane.showMessageDialog(this, outputFormat, System.getProperty("app.name"),
          JOptionPane.QUESTION_MESSAGE);
      format = outputFormat.getSelectedItem().toString();
    }

    // Translate
    for (File f : files) {
      createNewTab(f, format);
    }
    return files;
  }

  /**
   * Enables and disables buttons in the menu, depending on the current tabbed
   * pane content.
   */
  private void updateButtons() {
    GUITools.setEnabled(false, getJMenuBar(), BaseAction.FILE_SAVE_AS,
    // Action.TO_LATEX,
        BaseAction.FILE_CLOSE);
    TranslatorPanel<?> o = getCurrentlySelectedPanel();
    if (o != null) {
      o.updateButtons(getJMenuBar());
    }
  }

  /**
   * @return the currently selected TranslatorPanel from the {@link #tabbedPane}
   *         , or null if either no or no valid selection exists.
   */
  private TranslatorPanel<?> getCurrentlySelectedPanel() {
    if ((tabbedPane == null) || (tabbedPane.getSelectedIndex() < 0)) {
      return null;
    }
    Object o = ((JTabbedPane) tabbedPane).getSelectedComponent();
    if ((o == null) || !(o instanceof TranslatorPanel<?>)) {
      return null;
    }
    return ((TranslatorPanel<?>) o);
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#saveFile()
   */
  public File saveFile() {
    TranslatorPanel<?> o = getCurrentlySelectedPanel();
    if (o != null) {
      return o.saveToFile();
    }
    return null;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.awt.event.KeyListener#keyPressed(java.awt.event.KeyEvent)
   */
  public void keyPressed(KeyEvent e) {
    // Preferences for the "input file"
    PreferencesPanel.setProperty(prefsIO, e.getSource(), true);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.awt.event.KeyListener#keyReleased(java.awt.event.KeyEvent)
   */
  public void keyReleased(KeyEvent e) {
    // Preferences for the "input file"
    PreferencesPanel.setProperty(prefsIO, e.getSource(), true);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.awt.event.KeyListener#keyTyped(java.awt.event.KeyEvent)
   */
  public void keyTyped(KeyEvent e) {
    // Preferences for the "input file"
    PreferencesPanel.setProperty(prefsIO, e.getSource(), true);
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.awt.event.ItemListener#itemStateChanged(java.awt.event.ItemEvent)
   */
  public void itemStateChanged(ItemEvent e) {
    // Preferences for the "output format"
    PreferencesPanel.setProperty(prefsIO, e.getSource(), true);
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#additionalFileMenuItems()
   */
  @Override
  protected JMenuItem[] additionalFileMenuItems() {
    return new JMenuItem[] {
    /*
     * GUITools.createJMenuItem(this, Action.TO_LATEX,
     * UIManager.getIcon("ICON_LATEX_16"), KeyStroke .getKeyStroke('E',
     * InputEvent.CTRL_DOWN_MASK), 'E', false),
     */
    GUITools.createJMenuItem(this, Action.DOWNLOAD_PATHWAY, UIManager.getIcon("ICON_GEAR_16"),
        KeyStroke.getKeyStroke('D', InputEvent.CTRL_DOWN_MASK), 'D', true) };
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#closeFile()
   */
  public boolean closeFile() {
    if (tabbedPane.getSelectedIndex() < 0) {
      return false;
    }
    return closeTab(tabbedPane.getSelectedIndex());
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#createMainComponent()
   */
  protected Component createMainComponent() {
    // If you encounter an exception here PUT THE RESOURCES FOLDER ON YOUR CLASS
    // PATH!
    //TODO: create Icon
    ImageIcon logo = new ImageIcon(SBVCUI.class.getResource("img/logo_watermark.png"));

    // Crop animated loading bar from image.
    // logo.setImage(ImageTools.cropImage(logo.getImage(), 0, 0,
    // logo.getIconWidth(), logo.getIconHeight()-30));

    // Create the tabbed pane, with the SBVC logo.
    tabbedPane = new JTabbedLogoPane(logo);
    // Change active buttons, based on selection.
    tabbedPane.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
        updateButtons();
      }
    });
    return tabbedPane;
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#exit()
   */
  public void exit() {
    // Close all tab. If user want's to save a tab first, cancel the closing
    // process.
    while (tabbedPane.getTabCount() > 0) {
      if (!closeTab(0)) {
        return;
      }
    }

    // Close the app and save caches.
    setVisible(false);
    try {
      Translator.saveCache();

      SBProperties props = new SBProperties();
      { // Save SBVCIOOptions
        File f = getInputFile(toolBar);
        if (f != null && SBVCIOOptions.INPUT.getRange().isInRange(f, props)) {
          props.put(SBVCIOOptions.INPUT, f);
        }
        props.put(SBVCIOOptions.FORMAT, getOutputFileFormat(toolBar));
        SBPreferences.saveProperties(SBVCIOOptions.class, props);
      }
      props.clear();

      { // Save GUIOptions
        if (openDir != null && openDir.length() > 1) {
          props.put(GUIOptions.OPEN_DIR, openDir);
        }
        if (saveDir != null && saveDir.length() > 1) {
          props.put(GUIOptions.SAVE_DIR, saveDir);
        }
        if (props.size() > 0) {
          SBPreferences.saveProperties(GUIOptions.class, props);
        }
      }

    } catch (BackingStoreException exc) {
      exc.printStackTrace();
      // Unimportant error... don't bother the user here.
      // GUITools.showErrorMessage(this, exc);
    }
    dispose();
    System.exit(0);
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#getURLAboutMessage()
   */
  public URL getURLAboutMessage() {
    // return getClass().getResource("../html/about.html");
    // "../" does not work inside a jar.
    return SBVC.class.getResource("html/about.html");
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#getURLLicense()
   */
  public URL getURLLicense() {
    // return getClass().getResource("../html/license.html");
    // "../" does not work inside a jar.
    return SBVC.class.getResource("html/license.html");
  }

  /*
   * (non-Javadoc)
   * 
   * @see de.zbit.gui.BaseFrame#getURLOnlineHelp()
   */
  public URL getURLOnlineHelp() {
    // return getClass().getResource("../html/help.html");
    // "../" does not work inside a jar.
    return SBVC.class.getResource("html/help.html");
  }

  /* (non-Javadoc)
   * @see de.zbit.gui.BaseFrame#openFile(de.zbit.sbml.io.OpenedFile<org.sbml.jsbml.SBMLDocument>[])
   */
  @Override
  protected OpenedFile<SBMLDocument>[] openFile(OpenedFile<SBMLDocument>... files) {
    // TODO Auto-generated method stub
    return null;
  }

  /* (non-Javadoc)
   * @see de.zbit.gui.BaseFrame#saveFileAs()
   */
  @Override
  public File saveFileAs() {
    // TODO Auto-generated method stub
    return null;
  }

}
