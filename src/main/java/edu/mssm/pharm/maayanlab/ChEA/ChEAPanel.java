package edu.mssm.pharm.maayanlab.ChEA;

import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import edu.mssm.pharm.maayanlab.common.core.FileUtils;
import edu.mssm.pharm.maayanlab.common.core.SettingsChanger;
import edu.mssm.pharm.maayanlab.common.swing.FileDrop;
import edu.mssm.pharm.maayanlab.common.swing.UIUtils;

public class ChEAPanel extends JPanel {

	private static final long serialVersionUID = -1006178583380192229L;
	
	static Logger log = Logger.getLogger(ChEAPanel.class.getSimpleName());
	
	// JPanels
	private JPanel panel;
	
	// UI elements
	private JFileChooser openChooser, saveChooser;
	private JTextField openPath, savePath;
	private JTextArea inputTextArea, outputTextArea;
	private JButton openButton;
	private JComboBox dbCombo, bgCombo, sortByCombo;
	private JTextField selectTopText;
	
	// Output
	private String output;
	
	public static void main(String[] args) {

		if (args.length == 0) {			
			// Schedule a job for the EDT
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					createAndShowGUI();
				}
			});
		}
		// else{
		// 	ChEA.main(args);
		// }
	}
	
	private static void createAndShowGUI() {
		// Try to use Nimbus look and feel
		try {            
            UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
        } catch (Exception e) {
           log.warning("Nimbus: " + e);
        }
        
        // Create and set up the window
        JFrame appFrame = new JFrame("ChEA - ChIP Enrichment Analysis");
        appFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        // Add content to the window
        ChEAPanel appPanel = new ChEAPanel();
        appFrame.setContentPane(appPanel);
        
        // Display the window
        appFrame.setResizable(false);
        appFrame.pack();
        appFrame.setVisible(true);
	}
	
	public ChEAPanel() {
		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		if (!Boolean.getBoolean("verbose"))
            log.setLevel(Level.WARNING);
		
		// Attach instance to variable so nested classes can reference it
		panel = this;
		
		// File choosers
		openChooser = new JFileChooser(System.getProperty("user.dir"));
		openChooser.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				File file = openChooser.getSelectedFile();
				if (file.canRead() && e.getActionCommand().equals(JFileChooser.APPROVE_SELECTION))
					setupIO(file);
			}
		});
		saveChooser = new JFileChooser(System.getProperty("user.dir"));
		saveChooser.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				File file = saveChooser.getSelectedFile();
				if (file != null && e.getActionCommand().equals(JFileChooser.APPROVE_SELECTION)) {
					if (!file.getName().endsWith(".csv")) {
						file = new File(file.getAbsolutePath() + ".csv");
						saveChooser.setSelectedFile(file);
					}
					
					savePath.setText(file.getAbsolutePath());
				}
			}
		});
		
		// Select input/output file button
		JButton openFileButton = new JButton("Input Genes");
		openFileButton.setPreferredSize(new Dimension(300, 30));
		openFileButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				openChooser.showOpenDialog(panel);
			}
		});
		JButton saveFileButton = new JButton("Output TFs");
		saveFileButton.setPreferredSize(new Dimension(300, 30));
		saveFileButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveChooser.showSaveDialog(panel);
			}
		});
		
		// Text Fields
		openPath = new JTextField();
		savePath = new JTextField();
		
		// File Drop
		new FileDrop(openPath, new FileDrop.Listener() {
			public void filesDropped(File[] files) {
				if (files[0].canRead()) {
					setupIO(files[0]);
					openChooser.setSelectedFile(files[0]);
				}
			}
		});
		
		// Scroll panes
		inputTextArea = new JTextArea(20, 20);
		JScrollPane inputTextPane = new JScrollPane(inputTextArea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		outputTextArea = new JTextArea(20, 20);
		JScrollPane outputTextPane = new JScrollPane(outputTextArea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		
		// File Drop
		new FileDrop(inputTextArea, new FileDrop.Listener() {
			public void filesDropped(File[] files) {
				if (files[0].canRead()) {
					setupIO(files[0]);
					openChooser.setSelectedFile(files[0]);
				}
			}
		});
		
		// Open results
		openButton = new JButton("View Results");
		openButton.setEnabled(false);
		openButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					Desktop.getDesktop().open(new File(output));
				} catch (Exception e1) {
					JOptionPane.showMessageDialog(panel, "Unable to open " + output, "Unable to open file", JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		
		// Start button
		JButton runButton = new JButton("Find TFs");
		runButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				output = savePath.getText();
				ArrayList<String> inputList = UIUtils.getTextAreaText(inputTextArea); 
				
				try {
					if (!output.equals("") && FileUtils.validateList(inputList)) {
						ChEA chea = new ChEA();
						
						setSettings(chea);
						chea.run(inputList);
						UIUtils.setTextAreaText(outputTextArea, chea.getTopRankedList(Integer.parseInt(selectTopText.getText())));
						chea.writeFile(output);
						enableOutput(output);
					}
					else {
						JOptionPane.showMessageDialog(panel, "No save location specified.", "No Save Location", JOptionPane.WARNING_MESSAGE);
					}
				} catch (ParseException e1) {
					if (e1.getErrorOffset() == -1)
						JOptionPane.showMessageDialog(panel, "Input list is empty.", "Invalid Input", JOptionPane.WARNING_MESSAGE);
					else
						JOptionPane.showMessageDialog(panel, e1.getMessage() + " at line " + (e1.getErrorOffset() + 1) +" is not a valid Entrez Gene Symbol.", "Invalid Input", JOptionPane.WARNING_MESSAGE);
				}
			}
		});
		
		// Advanced Settings:
		// Select between transcription factor database
		JLabel dbLabel = new JLabel("Use transcription factor database derived from");
		String[] databases = {"Chip-X", "PWMs", "Genome Browser"};
		dbCombo = new JComboBox(databases);
		dbCombo.setSelectedIndex(0);
		dbCombo.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (dbCombo.getSelectedIndex() == 2) {
					bgCombo.setSelectedIndex(2);
					bgCombo.setEnabled(false);
				}
				else {
					bgCombo.setSelectedIndex(0);
					bgCombo.setEnabled(true);
				}
			}
		});
		JPanel dbBox = new JPanel();
		dbBox.add(dbLabel);
		dbBox.add(dbCombo);
		
		// Select between aggregate, mouse, and human backgrounds to use
		JLabel bgLabel1 = new JLabel("Use");
		String[] backgrounds = {"mouse + human", "mouse", "human"};
		bgCombo = new JComboBox(backgrounds);
		bgCombo.setSelectedIndex(0);
		JLabel bgLabel2 = new JLabel("as the background organism");
		JPanel bgBox = new JPanel();
		bgBox.add(bgLabel1);
		bgBox.add(bgCombo);
		bgBox.add(bgLabel2);
		
		// Sort results by p-value, combined score, or rank
		JLabel sortByLabel = new JLabel("Sort by");
		String[] options = {"p-value", "rank", "combined score"};
		sortByCombo = new JComboBox(options);
		sortByCombo.setSelectedIndex(2);
		JPanel sortByBox = new JPanel();
		sortByBox.add(sortByLabel);
		sortByBox.add(sortByCombo);
		
		// Select top x number of TFs for next step
		JLabel selectTopLabel1 = new JLabel("Select top");
		selectTopText = new JTextField("10");
		JLabel selectTopLabel2 = new JLabel("transcription factors");
		JPanel selectTopBox = new JPanel();
		selectTopBox.add(selectTopLabel1);
		selectTopBox.add(selectTopText);
		selectTopBox.add(selectTopLabel2);
		
		// Input and output box
		JPanel ioBox = new JPanel();
		ioBox.setLayout(new GridLayout(2,2));
		ioBox.add(openFileButton);
		ioBox.add(saveFileButton);
		ioBox.add(openPath);	
		ioBox.add(savePath);
		
		// Panes
		JPanel textPanesBox = new JPanel();
		textPanesBox.setLayout(new BoxLayout(textPanesBox, BoxLayout.LINE_AXIS));
		textPanesBox.add(inputTextPane);
		textPanesBox.add(outputTextPane);
		
		// Button box
		JPanel buttonBox = new JPanel();
		buttonBox.setLayout(new BoxLayout(buttonBox, BoxLayout.LINE_AXIS));
		buttonBox.add(runButton);
		buttonBox.add(openButton);
		
		// Advanced settings box
		JPanel advancedSettingsBox = new JPanel();
		advancedSettingsBox.setLayout(new BoxLayout(advancedSettingsBox, BoxLayout.PAGE_AXIS));
		advancedSettingsBox.setBorder(BorderFactory.createTitledBorder("Advanced Settings"));
		advancedSettingsBox.add(dbBox);
		advancedSettingsBox.add(bgBox);
		advancedSettingsBox.add(sortByBox);
		advancedSettingsBox.add(selectTopBox);
		
		// Add all the panels together
		this.add(ioBox);
		this.add(textPanesBox);
		this.add(Box.createRigidArea(new Dimension(0,10)));
		this.add(buttonBox);
		this.add(advancedSettingsBox);
	}
	
	// Use interface to set settings of applications
	public void setSettings(SettingsChanger changer) {
		switch (dbCombo.getSelectedIndex()) {
		case 0: changer.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.CHIPX); break;
		case 1: changer.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM); break;
		case 2: changer.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM_GB); break;
		}
		switch (bgCombo.getSelectedIndex()) {
		case 0: changer.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.BOTH); break;
		case 1: changer.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.MOUSE_ONLY); break;
		case 2: changer.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.HUMAN_ONLY); break;
		}
		switch (sortByCombo.getSelectedIndex()) {
		case 0: changer.setSetting(ChEA.SORT_BY, ChEA.PVALUE); break;
		case 1: changer.setSetting(ChEA.SORT_BY, ChEA.RANK); break;
		case 2: changer.setSetting(ChEA.SORT_BY, ChEA.COMBINED_SCORE); break;
		}
		changer.setSetting("number_of_top_TFs", selectTopText.getText());
	}
	
	// public accessor to set the input text area
	public void setInputTextArea(Collection<String> list) {
		UIUtils.setTextAreaText(inputTextArea, list);
	}
	
	// public accessor to set the output text area
	public void setOutputTextArea(Collection<String> list) {
		UIUtils.setTextAreaText(outputTextArea, list);
	}
	
	private void setupIO(File inputFile) {
		openPath.setText(inputFile.getAbsolutePath());
		UIUtils.setTextAreaText(inputTextArea, FileUtils.readFile(inputFile));
		
		File outputFile = new File(System.getProperty("user.dir"), FileUtils.stripFileExtension(inputFile.getName()) + ".results_tf.csv");
		saveChooser.setSelectedFile(outputFile);
		savePath.setText(outputFile.getAbsolutePath());
	}
	
	// Check if okay to enable view results button
	public void enableOutput(String output) {
		savePath.setText(output);
		this.output = output;
		if (Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.OPEN))
			openButton.setEnabled(true);
	}
	
}