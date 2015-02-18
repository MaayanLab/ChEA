package edu.mssm.pharm.maayanlab.ChEA;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.ListIterator;

import pal.statistics.FisherExact;
import edu.mssm.pharm.maayanlab.FileUtils;
import edu.mssm.pharm.maayanlab.SetOps;
import edu.mssm.pharm.maayanlab.Settings;
import edu.mssm.pharm.maayanlab.SettingsChanger;

public class ChEA implements SettingsChanger {

	private LinkedList<TranscriptionFactor> transcriptionFactors;
	
	protected final static String CHEA_BACKGROUND = "res/chea_background.csv";
	protected final static String TRANSFAC_BACKGROUND = "res/transfac_background.csv";
	protected final static String PWM_GB_BACKGROUND = "res/PWM-GB.csv";

	// Paths to the background rank files
	private final String MOUSE_CHEA_RANKS = "res/mouse_ChEA_ranks.txt";
	private final String HUMAN_CHEA_RANKS = "res/human_ChEA_ranks.txt";
	private final String COMBINED_CHEA_RANKS = "res/combined_ChEA_ranks.txt";
	private final String MOUSE_TRANSFAC_RANKS = "res/mouse_TRANSFAC_ranks.txt";
	private final String HUMAN_TRANSFAC_RANKS = "res/human_TRANSFAC_ranks.txt";
	private final String COMBINED_TRANSFAC_RANKS = "res/combined_TRANSFAC_ranks.txt";
	private final String PWM_GB_RANKS = "res/PWM-GB_ranks.txt";
			
	// Output header
	protected final String HEADER = "TF,Target/Input,Targets/Database,Fraction/Input,Fraction/Database,Difference,P-value,Z-score,Combined Score,Genes";
	
	// Default settings
	private final Settings settings = new Settings() {
		{
			// String: rank the TFs by the Fisher Exact test's p-value, rank against the background of random genes, or combined score of the two. [combined score/p-value/rank]
			set(ChEA.SORT_BY, ChEA.COMBINED_SCORE);
			// String: the organisms included in the transcription factor background database for enrichment analysis. [mouse/human/both]
			set(ChEA.INCLUDED_ORGANISMS, ChEA.BOTH);
			// String: the source of the transcription factor background database used for enrichment analysis. [ChIP-X/PWM-JT/PWM-GB]
			set(ChEA.BACKGROUND_DATABASE, ChEA.CHIPX);
		}
	};
	
	// Settings variables
	public final static String SORT_BY = "sort transcription factors by";
	public final static String COMBINED_SCORE = "combined score";
	public final static String PVALUE = "p-value";
	public final static String RANK = "rank";
	
	public final static String INCLUDED_ORGANISMS = "included organisms in the background database";
	public final static String MOUSE_ONLY = "mouse";
	public final static String HUMAN_ONLY = "human";
	public final static String BOTH = "both";
	
	public final static String BACKGROUND_DATABASE = "TF-target gene background database used for enrichment";
	public final static String CHIPX = "ChIP-X";
	public final static String PWM = "PWM-JT";
	public final static String PWM_GB = "PWM-GB";
	
	public static void main(String[] args) {
		if (args.length == 2) {
			ChEA chea = new ChEA();
			chea.run(FileUtils.readFile(args[0]));
			chea.writeFile(args[1]);
		}
		else if (args.length == 3) {
			ChEA chea = new ChEA();
			chea.run(args[0], args[1]);
			chea.writeFile(args[2]);
		}
		else
			System.err.println("Usage: java -jar chea.jar [background] genelist output");
	}
	
	// By default, load settings from file
	public ChEA() {
		settings.loadSettings();
	}
	
	// Load external settings, primarily for use with X2K
	public ChEA(Settings externalSettings) {
		settings.loadSettings(externalSettings);
	}
	
	@Override
	// Used for other methods to set settings
	public void setSetting(String key, String value) {
		settings.set(key, value);
	}
	
	// Used for other methods to get settings
	public String getSetting(String key) {
		return settings.get(key);
	}
	
	// Run for file names
	public void run(String background, String geneList) {
		ArrayList<String> inputList = FileUtils.readFile(geneList);
		
		try {
			if (FileUtils.validateList(inputList))
				readBackground(FileUtils.readFile(background));
				computeEnrichment(inputList);
		} catch (ParseException e) {
			if (e.getErrorOffset() == -1)
				System.err.println("Invalid input: " + "Input list is empty.");
			else
				System.err.println("Invalid input: " + e.getMessage() + " at line " + (e.getErrorOffset() + 1) +" is not a valid Entrez Gene Symbol.");
			System.exit(-1);	
		}
	}
	
	// Run for calling from other methods and pass in collection
	public void run(Collection<String> genelist) {
		if (settings.get(BACKGROUND_DATABASE).equals(PWM))
			readBackground(FileUtils.readResource(TRANSFAC_BACKGROUND));
		else if (settings.get(BACKGROUND_DATABASE).equals(PWM_GB)) {
			setSetting(INCLUDED_ORGANISMS, HUMAN_ONLY);
			readBackground(FileUtils.readResource(PWM_GB_BACKGROUND));			
		}
		else			
			readBackground(FileUtils.readResource(CHEA_BACKGROUND));
		computeEnrichment(genelist);
	}

	public void writeFile(String filename) {
		FileUtils.writeFile(filename, HEADER, transcriptionFactors);
	}
	
	// Used to pass information about transcription factors
	public Collection<TranscriptionFactor> getTopRanked(int ranks) {
		LinkedHashSet<TranscriptionFactor> topRanked = new LinkedHashSet<TranscriptionFactor>(ranks);
		
		Iterator<TranscriptionFactor> itr = transcriptionFactors.iterator();
		while (itr.hasNext() && topRanked.size() < ranks)
			topRanked.add(itr.next());
		
		return topRanked;
	}
	
	// Used to print lists
	public Collection<String> getTopRankedList(int ranks) {
		LinkedHashSet<String> topRanked = new LinkedHashSet<String>(ranks);				
		
		Iterator<TranscriptionFactor> itr = transcriptionFactors.iterator();
		while (itr.hasNext() && topRanked.size() < ranks)
			topRanked.add(itr.next().getSimpleName());
		
		return topRanked;
	}
	
	public Collection<TranscriptionFactor> getRankedList() {
		return transcriptionFactors;
	}
	
	@Deprecated
	// For internal use to get ranked list of tf-experiment names
	protected Collection<TranscriptionFactor> getTopRankedList() {
		return transcriptionFactors;
	}
	
	protected void readBackground(ArrayList<String> background) {
		HashMap<String, TranscriptionFactor> tfMap = new HashMap<String, TranscriptionFactor>();
		
		for (String record : background) {
			// Split record line
			String[] splitRecord = record.toUpperCase().split(",");
			String simpleName = splitRecord[1];
			String name = splitRecord[2];			
			String target = splitRecord[3];
			String species = splitRecord[7];
		
			// Ignore entries that are not in the included organisms list
			if (settings.get(INCLUDED_ORGANISMS).equals(MOUSE_ONLY) && !species.equals("MOUSE"))
				continue;
			else if (settings.get(INCLUDED_ORGANISMS).equals(HUMAN_ONLY) && !species.equals("HUMAN"))
				continue;
			
			if (tfMap.containsKey(name)) {
				TranscriptionFactor tf = tfMap.get(name);
				tf.addTarget(target);
				tf.setSpecies(species);
			}
			else {
				tfMap.put(name, new TranscriptionFactor(name, simpleName, species, target));
			}
		}
		
		// read correction ranks
		ArrayList<String> ranks;
		if (settings.get(BACKGROUND_DATABASE).equals(PWM)) {
			if (settings.get(INCLUDED_ORGANISMS).equals(MOUSE_ONLY))
				ranks = FileUtils.readResource(MOUSE_TRANSFAC_RANKS);
			else if (settings.get(INCLUDED_ORGANISMS).equals(HUMAN_ONLY))
				ranks = FileUtils.readResource(HUMAN_TRANSFAC_RANKS);
			else
				ranks = FileUtils.readResource(COMBINED_TRANSFAC_RANKS);			
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(PWM_GB)) {
			ranks = FileUtils.readResource(PWM_GB_RANKS);
		}
		else {
			if (settings.get(INCLUDED_ORGANISMS).equals(MOUSE_ONLY))
				ranks = FileUtils.readResource(MOUSE_CHEA_RANKS);
			else if (settings.get(INCLUDED_ORGANISMS).equals(HUMAN_ONLY))
				ranks = FileUtils.readResource(HUMAN_CHEA_RANKS);
			else
				ranks = FileUtils.readResource(COMBINED_CHEA_RANKS);
		}
		for (String rank : ranks) {
			String[] split_rank = rank.split("\\s");
			tfMap.get(split_rank[0]).setRankStats(Double.parseDouble(split_rank[1]), Double.parseDouble(split_rank[2]));
		}
		
		transcriptionFactors = new LinkedList<TranscriptionFactor>(tfMap.values());
	}
	
	protected void computeEnrichment(Collection<String> genes) {

		// filter genes from input list that are not associated with an upstream transcription factor
		HashSet<String> geneInputSet = new HashSet<String>();
		for (String gene : genes)
			geneInputSet.add(gene.toUpperCase());
		HashSet<String> geneBgSet = new HashSet<String>();
		for (TranscriptionFactor tf : transcriptionFactors)
			geneBgSet.addAll(tf.getTargets());
		geneInputSet.retainAll(geneBgSet);
		
		ListIterator<TranscriptionFactor> tfIterator = transcriptionFactors.listIterator();
				
		while (tfIterator.hasNext()) {
			TranscriptionFactor currentFactor = tfIterator.next();
			
			HashSet<String> targetBgGenes = currentFactor.getTargets();
			// Target input genes is the intersection of target background genes and input genes
			HashSet<String> targetInputGenes = SetOps.intersection(targetBgGenes, geneInputSet);
					
			double numOfTargetBgGenes = targetBgGenes.size();
			double totalBgGenes = geneBgSet.size();
			double totalInputGenes = geneInputSet.size();
			double numOfTargetInputGenes = targetInputGenes.size();
			
			if (targetInputGenes.size() > 0) {
				
				FisherExact fisherTest = new FisherExact(targetInputGenes.size()
						+ (geneInputSet.size() - targetInputGenes.size())
						+ targetBgGenes.size()
						+ (geneBgSet.size() - targetBgGenes.size()));
				
				double pvalue = fisherTest.getRightTailedP(targetInputGenes.size(),
						(geneInputSet.size() - targetInputGenes.size()), targetBgGenes.size(), 
						(geneBgSet.size() - targetBgGenes.size()));

				currentFactor.setEnrichedTargets(targetInputGenes);
				currentFactor.setFractionOfTargetsInInput(numOfTargetInputGenes/totalInputGenes);
				currentFactor.setFractionOfTargetsInBackground(numOfTargetBgGenes/totalBgGenes);
				currentFactor.setPValue(pvalue);
			}
			else {
				tfIterator.remove();
			}
		}
		
		// First, sort by p-value
		Collections.sort(transcriptionFactors);
		
		// Count current rank and compute z-score
		int counter = 1;
		for (TranscriptionFactor tf : transcriptionFactors) {
			tf.computeScore(counter);
			counter++;
		}
		
		if (settings.get(SORT_BY).equals(COMBINED_SCORE)) {
			// Sort by combined score
			Collections.sort(transcriptionFactors, new Comparator<TranscriptionFactor>() {
				@Override
				public int compare(TranscriptionFactor o1, TranscriptionFactor o2) {
					if (o1.getCombinedScore() < o2.getCombinedScore())				
						return 1;
					else if (o1.getCombinedScore() > o2.getCombinedScore())
						return -1;
					else
						return 0;
				}
			});
		}
		else if (settings.get(SORT_BY).equals(RANK)) {
			// Sort by z-score
			Collections.sort(transcriptionFactors, new Comparator<TranscriptionFactor>() {
				@Override
				public int compare(TranscriptionFactor o1, TranscriptionFactor o2) {
					if (o1.getZScore() > o2.getZScore())				
						return 1;
					else if (o1.getZScore() < o2.getZScore())
						return -1;
					else
						return 0;
				}
			});
		}
	}
}