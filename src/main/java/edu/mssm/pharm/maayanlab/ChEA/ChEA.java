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
import java.util.Set;

import pal.statistics.FisherExact;
import edu.mssm.pharm.maayanlab.common.core.FileUtils;
import edu.mssm.pharm.maayanlab.common.core.Settings;
import edu.mssm.pharm.maayanlab.common.core.SettingsChanger;
import edu.mssm.pharm.maayanlab.common.math.SetOps;

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
			
	// New
	protected static final String CHEA_2015_BACKGROUND = "res/ChEA_2015_background.txt";
	private final String MOUSE_CHEA_2015_RANKS = "res/ChEA_ranks_mouse.tsv";
	private final String HUMAN_CHEA_2015_RANKS = "res/ChEA_ranks_human.tsv";
	private final String COMBINED_CHEA_2015_RANKS = "res/ChEA_ranks.tsv";


	protected static final String TRANS_JASP_BACKGROUND = "res/TRANSFAC_and_JASPAR_PWMs_background.tsv";
	private final String MOUSE_TRANS_JASP_RANKS = "res/Transfac_and_Jaspar_PWMs_ranks_mouse.tsv";
	private final String HUMAN_TRANS_JASP_RANKS = "res/Transfac_and_Jaspar_PWMs_ranks_human.tsv";
	private final String COMBINED_TRANS_JASP_RANKS = "res/Transfac_and_Jaspar_PWMs_ranks.tsv";

	protected static final String CONSENSUS_BACKGROUND = "res/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X_background.txt";
	private final String CONSENSUS_RANKS = "res/ENCODE_and_ChEA_Consensus_ranks.tsv";

	protected static final String ENCODE_2015_BACKGROUND = "res/ENCODE_TF_ChIP-seq_2015_background.txt";
	private final String MOUSE_ENCODE_2015_RANKS = "res/ENCODE_TF_ChIP-seq_2015_ranks_mouse.tsv";
	private final String HUMAN_ENCODE_2015_RANKS = "res/ENCODE_TF_ChIP-seq_2015_ranks_human.tsv";
	private final String COMBINED_ENCODE_2015_RANKS = "res/ENCODE_TF_ChIP-seq_2015_ranks.tsv";


	protected final static String ARCHS4_BACKGROUND = "res/ARCHS4_TFs_Coexp.csv";
	private final String HUMAN_ARCHS4_RANKS = "res/ARCHS4_TFs_Coexp_ranks.txt";

	protected final static String ENRICHR_BACKGROUND = "res/Enrichr-Submissions-TF-Gene_Coocurrence.tsv";
	private final String COMBINED_ENRICHR_RANKS = "res/Enrichr_Submissions_TF-Gene_Coocurrence_ranks.txt";

	protected final static String CHEA_2016_BACKGROUND = "res/CHEA-2016_Both_TF.tsv";
	protected final static String CHEA_2016_RANKS = "";
	
	protected final static String CREEDS_BACKGROUND = "res/CREEDS_TF.tsv";
	protected final static String CREEDS_RANKS = "";
	
	// Output header
	protected final String HEADER = "TF,Target/Input,Targets/Database,Fraction/Input,Fraction/Database,Difference,P-value,Z-score,Combined Score,Genes";
	
	// Default settings
	private final Settings settings = new Settings() {
		{
			// String: rank the TFs by the Fisher Exact test's p-value, rank against the background of random genes, or combined score of the two. [combined_score/pvalue/rank]
			set(ChEA.SORT_BY, ChEA.COMBINED_SCORE);
			// String: the organisms included in the transcription factor background database for enrichment analysis. [mouse_only/human_only/both]
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
	public static final String ARCHS4 = "ARCHS4 TFs Coexp";
	public static final String CHEA_2015 = "ChEA 2015";
	public static final String CHEA_2016 = "ChEA 2016";
	public static final String CONSENSUS = "ChEA & ENCODE Consensus";
	public static final String CREEDS = "CREEDS";
	public static final String ENCODE_2015 = "ENCODE 2015";
	public static final String ENRICHR = "Enrichr Submissions TF-Gene Coocurrence";
	public static final String TRANS_JASP = "Transfac and Jaspar";
	
	
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
		else if (settings.get(BACKGROUND_DATABASE).equals(CHEA_2015)) {
			readBackground(FileUtils.readResource(CHEA_2015_BACKGROUND));
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(TRANS_JASP)) {
			readBackground(FileUtils.readResource(TRANS_JASP_BACKGROUND));
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(CONSENSUS)) {
			readBackground(FileUtils.readResource(CONSENSUS_BACKGROUND));
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(ENCODE_2015)) {
			readBackground(FileUtils.readResource(ENCODE_2015_BACKGROUND));
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(ARCHS4)) {
			setSetting(INCLUDED_ORGANISMS, HUMAN_ONLY);
			readBackground(FileUtils.readResource(ARCHS4_BACKGROUND));
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(ENRICHR)) {
			setSetting(INCLUDED_ORGANISMS, BOTH);
			readBackground(FileUtils.readResource(ENRICHR_BACKGROUND));
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(CHEA_2016)) {
			readBackground(FileUtils.readResource(CHEA_2016_BACKGROUND));
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(CREEDS)) {
			readBackground(FileUtils.readResource(CREEDS_BACKGROUND));
		}
		else {
			if(settings.get(BACKGROUND_DATABASE) == null) {
				System.err.println("WARN: BACKGROUND_DATABASE (null) couldn't be detected!");
			} else {
				System.err.println("WARN: BACKGROUND_DATABASE (" + settings.get(BACKGROUND_DATABASE) + ") couldn't be detected!");
			}
			readBackground(FileUtils.readResource(CHEA_BACKGROUND));
		}
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
			String splitter;
			if (record.contains("\t")) // TODO: use extension instead
				splitter = "\t";
			else
				splitter = ",";

			String[] splitRecord = record.split(splitter);
			String simpleName = splitRecord[1].toUpperCase();
			String name = splitRecord[2].toUpperCase();
			String target = splitRecord[3].toUpperCase();
			String species = splitRecord[7].toUpperCase();
			String meta = "{}";
			if(splitRecord.length >= 10)
				meta = splitRecord[9];

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
				tfMap.put(name, new TranscriptionFactor(name, simpleName, species, target, meta));
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
		else if (settings.get(BACKGROUND_DATABASE).equals(CHEA_2015)) {
			if (settings.get(INCLUDED_ORGANISMS).equals(MOUSE_ONLY))
				ranks = FileUtils.readResource(MOUSE_CHEA_2015_RANKS);
			else if (settings.get(INCLUDED_ORGANISMS).equals(HUMAN_ONLY))
				ranks = FileUtils.readResource(HUMAN_CHEA_2015_RANKS);
			else
				ranks = FileUtils.readResource(COMBINED_CHEA_2015_RANKS);			
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(TRANS_JASP)) {
			if (settings.get(INCLUDED_ORGANISMS).equals(MOUSE_ONLY))
				ranks = FileUtils.readResource(MOUSE_TRANS_JASP_RANKS);
			else if (settings.get(INCLUDED_ORGANISMS).equals(HUMAN_ONLY))
				ranks = FileUtils.readResource(HUMAN_TRANS_JASP_RANKS);
			else
				ranks = FileUtils.readResource(COMBINED_TRANS_JASP_RANKS);			
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(CONSENSUS)) {
			ranks = FileUtils.readResource(CONSENSUS_RANKS);			
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(ENCODE_2015)) {
			if (settings.get(INCLUDED_ORGANISMS).equals(MOUSE_ONLY))
				ranks = FileUtils.readResource(MOUSE_ENCODE_2015_RANKS);
			else if (settings.get(INCLUDED_ORGANISMS).equals(HUMAN_ONLY))
				ranks = FileUtils.readResource(HUMAN_ENCODE_2015_RANKS);
			else
				ranks = FileUtils.readResource(COMBINED_ENCODE_2015_RANKS);			
		}		
		else if (settings.get(BACKGROUND_DATABASE).equals(ARCHS4)) {
			ranks = FileUtils.readResource(HUMAN_ARCHS4_RANKS);		
		}	
		else if (settings.get(BACKGROUND_DATABASE).equals(ENRICHR)) {
			ranks = FileUtils.readResource(COMBINED_ENRICHR_RANKS);		
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(CHEA_2016)) {
			ranks = FileUtils.readResource(CHEA_2016_RANKS);		
		}
		else if (settings.get(BACKGROUND_DATABASE).equals(CREEDS)) {
			ranks = FileUtils.readResource(CREEDS_RANKS);		
		}
		else {
			if (settings.get(INCLUDED_ORGANISMS).equals(MOUSE_ONLY))
				ranks = FileUtils.readResource(MOUSE_CHEA_RANKS);
			else if (settings.get(INCLUDED_ORGANISMS).equals(HUMAN_ONLY))
				ranks = FileUtils.readResource(HUMAN_CHEA_RANKS);
			else
				ranks = FileUtils.readResource(COMBINED_CHEA_RANKS);
		}
		// Disabled Ranking for now
		// for (String rank : ranks) {
		// 	String[] split_rank = rank.toUpperCase().split("\\t");
		// 	TranscriptionFactor tf = tfMap.get(split_rank[0]);
		// 	if(tf == null)
		// 		throw new NullPointerException("setRankStats: "+rank);
		// 	tf.setRankStats(Double.parseDouble(split_rank[1]), Double.parseDouble(split_rank[2]));
		// }
		
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
			
			Set<String> targetBgGenes = currentFactor.getTargets();
			// Target input genes is the intersection of target background genes and input genes
			Set<String> targetInputGenes = SetOps.intersection(targetBgGenes, geneInputSet);
					
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
		
		// Disabled Ranking for now
		// Count current rank and compute z-score
		// int counter = 1;
		// for (TranscriptionFactor tf : transcriptionFactors) {
		// 	tf.computeScore(counter);
		// 	counter++;
		// }
		
		if (settings.get(SORT_BY).equals(COMBINED_SCORE)) {
			// Sort by combined score
			Collections.sort(transcriptionFactors, new Comparator<TranscriptionFactor>() {
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