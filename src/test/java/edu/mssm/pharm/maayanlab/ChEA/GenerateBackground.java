package edu.mssm.pharm.maayanlab.ChEA;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import edu.mssm.pharm.maayanlab.common.core.FileUtils;
import edu.mssm.pharm.maayanlab.common.math.Statistics;

public class GenerateBackground {
	
	private final static String humanGeneFile = "src/test/resources/human_genes.txt";
	private final static String mouseGeneFile = "src/test/resources/mouse_genes.txt";
	private final static String combinedGeneFile = "src/test/resources/combined_genes.txt";
	
	private final static int REPS = 1000;
	private final static int LENGTH = 300;
	
	private final static Random rng = new Random();
	
	public static void main(String[] args) {
//		generateHumanChEABackground();
//		generateMouseChEABackground();
//		generateCombinedChEABackground();
//		generateHumanTRANSFACBackground();
//		generateMouseTRANSFACBackground();
//		generateCombinedTRANSFACBackground();
		generatePWM_GBBackground();
	}
	
	private static Collection<String> generateRandomSample(ArrayList<String> list, int samples) {
		HashSet<String> sampleList = new HashSet<String>();
		while (sampleList.size() < samples)
			sampleList.add(list.get(rng.nextInt(list.size())));
		return sampleList;
	}
	
	private static ArrayList<String> generateOutputRanks(HashMap<String, ArrayList<Integer>> ranks) {
		ArrayList<String> output = new ArrayList<String>();
		for (String tf : ranks.keySet()) {
			double mean = Statistics.findMean(ranks.get(tf));
			double sd = Statistics.findStandardDeviation(ranks.get(tf), mean);
			output.add(tf + "\t" + mean + "\t" + sd);
		}
		
		return output;
	}
	
	private static void generateCombinedChEABackground() {
		long startTime = System.currentTimeMillis();
		
		ArrayList<String> humanGenes = FileUtils.readFile(combinedGeneFile);
		HashMap<String, ArrayList<Integer>> ranks = new HashMap<String, ArrayList<Integer>>();
		
		for (int i = 0; i < REPS; i++) {
			ChEA app = new ChEA();
			app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
			app.run(generateRandomSample(humanGenes, LENGTH));
			
			int counter = 1;
			for (TranscriptionFactor tf : app.getRankedList()) {
				if (!ranks.containsKey(tf.getName()))
					ranks.put(tf.getName(), new ArrayList<Integer>());
				ranks.get(tf.getName()).add(counter);
				counter++;
			}

			System.out.println("Run: " + i + " (" + (counter-1) + ")");
		}
		
		FileUtils.writeFile("combined_ChEA_ranks.txt", generateOutputRanks(ranks));
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: " + (endTime - startTime)/1000.0 + " seconds");
	}
	
	private static void generateHumanChEABackground() {
		long startTime = System.currentTimeMillis();
		
		ArrayList<String> humanGenes = FileUtils.readFile(humanGeneFile);
		HashMap<String, ArrayList<Integer>> ranks = new HashMap<String, ArrayList<Integer>>();
		
		for (int i = 0; i < REPS; i++) {
			ChEA app = new ChEA();
			app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
			app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.HUMAN_ONLY);
			app.run(generateRandomSample(humanGenes, LENGTH));
			
			int counter = 1;
			for (TranscriptionFactor tf : app.getRankedList()) {
				if (!ranks.containsKey(tf.getName()))
					ranks.put(tf.getName(), new ArrayList<Integer>());
				ranks.get(tf.getName()).add(counter);
				counter++;
			}

			System.out.println("Run: " + i + " (" + (counter-1) + ")");
		}
		
		FileUtils.writeFile("human_ChEA_ranks.txt", generateOutputRanks(ranks));
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: " + (endTime - startTime)/1000.0 + " seconds");
	}
	
	private static void generateMouseChEABackground() {
		long startTime = System.currentTimeMillis();
		
		ArrayList<String> mouseGenes = FileUtils.readFile(mouseGeneFile);
		HashMap<String, ArrayList<Integer>> ranks = new HashMap<String, ArrayList<Integer>>();
		
		for (int i = 0; i < REPS; i++) {
			ChEA app = new ChEA();
			app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
			app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.MOUSE_ONLY);
			app.run(generateRandomSample(mouseGenes, LENGTH));
			
			int counter = 1;
			for (TranscriptionFactor tf : app.getRankedList()) {
				if (!ranks.containsKey(tf.getName()))
					ranks.put(tf.getName(), new ArrayList<Integer>());
				ranks.get(tf.getName()).add(counter);
				counter++;
			}

			System.out.println("Run: " + i + " (" + (counter-1) + ")");
		}
		
		FileUtils.writeFile("mouse_ChEA_ranks.txt", generateOutputRanks(ranks));
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: " + (endTime - startTime)/1000.0 + " seconds");
	}
	
	private static void generateCombinedTRANSFACBackground() {
		long startTime = System.currentTimeMillis();
		
		ArrayList<String> humanGenes = FileUtils.readFile(humanGeneFile);
		HashMap<String, ArrayList<Integer>> ranks = new HashMap<String, ArrayList<Integer>>();
		
		for (int i = 0; i < REPS; i++) {
			ChEA app = new ChEA();
			app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
			app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM);
			app.run(generateRandomSample(humanGenes, LENGTH));
			
			int counter = 1;
			for (TranscriptionFactor tf : app.getRankedList()) {
				if (!ranks.containsKey(tf.getName()))
					ranks.put(tf.getName(), new ArrayList<Integer>());
				ranks.get(tf.getName()).add(counter);
				counter++;
			}

			System.out.println("Run: " + i + " (" + (counter-1) + ")");
		}
		
		FileUtils.writeFile("combined_TRANSFAC_ranks.txt", generateOutputRanks(ranks));
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: " + (endTime - startTime)/1000.0 + " seconds");
	}
	
	private static void generateHumanTRANSFACBackground() {
		long startTime = System.currentTimeMillis();
		
		ArrayList<String> humanGenes = FileUtils.readFile(humanGeneFile);
		HashMap<String, ArrayList<Integer>> ranks = new HashMap<String, ArrayList<Integer>>();
		
		for (int i = 0; i < REPS; i++) {
			ChEA app = new ChEA();
			app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
			app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM);
			app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.HUMAN_ONLY);
			app.run(generateRandomSample(humanGenes, LENGTH));
			
			int counter = 1;
			for (TranscriptionFactor tf : app.getRankedList()) {
				if (!ranks.containsKey(tf.getName()))
					ranks.put(tf.getName(), new ArrayList<Integer>());
				ranks.get(tf.getName()).add(counter);
				counter++;
			}

			System.out.println("Run: " + i + " (" + (counter-1) + ")");
		}
		
		FileUtils.writeFile("human_TRANSFAC_ranks.txt", generateOutputRanks(ranks));
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: " + (endTime - startTime)/1000.0 + " seconds");
	}
	
	private static void generateMouseTRANSFACBackground() {
		long startTime = System.currentTimeMillis();
		
		ArrayList<String> mouseGenes = FileUtils.readFile(mouseGeneFile);
		HashMap<String, ArrayList<Integer>> ranks = new HashMap<String, ArrayList<Integer>>();
		
		for (int i = 0; i < REPS; i++) {
			ChEA app = new ChEA();
			app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
			app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM);
			app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.MOUSE_ONLY);
			app.run(generateRandomSample(mouseGenes, LENGTH));
			
			int counter = 1;
			for (TranscriptionFactor tf : app.getRankedList()) {
				if (!ranks.containsKey(tf.getName()))
					ranks.put(tf.getName(), new ArrayList<Integer>());
				ranks.get(tf.getName()).add(counter);
				counter++;
			}

			System.out.println("Run: " + i + " (" + (counter-1) + ")");
		}
		
		FileUtils.writeFile("mouse_TRANSFAC_ranks.txt", generateOutputRanks(ranks));
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: " + (endTime - startTime)/1000.0 + " seconds");
	}
	
	private static void generatePWM_GBBackground() {
		long startTime = System.currentTimeMillis();
		
		ArrayList<String> humanGenes = FileUtils.readFile(humanGeneFile);
		HashMap<String, ArrayList<Integer>> ranks = new HashMap<String, ArrayList<Integer>>();
		
		for (int i = 0; i < REPS; i++) {
			ChEA app = new ChEA();
			app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
			app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM_GB);
			app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.HUMAN_ONLY);
			app.run(generateRandomSample(humanGenes, LENGTH));
			
			int counter = 1;
			for (TranscriptionFactor tf : app.getRankedList()) {
				if (!ranks.containsKey(tf.getName()))
					ranks.put(tf.getName(), new ArrayList<Integer>());
				ranks.get(tf.getName()).add(counter);
				counter++;
			}

			System.out.println("Run: " + i + " (" + (counter-1) + ")");
		}
		
		FileUtils.writeFile("human_PWM_GB_ranks.txt", generateOutputRanks(ranks));
		
		long endTime = System.currentTimeMillis();
		System.out.println("Elapsed time: " + (endTime - startTime)/1000.0 + " seconds");
	}
	
}
