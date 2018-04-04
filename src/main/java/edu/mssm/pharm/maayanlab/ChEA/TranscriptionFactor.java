package edu.mssm.pharm.maayanlab.ChEA;

import java.util.HashSet;
import java.util.Set;

import com.google.gson.annotations.Expose;
public class TranscriptionFactor implements Comparable<Object> {
	@Expose
	private String name;	// TF-ID
	@Expose
	private String simpleName;	// TF
	private int species;	// Species that transcription factor occurs in
	private Set<String> targets = new HashSet<String>();
	
	public final static int MOUSE = 1;
	public final static int HUMAN = 2;
	public final static int BOTH = 3; 
	
	private double mean;
	private double standardDeviation;
	
	private double fractionOfTargetsInInput;
	private double fractionOfTargetsInBackground;
	@Expose
	private double pvalue;
	@Expose
	private double zscore;
	@Expose
	private double combinedScore;
	@Expose
	private Set<String> enrichedTargets;
	
	public TranscriptionFactor(
			String name,
			String simpleName,
			String species,
			String target) {
		
		this.name = name;
		this.simpleName = simpleName;
		this.species = (species.equals("MOUSE")) ? TranscriptionFactor.MOUSE : TranscriptionFactor.HUMAN;
		targets.add(target);
	}

	// Get name with experiment id still attached
	public String getName() {
		return this.name;
	}
	
	// Get only TF name
	public String getSimpleName() {
		return this.simpleName;
	}
	
	public int getSpecies() {
		return this.species;
	}
	
	public void setSpecies(String species) {
		this.species |= (species.equals("MOUSE")) ? TranscriptionFactor.MOUSE : TranscriptionFactor.HUMAN;
	}
	
	public void addTarget(String target) {
		targets.add(target);
	}
	
	public Set<String> getTargets() {
		return targets;
	}
	
	protected void setTargets(Set<String> targets) {
		this.targets = targets;
	}
	
	public void setRankStats(double mean, double standardDeviation) {
		this.mean = mean;
		this.standardDeviation = standardDeviation;
	}
	
	public Set<String> getEnrichedTargets() {
		return enrichedTargets;
	}

	public void setEnrichedTargets(Set<String> enrichedTargets) {
		this.enrichedTargets = enrichedTargets;
	}
	
	public void setFractionOfTargetsInInput(double fractionOfTargetsInInput) {
		this.fractionOfTargetsInInput = fractionOfTargetsInInput;
	}
	
	public void setFractionOfTargetsInBackground(double fractionOfTargetsInBackground) {
		this.fractionOfTargetsInBackground = fractionOfTargetsInBackground;
	}
	
	public double getPValue() {
		return this.pvalue;
	}
	
	public void setPValue(double pvalue) {
		this.pvalue = pvalue;
	}
	
	public double getZScore() {
		return this.zscore;
	}
	
	public double getCombinedScore() {
		return this.combinedScore;
	}
	
	/*
	1. name of the transcription factor
	2. number of targets in the input gene-list
	3. number of genes that are targeted by the transcription factor
	4. the fraction of genes targeted compared to total number of genes in gene-list
	5. the fraction of genes targeted compared to total number of genes in background
	6. difference between the background fraction and the gene-list fraction
	7. p-value computed using the Fisher Test
	8. rank computed using z-test
	9. combined score computed from p-value and rank
	10. list of genes separated by a semi-colon
	 */	
	@Override
	public String toString() {
		StringBuilder outputString = new StringBuilder();
		outputString.append(name).append(",");
		outputString.append(enrichedTargets.size()).append(",");
		outputString.append(targets.size()).append(",");
		outputString.append(fractionOfTargetsInInput).append(",");
		outputString.append(fractionOfTargetsInBackground).append(",");
		outputString.append((fractionOfTargetsInInput - fractionOfTargetsInBackground)).append(",");
		outputString.append(pvalue).append(",").append(zscore).append(",").append(combinedScore).append(",");
		
		boolean firstTarget = true;
		for (String enrichedTarget : enrichedTargets) {
			if (firstTarget) {
				outputString.append(enrichedTarget);
				firstTarget = false;
			}
			else
				outputString.append(";").append(enrichedTarget);
		}
		
		return outputString.toString();
	}

	@Override
	public int compareTo(Object o) {
		if (this.pvalue > ((TranscriptionFactor) o).pvalue)
			return 1;
		else if (this.pvalue < ((TranscriptionFactor) o).pvalue)
			return -1;
		else
			return 0;
	}
	
	// Collapse TFs with the same name to the same one
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((simpleName == null) ? 0 : simpleName.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		TranscriptionFactor other = (TranscriptionFactor) obj;
		if (simpleName == null) {
			if (other.simpleName != null)
				return false;
		} else if (!simpleName.equals(other.simpleName))
			return false;
		return true;
	}

	public void computeScore(int currentRank) {
		if (mean == 0 && standardDeviation == 0)
			zscore = 0;
		else
			zscore = (currentRank - mean)/standardDeviation;
		combinedScore = Math.log(pvalue)*zscore;
	}
}
