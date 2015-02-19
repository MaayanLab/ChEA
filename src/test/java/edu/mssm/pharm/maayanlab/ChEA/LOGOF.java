package edu.mssm.pharm.maayanlab.ChEA;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;

import edu.mssm.pharm.maayanlab.common.core.FileUtils;

public class LOGOF {
	
	public static void main(String[] args) {
		LOGOF test = new LOGOF();
	}
	
	public LOGOF() {
		ArrayList<String> lines = FileUtils.readResource("logof.txt");
		HashMap<String, LinkedList<String>> experiments = new HashMap<String, LinkedList<String>>();
		for (String line : lines) {
			String[] splitLine = line.split("\t");
			if (splitLine[6].equals("LOF") && !splitLine[0].equals("sourceName")) { 
				String tfName = splitLine[0];
				String targetName = splitLine[2];
				if (!experiments.containsKey(tfName))
					experiments.put(tfName, new LinkedList<String>());			
				experiments.get(tfName).add(targetName);
			}
		}
		
		int tiedCount = 0, cheaCount = 0, pwmCount = 0, gbCount = 0, wtfCount = 0;
		for (String key : experiments.keySet()) {
			int score = scoreExperiments(key, experiments.get(key));
			if (score == 0)
				tiedCount++;
			else if (score == 1)
				cheaCount++;
			else if (score == 2)
				pwmCount++;
			else if (score == 3)
				gbCount++;
			else
				wtfCount++;
		}
		System.out.println(tiedCount);
		System.out.println(cheaCount);
		System.out.println(pwmCount);
		System.out.println(gbCount);
		System.out.println(wtfCount);
	}
	
	// 0 = tied
	// 1 = ChEA better
	// 2 = PWMs better
	// 3 = GB better
	// -1 = wtf?
	private int scoreExperiments(String tfName, Collection<String> geneList) {
		int cheaRank = 20, pwmRank = 20, gbRank = 20;
		ChEA chea = new ChEA();
//		chea.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.CHIPX);
//		chea.run(geneList);
//		ArrayList<String> cheaList = new ArrayList<String>(chea.getTopRankedList(20));		
//		for (int i = 0; i < 20; i++) {
//			if (tfName.equalsIgnoreCase(cheaList.get(i))) {
//				cheaRank = i;
//				break;
//			}
//		}
		chea.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM);
		chea.run(geneList);
		ArrayList<String> pwmList = new ArrayList<String>(chea.getTopRankedList(20));
		for (int i = 0; i < 20; i++) {
			if (tfName.equalsIgnoreCase(pwmList.get(i))) {
				pwmRank = i;
				break;
			}
		}
		
		chea.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM_GB);
		chea.run(geneList);
		ArrayList<String> gbList = new ArrayList<String>(chea.getTopRankedList(20));
		for (int i = 0; i < 20; i++) {
			if (tfName.equalsIgnoreCase(gbList.get(i))) {
				gbRank = i;
				break;
			}
		}
		
		if (cheaRank < pwmRank && cheaRank < gbRank)
			return 1;
		else if (pwmRank < cheaRank && pwmRank < gbRank)
			return 2;
		else if (gbRank < pwmRank && gbRank < cheaRank)
			return 3;
		else if (gbRank == pwmRank && pwmRank == cheaRank)
			return 0;
		else
			return -1;
	}
	
}
