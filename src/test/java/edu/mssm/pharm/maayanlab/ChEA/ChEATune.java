package edu.mssm.pharm.maayanlab.ChEA;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

import edu.mssm.pharm.maayanlab.common.core.FileUtils;

public class ChEATune {

	public static void findScoreMean() {
		ChEA app;
		
		ArrayList<String> backgroundFile = FileUtils.readResource(ChEA.CHEA_BACKGROUND);
		
		// 1000 trials
		for (int i = 0; i < 1000; i++) {
			app = new ChEA();
			app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.CHIPX);
			app.readBackground(backgroundFile);
			Collection<TranscriptionFactor> tfs = app.getRankedList();
			Iterator<TranscriptionFactor> tfIterator = tfs.iterator();
			
			while (tfIterator.hasNext()) {
				TranscriptionFactor tf = tfIterator.next();
				tf.setTargets(pickSample(tf.getTargets(), 200));
			}
			
//			app.computeEnrichment(genes);
		}
	}
	
	private static Set<String> pickSample(Set<String> population, int samplesNeeded) {
		HashSet<String> sampledPopulation = new HashSet<String>();
		Random r = new Random();
		
		Iterator<String> popIterator = population.iterator();
		int samplesLeft = population.size();
		
		while (samplesNeeded > 0) {
			int randomNumber = r.nextInt(samplesLeft--);
			if (randomNumber < samplesNeeded) {
				sampledPopulation.add(popIterator.next());
				samplesNeeded--;
			}
		}
		
		return sampledPopulation;
	}
	
}
