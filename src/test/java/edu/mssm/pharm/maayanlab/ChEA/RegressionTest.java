package edu.mssm.pharm.maayanlab.ChEA;

import java.util.Collection;
import java.util.Iterator;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import edu.mssm.pharm.maayanlab.common.core.FileUtils;

public class RegressionTest extends TestCase {

	private ChEA app;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		app = new ChEA();
	}

	/**
	 * @return the suite of tests being tested
	 */
	public static Test suite()
	{
		return new TestSuite( RegressionTest.class );
	}

//	public void testChEA() {
//		app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.CHIPX);
//		app.setSetting(ChEA.SORT_BY, ChEA.COMBINED_SCORE);
//		app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.BOTH);
//		app.run(FileUtils.readResource("test_list.txt"));
//		
//		assertEquivalentOutput("testChEA_results.csv");
//	}
//	
//	public void testPWMs() {
//		app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM);
//		app.setSetting(ChEA.SORT_BY, ChEA.COMBINED_SCORE);
//		app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.BOTH);
//		app.run(FileUtils.readResource("test_list.txt"));
//		
//		assertEquivalentOutput("testPWMs_results.csv");
//	}
//	
	public void testMouseOnly() {
		app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.CHIPX);
		app.setSetting(ChEA.SORT_BY, ChEA.COMBINED_SCORE);
		app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.MOUSE_ONLY);
		app.run(FileUtils.readResource("test_list.txt"));
		
		for (TranscriptionFactor tf : app.getRankedList())
			assertEquals(tf.getSpecies(), TranscriptionFactor.MOUSE);
	}
	
	public void testHumanOnly() {
		app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.PWM);
		app.setSetting(ChEA.SORT_BY, ChEA.COMBINED_SCORE);
		app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.HUMAN_ONLY);
		app.run(FileUtils.readResource("test_list.txt"));
		
		for (TranscriptionFactor tf : app.getRankedList())
			assertEquals(tf.getSpecies(), TranscriptionFactor.HUMAN);
	}
	
	public void testSorting() {
		app.setSetting(ChEA.BACKGROUND_DATABASE, ChEA.CHIPX);
		app.setSetting(ChEA.SORT_BY, ChEA.PVALUE);
		app.setSetting(ChEA.INCLUDED_ORGANISMS, ChEA.BOTH);
		app.run(FileUtils.readResource("test_list.txt"));
		
		double score = 0;
		for (TranscriptionFactor tf : app.getRankedList()) {
			assertTrue(score <= tf.getPValue());
			score = tf.getPValue();
		}
	}
	
	private void assertEquivalentOutput(String expectedFile) {
		Collection<TranscriptionFactor> tfs = app.getRankedList();
		Iterator<TranscriptionFactor> tf = tfs.iterator();
		Collection<String> testResults = FileUtils.readResource(expectedFile);
		Iterator<String> result = testResults.iterator();
		
		assertEquals(testResults.size(), tfs.size()+1);
		assertEquals(result.next(), app.HEADER);
		
		while (tf.hasNext())
			assertEquals(tf.next().toString(), result.next());
	}

}
