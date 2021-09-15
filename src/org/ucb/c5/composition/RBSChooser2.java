package org.ucb.c5.composition;
import java.lang.reflect.Array;
import java.util.*;

import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.sequtils.CalcEditDistance;
import org.ucb.c5.sequtils.HairpinCounter;
import org.ucb.c5.sequtils.RevComp;
import org.ucb.c5.sequtils.Translate;
import org.ucb.c5.utils.FileUtils;

/**
 * Second generation RBSChooser algorithm
 *
 * Employs a list of genes and their associated ribosome binding sites for
 * highly-expressed proteins in E. coli.
 *
 * @author J. Christopher Anderson
 * @author Milo
 */
public class RBSChooser2 {

    private List<RBSOption> rbss;
    private CalcEditDistance CalcEditDist;
    private HairpinCounter HPCounter;
    private RevComp Reverse;
    private Translate Trans;

    public void initiate() throws Exception {
        //initiate functions and variables
        rbss = new ArrayList<>();
        CalcEditDist = new CalcEditDistance();
        CalcEditDist.initiate();
        HPCounter = new HairpinCounter();
        HPCounter.initiate();
        Reverse = new RevComp();
        Reverse.initiate();
        Trans = new Translate();
        Trans.initiate();
        //Read the data file
        String data = FileUtils.readResourceFile("composition/data/coli_genes.txt");
        String options = FileUtils.readResourceFile("composition/data/rbs_options.txt");
        HashMap<String, String> geneList = new HashMap<String, String>();
        //Scan through each line, add names of genes mapped to cds
        String[] lines = data.split("\\r|\\r?\\n");
        for (int i = 0; i < lines.length; i++) {
            String line = lines[i];
            if(line.isEmpty()) {continue;}
            String[] tabs = line.split("\t");
            geneList.put(tabs[1], tabs[6]);
        }
        //scan through each line, add RBSOptions to rbss using information from hash map and options
        String[] linesR = options.split("\\r|\\r?\\n");
        for (int i = 0; i < linesR.length; i++) {
            String line = linesR[i];
            if(line.isEmpty()) {continue;}
            String[] tabs = line.split("\t");
            String name = tabs[0];
            String description = "RBS from " + tabs[0] + " gene";
            String RBS = tabs[1];
            String CDS = geneList.get(tabs[0]);
            String first6AA = Trans.run(CDS.substring(0, 30));
            RBSOption entry = new RBSOption(name, description, RBS, CDS, first6AA);
            rbss.add(entry);
        }
    }
    //Questions: why is it giving invalid file err? where do we get the description from?

    /**
     * Provided an ORF of sequence 'cds', this computes the best ribosome
     * binding site to use from a list of options.
     *
     * It also permits a list of options to exclude.
     *
     * @param cds The DNA sequence, ie ATGCATGAT...
     * @param ignores The list of RBS's to exclude
     * @return
     * @throws Exception
    .     */
    public RBSOption run(String cds, Set<RBSOption> ignores) throws Exception {
        //TODO:  Fill in

        return null;
    }


    public static void main(String[] args) throws Exception {
        //Create an example
        String cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATTATCACCATCGTTATATGACACAGCTTGGGTGGCTATGGTGCCGGAGCGTAGTTCTTCTCAA";

        //Initiate the chooser
        RBSChooser2 chooser = new RBSChooser2();
        chooser.initiate();
        /*
        //Make the first choice with an empty Set of ignores
        Set<RBSOption> ignores = new HashSet<>();
        RBSOption selected1 = chooser.run(cds, ignores);

        //Add the first selection to the list of things to ignore
        ignores.add(selected1);

        //Choose again with an ignore added
        RBSOption selected2 = chooser.run(cds, ignores);

        //Print out the two options, which should be different
        System.out.println("CDS starts with:");
        System.out.println(cds.substring(0, 18));
        System.out.println();
        System.out.println("Selected1:\n");
        System.out.println(selected1.toString());
        System.out.println();
        System.out.println("Selected2:\n");
        System.out.println(selected2.toString()); */
    }
}
