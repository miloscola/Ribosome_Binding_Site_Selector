package org.ucb.c5.composition;
import java.lang.reflect.Array;
import java.util.*;

import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.sequtils.*;
import org.ucb.c5.utils.FileUtils;

/**
 * Second generation RBSChooser algorithm
 *
 * Employs a list of genes and their associated ribosome binding sites for
 * highly-expressed proteins in E. coli.
 *
 * @author Milo
 */
public class RBSChooser2 {

    private List<RBSOption> rbss;
    private CalcEditDistance CalcEditDist;
    private HairpinCounter HPCounter;
    private RevComp Reverse;
    private Translate Trans;

    public void initiate() throws Exception {
        // Initiate functions and variables
        rbss = new ArrayList<>();
        CalcEditDist = new CalcEditDistance();
        CalcEditDist.initiate();
        HPCounter = new HairpinCounter();
        HPCounter.initiate();
        Reverse = new RevComp();
        Reverse.initiate();
        Trans = new Translate();
        Trans.initiate();

        // Read data file
        String data = FileUtils.readResourceFile("composition/data/coli_genes.txt");
        String options = FileUtils.readResourceFile("composition/data/rbs_options.txt");
        HashMap<String, String> geneList = new HashMap<>();

        // Scan through each line, add names of genes mapped to cds
        String[] lines = data.split("\\r|\\r?\\n");
        for (String line : lines) {
            if (line.isEmpty()) continue;
            String[] tabs = line.split("\t");
            geneList.put(tabs[1], tabs[6]);
        }

        // Scan through each line, add RBSOptions to rbss
        String[] linesR = options.split("\\r|\\r?\\n");
        for (String line : linesR) {
            if (line.isEmpty()) continue;
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

    public RBSOption run(String cds, Set<RBSOption> ignores) throws Exception {
        if (cds == null || cds.isEmpty()) {
            throw new IllegalArgumentException("Coding sequence is empty or null");
        }

        String first6AA = Trans.run(cds.substring(0, 30));
        List<RBSOption> curbest = new ArrayList<>();

        // Find the list of options with the lowest number of hairpins
        double minHairpins = Double.MAX_VALUE;
        for (RBSOption rbsOption : rbss) {
            if (ignores.contains(rbsOption)) continue;

            // Compute hairpin count for current RBS + CDS combination
            String rbsCdsCombo = rbsOption.getRbs() + cds;
            double hairpinCount = HPCounter.run(rbsCdsCombo);

            if (hairpinCount < minHairpins) {
                // Clear and update if new minimum found
                curbest.clear();
                curbest.add(rbsOption);
                minHairpins = hairpinCount;
            } else if (hairpinCount == minHairpins) {
                // Add if equal to current minimum
                curbest.add(rbsOption);
            }
        }

        if (curbest.isEmpty()) {
            throw new Exception("No viable RBS options. Try increasing the options in rbss or removing some ignores.");
        }

        // Select the best option based on amino acid sequence
        RBSOption best = curbest.get(0);
        for (int i = 1; i < curbest.size(); i++) {
            RBSOption current = curbest.get(i);
            int newDist = CalcEditDist.run(first6AA, current.getFirst6aas());
            int oldDist = CalcEditDist.run(first6AA, best.getFirst6aas());
            if (newDist < oldDist) {
                best = current;
            }
        }

        return best;
    }



    public static void main(String[] args) throws Exception {
        //Create an example
        String cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATTATCACCATCGTTATATGACACAGCTTGGGTGGCTATGGTGCCGGAGCGTAGTTCTTCTCAA";

        //Initiate the chooser
        RBSChooser2 chooser = new RBSChooser2();
        chooser.initiate();

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
        System.out.println(selected2.toString());
    }
}
