package org.ucb.c5.composition;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.ucb.c5.composition.model.RBSOption;
import org.ucb.c5.utils.FileUtils;

/**
 * Second generation RBSChooser algorithm
 *
 * Employs a list of genes and their associated ribosome binding sites for
 * highly-expressed proteins in E. coli.
 *
 * @author J. Christopher Anderson
 */
public class RBSChooser2 {

    private List<RBSOption> rbss;
    //TODO:  Fill in

    public void initiate() throws Exception {
        //Read the data file
        String data = FileUtils.readResourceFile("composition/data/coli_genes.txt");

        //TODO:  Fill in
    }

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
     */
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
