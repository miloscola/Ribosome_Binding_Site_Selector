package org.ucb.c5.tests;

import org.junit.BeforeClass;
import org.junit.Test;
import org.ucb.c5.composition.RBSChooser2;
import org.ucb.c5.composition.model.RBSOption;

import java.util.HashSet;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/** @author Milo */
public class TestRBSChooser2 {

    private static RBSChooser2 chooser;

    @BeforeClass
    public static void setUpClass() throws Exception {
        chooser = new RBSChooser2();
        chooser.initiate();
    }

    //tests an example with no hairpins
    @Test
    public void testCds1() throws Exception {
        String cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATTATCACCATCGTTATATGACACAGCTTGGGTGGCTATGGTGCCGGAGCGTAGTTCTTCTCAA";
        Set<RBSOption> ignores = new HashSet<>();
        RBSOption selected1 = chooser.run(cds, ignores);
        assertEquals(selected1.getName(), "yncE");
    }

    @Test
    public void testIgnores() throws Exception {
        String cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATTATCACCATCGTTATATGACACAGCTTGGGTGGCTATGGTGCCGGAGCGTAGTTCTTCTCAA";
        Set<RBSOption> ignores = new HashSet<>();
        RBSOption selected1 = chooser.run(cds, ignores);
        ignores.add(selected1);
        RBSOption selected2 = chooser.run(cds, ignores);
        assertEquals(selected2.getName(), "ppiB");
    }

    @Test
    public void testCds2() throws Exception {
        String cds = "ATGTGTGGCATACGCTAATACTTCATAATACTTGGGGGGATACTACTAGCCACTGAAACTGGCCTTGAATGTCGTGACAAGCCTGGTCGTCGG";
        Set<RBSOption> ignores = new HashSet<>();
        RBSOption selected1 = chooser.run(cds, ignores);
        assertEquals(selected1.getName(), "moaB");
    }

    @Test
    public void testCds3() throws Exception {
        String cds = "ATGCAATTCTAGGACCCCTCCTGACCGCCTTATGCTCCGCTTCGGGCCGTACTGATCTGGGTAAGTAAATTTTAAATGGGTAGTGGTACATGT";
        Set<RBSOption> ignores = new HashSet<>();
        RBSOption selected1 = chooser.run(cds, ignores);
        assertEquals(selected1.getName(), "dksA");
    }

}
