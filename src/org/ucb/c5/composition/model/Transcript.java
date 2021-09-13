package org.ucb.c5.composition.model;

/**
 * Encodes a monocistronic mRNA from an RBS and a coding sequence
 * 
 * @author J. Christopher Anderson
 */
public class Transcript {
    private final RBSOption rbs;
    private final String peptide;
    private final String[] codons;

    public Transcript(RBSOption rbs, String peptide, String[] codons) {
        this.rbs = rbs;
        this.peptide = peptide;
        this.codons = codons;
    }

    public RBSOption getRbs() {
        return rbs;
    }

    public String getPeptide() {
        return peptide;
    }

    public String[] getCodons() {
        return codons;
    }
    
    public String toSeq() throws Exception {
        //Construct the mRNA sequence one codon at a time
        StringBuilder out = new StringBuilder();
        for(int i=0; i<peptide.length(); i++) {
            out.append(codons[i]);
        }
        
        //Add a stop codon to the coding sequence
        out.append("TAA");
            
        //To visually distinguish, make the RBS lowercase and the CDS uppercase
        return rbs.getRbs().toLowerCase() + out.toString().toUpperCase();
    }
}
