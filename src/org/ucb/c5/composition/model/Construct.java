package org.ucb.c5.composition.model;

import java.util.List;

/**
 * Encodes a genetic construct described in terms of a single operon
 * encoding multiple cocistronic transcripts
 * 
 * @author J. Christopher Anderson
 */
public class Construct {
    private final List<Transcript> mRNAs;
    private final String promoter;
    private final String terminator;

    public Construct(List<Transcript> mRNAs, String promoter, String terminator) {
        this.mRNAs = mRNAs;
        this.promoter = promoter;
        this.terminator = terminator;
    }
    
    public String toSeq() throws Exception {
        StringBuilder out = new StringBuilder();
        out.append(promoter);
        for(Transcript mrna : mRNAs) {
            out.append(mrna.toSeq());
        }
        out.append(terminator);
        return out.toString();
    }
}
