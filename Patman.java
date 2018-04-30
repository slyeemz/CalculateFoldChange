/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CalculateFoldChange;

/**
 *
 * @author twh14ura
 */
public class Patman {
    
    String sequence;
    String chr;
    String strand;
    int abundance;
    int start;
    int end;
    
    
    public Patman(String sequence, String chr, String strand, int abundance, int start, int end){
        this.sequence = sequence;
        this.chr = chr;
        this.strand = strand;
        this.abundance = abundance;
        this.start = start;
        this.end = end;
    }
    
    public String getSequence(){
        return sequence;
    }
    
    public String getChr(){
        return chr;
    }
    
    public String getStrand(){
        return strand;
    }
    
    public int getAbund(){
        return abundance;
    }
    
    public int getStart(){
        return start;
    }
    
    public int getEnd(){
        return end;
    }
}
