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
public class MirBase {
    
    String mature;
    String precursor;
    String mature_mirID;
    String precursor_mirID;
    String chr;
    String strand;
    int preStart;
    int preEnd;
    
    
    public MirBase(String mature, String precursor, String mature_mirID, String precursor_mirID, String chr, String strand, int preStart, int preEnd){
        this.mature = mature;
        this.precursor = precursor;
        this.mature_mirID = mature_mirID;
        this.precursor_mirID = precursor_mirID;
        this.chr = chr;
        this.strand = strand;
        this.preStart = preStart;
        this.preEnd = preEnd;
    }
    
    public String getMature(){
        return mature;
    }
    
    public String getPrecursor(){
        return precursor;
    }
    
    public String getMature_mirID(){
        return mature_mirID;
    }
    
    public String getPrecursor_mirID(){
        return precursor_mirID;
    }
    
    public String getChr(){
        return chr;
    }
    
    public String getStrand(){
        return strand;
    }
    
    
    public int getStart(){
        return preStart;
    }
    
    public int getEnd(){
        return preEnd;
    }
}
