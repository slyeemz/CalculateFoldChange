/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CalculateFoldChange;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;
import java.util.stream.Stream;

/**
 *
 * @author twh14ura
 */
public class CalculateFoldChange {
    private final File runDir; // also output directory 
    private final File wtsrnaomeFile;
    private final File dicerFile;
    private final File paresnipFile;
    private final File mircat2DefFile;
    private final File mircat2LooseFile;
    private final File mirbaseFile;
    private final File mirbaseStemloopFile;
    private int total;  // dummy variable to hold some values
    private final int srna_redundant_size;
    private final int dicer_redundant_size;
    //private final int dpredicted_mirnas_size; // size of miRCat2 miRNA candidates using default parameters
    //private final int lpredicted_mirnas_size; // size of miRCat2 miRNA candidates using loose parameters
    private boolean known; // if sRNA is known miRNA
    private boolean isomir; // if sRNA is isomir
    private String mirID; // the mirID for sRNA if it was known, isomir, or found in mature miRNA precursor 
    private final HashMap<String,Integer> srnas; // WT srnas
    private final HashMap<String,Integer> dicers;
    private final HashMap<String,String> miRBasemirnas; // mature miRNAs
    private final HashMap<String,String> miRBase_stemloop;    // mature miRNAs precursors
    private final HashMap<String,Integer> sRNAwTarget;    // PAREsnip target predictions 
    private final HashSet<String> default_predicted_mirnas; // miRCat2 miRNA candidates using default parameters
    private final HashSet<String> loose_predicted_mirnas; // miRCat2 miRNA candidates using loose parameters
    private final HashSet<String> miRBase_sRNAs; // to store all sRNA that mapped to mature miRNA or precursor
    private final HashSet<String> DEsrnas; // Differentially expressed sRNAs
    //private int iteration;
    private final int[] truth_table;
    private int intersectionCount;
    private int fc_up_wt;
    private int fc_up_dc;
    private int miRBaseCountWT;
    private int miRBaseCountDcl1;
    private int miRBaseCountWTandDcl1;
    private HashSet<Patman> wtPatman;
    private HashSet<Patman> dcPatman;
        
    
    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     */
    public static void main(String[] args) throws Exception{
        // TODO code application logic here  
        /*
        String runDir = "C:/Users/twh14ura/Documents/NetBeansProjects/CalculateFoldChange/WTA";
        String wt_sRNA_file =  runDir + "/WTA.wtPatman.fa";
        String paresnipFile = runDir + "/WTA_COL_A_Category012_PAREsnipGUIout.csv";
        //String paresnipFile = runDir + "/paresnip2_result_conservedsRNAs/WTA_COL_A.csv";        
        String mircat2DefaultFile = runDir + "/WTA_miRCat2_defaultParameter_predicted_miRNA.fa";
        String mircat2LooseFile = runDir + "/WTA_miRCat2_looseParameter_predicted_miRNA.fa";
        
        /*
        String runDir = "C:\\Users\\twh14ura\\Documents\\NetBeansProjects\\CalculateFoldChange\\WTB_DcB";
        String wt_sRNA_file =  runDir + "\\WTB.wtPatman.fa";
        String paresnipFile = runDir + "\\WTB_COL_B_Category012_PAREsnipGUIout.csv";
        String mircat2DefaultFile = runDir + "\\WTB_miRCat2_defaultParameter_predicted_miRNA.fa";
        String mircat2LooseFile = runDir + "\\WTB_miRCat2_looseParameter_predicted_miRNA.fa";
        
        /*
        String runDir = "C:/Users/twh14ura/Documents/NetBeansProjects/CalculateFoldChange/WTC_DcC";
        String wt_sRNA_file =  runDir + "/WTC.wtPatman.fa";
        String paresnipFile = runDir + "/WTC_COL_C_Category012_PAREsnipGUIout.csv";
        String mircat2DefaultFile = runDir + "/WTC_miRCat2_defaultParameter_predicted_miRNA.fa";
        String mircat2LooseFile = runDir + "/WTC_miRCat2_looseParameter_predicted_miRNA.fa";
        */
        
        //String Dcl1_file = "C:/Users/twh14ura/Documents/NetBeansProjects/CalculateFoldChange/Dcl1_data/Dcl1A.fa";
        //String Dcl1_file = "C:/Users/twh14ura/Documents/NetBeansProjects/CalculateFoldChange/Dcl1_data/Dcl1B.fa";
        //String Dcl1_file = "C:/Users/twh14ura/Documents/NetBeansProjects/CalculateFoldChange/Dcl1_data/Dcl1C.fa";
        
        if (args.length != 6) {
            String exceptionStr = "";
            exceptionStr += "Expected 6 arguments:\n";
            exceptionStr += "args[0]: runDir\n";
            exceptionStr += "args[1]: sRNAome\n";
            exceptionStr += "args[2]: Dicer file\n";
            exceptionStr += "args[3]: PAREsnip file\n";              
            exceptionStr += "args[4]: miRCat2 default file\n";
            exceptionStr += "args[5]: miRCat2 loose file\n";

            System.err.println(exceptionStr);
            throw new Exception(exceptionStr);
        } else {
            //if(args.length == 6)
            CalculateFoldChange CFC = new CalculateFoldChange(args[0], new File(args[1]), new File(args[2]), new File(args[3]), new File(args[4]), new File(args[5]));
        
            //CalculateFoldChange CFC = new CalculateFoldChange(runDir ,new File(wt_sRNA_file), new File(Dcl1_file),new File(paresnipFile),new File(mircat2DefaultFile),new File(mircat2LooseFile));
            CFC.start(0);
            
            args[2] = "/gpfs/scratch/twh14ura/Data/sRNA/DCL1/Aligned/DCL1B.fa.patman.fa";
            CFC = new CalculateFoldChange(args[0], new File(args[1]), new File(args[2]), new File(args[3]), new File(args[4]), new File(args[5]));
            CFC.start(1);
            
            args[2] = "/gpfs/scratch/twh14ura/Data/sRNA/DCL1/Aligned/DCL1C.fa.patman.fa";
            CFC = new CalculateFoldChange(args[0], new File(args[1]), new File(args[2]), new File(args[3]), new File(args[4]), new File(args[5]));
            CFC.start(2);
        }
    }
    
    public void start(int i) throws FileNotFoundException , IOException {
        
        File outDir = new File(this.runDir.getAbsolutePath() + "/output" + i);
        outDir.mkdirs();
              
        //File infoFile = new File(outDir.getAbsolutePath() + "/information.txt");
        try (PrintWriter writer = new PrintWriter(new File(outDir.getAbsolutePath() + "/information.txt"))) {
            writer.println("This output table shows:");
            writer.println("WT abundance, Norm WT abund, Dicer abundance, Norm Dicer abund, FC, Log2(FC) , DE, miRBase, isomir, mirID, Target, Default miRCat2, Loose miRCat2");
            writer.println();
            writer.println("sRNAome file: " + this.wtsrnaomeFile.getAbsolutePath());
            writer.println("Dicer file: " + this.dicerFile.getAbsolutePath());
            writer.println("PAREsnip file: " + this.paresnipFile.getAbsolutePath());
            writer.println("miRCat2 default param file: " + this.mircat2DefFile.getAbsolutePath());
            writer.println("miRCat2 loose param file: " + this.mircat2LooseFile.getAbsolutePath());
            writer.println("miRBase file: " + this.mirbaseFile.getAbsolutePath());
            writer.println("miRBase strm-loop file: " + this.mirbaseStemloopFile.getAbsolutePath());
            writer.println("output dir: " + outDir.getAbsolutePath());
            writer.println("summary file: " + outDir.getAbsolutePath()+"/summary.txt");
        }
        
        outputsRNAvsDicer(outDir);
    }
    
    public CalculateFoldChange(String runDir, File wt_sRNA_file, File dicer_file, File paresnipFile, File mircat2DefaultFile, File mircat2LooseFile) throws FileNotFoundException {
        this.wtsrnaomeFile = wt_sRNA_file;
        this.dicerFile = dicer_file;
        this.paresnipFile = paresnipFile;
        this.mircat2DefFile = mircat2DefaultFile;
        this.mircat2LooseFile = mircat2LooseFile;
        this.mirbaseFile = new File( "/gpfs/scratch/twh14ura/Data/miRBase/mature.fa");
        this.mirbaseStemloopFile = new File( "/gpfs/scratch/twh14ura/Data/miRBase/miRBase_stem-loop.fa");
        //        /scratch/Salma/CalculateFoldChange/
        //      C:/Users/twh14ura/Documents/NetBeansProjects/CalculateFoldChange/
        this.runDir =  new File(runDir);
        this.srnas = getSmallRNAs(wt_sRNA_file);
        this.srna_redundant_size = this.total;
        this.dicers = getSmallRNAs(dicer_file);
        this.dicer_redundant_size = this.total;
        this.miRBasemirnas = readFileMirBase(this.mirbaseFile);
        this.miRBase_stemloop = readFileMirBase(mirbaseStemloopFile);
        //PAREsnip or PAREsnip2
        this.sRNAwTarget = getSmallRNAsFromPAREsnip2File(paresnipFile);// getSmallRNAsFromPAREsnipFileGUI(paresnipFile, this.srnas);
        this.default_predicted_mirnas = readFile(mircat2DefaultFile);
        //this.dpredicted_mirnas_size = this.total;
        this.loose_predicted_mirnas = readFile(mircat2LooseFile);
        //this.lpredicted_mirnas_size = this.total;
        this.miRBase_sRNAs = new HashSet<>();
        this.truth_table = new int[64];
        this.DEsrnas = new HashSet();
    }

    public void outputsRNAvsDicer(File outDir)  throws FileNotFoundException , IOException { //(HashMap<String, Integer> sRNAs, HashMap<String, Integer> dicers, HashMap<String, Integer> miRBase) throws FileNotFoundException {
//        int fc_up_wt = 0;
//        int fc_up_dc = 0;
//        int intersectionCount = 0;
//        int miRBaseCountWT = 0;
//        int miRBaseCountDcl1 = 0;
//        int miRBaseCountWTandDcl1 = 0;
        boolean DE;      // Differential expression is significant or not
        try (PrintWriter writer = new PrintWriter(outDir + "/output_" + this.wtsrnaomeFile.getName() + "_vs_" + this.dicerFile.getName() + "_FC_table.csv")) {
            writer.println("sRNAs , Length , WT abundance , Norm WT abund(MTC) , Dicer abundance , Norm Dicer abund(MTC) , Log2(FC) , DE , miRBase , isomir , mirID , Target , Category , Loose miRCat2"); // , Secondary structure");
            // norm abundance = sRNA abundance/total number of sRNAs x 1,000,000
            Integer wt;
            Integer da ;
            double norm_sa;
            double norm_da;
            //double fc;
            double log2fc;
            float mtc = (this.srna_redundant_size + this.dicer_redundant_size) / 2 ; // median of total count
            System.out.println("srna_redundant_size = " + this.srna_redundant_size);
            System.out.println("dicer_redundant_size = " + this.dicer_redundant_size);
            System.out.println("MTC = " + mtc);
            
            for(Patman p : this.wtPatman){
                String srna = p.getSequence();
//                wt = this.srnas.get(srna);
                wt = p.getAbund();
                da = this.dicers.remove(srna);  // Dicer abundance
                
                miRBase(p);
                if(this.known){
                    miRBaseCountWT++;
                    this.miRBase_sRNAs.add(srna);
                }
                if(this.isomir){
                    this.miRBase_sRNAs.add(srna);
                }
                if(da == null ){
                    da = 0;
                }else{
                    intersectionCount++;
                    if(this.known){
                        miRBaseCountDcl1++;
                        miRBaseCountWTandDcl1++;
                    }
                }
                
                norm_sa = ((double)(wt)/this.srna_redundant_size) * mtc;//5173806.0; // MTC (median of total count): median of WT triplicate
                norm_da = ((double)(da)/this.dicer_redundant_size) *  mtc; // 12296993.0; // median of Dcl1 triplicate

                   //fc = norm_da/norm_sa;
                log2fc = (log2(norm_da+10)-log2(norm_sa+10)); // add offset after notmalization as in C'seq thesis
                if(log2fc >= 1){
                    DE = true;
                    fc_up_dc++;
                    this.DEsrnas.add(srna);
                }else if(log2fc <= -1){
                    DE = true;
                    fc_up_wt++;
                    this.DEsrnas.add(srna);
                }else{
                    DE = false;
                }
              
                writer.println(srna + "," + srna.length() + "," + wt + "," + norm_sa + "," + da + "," + norm_da + "," + log2fc + "," + DE + "," + this.known + "," + this.isomir + "," + this.mirID + "," + hasTarget(srna) + "," + this.sRNAwTarget.get(srna) + "," + loosemiRcat2(srna)); // + " , " + (this.secondaryLfold.get(srna)) );
                
                //}
            }
            
            for(String srna : this.dicers.keySet()){
                miRBase(srna);
                if(this.known){
                    miRBaseCountDcl1++;
                    this.miRBase_sRNAs.add(srna);
                }
                if(this.isomir){
                    this.miRBase_sRNAs.add(srna);
                }
                wt = 0;
                da = this.dicers.get(srna);
                
                norm_sa = ((double)(wt)/this.srna_redundant_size)  * mtc; //*5173806.0; // MTC (median of total count): median of WT triplicate
                norm_da= ((double)(da)/this.dicer_redundant_size) *  mtc; // 12296993.0; // median of Dcl1 triplicate
                
                //fc = norm_da/norm_sa;
                log2fc = (log2(norm_da+10)-log2(norm_sa+10));
                if(log2fc >= 1){
                    DE = true;
                    fc_up_dc++;
                    this.DEsrnas.add(srna);
                }else if(log2fc <= -1){
                    DE = true;
                    fc_up_wt++;
                    this.DEsrnas.add(srna);
                }else{
                    DE = false;
                }

                writer.println(srna + "," + srna.length() + "," + wt + "," + norm_sa + "," + da + "," + norm_da +  "," + log2fc + "," + DE + "," + this.known + "," + this.isomir + "," + this.mirID + "," + false + "," + null  + "," + false); //  + " , " + (this.secondaryLfold.get(srna)) );
            }
            System.out.println("Generated: " + outDir + "/output_" + this.wtsrnaomeFile.getName() + "_vs_" + this.dicerFile.getName() + "_FC_table.csv");
        } // , Secondary structure"); // , Secondary structure"); // , Secondary structure"); // , Secondary structure"); // , Secondary structure"); // , Secondary structure"); // , Secondary structure"); // , Secondary structure");

        populateTruthTable();
        outputSummary(outDir);
    }  
    

    public void outputSummary(File outDir) throws FileNotFoundException{
        try (PrintWriter writer = new PrintWriter(outDir + "/summary.txt")) {
            
            writer.println("Total WT sRNAs: "   + this.srna_redundant_size + " , unique: " + this.srnas.size());
            writer.println("Total Dicer sRNAs: " + this.dicer_redundant_size + " , unique: " + this.dicers.size()); //+ totaldicer  + " ,, "
            writer.println("Common sRNAs between WT and Dcl: " + this.intersectionCount);
            writer.println("Up regulated in WT: " + this.fc_up_wt);
            writer.println("Up regulated in Dicer: " + this.fc_up_dc);
            writer.println("Known miRNAs in WT: " + this.miRBaseCountWT);
            writer.println("Known miRNAs in Dicer: " + this.miRBaseCountDcl1);
            writer.println("Common known miRNAs in WT and Dicer: " + this.miRBaseCountWTandDcl1);
            writer.println("Known and isomiR miRNAs in WT and Dicer: " + this.miRBase_sRNAs.size());
            writer.println("Functional sRNAs: " + this.sRNAwTarget.size());
            writer.println("Predicted miRNAs (default): " + this.default_predicted_mirnas.size());
            writer.println("Predicted miRNAs (loose): " + this.loose_predicted_mirnas.size());
            writer.println("===============================");
            
  //retu(sRNAome) + parseBoolean(dicer) * 2 + parseBoolean(miRCat) * 4 + parseBoolean(pareSnip) * 8 + parseBoolean(DE) * 16 + parseBoolean(miRBase) * 32;
    
        
            writer.println("| miRBase | DE | PAREsnip | miRCat | Dicer | sRNAome | N ");
            writer.println("--------------------------------------------");
            
            for (int i = 0; i < 64; i++) {
                String binaryStr = Integer.toBinaryString(i);
                while (binaryStr.length() < 6) {
                    binaryStr = "0" + binaryStr;
                }
                
                writer.println(String.format("|   %s    |    %s     |   %s    |    %s     |    %s    |    %s    | %d ", binaryStr.charAt(0), binaryStr.charAt(1), binaryStr.charAt(2), binaryStr.charAt(3), binaryStr.charAt(4), binaryStr.charAt(5), this.truth_table[i]));
                
            }
        }
    }
    public void miRBase(Patman p){
        String s = p.getSequence();
        this.isomir = false;
        this.known = false;
        this.mirID = "Unknown";
        boolean x = true;
        if(this.miRBasemirnas.containsKey(s)){
            this.known = true;
            mirID = this.miRBasemirnas.get(s);
        }else{
            for(String k : this.miRBasemirnas.keySet()){
                String v = this.miRBasemirnas.get(k);
                if(k.contains(s)){
                    this.isomir = true;
                    mirID = v;
                    x = false;
                    break;
                }
            }
            if(x){
                this.miRBase_stemloop.forEach((k2,v2) -> {
                    if(k2.contains(s)){
                        this.mirID = v2;
                        this.miRBasemirnas.forEach((k,v)-> {
                            if(s.contains(k) && v2.startsWith(v)){
                                this.isomir = true;
                                this.mirID = v;
                                return;
                            }else if(  (   s.contains(k.substring(2)) || s.contains(k.substring(0,k.length()-2))   ) && v2.startsWith(v)  ){
                                this.isomir = true;
                                this.mirID = v;
                                return;
                            }else if(  (   k.contains(s.substring(2)) || k.contains(s.substring(0,s.length()-2))   ) && v2.startsWith(v)  ){
                                this.isomir = true;
                                this.mirID = v;
                                return;
                            }
                        });
                        return;
                    }
                });
            }
        }
    }
    
    public boolean hasTarget(String s){
        return this.sRNAwTarget.containsKey(s);
    }
    
//    public boolean defaultmiRcat2(String seq){
//        return this.default_predicted_mirnas.contains(seq);
//    }
    
    public boolean loosemiRcat2(String s){
        return this.loose_predicted_mirnas.contains(s);
    }

    public static double log2(double n){
        return (Math.log(n) / Math.log(2));
    }

    public void readPatamnFile(String patmanFileName) throws IOException{
//        HashMap<String,Integer> sequences = new HashMap<>();
        HashSet<Patman> patman = new HashSet();
        //read file into stream, try-with-resources
//        try (Stream<String> stream = Files.lines(Paths.get(patmanFileName))) {
//            stream.forEach(line -> { 
//                if(!line.isEmpty()){
//                    if(line.split("\t").length>1){
//                        String sa = line.split("\t")[1].trim(); //sRNA (and abundance)
//                        //String sa = "gra(2)";
//                        String s = "";
//                        int a = 0;
//                        if (sequences.containsKey(s)) {
//                            a = Integer.parseInt(sequences.get(s).split(",")[0]);
//                        }
//                        if (sa.contains("(") && sa.contains(")") && sa.endsWith(")")){
//                            s = sa.split("\\(|\\)")[0].toUpperCase().replaceAll("U", "T").trim();
//                            try{
//                                a = Integer.parseInt(sa.split("\\(|\\)")[1].trim());
//                            } catch (NumberFormatException ex){}
//                        }else{
//                            s = sa.toUpperCase().replaceAll("U", "T");
//                            a++;
//                        }
//                        sequences.put(s, a);
//                    }
//                }
//            });
//        }
        try (Stream<String> stream = Files.lines(Paths.get(patmanFileName))) {
            stream.forEach(line -> { 
                if(!line.isEmpty()){
                    if(line.split("\t").length>1){
                        String sa = line.split("\t")[1].trim(); //sRNA (and abundance)
                        //String sa = "gra(2)";
                        String seq = "";
                        String chr = "";
                        String strand = "";
                        int a = 0;
                        int s = 0;
                        int e = 0;
//                        if (sequences.containsKey(seq)) {
//                            a = Integer.parseInt(sequences.get(seq).split(",")[0]);
//                        }
                        if (sa.contains("(") && sa.contains(")") && sa.endsWith(")")){
                            seq = sa.split("\\(|\\)")[0].toUpperCase().replaceAll("U", "T").trim();
                            chr = line.split("\t")[0].trim();
                            strand = line.split("\t")[4].trim();
                            try{
                                a = Integer.parseInt(sa.split("\\(|\\)")[1].trim());
                                s = Integer.parseInt(line.split("\t")[2].trim());
                                e = Integer.parseInt(line.split("\t")[3].trim());
                            } catch (NumberFormatException ex){}
                        }else{
                            System.out.println("Wrong formatting for Patman file.");
                        }
                        //sequences.put(seq, a+","+line);
                        patman.add(new Patman(seq,chr,strand,a,s,e));
                    }
                }
            });
        }
    }   
    
    public final HashMap getSmallRNAs(File sample) throws FileNotFoundException {
        HashMap<String,Integer> smallRNAs = new HashMap<>();
        this.total = 0;
        
        try (Scanner in = new Scanner(sample)) {
            String header = "";
            while (in.hasNextLine()) {
                String line = in.nextLine().trim();
                if (!line.startsWith(">") && !line.isEmpty()) {
                    line = line.toUpperCase().replaceAll("U", "T");
                    int abundance = 0;
                    if (smallRNAs.containsKey(line)) {
                        abundance = smallRNAs.get(line);
                    }
                    if (header.contains("(") && header.contains(")") && header.endsWith(")")) {
                        abundance += Integer.parseInt(header.split("[()]")[1]);
                        this.total += abundance;
                    } else {
                        abundance++;
                        this.total++;
                    }
                    smallRNAs.put(line, abundance);
                } else {
                    header = line.substring(1).trim();
                }
            }
        }
        return smallRNAs;
    }
    
    private HashSet readFile(File sample) throws FileNotFoundException {
        //this.total = 0;
        HashSet<String> sequences = new HashSet<>();
        try(Scanner in = new Scanner(sample)){
            while (in.hasNextLine()) {
                String line = in.nextLine().trim();
                if (!line.startsWith(">") && !line.isEmpty()) {
                    line = line.toUpperCase().replaceAll("U", "T");
                    sequences.add(line);
                    //this.total++;
                } 
            }
            in.close();
        }
        return sequences;
    }
    
    public static HashMap<String,String> readFileMirBase(File sample) throws FileNotFoundException {
        HashMap<String,String> sequences = new HashMap<>();
        try(Scanner in = new Scanner(sample)){
            String line;
            String mir;
            while (in.hasNextLine()) {
                line = in.nextLine().trim();
                if (line.startsWith(">")) {
                    mir = line.substring(1, line.indexOf(" ")).trim();
                    line = in.nextLine().trim();
                    line = line.toUpperCase().replaceAll("U", "T");
                    sequences.put(line,mir);
                } 
            }
            in.close();
        }
        return sequences;
    }
 
    public static HashMap<String,Integer> getSmallRNAsFromPAREsnip2File(File paresnipFile) throws FileNotFoundException {
        HashMap<String, Integer> srnas = new HashMap<>();
        try(Scanner in = new Scanner(paresnipFile)){
            int counter = 0;
            while (in.hasNextLine()) {
                String line = in.nextLine().trim();
                if (!line.startsWith("Record")) {

                    if (!line.isEmpty() && counter % 3 == 0) {
                        String[] fields = line.split("\"");
                        String srna = fields[1];
                        srna = srna.substring(srna.indexOf(">")+1);
                        srna = srna.toUpperCase().replaceAll("U", "T").replaceAll("-", "");
                        srna = srna.trim();
                        //int paresnipAbundance = Integer.parseInt(line.split(",")[7].replaceAll("[\"]", ""));
                        
                        in.nextLine();
                        line = in.nextLine();
                        counter += 2;

                        int cat = Integer.parseInt(fields[4].split(",")[1]);
                        if (srnas.containsKey(srna)) {
                            if (cat < srnas.get(srna)) {
                                srnas.put(srna, cat);
                            }
                        } else {
                            srnas.put(srna, cat);
                        }
                    }
                    if (!line.isEmpty()) {
                        counter++;
                    }
                }
            }
            //srnas.keySet().forEach( (srna,v) -> justSRNAs.put(srna,v));
            in.close();
        }
        return srnas;
    }
    
        private void clearTruthTable() {
        for (int i = 0; i < 64; i++) {
            this.truth_table[i] = 0;
        }
    }

    public static int parseBoolean(boolean b) {
        return b ? 1 : 0;
    }

    private int getIndex(boolean sRNAome, boolean dicer,  boolean miRCat, boolean pareSnip, boolean DE, boolean miRBase) {
        return parseBoolean(sRNAome) + parseBoolean(dicer) * 2 + parseBoolean(miRCat) * 4 + parseBoolean(pareSnip) * 8 + parseBoolean(DE) * 16 + parseBoolean(miRBase) * 32;
    }

    public void populateTruthTable() {

//        this.sRNAomeSize = sRNAome_sRNAs.size();
//        this.mirbaseSize = miRBase_sRNAs.size();
//        this.mircatSize = miRCat_sRNAs.size();
//        this.paresnipSize = paresnip_sRNAs.size();
        
        ArrayList<Set<String>> list = new ArrayList<>();
        list.add(this.srnas.keySet());
        list.add(this.dicers.keySet());
        list.add(this.DEsrnas);
        list.add(this.miRBase_sRNAs);
        list.add(this.sRNAwTarget.keySet());
        list.add(this.loose_predicted_mirnas);
        HashSet<String> all_sRNAs = getUnion(list);

        clearTruthTable();
        for (String sRNA : all_sRNAs) {
            boolean inSmallRNAome = this.srnas.containsKey(sRNA);
            boolean inDicer = this.dicers.containsKey(sRNA);
            boolean DE = this.DEsrnas.contains(sRNA);
            boolean inMiRBase = this.miRBase_sRNAs.contains(sRNA);
            boolean inPAREsnip = this.sRNAwTarget.containsKey(sRNA);
            boolean inMiRCat = this.loose_predicted_mirnas.contains(sRNA);
            int index = getIndex(inSmallRNAome, inDicer, DE, inMiRCat, inPAREsnip, inMiRBase);
            this.truth_table[index]++;
        }
    }

    public static HashSet<String> getUnion(ArrayList<Set<String>> lists) {
        HashSet<String> unionSet = new HashSet<>();
        lists.forEach((list) -> {
            list.stream().filter((item) -> (!unionSet.contains(item))).forEachOrdered((item) -> {
                unionSet.add(item);
            });
        });
        return unionSet;

    }
    
    public static HashMap<String,Integer> getSmallRNAsFromPAREsnipFileGUI(File paresnipFile) throws FileNotFoundException {
        //HashSet justSRNAs = new HashSet<>();
        HashMap<String, Integer> srnas = new HashMap<>();
        try(Scanner in = new Scanner(paresnipFile)){
            int counter = 0;
            while (in.hasNextLine()) {
                String line = in.nextLine().trim();
                if (!line.startsWith("Select")) {

                    if (!line.isEmpty() && counter % 3 == 0) {
                        String[] fields = line.split("\"");
                        String srna = fields[17];
                        srna = srna.substring(srna.indexOf("5'") + 2, srna.lastIndexOf("3'"));
                        srna = srna.toUpperCase().replaceAll("U", "T").replaceAll("-", "");
                        srna = srna.trim();

                        in.nextLine();
                        line = in.nextLine();
                        counter += 2;
                        
                        /*int paresnipAbundance = Integer.parseInt(line.split("\"")[6].replaceAll("[\"]", ""));
                        // is the srna in the srnaome?
                        if (!srnaome.containsKey(srna) || paresnipAbundance != srnaome.get(srna)) {
                            String bestMatch = "";
                            int bestMatchDiff = 0;

                            // loop through srnas in srnaome and find what we think is the sequence of intestest
                            for (String srna : srnaome.keySet()) {
                                if (srna.contains(srna) && srna.startsWith(srna) && paresnipAbundance == srnaome.get(srna)) {
                                    int sizeDiff = srna.length() - srna.length();
                                    if (sizeDiff < bestMatchDiff || bestMatch.isEmpty()) {
                                        bestMatch = srna;
                                        bestMatchDiff = sizeDiff;
                                    }
                                }
                            }
                            if (!bestMatch.isEmpty()) {
                                srna = bestMatch;
                            } else {
                                if (srnaome.containsKey(srna)) {
                                    //System.out.println(srna.length());
                                    //System.out.println("it was in the srnaome but not same abundance");
                                }
                                //System.out.println(srna.length());
                                //System.out.println("looking for: " + paresnipAbundance);
                                srna="";
                                //throw new Exception("Could not find sRNA from interaction duplex in sRNA file!");
                            }
                        }*/
                        int cat = Integer.parseInt(fields[5].replaceAll("[\"]", ""));
                        if (srnas.containsKey(srna)) {
                            if (cat < srnas.get(srna)) {
                                srnas.put(srna, cat);
                            }
                        } else {
                            srnas.put(srna, cat);
                        }
                    }
                    if (!line.isEmpty()) {
                        counter++;
                    }
                }
            }
            //srnas.keySet().forEach( srna -> justSRNAs.add(srna));
            in.close();
        }
        return srnas;
    }
   
    /*
    
    public String mapReads(File pattern, File database) throws IOException {
        String alignedFile = this.runDir.getAbsolutePath() + File.separator + pattern.getName() + ".wtPatman";

        File patman_out = new File(alignedFile);
        String[] command = new String[]{"wtPatman" , "-e" , Integer.toString(0), "-P",  pattern.getAbsolutePath() , "-D", database.getAbsolutePath() , "-o", patman_out.getAbsolutePath()};

        //File tmpDir = new File("tmp");
        
        Process process = Runtime.getRuntime().exec(command);

        return alignedFile;
    }
    */
    
}
        
        
        //miRBAse annotation for sRNAs WT
        /*for(String srna : this.srnas.keySet()){
            this.miRBasemirnas.forEach((srna,v)->{
            if(srna.contains(srna)){
                sRNAs_annotation.put(srna,v);
            }else if(srna.contains(srna)){
                this.miRBase_stemloop.forEach((k1,v1)->{
                    if(k1.contains(srna) && v.equals(v1)){
                        sRNAs_annotation.put(srna,v1);
                    }
                });
            }else{
                this.miRBase_stemloop.forEach((k1,v1)->{
                    if(k1.contains(srna)){
                        isomirs.put(srna,v1);
                    }
                });
            }
            });
        }*/



        //miRBAse annotation for sRNAs Dicer
        /*this.dicers.forEach((srna,z) ->{
            this.miRBasemirnas.forEach((srna,v)->{
                if(srna.contains(srna)){
                    sRNAs_annotation.put(srna,v);
                }else if(srna.contains(srna)){
                    this.miRBase_stemloop.forEach((k1,v1)->{
                        if(k1.contains(srna) && v.equals(v1)){
                            sRNAs_annotation.put(srna,v1);
                        }
                    });
                }else{
                    this.miRBase_stemloop.forEach((k1,v1)->{
                        if(k1.contains(srna)){
                            isomirs.put(srna,v1);
                        }
                    });
                }
            });
        });*/
        
        //----------
        /*public String[] RNALfold(String seq) throws IOException {
        String[] secondarymfe = {"","1000","1000"};
        this.secondaryLfold.forEach((srna,v) -> {
            if(srna.equals(seq)){
                secondarymfe[0] = v.substring(0, v.indexOf(" "));
                secondarymfe[1] = v.substring(v.indexOf(" ")).split("[()]")[1].trim();
                secondarymfe[2] = v.substring(v.lastIndexOf(")")+1).trim();
            }
        });
        return secondarymfe;
    }*/
    
        /*Runtime rt = Runtime.getRuntime();
        Process pr = rt.exec("C:\\release_MainApp\\ExeFiles\\win\\RNALfold.exe", null, new File("C:\\release_MainApp\\ExeFiles\\win"));
        OutputStreamWriter osw = new OutputStreamWriter(pr.getOutputStream());
        //System.out.println("Sequences is : "+seq);
        osw.write(seq);
        osw.close();
        
        String[] secondarymfe = {"","10"};
        try(BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));){
        //String[] fields = new String[2];
        
        
        String line=null;
        
        while((line=input.readLine()) != null) {
            if (line.startsWith(".") || line.startsWith("(.") || line.startsWith("((") || line.startsWith("()")) {
                if (!line.isEmpty()) {
                    //System.out.println("\n"+line);
                    secondarymfe[0] = line.substring(0, line.indexOf(" "));
                    secondarymfe[1] = line.substring(line.indexOf(" ")).split("[()]")[1].trim(); 
                    //break;
                }
            }
        }
        input.close();
        
        }*/
        //pr.getInputStream().close();

    /* format of RNALfold output:
        .((((........)))) ( -0.40)    7
        AAAACAGGUGGAUUCUCUCUCGC
        ( -0.40)
    
    public HashMap<String,String> readSecondaryFile(File secondaryFile){
        HashMap<String,String> secondaryLfold = new HashMap<String,String>();
        
        try (BufferedReader br = new BufferedReader(new FileReader(secondaryFile))) {
            String line;
            String secondary = "";
            String srna = "";
            while ((line = br.readLine()) != null) {
                if (line.startsWith(".") || line.startsWith("(.") || line.startsWith("((") || line.startsWith("()")){
                    secondary += "\""+line+"\"";
                    while((line = br.readLine()) != null){
                        if (line.startsWith(".") || line.startsWith("(.") || line.startsWith("((") || line.startsWith("()")) {
                            secondary += "\""+line+"\"";
                        }else { break;}
                    }
                }
                if (line.startsWith("A") || line.startsWith("U") || line.startsWith("G") || line.startsWith("C")) {
                    srna = line.toUpperCase().replaceAll("U", "T").trim();
                    secondaryLfold.put(srna, secondary);
                    //System.out.println(secondaryLfold.get(srna));
                    secondary = "";
                    br.readLine();
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
	}
        return secondaryLfold;
    }
    */


//            this.srnas.forEach((srna,a)->{
//                
//                if(srna.contains(srna)){
//                    known_mirnas.put(srna, v);
//                }else{
//                    this.miRBase_stemloop.forEach((k2,v2)->{
//                        if(k2.contains(srna)){
//                            if(srna.contains(srna) && v.equals(v2)){
//                                known_mirnas.put(srna,v);
//                            }else{
//                                isomirs.put(srna,v2);
//                            }
//                        }
//                    });
//                }
//            });
//            this.dicers.forEach((srna,a)->{
//                if(srna.contains(srna)){
//                    known_mirnas.put(srna, v);
//                }else{
//                    this.miRBase_stemloop.forEach((k2,v2)->{
//                        if(k2.contains(srna)){
//                            if(srna.contains(srna) && v.equals(v2)){
//                                known_mirnas.put(srna,v);
//                            }else{
//                                isomirs.put(srna,v2);
//                            }
//                        }
//                    });
//                }
//            });



/*for(String srna : srnas.keySet()){
            if(sRNAs_annotation.containsKey(srna)){
                mb = true;
                isomir = false;
                mirID = sRNAs_annotation.get(srna);
                miRBaseCount++;
            }else {
                mb = false;
                mirID = "Unknown";
                if(isomirs.containsKey(srna)){
                    isomir = true;
                    mirID = isomirs.get(srna);
                }
            }
            // To check if sRNA has target
            if(sRNAwTarget.contains(srna)){
                hasTarget = true;
            }
            // To check if miRCat2 predicted this sRNA
            if(this.default_predicted_mirnas.contains(srna)){
                defPredicted = true;
            }
            if(this.loose_predicted_mirnas.contains(srna)){
                loosePredicted = true;
            }
           
            Integer da = this.dicers.remove(srna);  // Dicer abundance
            if(da == null ){
                da = 0;
                DE = false;
                writer.println(srna + " , " + srna.length() + " , " + this.srnas.get(srna) + " , " + (double)(this.srnas.get(srna))/this.srna_redundant_size*1000000.0 + " , " + da + " , " + 0 + " , " + 0 + " , " + 1000 + " , " + DE + " , " + mb + " , " + isomir + " , " + mirID + " , " + hasTarget + " , " + defPredicted + " , " + loosePredicted); //  + " , " + (this.secondaryLfold.get(srna)) );
            }else{
                //if(sRNAs.get(srna)>1 && sRNAs.get(srna)<da){    
                float fc = (float)da/(float)(this.srnas.get(srna));
                double logfc = (log2(da)-log2(this.srnas.get(srna)));
                if(logfc < -1 && logfc>-1000){
                    DE = true;
                    significantFC++;
                }
                else 
                    DE = false;
                writer.println(srna + " , " + srna.length() + " , " + this.srnas.get(srna) + " , " + (double)(this.srnas.get(srna))/this.srna_redundant_size*1000000.0 + " , " + da + " , " + (double)da/this.dicer_redundant_size*1000000.0 + " , " + fc + " , " + logfc + " , " + DE + " , " + mb + " , " + isomir + " , " + mirID + " , " + hasTarget + " , " + defPredicted + " , " + loosePredicted); // + " , " + (this.secondaryLfold.get(srna)) );
                intersectionCount++;
            }
        } */

            /*if(sRNAs_annotation.containsKey(srna)){
                mb = true;
                isomir = false;
                mirID = sRNAs_annotation.get(srna);
                miRBaseCount++;
            }else {
                mb = false;
                mirID = "Unknown";
                if(isomirs.containsKey(srna)){
                    isomir = true;
                    mirID = isomirs.get(srna);
                }
            }*/




/*
        //HashMap<String,String> sRNAs_annotation = new HashMap<>(); 
        //HashMap<String,String> known_mirnas = new HashMap<>();
        //HashMap<String,String> isomirs = new HashMap<>(); 
        int significantFC = 0;
        int intersectionCount = 0;
        int miRBaseCount = 0;
//        String mirID;   // if sRNA is known
        //boolean mb; // known/uknown in miRBase
        boolean DE;      // Differential expression is significant or not
        //boolean isomir = false; // if sRNA is an isomir
//        boolean hasTarget =  false; // sRNA has target
//        boolean defPredicted = false; // sRNA is predicted by miRCat2 default param
//        boolean loosePredicted = false; // sRNA is predicted by miRCat2 loose param
 */