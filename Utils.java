/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CalculateFoldChange;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

/**
 *
 * @author ujy06jau
 */
public class Utils {

    public static final double NO_SCORE = -1000;

    public static void copyDir(File srcDir, File tgtDir) throws IOException {
        if (!srcDir.exists()) {
            return;
        }
        Files.copy(srcDir.toPath(), tgtDir.toPath(), StandardCopyOption.REPLACE_EXISTING);
        while (!tgtDir.exists());
        for (File f : srcDir.listFiles()) {
            if (f.isDirectory()) {
                File newDir = new File(tgtDir + "/" + f.getName());
                copyDir(f, newDir);
                while (!newDir.exists());
            } else {
                File newFile = new File(tgtDir + "/" + f.getName());
                Files.copy(f.toPath(), newFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                while (!newFile.exists());
            }
        }
    }

    public static void deleteDir(File dir) {

        if (!dir.exists()) {
            return;
        }
        for (File f : dir.listFiles()) {
            if (f.isDirectory()) {
                deleteDir(f);
            } else {
                f.delete();
                while (f.exists());
                System.out.println("deleted: " + f.getAbsolutePath());
            }
        }
        dir.delete();
        while (dir.exists());
        System.out.println("deleted: " + dir.getAbsolutePath());

    }

    public static int runExe(String[] arguments, File input, File output, File error, File workingDirectory) throws IOException, InterruptedException {

        ProcessBuilder pb = new ProcessBuilder(arguments);
        if (workingDirectory != null) {
            pb.directory(workingDirectory);
        }
        if (input == null) {
            pb.redirectInput(ProcessBuilder.Redirect.INHERIT);
        } else {
            pb.redirectInput(input);
        }
        if (output == null) {
            pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
        } else {
            pb.redirectOutput(output);
        }
        if (error == null) {
            pb.redirectError(ProcessBuilder.Redirect.INHERIT);
        } else {
            pb.redirectError(error);
        }
        Process process = pb.start();
        process.waitFor();
        return process.exitValue();
    }

    // write message string to end of file
    public static void write2File(File file, String str) throws IOException {
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(file, true)));
        out.println(str);
        out.close();
    }

    public static double doubleRandomInclusive(double max, double min) {
        double r = Math.random();
        if (r < 0.5) {
            return ((1 - Math.random()) * (max - min) + min);
        }
        return (Math.random() * (max - min) + min);
    }
}
