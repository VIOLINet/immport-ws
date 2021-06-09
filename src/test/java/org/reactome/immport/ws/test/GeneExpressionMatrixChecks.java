package org.reactome.immport.ws.test;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Test;

/**
 * Some Java methods used to look at some gene expression files.
 * @author wug
 *
 */
public class GeneExpressionMatrixChecks {

    public GeneExpressionMatrixChecks() {
    }
    
    @Test
    public void checkMatrixFile() throws Exception {
        String dir = "/Volumes/ssd/docs/Immport_VO_project_Nasim";
        String allGeneMatrixFile = dir + File.separator + "all_expr_all_genes_df_final.csv";
        // Get all genes
        Set<String> allGenes = new HashSet<>();
        List<String> headers = new ArrayList<>();
        Files.lines(Paths.get(allGeneMatrixFile))
             .forEach(line -> {
                 if (headers.size() == 0) {
                     headers.addAll(Arrays.asList(line.split(",")));
                     return;
                 }
                 String[] tokens = line.split(",");
                 if (tokens[0].startsWith("\"LOC"))
                     return;
                 if (tokens[0].equals("\"CXCL8\"")) {
                     System.out.println(line);
                     System.out.println("Total samples: " + (tokens.length - 1));
                     for (int i = 0; i < tokens.length; i++) 
                         System.out.println(headers.get(i) + "\t" + tokens[i]);
                 }
                 allGenes.add(tokens[0]);
             });
        System.out.println("Total genes: " + allGenes.size());
//        allGenes.stream().sorted().forEach(System.out::println);
    }
    
}
