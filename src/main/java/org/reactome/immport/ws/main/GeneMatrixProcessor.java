package org.reactome.immport.ws.main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.junit.Test;

/**
 * This class is used to do some further pre-processes to the grouped gene expression matrix file.
 * @author wug
 *
 */
public class GeneMatrixProcessor {
    private final String GENE_NAME_FILE = "/Volumes/ssd/datasets/genenames/human_gene_names_051121.txt";
    private final String MATRIX_FILE = "/Volumes/ssd/docs/Immport_VO_project_Nasim/all_expr_all_genes_df_final.csv";
    
    public GeneMatrixProcessor() {
    }
    
    /**
     * For the final matrix, another filtering is conducted to keep genes only for ones approaved by HUGO.
     * @throws IOException
     */
    @Test
    public void filterMatrixToStandarGeneSymols() throws IOException {
        String dir = "/Volumes/ssd/docs/Immport_VO_project_Nasim/";
        String srcFile = dir + "all_expr_all_genes_df_final_mapped_merged_051221.csv";
        String targetFile = dir + "all_expr_all_genes_df_final_mapped_merged_approved_genes_051321.csv";
        
        Set<String> approvedGenes = loadApprovedGenes();
        System.out.println("Total approved genes: " + approvedGenes.size());
        // Just for comparison
        Map<String, Set<String>> geneToSyns = loadGeneToSynonyms(false);
        System.out.println("Total genes in the genetoSyns map: " + geneToSyns.size());
        
        BufferedReader br = new BufferedReader(new FileReader(srcFile));
        PrintWriter pw = new PrintWriter(targetFile);
        String line = br.readLine();
        pw.println(line); // The header
        Set<String> genesInMatrix = new HashSet<>();
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split(",");
            genesInMatrix.add(tokens[0]);
            if (approvedGenes.contains(tokens[0]))
                pw.println(line);
        }
        br.close();
        pw.close();
        
        System.out.println("Total genes in the merged matrix file: " + genesInMatrix.size());
        Set<String> shared = new HashSet<>(genesInMatrix);
        shared.retainAll(approvedGenes);
        System.out.println("Shared with approved genes: " + shared.size());
        shared = new HashSet<>(genesInMatrix);
        shared.retainAll(geneToSyns.keySet());
        System.out.println("Shared with the geneToSyns map: " + shared.size());
    }
    
    /**
     * Because of the different platforms used in the matrix file, some of genes should be 
     * remapped and may need to merge together, e.g.: IL8 = CXCL8 (the correct name).
     * @throws IOException
     */
    @Test
    public void mapNameToGeneSymolInMatrix() throws IOException {
        String dir = "/Volumes/ssd/docs/Immport_VO_project_Nasim/";
        // Intermediate file with names mapped to latest gene symols
        String mappedFile = dir + "all_expr_all_genes_df_final_mapped_051221.csv";
        // Rows are merged if they have the same gene symols using median or mean?
        String mergedFile = dir + "all_expr_all_genes_df_final_mapped_merged_051221.csv";
        
        Map<String, Set<String>> symToName = loadGeneToSynonyms(false);
        Map<String, Set<String>> nameToSyms = switchKeyToValue(symToName);
        
        FileReader fr = new FileReader(MATRIX_FILE);
        BufferedReader br = new BufferedReader(fr);
        PrintWriter pr = new PrintWriter(mappedFile);
        int totalLines = 0;
        String line = br.readLine(); // Head line
        line = line.replaceAll("\"", "");
        pr.println(line);
        int counterOfMoreThanOne = 0;
        List<String> mappedSyms = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            totalLines ++;
            // Get rid of quotation marks
            line = line.replaceAll("\"", "");
            String[] tokens = line.split(",");
            String name = tokens[0];
            Set<String> syms = nameToSyms.get(name);
            if (syms == null || syms.size() > 1) {
                mappedSyms.add(name); // Use the original name directly
                pr.println(line); // Do nothing in these cases
                if (syms != null && syms.size() > 1) {
                    counterOfMoreThanOne ++;
                    System.out.println(counterOfMoreThanOne + ": " + name +
                                       " mappped to more than one symbols: " + String.join(",", syms));
                }
                continue;
            }
            // Replace the gene
            String sym = syms.stream().findAny().get();
            mappedSyms.add(sym);
            line = line.replace(tokens[0], sym); // Keep the quotation marks for consistent
            pr.println(line);
        }
        br.close();
        fr.close();
        pr.close();
        System.out.println("Total parsed lines: " + totalLines);
        // A sanity check to see how many lines should be merged.
        System.out.println("Total mapped symbols: " + mappedSyms.size());
        Map<String, Integer> mappedToLine = new HashMap<>();
        mappedSyms.forEach(sym -> {
            mappedToLine.compute(sym, (key, counter) -> {
               if (counter == null)
                   return 1;
               counter ++;
               return counter;
            });
        });
        Set<String> symsWithMoreThanOneLine = mappedToLine.keySet()
                                                          .stream()
                                                          .filter(sym -> mappedToLine.get(sym) > 1)
                                                          .collect(Collectors.toSet());
        System.out.println("Symbols having more than one line: " + symsWithMoreThanOneLine.size());
        
        // For genes having more than one line, we need to merge them via median
        mergeLinesWithSameNames(mappedFile,
                                mergedFile, 
                                symsWithMoreThanOneLine);
    }
    
    private void mergeLinesWithSameNames(String srcFile,
                                         String targetFile,
                                         Set<String> genesWithMultipleLines) throws IOException {
        FileReader fr = new FileReader(srcFile);
        BufferedReader br = new BufferedReader(fr);
        PrintWriter pr = new PrintWriter(targetFile);
        String line = br.readLine();
        int totalLine = 0;
        // Just copy. Add a new token
        pr.println("Gene" + line);
        Map<String, List<String>> geneToLines = new HashMap<>();
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split(",");
            if (genesWithMultipleLines.contains(tokens[0])) {
                String line1 = line;
                // keep it for the time being
                geneToLines.compute(tokens[0], (key, list) -> {
                    if (list == null)
                        list = new ArrayList<>();
                    list.add(line1);
                    return list;
                });
            }
            else {
                pr.println(line); // Copy directly
                totalLine ++;
            }
        }
        System.out.println("Total genes having more than one line: " + geneToLines.size());
        // Do merge for these genes
        for (String gene : geneToLines.keySet()) {
            List<String> lines = geneToLines.get(gene);
            String merged = mergeLines(gene, lines);
            pr.println(merged);
            totalLine ++;
        }
        pr.close();
        br.close();
        fr.close();
        System.out.println("Total lines after merging: " + totalLine);
    }
    
    /**
     * Merge multiple lines into one via median.
     * @param lines
     * @return
     */
    private String mergeLines(String gene, 
                              List<String> lines) {
        if (lines.size() > 2) {
            System.out.println(gene + " lines: " + lines.size());
        }
        List<List<Double>> allValues = new ArrayList<>();
        // Pull out the first line
        String line = lines.get(0);
        String[] tokens = line.split(",");
        for (int i = 1; i < tokens.length; i++) {
            List<Double> col = new ArrayList<>();
            allValues.add(col);
            if (!tokens[i].equals("NA")) {
                col.add(new Double(tokens[i]));
            }
        }
        for (int j = 1; j < lines.size(); j++) {
            line = lines.get(j);
            tokens = line.split(",");
            for (int i = 1; i < tokens.length; i++) {
                List<Double> col = allValues.get(i - 1);
                if (!tokens[i].equals("NA")) {
                    col.add(new Double(tokens[i]));
                }
            }
        }
        StringBuilder builder = new StringBuilder();
        builder.append(gene);
        for (List<Double> col : allValues) {
            builder.append(",");
            Optional<Double> value = mergeValues(gene, col);
            if (value.isPresent())
                builder.append(value.get());
            else
                builder.append("NA");
        }
        return builder.toString();
    }
    
    private Optional<Double> mergeValues(String gene, List<Double> values) {
        if (values.size() == 0)
            return Optional.empty();
        List<Double> sorted = values.stream().sorted().collect(Collectors.toList());
        int index = sorted.size() / 2;
        if (sorted.size() % 2 == 0) // This is an even number list. Get the mean of two middle values
            return Optional.of((sorted.get(index - 1) + sorted.get(index)) / 2.0d);
        else 
            return Optional.of(sorted.get(index));
    }
    
    @Test
    public void testMergeValues() {
        List<Double> values = Stream.of(1.0, 5.0, 2.0).collect(Collectors.toList());
        System.out.println(mergeValues(null, values).get());
        values = Stream.of(1.0, 20.0, 2.0, 5.0).collect(Collectors.toList());
        System.out.println(mergeValues(null, values).get());
    }
    
    /**
     * This check proves that genes used in the matrix file are unique. Nothing duplicated.
     * @throws IOException
     */
    @Test
    public void checkGeneNamesInMatrix() throws IOException {
        String fileName = MATRIX_FILE;
        boolean hasQuotations = true;
        // Check the merged, gene-name normalized file
        String dir = "/Volumes/ssd/docs/Immport_VO_project_Nasim/";
        String mergedFile = dir + "all_expr_all_genes_df_final_mapped_merged_051221.csv";
        fileName = mergedFile;
        
        String filteredFile = dir + "all_expr_all_genes_df_final_mapped_merged_approved_genes_051321.csv";
        fileName = filteredFile;
        hasQuotations = false;
        
        FileReader fr = new FileReader(fileName);
        BufferedReader br = new BufferedReader(fr);
        int totalLines = 0;
        Set<String> allGenes = new HashSet<>();
        String line = br.readLine(); // Head line
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split(",");
            String gene = null;
            if (hasQuotations)
                gene = tokens[0].substring(1, tokens[0].length() - 1);
            else
                gene = tokens[0];
            allGenes.add(gene);
            totalLines ++;
        }
        br.close();
        fr.close();
        System.out.println("Total genes: " + allGenes.size());
        Optional<String> gene = allGenes.stream().findAny();
        System.out.println("One gene: " + gene.get());
        System.out.println("Total lines: " + totalLines);
        Set<String> synsWithMoreThanOneGene = getSynonymsMappedToMoreThanOneGene();
        System.out.println("Total names mapped to more than one symbol: " + synsWithMoreThanOneGene.size());
        synsWithMoreThanOneGene.retainAll(allGenes);
        System.out.println("Used in the matrix file: " + synsWithMoreThanOneGene.size());
//        synsWithMoreThanOneGene.stream().sorted().forEach(System.out::println);
        // For these names that can be mapped to more than one gene symbol, try to use 
        // previous symbols only
        Map<String, Set<String>> geneToPreSymols = loadGeneToSynonyms(true);
//        Map<String, Set<String>> preSymolToGenes = switchKeyToValue(geneToPreSymols);
//        synsWithMoreThanOneGene.forEach(syn -> {
//            Set<String> preSyms = preSymolToGenes.get(syn);
//            if (preSyms == null) {
//                System.out.println(syn + " has no pre-symbol.");
//                return;
//            }
//            if (preSyms.size() > 1)
//                System.out.println(syn + "\t" + String.join(",", preSyms));
//            else
//                System.out.println(syn + " is good!");
//        });
        // Do a map
        Map<String, Set<String>> synToSym = loadGeneToSynonyms(false);
        List<String> mappedGenes = new ArrayList<>();
        for (String name : allGenes) {
            Set<String> syms = synToSym.get(name);
            if (syms == null || syms.size() > 1) 
                mappedGenes.add(name);
            else 
                mappedGenes.add(syms.stream().findAny().get());
        }
        System.out.println("Total mapped genes: " + mappedGenes.size());
        Set<String> mappedGeneSet = mappedGenes.stream().collect(Collectors.toSet());
        System.out.println("Removing duplications: " + mappedGeneSet.size());
        
        // The following should print out from the above two lines for MATRIX_FILE:
//        Total mapped genes: 41464
//        Removing duplications: 40878
    }
    
    @Test
    public void checkNALines() throws IOException {
        String fileName = "/Volumes/ssd/docs/Immport_VO_project_Nasim/all_expr_all_genes_df_final_mapped_merged_051221.csv";
        long total = Files.lines(Paths.get(fileName))
                         .map(line -> line.split(","))
                         .map(tokens -> String.join(",", Stream.of(tokens).skip(1).collect(Collectors.joining(","))))
                         .filter(line -> !line.contains("NA"))
                         .collect(Collectors.counting());
        System.out.println("Total lines without NA: " + total);
    }
    
    @Test
    public void checkGeneNames() throws IOException {
        Map<String, Set<String>> geneToSynonyms = loadGeneToSynonyms(false);
        System.out.println("Total collected genes: " + geneToSynonyms.size());
        // Check if a synonyms may be mapped to more than one gene
        Map<String, Set<String>> synToGene = switchKeyToValue(geneToSynonyms);
        System.out.println("Total synonyms: " + synToGene.size());
        synToGene.forEach((syn, genes) -> {
            if (genes.size() > 1)
                System.out.println(syn + "\t" + String.join(",", genes));
        });
    }
    
    private Set<String> getSynonymsMappedToMoreThanOneGene() throws IOException {
        Map<String, Set<String>> geneToSynonyms = loadGeneToSynonyms(false);
        // Check if a synonyms may be mapped to more than one gene
        Map<String, Set<String>> synToGene = switchKeyToValue(geneToSynonyms);
        return synToGene.keySet()
                        .stream()
                        .filter(syn -> synToGene.get(syn).size() > 1)
                        .collect(Collectors.toSet());
    }
    
    private Map<String, Set<String>> loadGeneToSynonyms(boolean needPreviousOnly) throws IOException {
        Map<String, Set<String>> geneToSynonyms = new HashMap<>();
        Files.lines(Paths.get(GENE_NAME_FILE))
             .skip(1)
             .map(line -> line.split("\t"))
             .filter(tokens -> tokens[3].equals("Approved")) // Not withdrawn
             .forEach(tokens -> {
                 if (tokens.length > 4)
                     pushNames(tokens[1], tokens[4], geneToSynonyms); // Previous symbols
                 if (!needPreviousOnly && tokens.length > 5)
                     pushNames(tokens[1], tokens[5], geneToSynonyms); // Alias symbols
             });
        return geneToSynonyms;
    }
    
    private Set<String> loadApprovedGenes() throws IOException {
        Set<String> genes = Files.lines(Paths.get(GENE_NAME_FILE))
                                 .skip(1)
                                 .map(line -> line.split("\t"))
                                 .filter(tokens -> tokens[3].equals("Approved"))
                                 .map(tokens -> tokens[1])
                                 .collect(Collectors.toSet());
        return genes;
    }
    
    private void pushNames(String gene,
                           String token, 
                           Map<String, Set<String>> geneToSynonyms) {
        if (token.length() == 0)
            return;
        String[] tokens = token.split(", ");
        for (String name : tokens) {
            geneToSynonyms.compute(gene, (key, set) -> {
                if (set == null)
                    set = new HashSet<>();
                set.add(name);
                return set;
            });
        }
    }
    
    private Map<String, Set<String>> switchKeyToValue(Map<String, Set<String>> keyToSet) {
        Map<String, Set<String>> rtn = new HashMap<>();
        keyToSet.forEach((key, set) -> {
            set.forEach(s -> {
                rtn.compute(s, (newKey, newSet) -> {
                    if (newSet == null)
                        newSet = new HashSet<>();
                    newSet.add(key);
                    return newSet;
                });
            });
        });
        return rtn;
    }
    
}
