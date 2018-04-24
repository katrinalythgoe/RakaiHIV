import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.RootedTree;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.lang3.ArrayUtils;
import sun.security.krb5.internal.crypto.Des;
import sun.security.ssl.SSLContextImpl;

import java.io.*;
import java.util.*;

/**
 * Created by jayna on 08/06/15.
 */
public class robustingCountingAnalysis {


    String treeFile;
    String outputFile;
    String name;
    DescriptiveStatistics treeLength_stat;
    List<RootedTree> trees;


    public robustingCountingAnalysis(String treefile, String name) {


        this.treeFile = treefile;
        this.outputFile = treefile.replace(".tre", "_output_N_S_over_time.csv");
        this.name = name;
        this.trees = readTrees(this.treeFile);
    }

    public void extract_N_and_S_overTime_raw() {

        try {
            NexusImporter importer = new NexusImporter(new FileReader(treeFile));

            //BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile));
            int tree_no = 1;
            while (importer.hasTree()) {

                RootedTree tree = (RootedTree) importer.importNextTree();

                Set<Node> nodes = tree.getNodes();
                nodes.remove(tree.getRootNode());

                for(Node n: nodes) {

                    double height = tree.getHeight(n);
                    double N_count = (Double)n.getAttribute("N");
                    double S_count = (Double)n.getAttribute("S");
                    String external_or_internal;
                    if(tree.isExternal(n)) {

                        external_or_internal = "external";

                    }
                    else {
                        external_or_internal = "internal";
                    }

                    System.out.println("tree"+tree_no+","+height+","+N_count+","+S_count+","+external_or_internal);
                }


            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void extract_N_and_S_overTime_per_bin(double mrsd, double binLength, String nodeType, String which_lineage) {


        File input = new File(treeFile);
        File parent_input = input.getParentFile();
        System.out.println(parent_input.toString());

        double [] binU;
        double [] binL;


        double maxRootHeight = 0;
        for(RootedTree t: trees) {

            if(t.getHeight(t.getRootNode()) > maxRootHeight) {

                maxRootHeight = t.getHeight(t.getRootNode());
            }

        }


        binL = new double[(int)Math.ceil(maxRootHeight/binLength)];
        binU = new double[(int)Math.ceil(maxRootHeight/binLength)];


        for(int b=0; b < binL.length; b++) {

            binL[b] = maxRootHeight-(b*binLength);
            binU[b] = maxRootHeight-((b+1)*binLength);

            System.out.println("bin"+(b+1)+": "+binL[b]+", "+ binU[b]);


        }



        double [][] N_counts_per_lineage_over_time = new double[binL.length][trees.size()];
        double [][] S_counts_per_lineage_over_time = new double[binL.length][trees.size()];
        double [][] cumN_counts_over_time = new double[binL.length][trees.size()];
        double [][] cumS_counts_over_time = new double[binL.length][trees.size()];
        double [][] N_rate_over_time = new double[binL.length][trees.size()];
        double [][] S_rate_over_time = new double[binL.length][trees.size()];
        double [][] cumN_rate_over_time = new double[binL.length][trees.size()];
        double [][] cumS_rate_over_time = new double[binL.length][trees.size()];
        double [][] N_rate_per_lineage_over_time = new double[binL.length][trees.size()];
        double [][] S_rate_per_lineage_over_time = new double[binL.length][trees.size()];
        double [][] cumN_rate_per_lineage_over_time = new double[binL.length][trees.size()];
        double [][] cumS_rate_per_lineage_over_time = new double[binL.length][trees.size()];

        double [][] uN_rate_over_time = new double[binL.length][trees.size()];
        double [][] uS_rate_over_time = new double[binL.length][trees.size()];
        double [][] uN_counts_over_time = new double[binL.length][trees.size()];
        double [][] uS_counts_over_time = new double[binL.length][trees.size()];
        double [][] uN_cum_rate_over_time = new double[binL.length][trees.size()];
        double [][] uS_cum_rate_over_time = new double[binL.length][trees.size()];
        double [][] uN_cum_counts_over_time = new double[binL.length][trees.size()];
        double [][] uS_cum_counts_over_time = new double[binL.length][trees.size()];
        double [][] dN_over_time = new double[binL.length][trees.size()];
        double [][] dS_over_time = new double[binL.length][trees.size()];
        double [][] dN_dS = new double[binL.length][trees.size()];
        double [][] cum_dN = new double[binL.length][trees.size()];
        double [][] cum_dS = new double[binL.length][trees.size()];

        //treeLength_stat = new DescriptiveStatistics();

        try {
            BufferedWriter writer_stats = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+name+"_"+nodeType+"_"+which_lineage+"_dN_dS_stats.txt"));
            writer_stats.write("State\tdN\tdS\tdn/dS\tN\tS\ttreeLength\tN_rate\tS_rate\n");


            int t_i = 0;

            treeLength_stat = new DescriptiveStatistics();




            for(RootedTree t: trees) {

                System.out.println(t_i);

                //just to get treelength right
                Set<Node> nodes = null;
                if(nodeType == "internal") {
                    nodes = t.getInternalNodes();

                    Set<Node> tempNodes = new HashSet<Node>();


                    if(which_lineage != null) {
                        for (Node n : nodes) {


                            if ((n.getAttributeMap().containsKey("lineage"))) {

                                String lineage = (String) n.getAttribute("lineage");


                                //Double trunk_weight = (Double)n.getAttribute("trunk_weight");

                                if (lineage.equals(which_lineage)) {
                                    tempNodes.add(n);
                                }
                            }

                        }
                        nodes.retainAll(tempNodes);
                    }
                    else{

                        for(Node n: nodes) {

                            if (n.getAttributeMap().containsKey("trunk_weight")){

                                Integer trunk_weight = (Integer) n.getAttribute("trunk_weight");

                                if(trunk_weight >= 10) {

                                    tempNodes.add(n);
                                }
                            }
                        }

                        if(tempNodes.size()>0) {
                            nodes.retainAll(tempNodes);
                        }
                    }

                }
                else if(nodeType == "external") {

                    nodes = t.getExternalNodes();

                }
                else{

                    nodes = t.getNodes();
                }

                double treeLength = 0;
                for(Node n: nodes) {

                    treeLength += t.getLength(n);
                }



                double[][] countsPerBin = new double[binL.length][];


                for(int b=0; b < binL.length; b++) {


                    countsPerBin[b] = getCountsPerBin(binL[b], binU[b], t, nodeType, which_lineage);
                    double dN;
                    double dS;
                    if(b>0) {

                        countsPerBin[b][10] = countsPerBin[b][0]+countsPerBin[b-1][10];
                        countsPerBin[b][11] = countsPerBin[b][1]+countsPerBin[b-1][11];
                        countsPerBin[b][12] = countsPerBin[b][2]+countsPerBin[b-1][12];
                        countsPerBin[b][13] = countsPerBin[b][3]+countsPerBin[b-1][13];
                        countsPerBin[b][14] = countsPerBin[b][4]+countsPerBin[b-1][14];
                        countsPerBin[b][15] = countsPerBin[b][5]+countsPerBin[b-1][15];
                        countsPerBin[b][16] = countsPerBin[b][6]+countsPerBin[b-1][16];
                        countsPerBin[b][17] = countsPerBin[b][7]+countsPerBin[b-1][17];
                        countsPerBin[b][18] = countsPerBin[b][8]+countsPerBin[b-1][18];
                        countsPerBin[b][19] = countsPerBin[b][9]+countsPerBin[b-1][19];

                        dN = countsPerBin[b][2]/countsPerBin[b][8];
                        dS = countsPerBin[b][3]/countsPerBin[b][9];

                        if(Double.isInfinite(dS)) {

                            dS = Double.NaN;
                        }
                        if(Double.isInfinite(dN)) {

                            dN = Double.NaN;
                        }

                        dN_dS[b][t_i] = dN/dS;


                        if(Double.isNaN(dN)) {

                            cum_dN[b][t_i] = cum_dN[b - 1][t_i];

                        }
                        else{

                            cum_dN[b][t_i] = dN + cum_dN[b - 1][t_i];

                        }

                        if(Double.isNaN(dS)) {

                            cum_dS[b][t_i] = cum_dS[b - 1][t_i];

                        }
                        else{

                            cum_dS[b][t_i] = dS + cum_dS[b - 1][t_i];

                        }

                        if(b == binL.length-1) {

                            double N = countsPerBin[b][12];
                            double S = countsPerBin[b][13];
                            //double treeLength = treeLength_stat.getValues()[t_i];
                            //System.out.println(treeLength_stat.getValues()[1]+","+treeLength_stat.getValues()[9]);

                            writer_stats.write(t_i+"\t"+cum_dN[b][t_i]+"\t"+cum_dS[b][t_i]+"\t"+(cum_dN[b][t_i]/cum_dS[b][t_i])+"\t"+N+"\t"
                                    +S+"\t"+treeLength+"\t"+(N/1212/treeLength)+"\t"+(S/1212/treeLength)+"\n");
                            writer_stats.flush();


                        }

                    }
                    else{
                        countsPerBin[b][10] = countsPerBin[b][0]; //cN_rate
                        countsPerBin[b][11] = countsPerBin[b][1]; //cS_rate
                        countsPerBin[b][12] = countsPerBin[b][2]; //cN_count
                        countsPerBin[b][13] = countsPerBin[b][3]; //cS_count
                        countsPerBin[b][14] = countsPerBin[b][4];
                        countsPerBin[b][15] = countsPerBin[b][5];
                        countsPerBin[b][16] = countsPerBin[b][6]; //uN_rate
                        countsPerBin[b][17] = countsPerBin[b][7]; //uS_rate
                        countsPerBin[b][18] = countsPerBin[b][8]; //uN_count
                        countsPerBin[b][19] = countsPerBin[b][9]; //uS_count

                        dN = countsPerBin[b][2]/countsPerBin[b][8];
                        dS = countsPerBin[b][3]/countsPerBin[b][9];

                        if(Double.isInfinite(dS)) {

                            dS = Double.NaN;
                        }
                        if(Double.isInfinite(dN)) {

                            dN = Double.NaN;
                        }


                        dN_dS[b][t_i] = dN/dS;

                        if(Double.isNaN(dN)) {
                            dN = 0.0;
                        }
                        if(Double.isNaN(dS)) {
                            dS = 0.0;
                        }

                        cum_dN[b][t_i] = dN;
                        cum_dS[b][t_i] = dS;

                    }


                    N_rate_over_time[b][t_i] = countsPerBin[b][0];
                    S_rate_over_time[b][t_i] = countsPerBin[b][1];
                    N_counts_per_lineage_over_time[b][t_i] = countsPerBin[b][2];
                    S_counts_per_lineage_over_time[b][t_i] = countsPerBin[b][3];
                    N_rate_per_lineage_over_time[b][t_i] = countsPerBin[b][4];
                    S_rate_per_lineage_over_time[b][t_i] = countsPerBin[b][5];
                    uN_rate_over_time[b][t_i] = countsPerBin[b][6];
                    uS_rate_over_time[b][t_i] = countsPerBin[b][7];
                    uN_counts_over_time[b][t_i] = countsPerBin[b][8];
                    uS_counts_over_time[b][t_i] = countsPerBin[b][9];

                    if(countsPerBin[b][8] > 0 && countsPerBin[b][9] > 0) {
                        dN_over_time[b][t_i] = countsPerBin[b][2] / countsPerBin[b][8];
                        dS_over_time[b][t_i] = countsPerBin[b][3] / countsPerBin[b][9];
                    }
                    else{
                        dN_over_time[b][t_i] = Double.NaN;
                        dS_over_time[b][t_i] = Double.NaN;
                    }

                    cumN_rate_over_time[b][t_i] = countsPerBin[b][10];
                    cumS_rate_over_time[b][t_i] = countsPerBin[b][11];
                    cumN_counts_over_time[b][t_i] = countsPerBin[b][12];
                    cumS_counts_over_time[b][t_i] = countsPerBin[b][13];
                    cumN_rate_per_lineage_over_time[b][t_i] = countsPerBin[b][14];
                    cumS_rate_per_lineage_over_time[b][t_i] = countsPerBin[b][15];
                    uN_cum_rate_over_time[b][t_i] = countsPerBin[b][16];
                    uS_cum_rate_over_time[b][t_i] = countsPerBin[b][17];
                    uN_cum_counts_over_time[b][t_i] = countsPerBin[b][18];
                    uS_cum_counts_over_time[b][t_i] = countsPerBin[b][19];



                }
                t_i++;

            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        Double median;
        Double mean;
        Double LPD;
        Double UPD;
        Double std;

        try {
            BufferedWriter writer1 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_N_S_counts_output.csv"));
            BufferedWriter writer2 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_N_S_rate_output.csv"));
            BufferedWriter writer3 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_N_S_cum_rate_output.csv"));
            BufferedWriter writer4 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_N_S_cum_counts_output.csv"));
            BufferedWriter writer5 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_N_S_rate_pl_output.csv"));
            BufferedWriter writer6 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_N_S_cum_rate_pl_output.csv"));
            BufferedWriter writer7 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_uN_uS_rate_output.csv"));
            BufferedWriter writer8 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_dN_dS_over_time_output.csv"));
            BufferedWriter writer9 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_cum_uN_uS_rate_output.csv"));
            BufferedWriter writer10 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_cum_uN_uS_counts_output.csv"));
            BufferedWriter writer11 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_uN_uS_counts_output.csv"));
            BufferedWriter writer12 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_dN_dS_ratio_counts_output.csv"));
            BufferedWriter writer13 = new BufferedWriter(new FileWriter(parent_input.toString()+"/"+nodeType+"_robustcounting_cum_dN_dS_over_over_time_output.csv"));


            writer1.write("Time,Mean_N,Median_N,LPD_N,HPD_N,std_N,Mean_S,Median_S,LPD_S,HPD_S,std_S\n");
            writer2.write("Time,Mean_N,Median_N,LPD_N,HPD_N,std_N,Mean_S,Median_S,LPD_S,HPD_S,std_S\n");
            writer3.write("Time,Mean_N,Median_N,LPD_N,HPD_N,std_N,Mean_S,Median_S,LPD_S,HPD_S,std_S\n");
            writer4.write("Time,Mean_N,Median_N,LPD_N,HPD_N,std_N,Mean_S,Median_S,LPD_S,HPD_S,std_S\n");
            writer5.write("Time,Mean_N,Median_N,LPD_N,HPD_N,std_N,Mean_S,Median_S,LPD_S,HPD_S,std_S\n");
            writer6.write("Time,Mean_N,Median_N,LPD_N,HPD_N,std_N,Mean_S,Median_S,LPD_S,HPD_S,std_S\n");
            writer7.write("Time,Mean_uN,Median_uN,LPD_uN,HPD_uN,std_uN,Mean_uS,Median_uS,LPD_uS,HPD_uS,std_uS\n");
            writer8.write("Time,Mean_dN,Median_dN,LPD_dN,HPD_dN,std_dN,Mean_dS,Median_dS,LPD_dS,HPD_dS,std_dS\n");
            writer9.write("Time,Mean_uN,Median_uN,LPD_uN,HPD_uN,std_uN,Mean_uS,Median_uS,LPD_uS,HPD_uS,std_uS\n");
            writer10.write("Time,Mean_uN,Median_uN,LPD_uN,HPD_uN,std_uN,Mean_uS,Median_uS,LPD_uS,HPD_uS,std_uS\n");
            writer11.write("Time,Mean_uN,Median_uN,LPD_uN,HPD_uN,std_uN,Mean_uS,Median_uS,LPD_uS,HPD_uS,std_uS\n");
            writer12.write("Time,Mean_dN_dS,Median_dN_dS,LPD_dN_dS,HPD_dN_dS,std_dN_dS\n");
            writer13.write("Time,Mean_cum_dN,Median_cum_dN,LPD_cum_dN,HPD_cum_dN,std_cum_dN,Mean_cum_dS,Median_cum_dS,LPD_cum_dS,HPD_cum_dS,std_cum_dS\n");



            for(int b=0; b < binL.length; b++) {

                double midpoint = mrsd - 0.5*(binL[b]+binU[b]);

                DescriptiveStatistics N_stats = new DescriptiveStatistics(processArray(N_counts_per_lineage_over_time[b]));
                DescriptiveStatistics S_stats = new DescriptiveStatistics(processArray(S_counts_per_lineage_over_time[b]));
                DescriptiveStatistics cumN_stats = new DescriptiveStatistics(processArray(cumN_counts_over_time[b]));
                DescriptiveStatistics cumS_stats = new DescriptiveStatistics(processArray(cumS_counts_over_time[b]));
                DescriptiveStatistics N_r_stats = new DescriptiveStatistics(processArray(N_rate_over_time[b]));
                DescriptiveStatistics S_r_stats = new DescriptiveStatistics(processArray(S_rate_over_time[b]));
                DescriptiveStatistics cumN_r_stats = new DescriptiveStatistics(processArray(cumN_rate_over_time[b]));
                DescriptiveStatistics cumS_r_stats = new DescriptiveStatistics(processArray(cumS_rate_over_time[b]));
                DescriptiveStatistics N_r_pl_stats = new DescriptiveStatistics(processArray(N_rate_per_lineage_over_time[b]));
                DescriptiveStatistics S_r_pl_stats = new DescriptiveStatistics(processArray(S_rate_per_lineage_over_time[b]));
                DescriptiveStatistics cumN_r_pl_stats = new DescriptiveStatistics(processArray(cumN_rate_per_lineage_over_time[b]));
                DescriptiveStatistics cumS_r_pl_stats = new DescriptiveStatistics(processArray(cumS_rate_per_lineage_over_time[b]));

                DescriptiveStatistics uN_r_stats = new DescriptiveStatistics(processArray(uN_rate_over_time[b]));
                DescriptiveStatistics uS_r_stats = new DescriptiveStatistics(processArray(uS_rate_over_time[b]));
                DescriptiveStatistics dN_r_stats = new DescriptiveStatistics(processArray(dN_over_time[b]));
                DescriptiveStatistics dS_r_stats = new DescriptiveStatistics(processArray(dS_over_time[b]));
                DescriptiveStatistics uN_counts = new DescriptiveStatistics(processArray(uN_counts_over_time[b]));
                DescriptiveStatistics uS_counts = new DescriptiveStatistics(processArray(uS_counts_over_time[b]));
                DescriptiveStatistics uN_cum_r_stats = new DescriptiveStatistics(processArray(uN_cum_rate_over_time[b]));
                DescriptiveStatistics uS_cum_r_stats = new DescriptiveStatistics(processArray(uS_cum_rate_over_time[b]));
                DescriptiveStatistics uN_cum_counts = new DescriptiveStatistics(processArray(uN_cum_counts_over_time[b]));
                DescriptiveStatistics uS_cum_counts = new DescriptiveStatistics(processArray(uS_cum_counts_over_time[b]));
                DescriptiveStatistics dN_dS_ratio = new DescriptiveStatistics(processArray(dN_dS[b]));
                DescriptiveStatistics cum_dN_stats = new DescriptiveStatistics(processArray(cum_dN[b]));
                DescriptiveStatistics cum_dS_stats = new DescriptiveStatistics(processArray(cum_dS[b]));



                //counts


                median = N_stats.getPercentile(50);
                mean = N_stats.getMean();
                LPD = N_stats.getPercentile(5);
                UPD = N_stats.getPercentile(95);
                std = N_stats.getStandardDeviation();

                String N_S_counts = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;
                //writer1.write(midpoint+","+mean+","+median+","+LQ+","+UQ+"\n");

                median = S_stats.getPercentile(50);
                mean = S_stats.getMean();
                LPD = S_stats.getPercentile(5);
                UPD = S_stats.getPercentile(95);
                std = S_stats.getStandardDeviation();


                N_S_counts += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";
                writer1.write(N_S_counts);
                writer1.flush();

                // rate

                median = N_r_stats.getPercentile(50);
                mean = N_r_stats.getMean();
                LPD = N_r_stats.getPercentile(5);
                UPD = N_r_stats.getPercentile(95);
                std = N_r_stats.getStandardDeviation();


                String N_S_rate = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = S_r_stats.getPercentile(50);
                mean = S_r_stats.getMean();
                LPD = S_r_stats.getPercentile(5);
                UPD = S_r_stats.getPercentile(95);
                std = S_r_stats.getStandardDeviation();

                N_S_rate += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer2.write(N_S_rate);
                writer2.flush();

                //cum rate
                median = cumN_r_stats.getPercentile(50);
                mean = cumN_r_stats.getMean();
                LPD = cumN_r_stats.getPercentile(5);
                UPD = cumN_r_stats.getPercentile(95);
                std = cumN_r_stats.getStandardDeviation();

                String N_S_r_cum = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = cumS_r_stats.getPercentile(50);
                mean = cumS_r_stats.getMean();
                LPD = cumS_r_stats.getPercentile(5);
                UPD = cumS_r_stats.getPercentile(95);
                std = cumS_r_stats.getStandardDeviation();

                N_S_r_cum += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer3.write(N_S_r_cum);
                writer3.flush();

                //cum counts
                median = cumN_stats.getPercentile(50);
                mean = cumN_stats.getMean();
                LPD = cumN_stats.getPercentile(5);
                UPD = cumN_stats.getPercentile(95);
                std = cumN_stats.getStandardDeviation();

                String N_S_cum = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = cumS_stats.getPercentile(50);
                mean = cumS_stats.getMean();
                LPD = cumS_stats.getPercentile(5);
                UPD = cumS_stats.getPercentile(95);
                std = cumS_stats.getStandardDeviation();

                N_S_cum += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer4.write(N_S_cum);
                writer4.flush();


                // rate  per lineage

                median = N_r_pl_stats.getPercentile(50);
                mean = N_r_pl_stats.getMean();
                LPD = N_r_pl_stats.getPercentile(5);
                UPD = N_r_pl_stats.getPercentile(95);
                std = N_r_pl_stats.getStandardDeviation();

                String N_S_rate_pl = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = S_r_pl_stats.getPercentile(50);
                mean = S_r_pl_stats.getMean();
                LPD = S_r_pl_stats.getPercentile(5);
                UPD = S_r_pl_stats.getPercentile(95);
                std = S_r_pl_stats.getStandardDeviation();


                N_S_rate_pl += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer5.write(N_S_rate_pl);
                writer5.flush();

                //cum rate per lineage
                median = cumN_r_pl_stats.getPercentile(50);
                mean = cumN_r_pl_stats.getMean();
                LPD = cumN_r_pl_stats.getPercentile(5);
                UPD = cumN_r_pl_stats.getPercentile(95);
                std = cumN_r_pl_stats.getStandardDeviation();

                String N_S_r_cum_pl = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = cumS_r_pl_stats.getPercentile(50);
                mean = cumS_r_pl_stats.getMean();
                LPD = cumS_r_pl_stats.getPercentile(5);
                UPD = cumS_r_pl_stats.getPercentile(95);
                std = cumS_r_pl_stats.getStandardDeviation();

                N_S_r_cum_pl += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer6.write(N_S_r_cum_pl);
                writer6.flush();


                //uN and uS rate

                median = uN_r_stats.getPercentile(50);
                mean = uN_r_stats.getMean();
                LPD = uN_r_stats.getPercentile(5);
                UPD = uN_r_stats.getPercentile(95);
                std = uN_r_stats.getStandardDeviation();

                String uN_uS_rate = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = uS_r_stats.getPercentile(50);
                mean = uS_r_stats.getMean();
                LPD = uS_r_stats.getPercentile(5);
                UPD = uS_r_stats.getPercentile(95);
                std = uS_r_stats.getStandardDeviation();

                uN_uS_rate += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer7.write(uN_uS_rate);
                writer7.flush();

                // uN and uS cum rate

                median = uN_cum_r_stats.getPercentile(50);
                mean = uN_cum_r_stats.getMean();
                LPD = uN_cum_r_stats.getPercentile(5);
                UPD = uN_cum_r_stats.getPercentile(95);
                std = uN_cum_r_stats.getStandardDeviation();

                String uN_uS_cum_rate = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = uS_cum_r_stats.getPercentile(50);
                mean = uS_cum_r_stats.getMean();
                LPD = uS_cum_r_stats.getPercentile(5);
                UPD = uS_cum_r_stats.getPercentile(95);
                std = uS_cum_r_stats.getStandardDeviation();

                uN_uS_cum_rate += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer9.write(uN_uS_cum_rate);
                writer9.flush();

                // uN and uS counts
                median = uN_counts.getPercentile(50);
                mean = uN_counts.getMean();
                LPD = uN_counts.getPercentile(5);
                UPD = uN_counts.getPercentile(95);
                std = uN_counts.getStandardDeviation();

                String uN_uS_counts = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = uS_counts.getPercentile(50);
                mean = uS_counts.getMean();
                LPD = uS_counts.getPercentile(5);
                UPD = uS_counts.getPercentile(95);
                std = uS_counts.getStandardDeviation();

                uN_uS_counts += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer11.write(uN_uS_counts);
                writer11.flush();

                // uN and uS cum counts
                median = uN_cum_counts.getPercentile(50);
                mean = uN_cum_counts.getMean();
                LPD = uN_cum_counts.getPercentile(5);
                UPD = uN_cum_counts.getPercentile(95);
                std = uN_cum_counts.getStandardDeviation();

                String uN_uS_cum_counts = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = uS_cum_counts.getPercentile(50);
                mean = uS_cum_counts.getMean();
                LPD = uS_cum_counts.getPercentile(5);
                UPD = uS_cum_counts.getPercentile(95);
                std = uS_cum_counts.getStandardDeviation();

                uN_uS_cum_counts += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer10.write(uN_uS_cum_counts);
                writer10.flush();


                //dN and dS
                median = dN_r_stats.getPercentile(50);
                mean = dN_r_stats.getMean();
                LPD = dN_r_stats.getPercentile(5);
                UPD = dN_r_stats.getPercentile(95);
                std = dN_r_stats.getStandardDeviation();

                String dN_dS_rate = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = dS_r_stats.getPercentile(50);
                mean = dS_r_stats.getMean();
                LPD = dS_r_stats.getPercentile(5);
                UPD = dS_r_stats.getPercentile(95);
                std = dS_r_stats.getStandardDeviation();

                dN_dS_rate += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer8.write(dN_dS_rate);
                writer8.flush();


                //dN/dS ratio

                median = dN_dS_ratio.getPercentile(50);
                mean = dN_dS_ratio.getMean();
                LPD = dN_dS_ratio.getPercentile(5);
                UPD = dN_dS_ratio.getPercentile(95);
                std = dN_dS_ratio.getStandardDeviation();

                String dN_dS_stats = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer12.write(dN_dS_stats);
                writer12.flush();


                median = cum_dN_stats.getPercentile(50);
                mean = cum_dN_stats.getMean();
                LPD = cum_dN_stats.getPercentile(5);
                UPD = cum_dN_stats.getPercentile(95);
                std = cum_dN_stats.getStandardDeviation();

                String cum_dN_dS_stats = midpoint+","+mean+","+median+","+LPD+","+UPD+","+std;

                median = cum_dS_stats.getPercentile(50);
                mean = cum_dS_stats.getMean();
                LPD = cum_dS_stats.getPercentile(5);
                UPD = cum_dS_stats.getPercentile(95);
                std = cum_dS_stats.getStandardDeviation();

                cum_dN_dS_stats += ","+mean+","+median+","+LPD+","+UPD+","+std+"\n";

                writer13.write(cum_dN_dS_stats);
                writer13.flush();



            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private List<RootedTree> readTrees(String file) {

        List<RootedTree> trees = new ArrayList<RootedTree>();
        int t_i = 0;
        try {
            NexusImporter importer = new NexusImporter(new FileReader(file));

            while (importer.hasTree()) {


                System.out.println(t_i);
                t_i++;
//                if(t_i == 500) {
//                    break;
//                }
                RootedTree tree = (RootedTree) importer.importNextTree();
                trees.add(tree);
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (ImportException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return trees;

    }

    private double [] processArray(double [] array) {

        Double [] a = ArrayUtils.toObject(array);
        List<Double> tempList = new ArrayList<Double>(Arrays.asList(a));

        while(tempList.contains(Double.POSITIVE_INFINITY)) {

            tempList.remove(Double.POSITIVE_INFINITY);
        }

        while(tempList.contains(Double.NaN)) {

            tempList.remove(Double.NaN);
        }



        double [] tempArray = new double[tempList.size()];
        for(int i=0; i < tempList.size(); i++) {

            tempArray[i] = (Double)tempList.get(i);
        }

        return tempArray;

    }

    private double[] getCountsPerBin(double binL, double binU, RootedTree tree, String nodeType, String which_lineage) {

        double [] counts = new double[20];
        Arrays.fill(counts, 0.0);


        Set<Node> nodes = null;
        if(nodeType == "internal") {
            nodes = tree.getInternalNodes();

            Set<Node> tempNodes = new HashSet<Node>();


            if(which_lineage != null) {
                for (Node n : nodes) {


                    if ((n.getAttributeMap().containsKey("lineage"))) {

                        String lineage = (String) n.getAttribute("lineage");


                        //Double trunk_weight = (Double)n.getAttribute("trunk_weight");

                        if (lineage.equals(which_lineage)) {
                            tempNodes.add(n);
                        }
                    }

                }
                nodes.retainAll(tempNodes);
            }
            else{

                for(Node n: nodes) {

                    if (n.getAttributeMap().containsKey("trunk_weight")){

                        Integer trunk_weight = (Integer) n.getAttribute("trunk_weight");

                        if(trunk_weight >= 10) {

                            tempNodes.add(n);
                        }
                    }
                }

                if(tempNodes.size()>0) {
                    nodes.retainAll(tempNodes);
                }
            }

        }
        else if(nodeType == "external") {

            nodes = tree.getExternalNodes();

        }
        else{

            nodes = tree.getNodes();
        }

        double binLength = binL - binU;
        int lineages = 1;
        //System.out.println(tree.getAttributeNames());
        int n_i = 1;
        double total_time = 0.0; // total branch length during the interval
        double treeLength = 0.0;

        //treeLength_stat = new DescriptiveStatistics();
        for(Node n: nodes) {

            //treeLength+= tree.getLength(n);


            if (!tree.isRoot(n) && tree.getParent(n)!=null) {
                double parentHeight = tree.getHeight(tree.getParent(n));
                double height = tree.getHeight(n);

                double N_count = 0.0;
                double S_count = 0.0;
                double uN_count = 0.0;
                double uS_count = 0.0;


                if(n.getAttributeNames().contains("N")) {
                    //if(NumberUtils.isNumber((String)n.getAttribute("N"))) {

                    try {

                        N_count = (Double) n.getAttribute("N");
                    }
                    catch(Exception e){
                        System.out.println("class: " + n.getAttribute("N").getClass() + " " + n.getAttribute("N").toString());
                        throw new RuntimeException(e);
                    }

                    if (Double.isNaN(N_count)) {

                        N_count = 0;
                    }
                    //}
                }
                if(n.getAttributeNames().contains("S")) {
                    try {
                        S_count = (Double) n.getAttribute("S");
                    }
                    catch(Exception e){
                        System.out.println("class: " + n.getAttribute("S").getClass() + " " + n.getAttribute("S").toString());
                        throw new RuntimeException(e);

                    }
                    if (Double.isNaN(S_count)) {

                        N_count = 0;
                    }
                    // }
                }
                if(n.getAttributeNames().contains("b_u_N")) {
                    //if(NumberUtils.isNumber((String)n.getAttribute("b_u_N"))) {
                    try {
                        uN_count = (Double) n.getAttribute("b_u_N");
                    }
                    catch(Exception e){
                        System.out.println("class: " + n.getAttribute("b_u_N").getClass() + " " + n.getAttribute("b_u_N").toString());
                        throw new RuntimeException(e);

                    }
                    if (Double.isNaN(uN_count)) {
                        uN_count = 0;
                    }
                    //}
                }

                if(n.getAttributeNames().contains("b_u_S")) {
                    //if(NumberUtils.isNumber((String)n.getAttribute("b_u_S"))) {
                    try {
                        uS_count = (Double) n.getAttribute("b_u_S");
                    }
                    catch(Exception e){
                        System.out.println("class: " + n.getAttribute("b_u_S").getClass() + " " + n.getAttribute("b_u_S").toString());
                        throw new RuntimeException(e);
                    }
                    if (Double.isNaN(uS_count)) {
                        uS_count = 0;
                    }
                    // }
                }

                double length = parentHeight - height;
                double prop = 0.0;


                if (length > 0.0) {

                    if (parentHeight < binL && height > binU && parentHeight > binU) {

                        counts[0] += N_count;
                        counts[1] += S_count;
                        counts[2] += N_count;
                        counts[3] += S_count;
                        counts[4] += N_count;
                        counts[5] += S_count;
                        counts[6] += uN_count;
                        counts[7] += uS_count;
                        counts[8] += uN_count;
                        counts[9] += uS_count;

                        lineages++;
                        //System.out.println("length1 " + prop + " count " + counts[0] + " height: " + tree.getHeight(n) + " parentheight: " + tree.getHeight(tree.getParent(n)) + " binL " + binL + " binU " + binU);

                        total_time += length;

                    }
                    if (parentHeight < binL && height < binU && parentHeight > binU) {


                        prop = (parentHeight-binU) / (length);

                        counts[0] += N_count * prop;
                        counts[1] += S_count * prop;
                        counts[2] += N_count * prop;
                        counts[3] += S_count * prop;
                        counts[4] += N_count * prop;
                        counts[5] += S_count * prop;
                        counts[6] += uN_count* prop;
                        counts[7] += uS_count* prop;
                        counts[8] += uN_count* prop;
                        counts[9] += uS_count* prop;

                        lineages++;

                        total_time += prop*length;
                        //System.out.println("node"+n_i+" length2 " + prop + " count " + counts[0] + " height: " + tree.getHeight(n) + " parentheight: " + tree.getHeight(tree.getParent(n)) + " binL " + binL + " binU " + binU);


                    }
                    if (parentHeight > binL && height < binU) {

                        prop = binLength / length;
                        counts[0] += N_count * prop;
                        counts[1] += S_count * prop;
                        counts[2] += N_count * prop;
                        counts[3] += S_count * prop;
                        counts[4] += N_count * prop;
                        counts[5] += S_count * prop;
                        counts[6] += uN_count* prop;
                        counts[7] += uS_count* prop;
                        counts[8] += uN_count* prop;
                        counts[9] += uS_count* prop;

                        lineages++;
                        //System.out.println("node"+n_i+" length3 " + prop + " count " + counts[0] + " height: " + tree.getHeight(n) + " parentheight: " + tree.getHeight(tree.getParent(n)) + " binL " + binL + " binU " + binU);

                        total_time += prop*length;

                    }
                    if (parentHeight > binL && height > binU && height < binL && parentHeight >= binU) {

                        prop = (binL-height) / length;
                        counts[0] += N_count * prop;
                        counts[1] += S_count * prop;
                        counts[2] += N_count * prop;
                        counts[3] += S_count * prop;
                        counts[4] += N_count * prop;
                        counts[5] += S_count * prop;
                        counts[6] += uN_count* prop;
                        counts[7] += uS_count* prop;
                        counts[8] += uN_count* prop;
                        counts[9] += uS_count* prop;

                        lineages++;
                        //System.out.println("node"+n_i+" length4 "+prop+" count "+counts[0]+" height: "+tree.getHeight(n)+" parentheight: "+tree.getHeight(tree.getParent(n))+" binL "+binL+" binU "+binU);
                        total_time += prop*length;
                    }

                    //System.out.println("length "+prop+" height: "+tree.getHeight(n)+" parentheight: "+tree.getHeight(tree.getParent(n))+" binL "+binL+" binU "+binU);

                    if(Double.isNaN(counts[0])) {
                        counts[0] = 0;
                    }
                    if(Double.isNaN(counts[1])) {
                        counts[1] = 0;
                    }

                    if(Double.isNaN(counts[2])) {
                        counts[2] = 0;

                    }
                    if(Double.isNaN(counts[3])) {
                        counts[3] = 0;
                    }
                    if(Double.isNaN(counts[4])) {
                        counts[4] = 0;
                    }
                    if(Double.isNaN(counts[5])) {
                        counts[5] = 0;
                    }
                    if(Double.isNaN(counts[6])) {
                        counts[6] = 0;
                    }
                    if(Double.isNaN(counts[7])) {
                        counts[7] = 0;
                    }
                    if(Double.isNaN(counts[8])) {
                        counts[8] = 0;
                    }
                    if(Double.isNaN(counts[9])) {
                        counts[9] = 0;
                    }

                    n_i++;
                    //System.out.println(lineages);

                }

            }

        }
        //System.out.println(treeLength);

        //treeLength_stat.addValue(treeLength);

        counts[0] /= ((total_time));
        counts[1] /= ((total_time));
        //counts[2] /= (lineages-1);
        //counts[3] /= (lineages-1);
        counts[4] /= ((total_time)*(lineages-1));
        counts[5] /= ((total_time)*(lineages-1));
        counts[6] /= ((total_time));
        counts[7] /= ((total_time));

        if(Double.isNaN(counts[0])) {
            counts[0] = 0;
        }
        if(Double.isNaN(counts[1])) {
            counts[1] = 0;
        }
        if(Double.isNaN(counts[2])) {
            counts[2] = 0;

        }
        if(Double.isNaN(counts[3])) {
            counts[3] = 0;
        }
        if(Double.isNaN(counts[4])) {
            counts[4] = 0;

        }
        if(Double.isNaN(counts[5])) {
            counts[5] = 0;
        }
        if(Double.isNaN(counts[6])) {
            counts[6] = 0;
        }
        if(Double.isNaN(counts[7])) {
            counts[7] = 0;
        }
        if(Double.isNaN(counts[8])) {
            counts[8] = 0;
        }
        if(Double.isNaN(counts[9])) {
            counts[9] = 0;
        }

        return counts;

    }

    public static void main(String args[]) {

//        robustingCountingAnalysis analysis1 = new robustingCountingAnalysis("/Users/jayna/Documents/Projects/AMC_HCV_DATA/BEAST_201605_randomsample/renaissance_counting/p4/p4_combined_renaissance_counting.trees");
//        analysis1.extract_N_and_S_overTime_per_bin(2013.4521, 0.5, "internal");
//        analysis1.extract_N_and_S_overTime_per_bin(2013.4521, 0.5, "external");
//        analysis1.extract_N_and_S_overTime_per_bin(2013.4521, 0.5, "all");
//
//        robustingCountingAnalysis analysis2 = new robustingCountingAnalysis("/Users/jayna/Documents/Projects/AMC_HCV_DATA/BEAST_201605_randomsample/renaissance_counting/p37/p37_combined_renaissance_counting.trees");
//        analysis2.extract_N_and_S_overTime_per_bin(2007.3507, 0.5, "internal");
//        analysis2.extract_N_and_S_overTime_per_bin(2007.3507, 0.5, "external");
//        analysis2.extract_N_and_S_overTime_per_bin(2007.3507, 0.5, "all");
//
//
//        robustingCountingAnalysis analysis3 = new robustingCountingAnalysis("/Users/jayna/Documents/Projects/AMC_HCV_DATA/BEAST_201605_randomsample/renaissance_counting/p53/p53_combined_renaissance_counting.trees");
//        analysis3.extract_N_and_S_overTime_per_bin(2008.347, 0.5, "internal");
//        analysis3.extract_N_and_S_overTime_per_bin(2008.347, 0.5, "external");
//        analysis3.extract_N_and_S_overTime_per_bin(2008.347, 0.5, "all");
//
//
//        robustingCountingAnalysis analysis = new robustingCountingAnalysis("/Users/jayna/Documents/Projects/AMC_HCV_DATA/BEAST_201605_randomsample/renaissance_counting/p61/p61_combined_renaissance_counting.trees");
//        analysis.extract_N_and_S_overTime_per_bin(2013.6438, 0.5, "internal");
//        analysis.extract_N_and_S_overTime_per_bin(2013.6438, 0.5, "external");
//        analysis.extract_N_and_S_overTime_per_bin(2013.6438, 0.5, "all");


//        robustingCountingAnalysis analysis = new robustingCountingAnalysis("/Users/jayna/Dropbox/alignments_and_tree_files/HA-yam/trunk_lineages.yam.HA.1000.trees");
//        analysis.extract_N_and_S_overTime_per_bin(2015.5, 0.5, "internal", "clade1" );

        //robustingCountingAnalysis analysis1 = new robustingCountingAnalysis("/Users/jayna/Dropbox/alignments_and_tree_files/HA-vic/2/no_cpvals.burned.vic.HA.chain1.trees");
//        robustingCountingAnalysis analysis1 = new robustingCountingAnalysis("/Users/jayna/Dropbox/alignments_and_tree_files/HA-vic/trunk_mutations.vic.HA_2.1000.trees");
//        analysis1.extract_N_and_S_overTime_per_bin(2015.3, 0.5, "external", null);
//
//        robustingCountingAnalysis analysis2 = new robustingCountingAnalysis("/Users/jayna/Dropbox/alignments_and_tree_files/HA-yam/trunk_lineages.yam.HA.1000.trees");
//        analysis2.extract_N_and_S_overTime_per_bin(2015.5, 0.5, "all", null);
//
//        robustingCountingAnalysis analysis3 = new robustingCountingAnalysis("/Users/jayna/Dropbox/alignments_and_tree_files/HA-yam/trunk_lineages.yam.HA.1000.trees");
//        analysis3.extract_N_and_S_overTime_per_bin(2015.5, 0.5, "external", null);



        robustingCountingAnalysis analysis1 = new robustingCountingAnalysis("/Users/jayna/Dropbox/ZiBRA_final2_dNdS/evolutionary_rates/genome_wide_dNdS/Pre-Am-ZIKV-lineage_output.trees", "PreAm-ZIKV");
        analysis1.extract_N_and_S_overTime_per_bin(2016.16, 10, "all", null);

//        robustingCountingAnalysis analysis2 = new robustingCountingAnalysis("/Users/jayna/Documents/Projects/DRCandEastAfrica/subA2/combined_A2_RC.trees");
//        analysis2.extract_N_and_S_overTime_per_bin(2011.12716, 0.5, "all", null);
//
//        robustingCountingAnalysis analysis3 = new robustingCountingAnalysis("/Users/jayna/Documents/Projects/DRCandEastAfrica/subD1/subD_skygrid_RC_combined.trees");
//        analysis3.extract_N_and_S_overTime_per_bin(2010.87991, 0.5, "all", null);


        //System.out.println("av treelength: "+analysis.treeLength_stat.getMean() +" "+analysis.treeLength_stat.getValues().length);


//        a2 = 2011.12716;
//        c2 = 2011.032877;
//        d1 = 2010.87991;
    }
}
