#!/bin/bash

head -n 34759 results/recount/tmp/leiden_rsweep.seed > results/recount/clustering/leiden_rsweep_pca_over50_pc3sd_tpm_log_nv100_90th_var_genes_k30_sim.tsv

for i in 0.007 0.0075 0.008 0.0085 0.009 0.0095 0.01 0.0105 0.011 0.0115 0.012 0.0125 0.013
do

    java -jar lib/NetworkClustering/RunNetworkClustering.jar -r $i -o results/recount/tmp/leiden_pca_over50_pc3sd_tpm_log_nv100_k30_sim_90th_var_genes.tmp results/recount/clustering/nn_sim_weights_recount_over50_pc3sd_tpm_log_nv100_90th_var_genes_k30_filtered.tsv 

    cut -f2,2 results/recount/tmp/leiden_pca_over50_pc3sd_tpm_log_nv100_k30_sim_90th_var_genes.tmp | paste results/recount/clustering/leiden_rsweep_pca_over50_pc3sd_tpm_log_nv100_90th_var_genes_k30_sim.tsv - > results/recount/tmp/leiden_rsweep_pca_over50_pc3sd_tpm_log_nv100_90th_var_genes_k30_sim.tsv

    mv results/recount/tmp/leiden_rsweep_pca_over50_pc3sd_tpm_log_nv100_90th_var_genes_k30_sim.tsv results/recount/clustering/leiden_rsweep_pca_over50_pc3sd_tpm_log_nv100_90th_var_genes_k30_sim.tsv
done
