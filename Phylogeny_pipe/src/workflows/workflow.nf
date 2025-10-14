#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {rename_fasta} from "${params.nfpath}/modules/module.nf"
include {distance_matrix_mash} from "${params.nfpath}/modules/module.nf"
include {distance_matrix_analysis} from "${params.nfpath}/modules/module.nf"

include {map_bowtie2} from "${params.nfpath}/modules/module.nf"
include {polish_pilon} from "${params.nfpath}/modules/module.nf"
include {duplicate_masker_repeatmasker} from "${params.nfpath}/modules/module.nf"
include {duplicate_masker_bedops} from "${params.nfpath}/modules/module.nf"
include {create_input_tab} from "${params.nfpath}/modules/module.nf"
include {core_snps_snippy} from "${params.nfpath}/modules/module.nf" 
include {vc_recombination_analysis_gubbins} from "${params.nfpath}/modules/module.nf"
include {vc_final_snp_matrix} from "${params.nfpath}/modules/module.nf"

include {ref_phylogeny} from "${params.nfpath}/modules/module.nf" 
include {create_input_tab_phylogeny} from "${params.nfpath}/modules/module.nf" 
include {core_snps_snippy_phylogeny} from "${params.nfpath}/modules/module.nf" 
include {snps_tree_iqtree} from "${params.nfpath}/modules/module.nf" 
include {rec_removal_clonalframeml} from "${params.nfpath}/modules/module.nf"
include {dating_treetime} from "${params.nfpath}/modules/module.nf"
include {recombination_analysis_gubbins} from "${params.nfpath}/modules/module.nf"

// workflow script
workflow bacteria_mash_clustering {
    take:
        ch_replicons

    main:

          //Step0- Grouping samples by replicon
          rename_fasta(ch_replicons)
          rename_fasta.out.fasta_renamed.set{ ch_replicons }
          ch_replicons.groupTuple(by: 0).set{ ch_replicons }

          //Step1- Compute MASH distance matrix
          if ( params.mash ) {
               distance_matrix_mash(ch_replicons)
               distance_matrix_mash.out.replicons_ch.set{ ch_replicons }
               distance_matrix_analysis(ch_replicons)
          }
}

workflow bacteria_variant_calling {
    take:
        ch_replicons

    main:

          //Step0- Preparation of reference genome
          if ( params.ref_preparation ) {
               map_bowtie2(ch_replicons)
               map_bowtie2.out.sorted_bam_files.set{ ch_replicons }
               polish_pilon(ch_replicons)
               polish_pilon.out.polished_reference.set{ ch_replicons }
               duplicate_masker_repeatmasker(ch_replicons)
               duplicate_masker_repeatmasker.out.replicons_ch.set{ ch_replicons }               
               duplicate_masker_bedops(ch_replicons)
               duplicate_masker_bedops.out.replicons_ch.set{ ch_replicons }
          }

          //Step1- Core SNPs calling
          if ( params.snippy ) {
               create_input_tab(ch_replicons)
               create_input_tab.out.input_tab.set{ ch_input_tab }
               core_snps_snippy(ch_input_tab)
               core_snps_snippy.out.clean_core_snps.set{ ch_core_snps }
               vc_recombination_analysis_gubbins(ch_core_snps)
               vc_recombination_analysis_gubbins.out.aln_without_rec.set{ ch_core_snps }
               vc_final_snp_matrix(ch_core_snps)
          }
}

workflow bacteria_phylogeny {
    take:
        ch_replicons

    main:
          
          //Step0- Preparation of reference genome
          if ( params.ref ) {
               ref_phylogeny(ch_replicons)
               ref_phylogeny.out.replicons_ch.set{ ch_replicons }
          }

          //Step1- Core SNPs calling
          if ( params.snippy ) {
               create_input_tab_phylogeny(ch_replicons)
               create_input_tab_phylogeny.out.input_tab.set{ ch_input_tab }
               core_snps_snippy_phylogeny(ch_input_tab)
               core_snps_snippy_phylogeny.out.core_snps_aln.set{ ch_core_snps }
          }

          //Step2- Core Tree construction on SNPS
          if ( params.iqtree ) {
               snps_tree_iqtree(ch_core_snps)
               snps_tree_iqtree.out.treefile.set{ ch_treefiles }
          }

          //Step3- Recombination correction on core tree
          if ( params.clonalframeml ) {
               rec_removal_clonalframeml(ch_treefiles)
               rec_removal_clonalframeml.out.tree_without_rec.set{ ch_treefiles }
          }

          //Step4- Phylogenetic dating
          if ( params.treetime ) {
               //dating_treetime(ch_treefiles)
               //dating_treetime.out.snippy_aln.set{ ch_treefiles }
          }

          //Step5- Gubbins recombination analysis
          if( params.gubbins ) {
               recombination_analysis_gubbins(ch_treefiles)
               recombination_analysis_gubbins.out.aln_without_rec.set{ ch_treefiles }
          }




}
