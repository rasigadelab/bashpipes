#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {annotate_prokka} from "${params.nfpath}/modules/module.nf"
include {pan_genome_panaroo} from "${params.nfpath}/modules/module.nf"
include {core_tree_iqtree} from "${params.nfpath}/modules/module.nf"
include {create_input_tab} from "${params.nfpath}/modules/module.nf"
include {core_snps_snippy} from "${params.nfpath}/modules/module.nf" 


// workflow script
workflow bacteria_phylogeny {
    take:
        ch_replicons

    main:
       
       //Step0- Assembly annotation
       if ( params.annotate ) {
            annotate_prokka(ch_replicons)
            annotate_prokka.out.annotation_file.set{ ch_replicons }
       }
        
        ch_replicons.groupTuple(by: 0).set{ ch_replicons }
         
       //Step1- Pan-genome analysis
       if ( params.panaroo ) {
            pan_genome_panaroo(ch_replicons)
            pan_genome_panaroo.out.core_genome_aln.set{ ch_core_aln }
            pan_genome_panaroo.out.core_genome_ref.set{ ch_core_ref }
       }

       //Step2- Fast Core Tree construction
       if ( params.iqtree_fast ) {
            core_tree_iqtree(ch_core_aln, "iqtree_after_panaroo")
       }
        
        //Step3- Core SNPs calling
        if ( params.snippy ) {
            create_input_tab(ch_replicons)
            create_input_tab.out.input_tab.set{ ch_input_tab }
            ch_core_ref.join(ch_input_tab).set{ ch_core_ref }
            core_snps_snippy(ch_core_ref)
            core_snps_snippy.out.core_snps_aln.set{ ch_core_snps }
            core_snps_snippy.out.core_full_aln.set{ ch_core_full }
        }

        //Step4- Core Tree construction on SNPS
        if ( params.iqtree_snps ) {
          //core_tree_iqtree(ch_core_full, "iqtree_after_snippy_full")
          //snps_tree_iqtree(ch_core_aln, "iqtree_after_snippy_snps")
        }


}