#!/usr/bin/env nextflow

params.output		= "$PWD/Piscator_out"
params.primers_csv 	= "primers.csv"
params.blast_RefSeq 	= "$PWD/M32703.fasta"
params.taxonomy_mapping = "taxonomy_mapping.txt"
params.cd_hit_threads  	= 4
params.clustalo_threads = 2
params.ref_fasta 	= "/opt/DB_REF/sel_SILVA.fasta"
cd_hit_clusters 	= Channel.from(0.97)






//sel_SILVA_132_SSURef_Nr99_tax_silva_full_align_trunc.fasta
// //primers_csv = file(params.primers_csv)
// // blast_RefSeq = Channel.value(params.blast_RefSeq)
// ref_fasta = file(ref_fasta)
// taxonomy_mapping = file(params.taxonomy_mapping)
// taxa_coverage_dir = "$output/taxa_coverage"
// hits_ext = '_'+ref_fasta.replace('.fasta','')+'_hits.txt'



output = params.output

Channel.fromPath(params.ref_fasta).into{ref1; ref2; ref3}

hostname = ( params.dockerIP != 'None' ) ? params.dockerIP : 'localhost'   

log.info """
==========================================================================

Piscator: A nexflow pipeline for rapid evaluation of metabarcoding primers
==========================================================================
params.output				= ${params.output}	       
params.primers_csv 		   	= ${params.primers_csv}
mysql DB 				= ${hostname}

"""

"""
params.blast_RefSeq 		   	= ${params.blast_RefSeq}    
params.taxonomy_mapping 	   	= ${params.taxonomy_mapping}
params.cd_hit_threads  		   	= ${params.cd_hit_threads}  
params.clustalo_threads        		= ${params.clustalo_threads}                                                                  
------------------------------------------------------------------------

"""




process Init_PrimerDB{
  
    echo true

    input:
	val hostname
	file params.primers_csv

    output:
	val 'success' into compiled_DB 

        
"""

    init_primerDB.py \
    ${hostname}

    load_primerDB.py \
    ${hostname} \
    -p ${params.primers_csv} 
    
""" 


}




process get_tsv{

    //echo true
    publishDir path: output, mode: 'copy'
    input:
        val hostname
	val return_msge from compiled_DB 

    output:
	stdout Fwd_Rev_pair 
        file ("tsv_files/*.tsv") into (data_tsv,primer_fasta, primers_physchem, primers_analyze) mode flatten
    script:
	tsv_dir="tsv_files"
    
"""
    mkdir ${tsv_dir}
    gen_tsv_primers.py ${hostname}  -o  ${tsv_dir} 
     
"""

}




Fwd_Rev_pair.splitCsv()
	    .flatten()
	    .collate(2)
	    .into{primers_pair; temp }




process Generate_fasta{
    
    publishDir path: "${output}/primers_fasta", mode: 'copy'

    input:
       file primer_tsv from primer_fasta.flatten()
    
    output:
	file primer_tsv into primers_fasta

   
"""

    generate_fasta.py  -p  ${primer_tsv}
      
"""

}




process Get_PhysProp {

    input:
        val hostname
	val primer from primers_physchem.flatten()
    
    output:
	val primer into physchem_data

    
"""
    get_physprop.py \
    ${hostname} \
   -p ${primer}

"""
}




process physchem_plots{

    publishDir path: "${output}/physchem_plots", mode: 'copy'
    input:
         val hostname
	 physchem_data.collect()

    output:
       file '*.pdf' into propplots
   
"""
    plot_physprop.py \
    ${hostname}

"""

}




process phase_pairs{

    echo true
    input:
	file "*" from data_tsv.collect()
	set fwd, rev from  primers_pair

    output:
	file '*_*.tsv' into primersearch_tsv

"""
 
   get_pairs.py \
   $fwd \
   $rev 
   
 
"""

}




process emboss_primersearch{

    input:
	file  ref_fasta from ref1
	file  primer_pair from primersearch_tsv

    output:
        file 'amplicons.rst' into amplicon_results mode flatten
    
    
"""
    primersearch \
    -seqall ${ref_fasta} \
    -infile ${primer_pair} \
    -mismatchpercent 15 \
    -verbose \
    -outfile amplicons.rst 

"""

}




process amplicon_stats{

    publishDir path: output, mode: 'copy'
        echo true
  
    input:
        val hostname
        file amplicons from amplicon_results

    output:
        val '*.tsv' into (amplicon_stats1, amplicon_stats2 ) mode flatten

"""

    amplicon_stats.py \
    ${hostname} \
    --amplicon $amplicons

"""

}



   
process  plot_ampliconData{

    publishDir path: "${output}/EmbossAmplicons", mode: 'move'
    input:
        val  stats from amplicon_stats1.collect()
        file ref_fasta from ref2

    output:
         file '*.pdf'  into Plots
        

"""   

   plot_ampliconData.py \
   ${hostname} \
   -r ${ref_fasta}

"""

}




process amplicon_fasta_gen{

    publishDir path: "${output}/emboss_fasta", mode: 'copy'
    
    input:
        val hostname
        val stats from amplicon_stats2.collect()
        file ref_fasta from ref3
    
    output:
        file "*.fasta" into amplicon_fasta mode flatten


"""

    seqret_amplicons.py \
    ${hostname} \
    -r ${ref_fasta}

"""

}




process cd_hit_est{

    maxForks 1

    publishDir path: "${output}/CD-hit", mode: 'copy'
  
    input:
        each cluster_perc from  cd_hit_clusters
        each amplicons from amplicon_fasta
   
    output:
        file "*cd_hits" into cd_hits
        file "*.cd_hits.clstr" into cdhit_clusters
   
    script:
	cdhit_clusters = amplicons.baseName


"""
    cd-hit-est \
    -i $amplicons \
    -c $cluster_perc \
    -T ${params.cd_hit_threads} \
    -d 0 \
    -r 0 \
    -p 1 \
    -g 1 \
    -o ${cdhit_clusters}_${cluster_perc}.cd_hits 
    
"""

}




process Analyse_clusters{
  
    publishDir path: output, mode: 'copy'  

    input:  
        file fasta_cluster from cdhit_clusters.collect()

    output:
       val  'CDhit_counts*' into hit_plots
 
   
""" 
 
    analyse_clusters.py
     

"""

}




process Analyse_Genetic_dist{

    errorStrategy 'ignore'
    
    publishDir path: output, mode: 'move'
    
    input:  
	file seq_hits from cd_hits.collect()
    
    output:
	file "PrimerVar_*" into distance_plots

   
""" 

    get_distances.py -t  $params.clustalo_threads

     
"""

}




// // process  BlastPrimers{

// // echo true  
  
// // input:
// //     file primer_fasta from primers_fasta.flatten()

    
// // """

// //     BlastPrimers.py -q $primer_fasta -s $params.blast_RefSeq  -d $primerDB

// // """    

// // }


// process Analyze_primers {
  
// publishDir path: output, mode: 'copy'
  
// input:
//       val  primer from primers_analyze.flatten() 
//       file ref_fasta
      
// output:
//      file '*_hits.txt' into primer_results1
//      file "*.ps" into analysis_charts
     
// """

//    analyze_primers.py -f  $ref_fasta -P $primer

// """

// }




// process Pair_hits{

// maxForks 1

// input:
//     set fwd, rev from  primers_raw
                                
// output:
//     stdout data_out  into hits
//     set file(fwd_hit),  file(rev_hit) into primer_hits1,
//                                primer_hits2

// script:
//     fwd_hit = fwd.replace('.tsv','')+hits_ext 
//     rev_hit = rev.replace('.tsv','')+hits_ext

    
// """

//     pairing_block.py -p  "$output/$fwd_hit" "$output/$rev_hit"
  
     
// """

// }


// process Get_amplicons_and_reads{

//   //publishDir path: output, mode: 'copy'

// input:
//      set fwd_hit, rev_hit from primer_hits1
     
// output:
//      file '*amplicons.fasta'  into amplicon_seqs
//      file '*reads.fasta'  into amplicon_reads


// """
   
//    get_amplicons_and_reads.py -f $ref_fasta -R 300 -d p -i $fwd_hit:$rev_hit

// """

// }




// process Analyse_amplicons{

// errorStrategy 'ignore'
// publishDir path: output, mode: 'copy'

// input:
//     file amplicons from amplicon_seqs

// output:
//     file '*.ps' into amplicon_histograms

    
// """

//    amplicons_histograms.py -f $amplicons

// """

// }




// process Read_taxonomy{
  
// errorStrategy 'ignore'
// publishDir path: output, mode: 'copy'

  
// input:
//    file reads from amplicon_reads.flatten()
//    file taxonomy_mapping
   
// output:
//     file "*reads_accuracy_report.txt" into reads_accuracy_reports
//     file "*reads_assignments.txt" into reads_assignments

// """

//    taxa_assignment_report.py -t $params.taxonomy_mapping -f $reads

// """

// }




// process Taxa_coverage{
 

// publishDir path: taxa_coverage_dir, mode: 'copy'

// input:
//     set fwd_hit,  rev_hit from primer_hits2
//     file taxonomy_mapping

    

// """

//    taxa_coverage.py -i $fwd_hit:$rev_hit --taxa_fp $taxonomy_mapping -p -d 7 -o $taxa_coverage_dir

// """    

// }


// workflow.onComplete {

//   println """

//     Pipeline execution summary
//     ---------------------------
//     Completed at: ${workflow.complete}
//     Duration    : ${workflow.duration}
//     Success     : ${workflow.success}
//     workDir     : ${workflow.workDir}
//     exit status : ${workflow.exitStatus}
//     Error report: ${workflow.errorReport ?: '-'}
//     """
// }
