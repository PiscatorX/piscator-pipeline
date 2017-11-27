#!/usr/bin/env nextflow

output = "$PWD/Piscator_out"
//params.primers_csv = "primer_select.csv"
params.primers_csv = "OptimusPrimerDB.csv"
params.db_name = "Piscator.db"
primerDB = "$output/$params.db_name"
primers_csv = file(params.primers_csv)
params.blast_RefSeq = "$PWD/M32703.fasta"
//params.ref_fasta = "SILVA_DB.fasta"
params.ref_fasta = "temp2.fasta"
params.taxonomy_mapping = ("taxonomy_mapping.txt")
// blast_RefSeq = Channel.value(params.blast_RefSeq)
ref_fasta = file(params.ref_fasta)
taxonomy_mapping = file(params.taxonomy_mapping)
taxa_coverage_dir = "$output/taxa_coverage"
cd_hit_clusters = Channel.from(0.97, 1.00)
//cd_hit_clusters = Channel.from(0.80, 0.90, 0.97, 1.00)
cd_hit_threads  =  2
clustalo_threads = 2
hits_ext = '_'+params.ref_fasta.replace('.fasta','')+'_hits.txt'




process Init_PrimerDB{
  
echo true

input:
    file primers_csv
       
output:
    val 'sucess' into DB

        
"""
    init_primerDB.py --primers-file $primers_csv 
    
""" 

}



process Load_PrimerDB{

echo true
  
input:
    val return_msge  from DB

output:
    val 'sucess' into compiled_DB

    
    
"""
    load_primerDB.py  --primers-file $primers_csv 
    
""" 

}



process Get_pairs{
  
publishDir path: output, mode: 'copy'
  
input:
   val return_msge from compiled_DB 

output:
    stdout fwd_rev_pair 
    file '*.tsv' into  primer_files
 
    
"""
     gen_tsv_primers.py 
   
"""

}



process Data_splitter{

//This process does nothing but create other channels

input:
   
    set fwd, rev  from fwd_rev_pair.splitCsv()
                                   .flatten()
                                   .map{ it }
                                   .collate(2)
output:
   set fwd, rev into primers_analyze,
                     primers_raw,
   		     primers_tsv,
   		     primers_files,
		     primers_physchem,
		     primer_data
		     
"""
  :

"""

}



process Generate_fasta{

publishDir path: output, mode: 'copy'

input:
   val primer_tsv from primer_files.flatten()
   
output:
   file '*.fasta' into primers_fasta

   
"""

     generate_fasta.py  -p  $primer_tsv
      
"""

}




process Get_PhysProp {

input:
    val primer from primers_physchem.flatten()
    
output:
    val primer into physchem_data
    
"""
   get_physprop.py -p $output/$primer

"""
}


save_data = physchem_data.collectFile()


process physchem_plots{

publishDir path: output, mode: 'copy'

input:
   val complete from save_data

   
"""

    plot_physprop.py

"""

}


process Get_Pairs{
  
input:
   set  fwd, rev from primers_tsv

output:
   file '*.tsv' into primersearch_tsv

   
"""

    get_pairs.py $output/$fwd $output/$rev 

"""

}


process emboss_primersearch{
  
input:
    file  ref_fasta
    file  primer_pair from primersearch_tsv

output:
    file 'amplicons.rst' into amplicon_results mode flatten
    
    
"""
   primersearch -seqall $ref_fasta -infile $primer_pair -mismatchpercent 15 -verbose  -outfile amplicons.rst 

"""

}


process Amplicon_stats{

publishDir path: output, mode: 'copy'
echo true
  
input:
    file amplicons from amplicon_results

output:
    val '*.tsv' into Amplicon_stats mode flatten

    
"""
    Amplicon_stats.py --amplicon $amplicons

"""

}


plotStats = Amplicon_stats.collectFile()
 	    
   
process  plot_ampliconData{
	    
publishDir path: output, mode: 'copy'

input:
     val stats from plotStats
     
output:
    file '*.pdf'  into Plots
    val  'success' into AmpliconPlots
    
"""   
   plot_ampliconData.py   -r  $ref_fasta

"""

}


process amplicon_fasta_gen{

publishDir path: output, mode: 'copy'  
  
input:
    val return_msge from AmpliconPlots
    file ref_fasta
    
output:
    file "*.fasta" into amplicon_fasta mode flatten
            
"""
   seqret_amplicons.py -r $ref_fasta

"""

}


process cd_hit_est{

maxForks 1

publishDir path: output, mode: 'copy'

  
input:
   each cluster_perc from  cd_hit_clusters
   each amplicons from amplicon_fasta
   
output:
   file "*cd_hits" into cd_hits
   file "*.cd_hits.clstr" into cdhit_clusters
   
script:
   cdhit_clusters = amplicons.baseName


"""
    cdhit-est -i $amplicons  -c $cluster_perc -T $cd_hit_threads -d 0 -r 0 -p 1 -g 1  -o \
    ${cdhit_clusters}_${cluster_perc}.cd_hits 
    
"""

}


seq_clusters = cdhit_clusters.collectFile()
                             .toList()



process Analyse_clusters{
  
publishDir path: output, mode: 'copy'  


input:  
   file fasta_cluster from seq_clusters


output:
   val  'CDhit_counts*' into hit_plots

   

""" 
 
    analyse_clusters.py
     

"""

}


hit_clusters = cd_hits.collectFile()
                             .toList()


 
process Analyse_genetic_distances{

echo true
  
publishDir path: output, mode: 'copy'  


input:  
    file seq_hits from hit_clusters


output:
   file "PrimerVar_*" into distance_plots

   
""" 

    get_distances.py
     
"""

}




// process  BlastPrimers{

// echo true  
  
// input:
//     file primer_fasta from primers_fasta.flatten()

    
// """

//     BlastPrimers.py -q $primer_fasta -s $params.blast_RefSeq  -d $primerDB

// """    

// }


process Analyze_primers {
  
publishDir path: output, mode: 'copy'
  
input:
      val primer from primers_analyze.flatten() 
      file ref_fasta
      
output:
     file '*_hits.txt' into primer_results1
     file "*.ps" into analysis_charts
     
"""

   analyze_primers.py -f  $ref_fasta -P $output/$primer

"""

}




process Pair_hits{

maxForks 1

input:
    set fwd, rev from  primers_raw
                                
output:
    stdout data_out  into hits
    set file(fwd_hit),  file(rev_hit) into primer_hits1,
                               primer_hits2

script:
    fwd_hit = fwd.replace('.tsv','')+hits_ext 
    rev_hit = rev.replace('.tsv','')+hits_ext

    
"""

    pairing_block.py -p  "$output/$fwd_hit" "$output/$rev_hit"
  
     
"""

}


process Get_amplicons_and_reads{

publishDir path: output, mode: 'copy'

input:
     set fwd_hit, rev_hit from primer_hits1
     
output:
     file '*amplicons.fasta'  into amplicon_seqs
     file '*reads.fasta'  into amplicon_reads


"""
   
   get_amplicons_and_reads.py -f $ref_fasta -R 300 -d p -i $fwd_hit:$rev_hit

"""

}




process Analyse_amplicons{

errorStrategy 'ignore'
publishDir path: output, mode: 'copy'

input:
    file amplicons from amplicon_seqs

output:
    file '*.ps' into amplicon_histograms

    
"""

   amplicons_histograms.py -f $amplicons

"""

}




process Read_taxonomy{
  
errorStrategy 'ignore'
publishDir path: output, mode: 'copy'

  
input:
   file reads from amplicon_reads.flatten()
   file taxonomy_mapping
   
output:
    file "*reads_accuracy_report.txt" into reads_accuracy_reports
    file "*reads_assignments.txt" into reads_assignments

"""

   taxa_assignment_report.py -t $params.taxonomy_mapping -f $reads

"""

}




process Taxa_coverage{
 

publishDir path: taxa_coverage_dir, mode: 'copy'

input:
    set fwd_hit,  rev_hit from primer_hits2
    file taxonomy_mapping
    

"""

   taxa_coverage.py -i $fwd_hit:$rev_hit --taxa_fp $taxonomy_mapping -p -d 7 -o $taxa_coverage_dir

"""    

}


workflow.onComplete {

  println """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}
