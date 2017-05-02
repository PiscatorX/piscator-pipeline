#!/usr/bin/env nextflow

output = "Piscator_results"
params.primers_csv = "primer_select.csv"
params.db_name = "Piscator.db"
params.primerDB = "$baseDir/$output/$params.db_name"
primers_csv = file(params.primers_csv)
params.blast_RefSeq = "$PWD/M32703.fasta"
// params.ref_fasta = "Ref_SeqX.fasta"
// params.ref_fasta = "temp.fasta"
// params.taxonomy_mapping = ("taxonomy_mapping.txt")
// blast_RefSeq = Channel.value(params.blast_RefSeq)
// ref_fasta = file(params.ref_fasta)
// get_amplicon_sql = file("get_amplicon.sql")
// taxonomy_mapping = file(params.taxonomy_mapping)
// taxa_coverage_dir = "$output/taxa_coverage"
// cd_hit_clusters = Channel.from(0.80, 0.90, 0.95, 0.97, 1.00)



process Init_PrimerDB{

publishDir path: output, mode: 'copy'

input:
       file primers_csv
       val  params.db_name  

output:
       file '*.db' into primers
       
"""
    Init_PrimerDB.py --primers-file $primers_csv --db-name $params.db_name
     
""" 

}


process gen_tsv_primers{
  
publishDir path: output, mode: 'copy'

  
input:
   val primerDB from primers 

output:
    stdout fwd_rev_pair 
    file '*.tsv' into  primer_files

    
"""
     gen_tsv_primers.py  --out $output --db $primerDB 
   
"""

}


process generate_fasta{

echo true

publishDir path: output, mode: 'copy'

input:
   val primer_tsv from primer_files.flatten()

   
output:
   file '*.fasta' into primers_fasta

   
"""

     generate_fasta.py  -p  $primer_tsv
      
"""
}


process  BlastPrimers{

echo true  
  
input:
    file primer_fasta from primers_fasta.flatten()

    
"""

    BlastPrimers.py -q $primer_fasta -s $params.blast_RefSeq  -d $params.primerDB

"""    

}


//
// process Data_splitter{
//   //This process does nothing but create other channels
//   input:
//     val 
//     set fwd, rev  from fwd_rev_pair.splitCsv()
//                                    .flatten()
//                                    .map{ it }
//                                    .collate(2)
// output:
//    set fwd, rev into primers_analyze,
//                      primers_raw,
//    		     primers_tsv,
//    		     primers_files
// """
//   :
// """
// }
// process get_PhysChem {
// echo true
// input:
//     val primer from primer_files.flatten()
// output:
//     val primer into physchem_data  
// """
//    echo get_PhysChem.py -p $primer -d  $params.primerDB
// """
// }
// save_data = Channel.create()
// physchem_data.collectFile()
//              .toList()
// 	     .into(save_data)
// process physchem_plots{
// input:
//    val complete from save_data
// """
//     physchem.R  -c $baseDir  
// """
// }
// primers_fasta.subscribe{ print it }
// process analyze_primers {
//   publishDir path: output, mode: 'copy'
// input:
//       val primer from primers_analyze.flatten() 
//       file ref_fasta
// output:
//      file '*_hits.txt' into primer_results1,
//                             primer_results2 
// """
//    analyze_primers.py -f  $ref_fasta -P $primer
// """
// }
// hits_ext = '_'+params.ref_fasta.replace('.fasta','')+'_hits.txt'
// (primer_data,
// primer_hits) = primers_raw.flatten()
//                       .map{it.replace('.tsv','')+hits_ext}
// 		      .collate(2)
//                       .separate(2){hits -> [hits, hits]}                   
// process get_amplicons_and_reads{
//   publishDir path: output, mode: 'copy'
// input:
//      set fwd, rev from primer_hits
//      val hits from primer_results1      
// output:
//      file '*amplicons.fasta'  into amplicon_seqs
//      file '*reads.fasta'  into amplicon_reads
// """ 
//    get_amplicons_and_reads.py -f $ref_fasta -i $fwd:$rev
// """
// }
// process get_pairs{
// input:
//    set  fwd, rev from primers_tsv
// output:
//    file '*.tsv' into primersearch_tsv
// """
//     get_pairs.py $fwd $rev 
// """
// }
// process analyse_amplicons{
//   publishDir path: output, mode: 'copy'
// errorStrategy 'ignore'
// input:
//     file amplicons from amplicon_seqs
// output:
//     file '*.ps' into amplicon_histograms
// """
//    amplicons_histograms.py -f $amplicons
// """
// }
// process read_taxonomy{
//   publishDir path: output, mode: 'copy'
// input:
//    file reads from amplicon_reads
//    file taxonomy_mapping
// output:
//     file "*reads_accuracy_report.txt" into reads_accuracy_reports
//     file "*reads_assignments.txt" into reads_assignments
// """
//    taxa_assignment_report.py -t $params.taxonomy_mapping -f $reads
// """
// }
// process taxa_coverage{
//   publishDir path: taxa_coverage_dir, mode: 'copy'
// input:
//     set fwd, rev from primer_data
//     val hits from primer_results2
//     file taxonomy_mapping
// """
//    taxa_coverage.py -i $fwd:$rev --taxa_fp $taxonomy_mapping -p -d 7 -o $taxa_coverage_dir
// """    
// }
// process emboss_primersearch{
// input:
//     file  ref_fasta
//     file  primer_pair from primersearch_tsv
// output:
//     file 'amplicons.rst' into amplicon_results mode flatten
// """
//   primersearch -seqall $ref_fasta -infile $primer_pair -mismatchpercent 15 -verbose  -outfile amplicons.rst 
// """
// }
// process Amplicon_stats{
// input:
//     file amplicons from amplicon_results
// output:
//     file '*.tsv' into Amplicon_stats mode flatten
//     val 1 into results
// """
//     Amplicon_stats.py --amplicon $amplicons
// """
// }
// process Init_ampliconData{
// input:
//     file ref_fasta
//     val  primerDB
// """
//    Init_ampliconData.py --ref-fasta $ref_fasta --db-name $primerDB 
// """
// }
// process  Compile_ampliconData{
// maxForks 1
// errorStrategy 'retry'
// maxErrors 5
// input:
//     val  primerDB 
//     file stats from Amplicon_stats
// output:
//     val  primerDB into compiledata
// """
//    Compile_ampliconData.py --tsv $stats --db $primerDB
// """
// }
// summarise_data = compiledata.collectFile()
// process amplicon_summary{
// input:
//      val complete from summarise_data
//      val primerDB
//      file get_amplicon_sql 
// output:
//      file 'amplicon_data.tsv' into amplicon_data
    
// """
//      sqlite3 $primerDB < $get_amplicon_sql
// """
// }
// process amplicon_fasta_gen{
// input:
//     file amplicons from amplicon_data
//     file ref_fasta
// output:
//      file "*.fasta" into amplicon_fasta mode flatten
// """
//    seqret_amplicons.py  $ref_fasta $amplicons
// """
// }
// process cd_hit_est{
// maxForks 1 
// publishDir path: output, mode: 'copy'
// input:
//    each cluster_perc from  cd_hit_clusters
//    each amplicons  from  amplicon_fasta
// output:
//    file "*cd_hits*" into cd_hits
//    stdout output into  pipe_out 
// script:
// cd_hit_clusters = amplicons.baseName
// """
//     cd-hit-est -i $amplicons  -c $cluster_perc -d 0 -r 0 -p 1 -g 1  -o \
//     ${cd_hit_clusters}_${cluster_perc}.cd_hits 
// """
// }
// verbose.subscribe{print it}
