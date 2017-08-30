#!/usr/bin/env nextflow

output = "$PWD/Piscator_results"
params.primers_csv = "primer_select.csv"
//params.primers_csv = "OptimusPrimerDB.csv"
params.db_name = "Piscator.db"
primerDB = "$output/$params.db_name"
primers_csv = file(params.primers_csv)
params.blast_RefSeq = "$PWD/M32703.fasta"
//params.ref_fasta = "SILVA_DB.fasta"
params.ref_fasta = "temp.fasta"
params.taxonomy_mapping = ("taxonomy_mapping.txt")
// blast_RefSeq = Channel.value(params.blast_RefSeq)
ref_fasta = file(params.ref_fasta)
get_amplicon_sql = file("get_amplicon.sql")
taxonomy_mapping = file(params.taxonomy_mapping)
taxa_coverage_dir = "$output/taxa_coverage"
cd_hit_clusters = Channel.from(0.80, 0.90, 0.97, 1.00)
hits_ext = '_'+params.ref_fasta.replace('.fasta','')+'_hits.txt'


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
     gen_tsv_primers.py  --output  $output  --db $primerDB 
   
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
		     primers_physchem
"""
  :

"""

}


process get_pairs{
  
input:
   set  fwd, rev from primers_tsv
output:
   file '*.tsv' into primersearch_tsv
   
"""
    get_pairs.py $fwd $rev 
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
input:
    file amplicons from amplicon_results
output:
    file '*.tsv' into Amplicon_stats mode flatten
    
"""
    Amplicon_stats.py --amplicon $amplicons
"""
}


process  Compile_ampliconData{

maxForks 1

input:
     file stats from Amplicon_stats
     val  primerDB
output:
    file  stats into compiledata
    
"""
   Init_ampliconData.py --ref-fasta $ref_fasta --db-name $primerDB 
   Compile_ampliconData.py --tsv $stats --db $primerDB  -o $output  -r  $ref_fasta

"""

}

summarise_data = compiledata.collectFile()

process amplicon_summary{

input:
     val complete from summarise_data
     val primerDB
     file get_amplicon_sql 
output:
     file 'amplicon_data.tsv' into amplicon_data
    
"""
     sqlite3 $primerDB < $get_amplicon_sql
"""

}

process amplicon_fasta_gen{

input:
    file amplicons from amplicon_data
    file ref_fasta
output:
    file "*.fasta" into amplicon_fasta mode flatten
     
        
"""
   seqret_amplicons.py  $ref_fasta $amplicons

"""

}


process cd_hit_est{
 
publishDir path: output, mode: 'copy'

input:
   each cluster_perc from  cd_hit_clusters
   each amplicons from amplicon_fasta
   
output:
   file "*cd_hits*" into cd_hits
   stdout output into cd_clusters
   
script:
   cd_hit_clusters = amplicons.baseName


"""
    cdhit-est -i $amplicons  -c $cluster_perc -d 0 -r 0 -p 1 -g 1  -o \
    ${cd_hit_clusters}_${cluster_perc}.cd_hits 
    
"""

}

generate_plots = cd_clusters.collectFile()

process Genplots{

input:  
   val generate_plots



"""
  
    Align_Dist.py -o $output
    GenPlots.py -w $output

"""

}


process  BlastPrimers{

echo true  
  
input:
    file primer_fasta from primers_fasta.flatten()

    
"""

    BlastPrimers.py -q $primer_fasta -s $params.blast_RefSeq  -d $primerDB

"""    

}


process get_PhysChem {
echo true
input:
    val primer from primers_physchem.flatten()
output:
    val primer into physchem_data
    
"""
   get_PhysChem.py -p $primer -d  $primerDB
"""
}

save_data = Channel.create()
physchem_data.collectFile()
             .toList()
	     .into(save_data)

process physchem_plots{
input:
   val complete from save_data
"""
    physchem.R  -c $baseDir -o $output
"""
}

process analyze_primers {
  publishDir path: output, mode: 'copy'
input:
      val primer from primers_analyze.flatten() 
      file ref_fasta
output:
     file '*_hits.txt' into primer_results1,
                            primer_results2
			    
     file "*.ps" into analysis_charts
     
"""
   analyze_primers.py -f  $ref_fasta -P $primer

"""

}


process pair_hits{

maxForks 1
errorStrategy 'ignore'

input:
   set fwd, rev from  primers_raw
                                
output:
    stdout data_out  into hits
    set fwd_hit,  rev_hit into primer_data,
                               primer_hits

script:
   fwd_hit = fwd.replace('.tsv','')+hits_ext 
   rev_hit = rev.replace('.tsv','')+hits_ext

"""
    pairing_block.py -p  $fwd_hit $rev_hit

"""

}


process get_amplicons_and_reads{

publishDir path: output, mode: 'copy'

input:
     set fwd, rev from primer_hits
     val hits from primer_results1      
output:
     file '*amplicons.fasta'  into amplicon_seqs
     file '*reads.fasta'  into amplicon_reads

"""
   
   get_amplicons_and_reads.py -f $ref_fasta -R 300 -d p -i $fwd:$rev

"""

}

process analyse_amplicons{
  publishDir path: output, mode: 'copy'
errorStrategy 'ignore'
input:
    file amplicons from amplicon_seqs
output:
    file '*.ps' into amplicon_histograms
"""
   amplicons_histograms.py -f $amplicons
"""
}

process read_taxonomy{
  
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

process taxa_coverage{
  publishDir path: taxa_coverage_dir, mode: 'copy'
input:
    set fwd, rev from primer_data
    val hits from primer_results2
    file taxonomy_mapping

"""

   taxa_coverage.py -i $fwd:$rev --taxa_fp $taxonomy_mapping -p -d 7 -o $taxa_coverage_dir

"""    

}

