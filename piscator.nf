#!/usr/bin/env nextflow

python_virtualenv       = workflow.projectDir+"/python_virtualenv"
//virtualenv env path for python modules
//this allows for new versions of python plotting modules
params.output		= "$PWD/piscator.Out"
params.primers_csv 	= "primers.csv4"
params.blast_RefSeq 	= "$PWD/M32703.fasta"
params.ref_fasta 	= "$PWD/SILVA.fasta"
params.taxonomy_mapping	= "$PWD/SILVA.taxonomy"
cd_hit_clusters 	= Channel.from(0.80, 0.97, 1.0)
ref_fasta               = Channel.value(params.ref_fasta)
taxonomy_mapping        = Channel.value(params.taxonomy_mapping)
ref_fasta_path          = file(params.ref_fasta)
output                  = params.output
primers_csv             = Channel.fromPath(params.primers_csv)
hostname                = 'localhost'
params.taxa_depth	= 6

//flags for debugging
params.analyze_primers  = true
params.analyse_amplicons= true
params.taxa_coverage	= true
taxa_depth              = Channel.value(params.taxa_depth)

if ( !ref_fasta.endsWith(".fasta")) {
    
    println "Fail! reference fasta sequence must have `.fasta`."
    
    exit(1)
}


hits_base = '_' + ref_fasta_path.getName().replace('.fasta','_hits.txt')


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
    errorStrategy 'retry'
    maxRetries 3
    input:
	val hostname
	file primers_csv
	val taxonomy_mapping
	
    output:
	val 'success' into compiled_DB  

    script:
	hostname = ( params.dockerIP != 'None' ) ? params.dockerIP : 'localhost'   

"""
    
    init_primerDB.py \
    ${hostname}

    load_primerDB.py \
    ${hostname} \
    -p ${primers_csv} \
    -t ${taxonomy_mapping}
    
""" 
 
}




process get_tsv{

    //echo true
    publishDir path: "${output}/tsv_files", mode: 'copy'
    
    input:
        val return_msge from compiled_DB 
        val hostname


    output:
	 stdout Fwd_Rev_pair 
         file ("*.tsv") into (pair_tsv, primer_fasta, primers_physchem, primers_analyze) mode flatten

"""

    gen_tsv_primers.py \
    ${hostname}  
    
"""

}




Fwd_Rev_pair.splitCsv()
	    .flatten()
	    .collate(2)
	    .into{primers_pair1;
	          primers_pair2}




process analyze_primers{
    
    
    publishDir path: "${output}/primers_analysis", mode: 'copy'
    
    input:
         file primer from primers_analyze
         val ref_fasta
	 
    output:
        file '*_hits.txt' into (primer_results1, primer_results2)
        file "*.ps" into analysis_charts
	
    when:
       params.analyze_primers == true
    
"""

    analyze_primers.py \
        -f  $ref_fasta \
        -P $primer

"""

}




process get_amplicons_and_reads{

    publishDir path: "${output}/amplicons", mode: 'copy'
    
    input:
	set fwd_tsv, rev_tsv  from primers_pair1
	file(hits) from primer_results1.collect()
        val  ref_fasta	

   output:
        file '*reads.fasta'  into reads
        file '*amplicons.fasta'  into amplicons1, amplicons2 
        set fwd_hit, rev_hit into primer_hits
	
    when:
        params.taxa_coverage == true

   script:  
       
       fwd_hit = fwd_tsv.replace('.tsv','') + hits_base
       rev_hit = rev_tsv.replace('.tsv','') + hits_base 


    
"""

    get_amplicons_and_reads.py \
	-i $fwd_hit:$rev_hit \
	-R 300 \
	-f $ref_fasta \
	-d p \
	-v

"""

}




process analyse_amplicons{

    errorStrategy 'ignore'
    publishDir path: "${output}/amplicon_plots", mode: 'copy'

    input:
        file amplicons from amplicons1

    output:
        file '*.ps' into amplicon_histograms

    when:
       params.analyse_amplicons == true
    
"""

   amplicons_histograms.py \
       -f $amplicons

"""

}




process taxa_coverage{
 
    
    publishDir path: "${output}/taxonomy_coverage", mode: 'copy'
    input:
        set fwd_hit, rev_hit from primer_hits
        file hits from primer_results2.collect()	
        val taxonomy_mapping
	val taxa_depth

    output:
	file("${pair_id}") into taxa_coverage_plots

    when:
        params.taxa_coverage == true

    script:
        pair_id =  fwd_hit.replace(hits_base,'')+'_'+rev_hit.replace(hits_base,'')        
     	 
"""

   taxa_coverage.py \
       -i ${fwd_hit}:${rev_hit} \
       --taxa_fp ${taxonomy_mapping} \
       -d ${taxa_depth} \
       -o ${pair_id} \
       -p 

"""    

}




process compile_taxa_coverage{

    //collects taxonomy coverage files into levels for comparison
    //will fail for some taxonomy levels
    errorStrategy 'ignore'
    publishDir path: "${output}/coverage_levels", mode: 'copy'
    input:
        file files from taxa_coverage_plots.collect()
	val taxa_depth

    output:
	file "*" into coverage
	file "taxonomy_reports/*hits_*_hits_coverage.txt" into  coverage_reports

 
"""

for level in `seq ${taxa_depth}` 

do
     mkdir -pv taxonomy_level_\${level}
     for  pdf  in  `ls  */*primers_coverage/*_level_\${level}.pdf`
         do
 	    pdf=`realpath \${pdf}`
 	    mv \${pdf}  taxonomy_level_\${level}
         done
     
done

mkdir -pv taxonomy_reports

mv -nv  */*/*.txt  taxonomy_reports

"""
}




process taxa_assignment{

    //errorStrategy 'ignore'
    publishDir path: "${output}/taxa_assignment_report", mode: 'copy'
    input:
	file seq from amplicons2
        val taxonomy_mapping
        val  ref_fasta
    output:
         file "*_amplicons_assignments.txt" into amplicon_assignments
	 file "*_amplicons_accuracy_report.txt"

    when:
       params.analyse_amplicons == true

    
"""	
    
    taxa_assignment_report.py \
	-t ${taxonomy_mapping} \
	-f ${seq} \
	-T ${ref_fasta} \
	-d 6 \
        -c 0.5
    

"""


}




process load_assignment_accuracy{

    input:
	file report from amplicon_assignments
	
    output:
         val 'success' into accuracy_reports


"""	

     load_taxa_report.py  localhost \
        -r ${report}

"""

}




process plot_accuracy{

    publishDir path: "${output}/assignment_accuracy", mode: 'copy'
    input:
	val accuracy from accuracy_reports.collect()
	
    output:
        file  "assignment_accuracy.*"

"""	

     source activate_env.sh
     plot_accuracy.py localhost


"""

}




process load_coverage{

    input:
	file coverage from coverage_reports.flatten()
	


"""	

  load_coverage.py localhost \
     -r ${coverage}

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
	 val data from physchem_data.collect()

    output:
       file '*.pdf' into propplots

    
"""    
   source activate_env.sh
   
   plot_physprop.py \
   ${hostname}

"""

}




process pair_tsv{

   
    input:
	file "*" from pair_tsv.collect()
	set fwd, rev from  primers_pair2

    output:
	file '*_*.tsv' into primersearch_tsv

    
"""
 
   get_pairs.py \
   $fwd \
   $rev 
   
 
"""

}




process emboss_primersearch{


    publishDir path: "${output}/primersearch", mode: 'copy'
    input:
	file primer_pair from primersearch_tsv
        val ref_fasta	

    output:
        file("${outfile}") into amplicon_results mode flatten

    script:
	outfile = primer_pair.getName().replace('tsv','rst') 
    
    
"""
    primersearch \
    -seqall ${ref_fasta} \
    -infile ${primer_pair} \
    -mismatchpercent 20 \
    -outfile ${outfile} \
    -verbose 

"""

}


process amplicon_stats{

   
    input:
        val hostname
        file amplicons from amplicon_results

    output:
        file amplicons into (amplicon_stats1, amplicon_stats2 ) mode flatten

"""

    amplicon_stats.py \
    ${hostname} \
    --amplicon $amplicons
      
"""

}




process  plot_ampliconData{

    
    publishDir path: "${output}/emboss_amplicons", mode: 'copy'
    input:
        val stats from amplicon_stats1.collect()
        val ref_fasta

    output:
         file '*' 
        

"""   
   source activate_env.sh

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
        val ref_fasta
    
    output:
        file "*.fasta" into amplicon_fasta mode flatten


"""

    seqret_amplicons.py \
    ${hostname} \
    -R ${ref_fasta}

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




process analyse_clusters{
  
    publishDir path: "${output}/analyse_clusters", mode: 'copy'  
    input:  
        file fasta_cluster from cdhit_clusters.collect()

    output:
        file  "CDhit_counts*" into hit_plots



""" 
 
    analyse_clusters.py
     

"""

}




process analyse_genetic_distances{

    errorStrategy 'ignore'
    publishDir path: "${output}/genetic_distances", mode: 'copy'
    input:  
	file seq_hits from cd_hits.collect()
    
    output:
	file "*.pdf" into distance_plots


"""    
    source activate_env.sh 

    get_distances.py \
    -t $params.clustalo_threads
     
"""

}
