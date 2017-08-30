.separator "\t"
.output 'amplicon_data.tsv'
select * from amplicons order by primer_id asc;
