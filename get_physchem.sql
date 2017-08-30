.headers on
.mode csv


.output 'Piscator_results/physchem1.csv'
select distinct primer_ID,
       Length,
       GC,
       Tm,
       TmProd from primerphyschem;
       
.output 'Piscator_results/physchem2.csv'

select distinct primer_ID,
       DeltaS,
       DeltaG,
       DeltaH from primerphyschem;


