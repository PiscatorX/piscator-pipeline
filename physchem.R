#!/usr/bin/Rscript

library(ggplot2)
library(reshape)
library(ggrepel)
library(optparse)
library(optparse)
library(stringi)


option_list = list(
make_option(c("-c", "--cwd"), type="character", default=NULL,
               help="dataset file name", metavar="cwd"),
make_option(c("-o", "--output"), type="character", default=NULL,
               help="dataset file name", metavar="out"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



setwd(opt$cwd)

output_dir = (opt$out)

print(getwd())

db_path = file.path(output_dir, "Piscator.db")

run_sql = paste("sqlite3 ",db_path," < get_physchem.sql")

system(run_sql)  

data_files = file.path(output_dir, c("physchem1.csv", "physchem2.csv"))

plot_files = gsub('csv','pdf',data_files)

n = length(data_files)
csv_data = list(n)
csv_df   = list(n)

for (i in 1:n) {
  csv_data[[i]] <- melt(read.csv(data_files[i]))
  ggplot(csv_data[[i]], aes(x=primer_ID, y = value, fill = primer_ID )) +
        geom_bar(stat='identity', width=0.9, position = "dodge") +
        xlab("Primer") +
        ylab(NULL) +
        facet_grid(variable ~ . , scales="free")  +
	theme_bw() +
        theme(axis.text.x = element_text(angle=45, h=1)) +
        scale_fill_discrete(name = "Primer") +
        geom_text_repel(aes(x=primer_ID,y=value,label=value), 
                size = 3,
                color = 'black',
                box.padding = unit(0, "lines"),
                point.padding = unit(0, "lines"),
                segment.color = 'grey50')
  
  ggsave(plot_files[[i]])      
  
  }
        
  

