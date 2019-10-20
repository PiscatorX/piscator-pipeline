

rm(list=ls())
cat("\f")
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(seqinr)
library(ape)
library(ips)
library(HDMD)
library(grid)
library(zoo)


primers <- read.csv("/home/drewx/Documents/piscator-pipeline/primers.csv4", strip.white = TRUE)
primer_blast <-  read.csv("/home/drewx/Documents/piscator-pipeline/primer.map/primer_blast/all_primer.tsv", strip.white = TRUE)

############################################################# PrimerBlast >>  primers ######################################################################################

primer_blast <- primer_blast %>%  
                dplyr::select(primer, sbjct_start, sbjct_end, score) %>%
                distinct()  %>%
                arrange(primer)
            
primers$loc <- 1:nrow(primers)  

primer_blast$loc <- NA
    
for (i in  1:nrow(primer_blast)){
        blast_record <- primer_blast[i,]
        query_id  <-blast_record$primer
        query_lc <- tolower( gsub("[+*]", "", query_id))
        for (j in 1:nrow(primers)){
            loc <- primers[j,]$loc
            primer_pair  <- primers[j,]
            primer_rev <-  gsub("[+*]", "", primer_pair$Rev_id)
            primer_rev_lc <- tolower(primer_rev)
            
            if (str_detect(query_lc, primer_rev_lc)){
                  #print(c(query_lc, primer_rev_lc))
                  primer_blast[i,]$loc <- loc
            }
            
            primer_fwd <-  gsub("[+*]", "", primer_pair$Fwd_id)
            primer_fwd_lc <- tolower(primer_fwd)
    
            if (str_detect(query_lc, primer_fwd_lc)){
                primer_blast[i,]$loc <- loc
            }
  }
}


##################################################################### Alignment ###################################################################


silva_aln <- read.alignment("/opt/DB_REF/SILVA/Align/SILVA_map.fasta", format  = "fasta")
silva_clean <- del.colgapsonly(silva_aln)

dim(silva_clean)

write.dna(silva_clean,file = "/opt/DB_REF/SILVA/Align/SILVA_ungapped.fasta", format = "fasta")

ref_seq <- read.FASTA("/home/drewx/Documents/piscator-pipeline/Z75578.fasta")

mapped_aln <- mafft(silva_clean, 
                    ref_seq,  
                    exec = "/usr/bin/mafft",
                    thread = 4)

write.dna(mapped_aln, file = "/opt/DB_REF/SILVA/Align/SILVA_mafft.fasta", format = "fasta")

degap_map <- function(alignment, ref_loc){
  
        print(alignment)
        print(dim(alignment))
        ref_seq <- as.character.DNAbin(alignment[ref_loc,])
        gaps <- (ref_seq != "-")
        dgpd_aln <- alignment[,gaps]
        print(dgpd_aln)
        print(dim(dgpd_aln))
        return(dgpd_aln)
  
}

dgpd_aln <- degap_map(mapped_aln,  4560)

write.FASTA(dgpd_aln, file = "/opt/DB_REF/SILVA/Align/SILVA_mafft_final.fasta")

aln_matrix  <- toupper(as.character.DNAbin(dgpd_aln))

Entropy <- MolecularEntropy(aln_matrix[-c(4560),], "DNA")

aln_entropy <- data.frame(position = 1:length(Entropy$H), entropy = Entropy$H)

aln_entropy$rollmeans <- rollmean(aln_entropy$entropy, 10, fill = 0.5 )


#################################################taxa_coverage#######################################################

taxa_coverage = read.csv("/home/drewx/Documents/piscator-pipeline/primer.map/primer_blast/taxa_coverage.tsv", sep = "\t",  stringsAsFactors = FALSE)

taxa_coverage$percent_passing <- round(taxa_coverage$percent_passing *100, 2)

primer_blast$cover <- NA

for (i in  1:nrow(primer_blast)){
  blast_record <- primer_blast[i,]
  query_id  <-blast_record$primer
  query_lc <- tolower( gsub("[+*]", "", query_id))
  for (j in 1:nrow(taxa_coverage)){
     cover <- taxa_coverage[j,]$percent_passing
     primer_pair  <- unlist(strsplit(taxa_coverage[j,]$primer_pair, "_"))
     primer_rev <-  gsub("[+*]", "", primer_pair[2])
     primer_rev_lc <- tolower(primer_rev)
     if (str_detect(query_lc, primer_rev_lc)){
        primer_blast[i,]$cover <- cover
     }
     primer_fwd <-  gsub("[+*]", "", primer_pair[1])
     primer_fwd_lc <- tolower(primer_fwd)
     if (str_detect(query_lc, primer_fwd_lc)){
       primer_blast[i,]$cover <- cover
     }
  }
}

primer_blast$loc <- primer_blast$loc + 1

normalize_vars <- function(data_vector, my_min = 0,  my_max = 1,  min_var = NA,  max_var = NA){
  
    print(c(my_min, my_max, min_var, max_var))
  
    if (is.na(min_var)){
        min_var <- min(data_vector)  
    } 
    if (is.na(max_var)){
        max_var <- max(data_vector)
    }
    print(c(my_min, my_max, min_var, max_var))
    a = (my_max - my_min)/(max_var - min_var)
    b = (my_max - (a * max_var))
    c = (my_min - (a * min_var))
    
    return(a * data_vector + b)
}


# length  0.5  1.5
# size    0.5  1
primer_blast$len <- normalize_vars(primer_blast$cover,  0.5, 1.25)
primer_blast$size  <- normalize_vars(primer_blast$cover,  0.5, 1)
################################################## Variable regions ###############################################################

z75578 <- read.csv("/home/drewx/Documents/piscator-pipeline/z75578_positions.csv")

z75578$label <- (z75578$Start_position + z75578$End_position)/2


################################################ combine plots ##################################################################

primer_blast[20,] <- primer_blast[primer_blast$primer == '1510r',][1,]
primer_blast[20,]$loc <- primer_blast[primer_blast$primer == '1380f',]$loc
primer_blast[20,]$cover <- taxa_coverage[taxa_coverage$primer_pair == "1380f_1510r",]$percent_passing

#   scale_color_manual(name = "slidding_mean",
#                      values = c("(-Inf,0.1]" = "red", "(0.1,0.5]" = "black","(0.5, Inf]" = "blue"))
# 


p1 <-ggplot(aln_entropy) + 
  geom_point(aes(x = position, 
                 y = entropy,
                 colour = cut(rollmeans, 
                              c(-Inf, 0.1, Inf))),
             size = 0.75) +
  theme_bw() + 
  theme(axis.title.y = element_text(size=12, 
                                    face="bold",
                                    margin = unit(c(0,0,0,-3), "cm")),
        axis.text=element_text(size=12),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin=unit(c(0,0,0,0), "cm")) +
  ylab("Normalised Entropy") +
  xlab("Position") +
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0.05)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits = c(0, 1900),
                     breaks = z75578$label,
                     labels = z75578$Region,
                     position = "top")


p1




p2 <- ggplot(primer_blast, aes(color = cover  )) + 
          geom_segment(aes(x = sbjct_start,  y =  loc, xend = sbjct_end, yend = loc),
                       arrow=arrow(angle = 30, 
                                   length = unit(primer_blast$len, "mm"), 
                                   type = "closed"), 
                       size = primer_blast$size) +
          geom_text_repel(aes(x = sbjct_start,
                              y = loc,  
                              label = primer), 
                              colour = "black",
                              nudge_y = 0.25,
                              size = 4) +
          theme_bw() +
          theme(panel.background = element_rect(fill = "transparent"),
                axis.title.x = element_text(size=12,
                                            face="bold",
                                            margin = unit(c(0.25,0,0,0), "cm")),
                legend.position="bottom",
                legend.text = element_text(size=10),
                legend.title = element_text(size = 12, face = "bold", vjust = 0.85 ),
                panel.grid.minor = element_blank(),
                axis.text=element_text(size=12),
                axis.ticks.y = element_blank(),
                panel.grid.major = element_blank(),
                axis.text.y = element_blank(),#t, r, b, l 
                plot.margin=unit(c(0,0,0,2), "cm")) +
          xlab("Position (bp)") +
          ylab(NULL) +
          scale_y_continuous(expand = c(0, 0), 
                             breaks = seq(2,10,2), 
                             limits=c(1, 12)) +
          scale_x_continuous(expand = c(0, 0), 
                             breaks = seq(0,1800,200), 
                             limits=c(0, 1900)) +
          guides(color = guide_colourbar(barwidth = 10, barheight = 1)) + 
          scale_colour_gradient(name = "Coverage (%)", low = "dark green",  high = "red",
                                limits = c(25,100))




  
p2


setwd("~/Documents/piscator-pipeline/")

grid.newpage()
pdf("primer_map.pdf", width = 7, height = 7)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

dev.off()

