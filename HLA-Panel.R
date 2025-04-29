##HLA-Probe Designing

##In terminal 
#Download human genoome annotatons and sequences (with haplotypes)

#curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf.gz
#&& gunzip gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf.gz



#curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.p13.genome.fa.gz && gunzip GRCh38.p13.genome.fa.gz

##get a tsv file containing all exons with you desired HLA-genelist

#cat gencode.v40.chr_patch_hapl_scaff.basic.annotation.gtf | grep -w 'HLA-DRA\|HLA-DRB1\|HLA-DRB3\|HLA-DRB4\|HLA-DRB5\|HLA-DQA1\|HLA-DQB1\|HLA-DPA1\|HLA-DPB1' > exon_all.tsv

#In R 
setwd('/Users/prashant/Documents/AMC/ProbePannel') ##Set path for a working directory:Replace with yours

targets=read.table('exon_all.tsv',header = T,sep = '\t') #Read in the tsv file

#extract useful information and rename columns
colnames(targets)=c("chr","database","feature","start","end","score",
                    "strand","frame","isoform")

targets=targets[targets$feature=="exon",]
targets$exon_num <- targets$isoform
targets$isoform <- gsub("^.*transcript_name ","",targets$isoform)
targets$isoform <- gsub("; exon_number.*$","",targets$isoform)
targets$exon_num <- gsub("^.*exon_number ","",targets$exon_num)
targets$exon_num <- gsub("; exon_id.*$","",targets$exon_num)
targets$exon_num <- as.numeric(targets$exon_num)

#fix the chromsome name and make sure they are consistent with chr names in 
#the fasta file
targets$chr <- gsub("chr","",targets$chr)

##Check the gene-names if all genes are captured-

targets$gene_name <- sub("^(([^-]+)-[^-]+)-.*", "\\1", targets$isoform)
levels(factor(targets$gene_name))

#take 1 probe / per kb, calculate how many probes we need to cover each exon
targets$length <- abs(targets$start - targets$end)
targets$probe_num <- (targets$length %/% 1000) + 1

#find the position of exons. 
#the first exon in a transcripts is called "start"
#the last exon in a transcripts is called "end"
#others are in the "middle"

targets$pos <- "middle"
for (i in unique(targets$isoform)) {
  x <- targets[targets$isoform == i,]
  x[1,]$pos <- "start"
  x[nrow(x),]$pos <- "end"
  targets[targets$isoform == i,] <- x
}


#For GRCH38.p13(included haplotypes etc) use-GRCh38.p13.genome.fa and load it in R-

library(Biostrings)

# Load the genome
genome <- readDNAStringSet('GRCh38.p13.genome.fa')

# Clean sequence names if needed (removes stuff after the first space)
names(genome) <- sub(" .*", "", names(genome))


# Get chromosome names from genome
genome_names <- names(genome)

# Create a mapping using grep ( match each chr name in targets)
chr_index_map <- sapply(unique(targets$chr), function(chr) {
  # Find genome sequence names containing the target chr
  match <- grep(chr, genome_names, value = FALSE, fixed = TRUE)
  if (length(match) > 0) {
    return(match[1])  # take the first match
  } else {
    return(NA)
  }
})

# Result: named vector of genome indices
names(chr_index_map) <- unique(targets$chr)

# Map each row's chr to its genome index
targets$chr_index <- chr_index_map[as.character(targets$chr)]

# Check how many unmatched
sum(is.na(targets$chr_index))  # should be small, ideally 0


##Function to extract sequences

extractSeq <- function(x) {
  chr <- as.integer(x["chr_index"])  
  
  if (is.na(chr)) return(NULL)
  
  start <- as.integer(x["start"])
  end   <- as.integer(x["end"])
  
  seq <- genome[[chr]][start:end]
  
  if (x["strand"] == "-") {
    seq <- reverseComplement(seq)
  }
  
  return(seq)
}


##Create exons.fa

if (file.exists("exons.fa")) file.remove("exons.fa")

for (i in 1:nrow(targets)) {
  x <- targets[i, ]
  seq <- extractSeq(x)
  if (!is.null(seq)) {
    seq <- DNAStringSet(seq)
    names(seq) <- paste0(x$isoform, "_exon_num_", x$exon_num)
    writeXStringSet(seq, "exons.fa", append = TRUE)
  }
}

### make up and down concatenations
origin_concatenation_pathway <- 'orig.fa'
up_concatenation_pathway <- 'up.fa'
down_concatenation_pathway <- 'dwn.fa'


#make the origin fasta (containing all exon >= 120bp) and up-concatenations
#(short exon concatenated to last exon)

#make the origin fasta (containing all exon >= 120bp) and up-concatenations (short exon concatenated to last exon)
for (i in unique(targets$isoform)) {
  x <- targets[targets$isoform == i,]
  x <- x[order(x$exon_num, decreasing = F),]
  seq.last01 <- DNAString("")
  seq.last02 <- DNAString("")
  for (j in 1:nrow(x)) {
    y <- x[j,]
    if(y$length < 120){
      seq.current <- extractSeq(y)
      seq.concatenation <- c(seq.last01, seq.current)
      
      if (length(seq.concatenation) < 120){ # if the up-concatenated exons is still < 120bp, 
        seq.concatenation <- c(seq.last02, seq.concatenation)
        # now, 3 exons are concatenated. write the concatenation in fasta file if it is >= 120bp
        if (length(seq.concatenation) >= 120){
          seq.write <- DNAStringSet(seq.concatenation)
          names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_up2", sep = "")
          writeXStringSet(seq.write, up_concatenation_pathway, append = T)
        } 
      } else {
        seq.write <- DNAStringSet(seq.concatenation)
        names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_up1", sep = "")
        writeXStringSet(seq.write, up_concatenation_pathway, append = T)
      }
      seq.last02 <- seq.last01
      seq.last01 <- seq.current
    } else {
      seq.o <- extractSeq(y)
      seq.write <- DNAStringSet(seq.o)
      names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_origin", sep = "")
      writeXStringSet(seq.write, origin_concatenation_pathway, append = T)
      seq.last02 <- seq.last01
      seq.last01 <- seq.o
    }
  }
}

#make down-concatenations (short exon concatenated to next exon)
for (i in unique(targets$isoform)) {
  x <- targets[targets$isoform == i,]
  x <- x[order(x$exon_num, decreasing = T),]
  seq.last01 <- DNAString("")
  seq.last02 <- DNAString("")
  for (j in 1:nrow(x)) {
    y <- x[j,]
    if(y$length < 120) {
      seq.current <- extractSeq(y)
      seq.concatenation <- c(seq.current,seq.last01)
      if(length(seq.concatenation) < 120){ # if the down-concatenated exons is still < 120bp,
        seq.concatenation <- c(seq.concatenation,seq.last02)
        if(length(seq.concatenation) >= 120){
          seq.write <- DNAStringSet(seq.concatenation)
          names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_down2", sep = "")
          writeXStringSet(seq.write, down_concatenation_pathway, append = T)
        } 
      } else {
        seq.write <- DNAStringSet(seq.concatenation)
        names(seq.write) <- paste(y$isoform, "_exon_num", y$exon_num, "_POS", y$pos, "_L", y$length, "_PN", y$probe_num, "_down1", sep = "")
        writeXStringSet(seq.write, down_concatenation_pathway, append = T)
      }
      seq.last02 <- seq.last01
      seq.last01 <-seq.current
    } else {
      seq.current <- extractSeq(y)
      seq.last02 <- seq.last01
      seq.last01 <- seq.current
    }
  }
}

#In terminal, merge all .fa files (i.e- up.fa, dwn.fa,orig.fa)
#cat *.fa > combined.fa
# put combined.fa to IDT website's design tool-
#https://www.idtdna.com/site/order/designtool/index/XGENDESIGN/XGEN.37F60E82163A47B9AF16A88E1EB1BDEB


##Read in IDT files (once you have recied files from them)

probe.seq <- readxl::read_excel('DO.xlsx')

#exon names
probe.seq$exon_name <- gsub("_POS.*$","",probe.seq$Chromosome)
#probe num, how many probes we need to cover this exon
probe.seq$probe_num <- gsub("^.*_PN","",probe.seq$Chromosome)
probe.seq$probe_num <- gsub("_down.*$","",probe.seq$probe_num)
probe.seq$probe_num <- gsub("_up.*$","",probe.seq$probe_num)
probe.seq$probe_num <- gsub("_origin.*$","",probe.seq$probe_num)
probe.seq$probe_num <- as.numeric(probe.seq$probe_num)

#concatenation
probe.seq$concatenation <- gsub("^.*down.*$","down",probe.seq$Chromosome)
probe.seq$concatenation <- gsub("^.*up.*$","up",probe.seq$concatenation)
probe.seq$concatenation <- gsub("^.*origin.*$","origin",probe.seq$concatenation)
#GC% similarity
probe.seq$drift <- abs(probe.seq$GC - mean(probe.seq$GC))
#length of exon
probe.seq$length <- gsub("^.*_L","",probe.seq$Chromosome)
probe.seq$length <- gsub("_PN.*$","",probe.seq$length)
probe.seq$length <- as.numeric(probe.seq$length)


targets$exon_name <- paste(targets$isoform, "_exon_num", targets$exon_num, 
                           sep = "")
probe.seq <- probe.seq[probe.seq$exon_name %in% targets$exon_name,]

if (sum(unique(targets$exon_name) %in% unique(probe.seq$exon_name)) / length(unique(targets$exon_name)) == 1) {
  print("All exons are covered")
} else {
  print("Following are exons not covered:")
  check_list <- unique(targets$exon_name)[!(unique(targets$exon_name) %in% unique(probe.seq$exon_name))]
  print(check_list)
}


probe.select <- data.frame()
probe.up <- dplyr::filter(probe.seq, concatenation == "up")
for (i in unique(probe.up$exon_name)) {
  x <- probe.up[probe.up$exon_name == i,]
  x <- x[order(x$Start),]
  probe.select <- rbind(probe.select,x[nrow(x),])
}

probe.down <- dplyr::filter(probe.seq, concatenation == "down")
for (i  in unique(probe.down$exon_name)) {
  x <- probe.down[probe.down$exon_name == i,]
  x <- x[order(x$Start),]
  probe.select <- rbind(probe.select,x[1,])
}

probe.origin <- dplyr::filter(probe.seq, concatenation == "origin")
for (i in unique(probe.origin$exon_name)) {
  x <- probe.origin[probe.origin$exon_name == i,]
  pn <- x[1,]$probe_num
  if(pn == 1){
    x <- x[order(x$drift),]
    probe.select <- rbind(probe.select, x[1,])
  } else if (pn > 1) {
    x <- x[order(x$Start),]
    len <- x$length[1]
    n <- len %/% 1000
    for (j in 0:n) {
      g1 <- j*1000
      g2 <- (j+1)*1000
      y <- x[x$Start>g1 & x$Stop<g2,]
      y <- y[order(y$drift),]
      probe.select <- rbind(probe.select, y[1,])
    }
  }
}

#remove duplicates
probe.select=probe.select[!duplicated(probe.select$Seq),]
#remove NA
probe.select=probe.select[!is.na(probe.select$Seq),]
dim(probe.select)

##Check if all the genes are covered
# Create 'gene_name' column by extracting before the second hyphen
probe.select$gene_name <- sub("^(([^-]+)-[^-]+)-.*", "\\1", probe.select$exon_name)
levels(factor(probe.select$gene_name))



#Get the probe_panel to order from IDT
write.csv(probe.select,file = 'HLA_panel.csv',row.names=T)


