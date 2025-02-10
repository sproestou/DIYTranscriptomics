##Author: Savana Hadjipanteli
##Date: February 5th-10th, 2025
##Course: DIY Transcriptomics

###Lab 4: Annotating gene expression data

#load modules
library(tidyverse) #for dataframe manipultation and ggplot
library(ensembldb) #to interface with ensembl
library(biomaRt) #to create marts and process

#WARNING: biomaRt/ensembldb mask certain dplyr functions, which then require additional path information to call. Ex. call dplyr::select() rather than select()
  #for this reason, I also unload these packages at the end of the script


##Task 1: 
#You’re a new grad student starting their first rotation in a viral pathogenesis lab. 
#A previous postdoc in the lab carried out RNA-seq on lung tissue collected from ferrets infected with influenza (ferrets are a great animal model for pathogenic respiratory viruses). 
#For your rotation, your PI asks you to analyze this dataset, and she is particularly interested in antiviral genes. 
#To begin this project, you must first find annotation data for ferrets (Mustela putorius furo). 
#To complete this task, you must use BiomaRt to locate this annotation data, and generate a dataframe that contains the following information:
#for each transcript in the ferret genome: transcript ID, start position, end position, gene name, gene description, entrez gene ID, and pfam domains.

#first pull ferret genome annotation and save to dataframe using biomart

#retrieve ferret gene annotaitons as a mart using:
#set mart to ENSEMBL_MART_ENSEMBL --> this is Ensembl Genes version 113 right now
#call dataset: mpfuro_gene_ensembl , which is: Ferret genes (MusPutFur1.0)
ferret <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl")

#this command just shows what is available in annotation, that way I can use the right names to call everything
ferret.atts <- listAttributes(ferret)

#pull data out of the mart (ferret, created above) and into a dataframe 

#must specify attributes to retrieve, which are the ones specified above, plus 
#note: I used transcript id version because it's more informative than just id, pfam is an id not a list of domains? but that matched the prompt description best
ferret.genes <- getBM(attributes=c('ensembl_transcript_id_version',
                                   'start_position', 
                                   'end_position',
                                   'external_gene_name',
                                   'description',
                                   'entrezgene_id',
                                   'pfam'),
                mart = ferret)


##Task 2:
# After showing your PI how you successfully used R to solve the annotation problem above, she gets excited and decides you might be able to help her with a related question. 
#Given her interest in antiviral genes, she asks if you could retreive the ferret promoter sequences (1kb upstream) for her 5 favorite antiviral genes, IFIT2, OAS2, IRF1, IFNAR1, and MX1. 
#She hopes to use these sequences to engineer some reporter constructs. To get started on this task, you may want to use RStudio to access the help documentation for the getSequence function that is part of the BiomaRt package.

#set variable with how much sequence upstream from the transcript we want
dist = 1000

ferret_promoter_genes <- c("IFIT2", "OAS2", "IRF1", "IFNAR1", "MX1")

#first filter transcript information for each gene to retrieve transcript names with different start sites for each gene
#since these appear to have repeats based on multiple pfam domains, I also set it to collapse these lists given the external gene name and transcript start position are the same...
ferret.select.genes <- ferret.genes %>% dplyr::select(-pfam) %>% dplyr::filter(external_gene_name %in% ferret_promoter_genes) %>% distinct(external_gene_name, start_position, .keep_all = TRUE) 

#then use getSequence to generate a list of sequences...
#check listFilters(ferret) for possible type arguments in getSequence

#retrieve sequences

#here I use specific transcript id version to call sequence location, since sometimes a single gene can have multiple transcripts with different start sites, I filtered for distinct transcript start coordinates and gene names (which in this case they all just had the same start coordinate)
#seqType is transcript_flank so that it just returns the regions up/downstream of transcript (not including UTR, which should be included in the marked start coordinate for the transcript) --> if we just wanted the nucleotide sequence we would call "cdna"
getSequence(id = ferret.select.genes[["ensembl_transcript_id_version"]], 
            type = "ensembl_transcript_id_version", 
            seqType = "transcript_flank", 
            upstream = dist,
            mart = ferret)

#this returns a single sequence of 1000 bp of nucleotide sequences upstream of the transcirption start site

##########Alternative attempt, but the results were not consistent 
#this just uses gene names, gets upstream transcript flanks (not including UTR)

#getSequence(id = c("IFIT2", "OAS2", "IRF1", "IFNAR1", "MX1"), 
#                  type = "external_gene_name", 
#                  seqType = "transcript_flank", 
#                  upstream = 1000,
#                 mart = ferret)
#for whatever reason this returns multiple sequences for MX1, which have different sequences, which doesn't make sense except that there are multiple transcripts with different start sites, but there are not, at least according to our annotation...may be an issue with the list they have of genes but I don't really want to go on a manhunt for the specific sequences we're pulling
#and using hgnc_symbol excludes MX1...?
##########

#unload packages
detach("package:biomaRt", unload=TRUE)
detach("package:ensembldb", unload=TRUE)
