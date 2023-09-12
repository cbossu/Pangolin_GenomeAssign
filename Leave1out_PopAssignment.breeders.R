
library(tidyverse)
library(rubias)
library(dplyr)
library(xlsx)

## Population assignment using rubias- leave one out of 88 steps: need to automate it

#Rubias is an R software program that probabilistically assigns individuals to specific reporting units (i.e. genetic clusters diagnosed in the breeding region). Originally used for fish stock identification, we've co-opted it for bird genoscapes! I've attached the two rubias files you'll need. The known breeder genotypes and the genotypes of the migrating birds. Note, since we don't know where the migrating birds are coming from, we consider them unknown and we are going to use rubias to estimate mixing proportions and probability of assignment to their breeding origin. Thus, we consider the unknown birds a potential mixture of individuals from distinct genetic clusters.


# get command line arguments.
comarg     <- commandArgs()
input   <- comarg[6]
#input <-1
TEST=paste("test",input,sep="")

# load the data
#datad<-read_delim("../test88/test1.SNPsubset.recode.vcf.ped.head",delim="\t") %>% rename(indiv=IID) %>% dplyr::select(-FID,-X3,-X4,-X5,-X6)
datad <- read_delim(paste(TEST,"SNPsubset100.recode.vcf.ped.head",sep="."),delim="\t")  %>% 
  dplyr::select(-FID,-X3,-X4,-X5,-X6)

t<-read.xlsx("Pangolin_tests_assertainment.xlsx",sheetIndex = 1) %>% filter(test==TEST)
meta<-read.xlsx("../../meta_data/WB Pangolin WGS Libraries_metadata_CG_09.20.2021.xlsx",sheetIndex = 1) %>% rename(indiv=Novogene_ID,collection=Location) %>%  dplyr::select(repunit,collection,indiv)

data<-datad %>% left_join(meta) %>% dplyr::select(repunit,collection,indiv,everything()) %>% filter(indiv!="34501LE")



mix<-data %>% filter(indiv %in% t$ind_removed) %>%  mutate(sample_type="mixture",repunit2="NA") %>% dplyr::select(-repunit) %>% rename(repunit=repunit2) %>%  dplyr::select(sample_type,repunit,everything())
str(mix)
mix[,-(1:4)][mix[,-(1:4)] == 0] <- NA
mix[,-(1:4)][mix[,-(1:4)] == "A"] <- "1"
mix[,-(1:4)][mix[,-(1:4)] == "C"] <- "2"
mix[,-(1:4)][mix[,-(1:4)] == "G"] <- "3"
mix[,-(1:4)][mix[,-(1:4)] == "T"] <- "4"

mix %>% write.table(paste(TEST,"mix.rubias_input.txt",sep="."),row.names = F,col.names = T,quote=F,sep="\t")


ref<-data %>% filter(indiv!= t$ind_removed) %>% 
  mutate(sample_type="reference") %>% dplyr::select(sample_type,repunit,everything())
str(ref)
ref[,-(1:4)][ref[,-(1:4)] == 0] <- NA
ref[,-(1:4)][ref[,-(1:4)] == "A"] <- "1"
ref[,-(1:4)][ref[,-(1:4)] == "C"] <- "2"
ref[,-(1:4)][ref[,-(1:4)] == "G"] <- "3"
ref[,-(1:4)][ref[,-(1:4)] == "T"] <- "4"

ref %>% write.table(paste(TEST,"reference.rubias_input.txt",sep="."),row.names = F,col.names = T,quote=F,sep="\t")

##Read in mixture file and remove any columns with NA in here and in reference
pang_mix<-read_delim(paste(TEST,"mix.rubias_input.txt",sep="."),delim="\t") %>% distinct()

na_col<-names(which(colSums(is.na(pang_mix[,-(1:4)]))>0)) 
write.table(na_col,paste(TEST,"NA_columns_removed.txt",sep="."),row.names = F,quote=F,sep="\t")
pang_mix_final <- pang_mix %>% dplyr::select(-c(na_col)) %>% mutate_if(is.numeric, as.character)



#Read in the reference file, and make sure all fields are correct type. 

pang_ref <- read_delim(paste(TEST,"reference.rubias_input.txt",sep="."),delim="\t") %>% distinct()
pang_ref %>% dplyr::select(repunit) %>% distinct()

pang_ref_final <- pang_ref %>% dplyr::select(-c(na_col)) %>% mutate_if(is.numeric, as.character)

#Run rubias
mix_estC <- infer_mixture(reference = pang_ref_final, mixture = pang_mix_final, gen_start_col = 5)


#The result comes back as a list of four tidy data frames:

#1. mixing_proportions: the mixing proportions. The column pi holds the estimated mixing proportion for each collection.

#2. indiv_posteriors: this holds, for each individual, the posterior means of group membership in each collection. Column PofZ holds those values. Column log_likelihood holds the log of the probability of the individuals genotype given it is from the collection. Also included are n_non_miss_loci and n_miss_loci which are the number of observed loci and the number of missing loci at the individual. A list column missing_loci contains vectors with the indices (and the names) of the loci that are missing in that individual. It also includes a column z_score which can be used to diagnose fish that don’t belong to any samples in the reference data base (see below).

#3. mix_prop_traces: MCMC traces of the mixing proportions for each collection. You will use these if you want to make density estimates of the posterior distribution of the mixing proportions or if you want to compute credible intervals.

#4. bootstrapped_proportions: This is NULL in the above example, but if we had chosen method = "PB" then this would be a tibble of bootstrap-corrected reporting unit mixing proportions.

#What we are most interested in is the individual posteriors. We can look at the number of individuals assigned to different genetic clusters and associate those individuals to the meta data we have on them: location (MAPS site) and timing (Month,Day,Year).


options(scipen = 999)

#Take top posterior probability
rep_indiv_estsC <- mix_estC$indiv_posteriors %>% group_by(mixture_collection, indiv, repunit) %>% summarise(rep_pofz = sum(PofZ)) %>% 
  dplyr::select(mixture_collection,indiv,repunit,rep_pofz, everything()) %>% arrange(-rep_pofz) %>% slice(1)

#print assignment
rep_indiv_estsC %>% write.table(paste(TEST,"repunit_pofz.s100.txt",sep="."),row.names=F,quote=F,sep="\t")
rep_indiv_estsC 


###Write out a snp table for supplement

map96b<-read_delim("~/Dropbox/BGP/Pangolin/Fluidigm_genotypes/old_files/Pangolin.inclHK1_8.geno_LCWG.FinalPanel.96.2mind.map",delim="\t",col_names = F)

seq<-read.xlsx("~/Dropbox/BGP/Pangolin/FST/design_assays/assays_submitted/Pangolin.Fluidigm_SNP_Type_Assays_Target_Import_NoRs_Template.chosen.xlsx",sheetIndex = 1) %>% rename(Assay_name=Fluidigm.SNP.Type.Assays.Target.Import.by.Sequence.STS1,Sequence= `NA.`) %>% filter(Assay_name %in% map96b$X1) %>% dplyr::select(Assay_name,Sequence)
str(seq)

seq_info<-read.xlsx("~/Dropbox/BGP/Pangolin/FST/design_assays/Pangolin_choice1_3.126.FST_rank.xlsx",sheetIndex = 1) %>% rename(Sequence=Seq)

seq %>% left_join(seq_info) %>% dplyr::select(Assay_name,LeftFlank,RightFlank,GC_content,BasesPresent,IUPAC,CHROM,Sequence) %>% write.table("~/Dropbox/BGP/Pangolin/manuscript/Revisions/Pangolin.TableS1.Assay_info.txt",row.names = F,quote=F,sep="\t")
