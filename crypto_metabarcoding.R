
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
# BiocManager::install("phyloseq")

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library('plyr')
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
# library(Manu)
library(RColorBrewer)
library(viridis)
library(scales)

setwd("C:/Users/dtshayma/OneDrive - Massey University/Documents/Massey_2025/crypto_code_metabarcoding")
otu_mat<- read_excel("C:/Users/dtshayma/OneDrive - Massey University/Documents/Massey_2025/crypto_code_metabarcoding/crypto_metabarcoding.xlsx", sheet = "otuAbundance")
tax_mat<- read_excel("C:/Users/dtshayma/OneDrive - Massey University/Documents/Massey_2025/crypto_code_metabarcoding/crypto_metabarcoding.xlsx", sheet = "otuTaxonomy")
samples_df<- read_excel("C:/Users/dtshayma/OneDrive - Massey University/Documents/Massey_2025/crypto_code_metabarcoding/crypto_metabarcoding.xlsx", sheet = "sampleMeta")

# define row names

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") 

tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample")

# transform into matrices and tax tables

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

crypto <- phyloseq(OTU, TAX, samples)
crypto

# Only include crypto sequences
crypto <- subset_taxa(crypto, Family %in% c("Cryptosporidiidae"))
crypto

# check data

sample_names(crypto)

rank_names(crypto)

sample_variables(crypto)

# add a colour palette for people with deuteranopia

cblind <- (c("#1D91C0", "#624B27", "#CB181D", "#F46D43", "#FAE093", "#A6CEE3", "#74C476", "#EF3B2C", "#000000", "#004949", "#009292", "#F9E211", "#BA2F00", "#425266", "#3D4928", "#008CEC"))

# change the plotting here for the top 50, but this time, can plot the data based on Subtype (or gene family)
neworder = c("Auckland_2010","Christchurch_2010","Hawke's Bay_2013","Waikato_2013","Wellington_2013",
             "Taranaki_2013","Auckland_2015","Auckland_2017","Blenheim_2017","Wellington_2018",
             "Routine Surveillance_2015","Routine Surveillance_2017","Routine Surveillance_2018")
sample_data(crypto)$Outbreak<-factor(sample_data(crypto)$Outbreak, levels = neworder)

# Rename levels
sample_data(crypto)$Outbreak <- revalue(sample_data(crypto)$Outbreak, c("Auckland_2010"="Auckland 2010",
                                                                        "Christchurch_2010"="Christchurch 2010",
                                                                        "Hawke's Bay_2013"="Hawke's Bay 2013",
                                                                        "Waikato_2013"="Waikato 2013",
                                                                        "Wellington_2013"="Wellington 2013",
                                                                        "Taranaki_2013"="Taranaki 2013",
                                                                        "Auckland_2015"="Auckland 2015",
                                                                        "Auckland_2017"="Auckland 2017",
                                                                        "Blenheim_2017"="Blenheim 2017",
                                                                        "Wellington_2018"="Wellington 2018",
                                                                        "Routine Surveillance_2015"="Routine Surveillance 2015",
                                                                        "Routine Surveillance_2017"="Routine Surveillance 2017",
                                                                        "Routine Surveillance_2018"="Routine Surveillance 2018"))

# View modified levels
levels(sample_data(crypto)$Outbreak)

# levels of subtypes
length(unique(tax_table(crypto)[,"Subtype"]))

# omit sequences with zero reads
# NB - none
cryptoNoZero <- prune_taxa(taxa_sums(crypto)>=1, crypto)
cryptoNoZero

# plotting diversity by outbreak using the Shannon and Simpson measures

richness<-plot_richness(cryptoNoZero, x = "NID", measures=c("Shannon", "Simpson"), color = "Outbreak") +
  scale_colour_manual(values = cblind) +
  labs(x = "Sample")
richness

ggsave("Suppl_3_richness_comparison.png", richness, 
        height=6, width=10)

richness_facet<-plot_richness(cryptoNoZero, x = "LibraryID", measures=c("Shannon"), color = "Outbreak") + scale_colour_manual(values = cblind) +
 facet_grid(~Outbreak,scales="free_x", space = "free_x")+
  theme(axis.text.x = element_blank())+
  theme(strip.text = element_blank()) +
  labs(x = "Sample")
richness_facet
 
ggsave("Suppl_4_richness_outbreaks.png", richness_facet, 
        height=6, width=10)

# a quick look at the data via an ordination method

ord.nmds.bray <- ordinate(cryptoNoZero, method="NMDS", distance="bray")

plot_ordination(cryptoNoZero, ord.nmds.bray, color="Outbreak", title="Bray NMDS") + scale_colour_manual(values = cblind) #+ theme_bw()

 # what about the top 50 taxa?
 
#  top50 <- names(sort(taxa_sums(cryptoNoZero), decreasing=TRUE))[1:50]
#  head(top50)
# # 
#  crypto.top50 <- transform_sample_counts(cryptoNoZero, function(OTU) OTU/sum(OTU))
#  crypto.top50 <- prune_taxa(top50, crypto.top50)
# # 
#  crypto.top50AG <- tax_glom(crypto.top50, "Subtype")

 top300 <- names(sort(taxa_sums(cryptoNoZero), decreasing=TRUE))[1:300]

 head(top300)

crypto.top300 <- transform_sample_counts(cryptoNoZero, function(OTU) OTU/sum(OTU))
crypto.top300 <- prune_taxa(top300, crypto.top300)

crypto.top300AG <- tax_glom(crypto.top300, "Subtype")

cblind7 <- (c("#74C476", "#009292", "#F9E211", "#BA2F00", "#425266", "#3D4928", "#008CEC"))

comp_plot_300 <- plot_bar(crypto.top300AG, x="Sanger_type", fill="Subtype_family") + 
  facet_wrap(~Outbreak, ncol = 3) + 
  scale_fill_manual(values = cblind) +
  labs(fill = "NGS\nsubtype\nfamily", x = "Sanger subtype family") +
  scale_x_discrete(labels=c("cuniculus" = "V", 
                            "erinacei" = "XIII",
                            "tyzzeri" = "IX"))+
  theme(axis.title.y = element_blank())

comp_plot_300

ggsave("Fig_S2_sanger_ngs_300.png", comp_plot_300, 
       height=8, width=10)

comp_plot_300_all_no_facet_cols<-plot_bar(crypto.top300AG, x = "Sanger_type", fill = "Subtype_family") +
  # facet_wrap(~Outbreak, ncol = 3) +
  scale_fill_manual(values = cblind) +
  labs(fill = "NGS\nsubtype\nfamily", x = "Sanger subtype family") +
  scale_x_discrete(labels = c("cuniculus" = "V", 
                              "erinacei" = "XIII",
                              "tyzzeri" = "IX")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(legend.text=element_text(size=16),legend.title=element_text(size=20)) +
  geom_bar(stat = "identity", position = "stack", color = NA)
comp_plot_300_all_no_facet_cols

ggsave("Fig_1_sanger_ngs_300_no_facet_cols_R1.png", comp_plot_300_all_no_facet_cols, 
       height=8, width=10)

 frig <- sample_names(crypto.top50)

length(unique(tax_table(crypto.top50)[,"Subtype"]))
 print(unique(tax_table(crypto.top50)[,"Subtype"]))
 
 neworder_st = c("IbA8G2",    
              "IbA9G2",    
              "IbA10G2",   
              "IbA11G2",   
              "IfA11G1",   
              "IfA12G1",   
              "IgA9",
              "IgA10",     
              "IgA11",     
              "IgA12",     
              "IgA13",     
              "IgA14",     
              "IgA15",     
              "IgA16",  
              "IgA17",
              "IgA18",     
              "IgA19",     
              "IgA20", 
              "IgA21",     
              "IIaA13G1R1",
              "IIaA14G1R1",
              "IIaA14G2R1",
              "IIaA15G2R1",
              "IIaA16G3R1",
              "IIaA17G3R1",
              "IIaA18G3R1",
              "IIaA18G4R1",
              "IIaA19G3R1",
              "IIaA19G4R1",
              "IIdA16G1",
              "IIdA17G1",  
              "IIdA20G1",  
              "IIdA23G1",  
              "IIdA24G1",  
              "VbA23",
              "VbA24",     
              "VbA25",     
              "IXbA5",     
              "IXbA6")
 length(neworder_st)
 
  pH_outbreak_wrap <- plot_heatmap(crypto.top50, method = "NMDS", distance = "bray", sample.label = "Outbreak", 
                                          sample.order = frig, taxa.label = "Subtype", taxa.order = "Subtype")+
             facet_grid(.~Outbreak, scales="free_x", space = "free_x", switch="both") +
   theme(axis.text.x = element_blank())+
   theme(strip.text = element_text(angle = 90))+
   theme(axis.ticks.length=unit(0.05,"cm"))+
   labs(fill = "Relative\nabundance")
 
 pH_outbreak_wrap
 
 ggsave("Suppl_S6_all_raw_50.png", pH_outbreak_wrap, 
        height=6.5, width=10)
 
my_col<-c(rep("#425266",3),rep("#3D4928",28),rep("#425266",16),rep("#3D4928",3))

pH_outbreak_trans_wrap_1 <- plot_heatmap(crypto.top50, 
                                         method = "NMDS", 
                                         distance = "bray", 
                                         sample.label = "Outbreak",
                                         
                                         sample.order = frig, 
                                         taxa.label = "Subtype", 
                                         taxa.order = "Subtype",
                                         trans = identity_trans()) + facet_grid(.~Outbreak, scales="free_x", space = "free_x", 
                                                                                switch="both") +
  theme(axis.text.y = element_text(color = my_col))+
  theme(axis.text.x = element_blank())+
  theme(strip.text = element_text(angle = 90))+
  theme(axis.ticks.length=unit(0.05,"cm"))+
  labs(fill = "Relative\nabundance")

pH_outbreak_trans_wrap_1

ggsave("Suppl_S7_all_scaled_col_50.png", pH_outbreak_trans_wrap_1, 
       height=6.5, width=10)

##

# BiocManager::install("microbiome")
library(microbiome)

s <- aggregate_taxa(crypto.top50, 'Subtype')

s


my_col<-c(rep("grey35",19),
          rep("black",15),
          rep("grey35",3),
          rep("black",3))


# Use neworder_st for the taxa order
pH_outbreak_s_bw <- plot_heatmap(
  s,
  method = "NMDS",
  distance = "bray",
  sample.label = "Outbreak",
  sample.order = frig,
  taxa.label = "Subtype",
  taxa.order = neworder_st,  # Apply the custom order for Subtype
  trans = identity_trans(),
  low = "white",
  high = "#000033"
) + 
  facet_grid(. ~ Outbreak, scales = "free_x", space = "free_x", switch = "both") +
  theme(axis.text.y = element_text(color = my_col)) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text = element_text(angle = 90)) +
  theme(axis.ticks.length = unit(0.05, "cm")) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 19.5) +
  geom_hline(yintercept = 34.5) +
  geom_hline(yintercept = 37.5) +
  geom_hline(yintercept = 39.5) +
  labs(fill = "Relative\nabundance",
       subtitle = "B")

pH_outbreak_s_bw

ggsave("Fig_2b_agg_scaled_bw_top50_R1.png", pH_outbreak_s_bw, 
       height=6.5, width=10)

##
length(unique(tax_table(crypto.top300AG)[,"Subtype"]))
print(unique(tax_table(crypto.top300AG)[,"Subtype"]))

neworder_300 = c("IbA8G2",    
                 "IbA9G2",     
                 "IbA9G3",   
                 "IbA10G1",
                 "IbA10G2",
                 "IbA11G1",  
                 "IbA11G2",  
                 "IbA12G2", 
                 "IfA8G1",
                 "IfA9G1",  
                 "IfA11G1",    
                 "IfA12G1",    
                 "IfA13G1",
                 "IfA14G1",
                 "IgA6",  
                 "IgA7",
                 "IgA8",  
                 "IgA9",     
                 "IgA10",
                 "IgA11",   
                 "IgA12",
                 "IgA13",     
                 "IgA14",   
                 "IgA15",     
                 "IgA16",   
                 "IgA16G1",    
                 "IgA17", 
                 "IgA18",   
                 "IgA19",    
                 "IgA19G1",
                 "IgA20",     
                 "IgA20G1",   
                 "IgA21",  
                 "IgA22",      
                 "IgA23",  
                 "IgA24",
                 "IgA25",
                 "IgA26",  # 38
                 "IIaA6G3R1",
                 "IIaA7G3R1",
                 "IIaA8G3R1",
                 "IIaA9G3R1",
                 "IIaA10G1R1",
                 "IIaA11G3R1",
                 "IIaA12R1",
                 "IIaA12G1R1", 
                 "IIaA12G2R1", 
                 "IIaA12G3R1",
                 "IIaA13G1R1",
                 "IIaA13G2R1",
                 "IIaA13G3R1",
                 "IIaA14G1R1",
                 "IIaA14G2R1",
                 "IIaA14G3R1",
                 "IIaA15G2R1", 
                 "IIaA15G3R1",
                 "IIaA16",    
                 "IIaA16G1R1", 
                 "IIaA16G2R1", 
                 "IIaA16G3R1",
                 "IIaA16G4R1", 
                 "IIaA17G2",
                 "IIaA17G2R1", 
                 "IIaA17G3R1",
                 "IIaA17G4R1",
                 "IIaA18G2R1", 
                 "IIaA18G3R1",
                 "IIaA18G4R1",
                 "IIaA18G5R1",
                 "IIaA19G2R1",
                 "IIaA19G3R1",
                 "IIaA19G4R1",
                 "IIaA20G3",   
                 "IIaA20G3R1", 
                 "IIaA20G4R1",
                 "IIaA21G4",   
                 "IIdA12",
                 "IIdA13",  
                 "IIdA16G1",  
                 "IIdA17G1",
                 "IIdA18G1",  
                 "IIdA19G1",  
                 "IIdA20G1",  
                 "IIdA21G1",   
                 "IIdA22G1",  
                 "IIdA22G3",   
                 "IIdA23G1",  
                 "IIdA23G2",
                 "IIdA24G1",
                 "IIdA24G2",
                 "IIdA25G1",  
                 "IIdA26G1",
                 "IIeA7G1",   # 56
                 "VbA17",
                 "VbA21", 
                 "VbA22",  
                 "VbA23",  
                 "VbA24",     
                 "VbA25",
                 "VbA26", # 7  
                 "IXbA5",
                 "IXbA6", # 2
                 "XIIIaA8R10",
                 "XIIIaA9R10", 
                 "XIIIaA10R10",
                 "XIIIaA11R10",
                 "XIIIaA12R10",
                 "XIIIaA13R10",
                 "XIIIaA14R5", 
                 "XIIIaA14R6",
                 "XIIIaA14R10",
                 "XIIIaA15R10",
                 "XIIIaA16R10") # 11
length(neworder_300)

s300 <- aggregate_taxa(crypto.top300AG, 'Subtype')

s300

my_col_300<-c(rep("grey35",38),
          rep("black",55),
          rep("grey35",7),
          rep("black",2),
          rep("grey35",11))


pH_outbreak_300_bw <- plot_heatmap(
  s300,
  method = "NMDS",
  distance = "bray",
  sample.label = "Outbreak",
  sample.order = frig,
  taxa.label = "Subtype",
  taxa.order = neworder_300,  # Apply the custom order for Subtype
  trans = identity_trans(),
  low = "white",
   high = "#000033") + 
  facet_grid(. ~ Outbreak, scales = "free_x", space = "free_x", switch = "both") +
  theme(axis.text.y = element_text(size = 10,color = my_col_300)) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text = element_text(angle = 90)) +
  theme(axis.ticks.length = unit(0.05, "cm")) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 38.5) +
  geom_hline(yintercept = 93.5) +
  geom_hline(yintercept = 100.5) +
  geom_hline(yintercept = 102.5) +
  geom_hline(yintercept = 113.5) +
    labs(fill = "Relative\nabundance")

pH_outbreak_300_bw

ggsave("Suppl_S5_agg_scaled_bw_300_R1.png", pH_outbreak_300_bw, 
       height=14, width=10)

## aggregate by Subtype_family
length(unique(tax_table(crypto.top300AG)[,"Subtype_family"]))
print(unique(tax_table(crypto.top300AG)[,"Subtype_family"]))

s300_sf <- aggregate_taxa(crypto.top300AG, 'Subtype_family')
s300_sf

neworder_300sf = c("Ib",
                 "If",       
                 "Ig",        
                 "IIa",        
                 "IId",
                 "IIe",
                 "V",      
                 "IX",       
                 "XIII") # 9
length(neworder_300sf)

pH_outbreak_300sf_bw <- plot_heatmap(
  s300_sf,
  method = "NMDS",
  distance = "bray",
  sample.label = "Outbreak",
  sample.order = frig,
  #taxa.label = "Subtype_family",
  taxa.order = neworder_300sf,  # Apply the custom order for Subtype
  trans = identity_trans(),
  low = "white",
  high = "#000033") + 
  facet_grid(. ~ Outbreak, scales = "free_x", space = "free_x", switch = "both") +
 # theme(axis.text.y = element_text(color = my_col_300)) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text = element_text(angle = 90)) +
   theme(axis.ticks.length = unit(0.05, "cm")) +
   geom_hline(yintercept = 0.5) +
   geom_hline(yintercept = 3.5) +
   geom_hline(yintercept = 6.5) +
   geom_hline(yintercept = 7.5) +
   geom_hline(yintercept = 8.5) +
   geom_hline(yintercept = 9.5) +
  labs(fill = "Relative\nabundance",
       y = "Subtype family",
       subtitle = "A")

pH_outbreak_300sf_bw

ggsave("Fig_2a_agg_scaled_bw_300sf.png", pH_outbreak_300sf_bw, 
       height=3.5, width=10)

