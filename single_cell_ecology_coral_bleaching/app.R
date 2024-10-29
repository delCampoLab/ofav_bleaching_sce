#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

#Libraries

library(shiny)
library(shiny)
library(html5)
library(tidyverse)
library(DT)
library("phyloseq")
library(stringr)
library("tidyverse")
library(randomcoloR)
library(data.table)
library("Seurat")
library(ggplot2)
library(shinydashboard)
library(gridExtra)
library(shinycssloaders)

#Creating main data objects to display
##From the foundational code

psf_16S <- readRDS("psf_16S.rds")
psf_18S <- readRDS("psf_18S.rds")
psf_16S@sam_data$Combine <- c("Healthy Control", "Healthy Control", "Healthy Control",  "Healthy Treatment", "Healthy Treatment", "Healthy Treatment", "Mid-Bleach-1 Control", "Mid-Bleach-1 Control", "Mid-Bleach-1 Control",   "Mid-Bleach-1 Treatment", "Mid-Bleach-1 Treatment", "Mid-Bleach-1 Treatment", "Mid-Bleach-2 Control", "Mid-Bleach-2 Control", "Mid-Bleach-2 Control",  "Mid-Bleach-2 Treatment", "Mid-Bleach-2 Treatment", "Mid-Bleach-2 Treatment",  "Bleached Control", "Bleached Control", "Bleached Control",  "Bleached Treatment", "Bleached Treatment", "Bleached Treatment")
psf_18S@sam_data$Combine <- c("Healthy Control", "Healthy Control", "Healthy Control",  "Healthy Treatment", "Healthy Treatment", "Healthy Treatment", "Mid-Bleach-1 Control", "Mid-Bleach-1 Control", "Mid-Bleach-1 Control",   "Mid-Bleach-1 Treatment", "Mid-Bleach-1 Treatment", "Mid-Bleach-1 Treatment", "Mid-Bleach-2 Control", "Mid-Bleach-2 Control", "Mid-Bleach-2 Control",  "Mid-Bleach-2 Treatment", "Mid-Bleach-2 Treatment", "Mid-Bleach-2 Treatment",  "Bleached Control", "Bleached Control", "Bleached Control",  "Bleached Treatment", "Bleached Treatment", "Bleached Treatment")


tax.clean <- psf_16S@tax_table
tax.clean <- data.frame(row.names = row.names(tax.clean),
                        Kingdom = str_replace(tax.clean[,1], "D_0__",""),
                        Phylum = str_replace(tax.clean[,2], "D_1__",""),
                        Class = str_replace(tax.clean[,3], "D_2__",""),
                        Order = str_replace(tax.clean[,4], "D_3__",""),
                        Family = str_replace(tax.clean[,5], "D_4__",""),
                        Genus = str_replace(tax.clean[,6], "D_5__",""),
                        Species = str_replace(tax.clean[,7], "D_6__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,2] == ""){kingdom <- paste(tax.clean[i,1],"_X", sep = "")
  tax.clean[i, 2] <- kingdom
  tax.clean[i, 3] <- paste(kingdom,"X", sep = "")
  tax.clean[i, 4] <- paste(kingdom,"XX", sep = "")
  tax.clean[i, 5] <- paste(kingdom,"XXX", sep = "")
  tax.clean[i, 6] <- paste(kingdom,"XXXX", sep = "")
  tax.clean[i, 7] <- paste(kingdom,"XXXX_sp.", sep = "")
  } else 
    if (tax.clean[i,3] == ""){phylum <- paste(tax.clean[i,2], "_X", sep = "")
    tax.clean[i, 3] <- phylum
    tax.clean[i, 4] <- paste(phylum,"X", sep = "")
    tax.clean[i, 5] <- paste(phylum,"XX", sep = "")
    tax.clean[i, 6] <- paste(phylum,"XXX", sep = "")
    tax.clean[i, 7] <- paste(phylum,"XXX_sp.", sep = "")
    } else 
      if (tax.clean[i,4] == ""){class <- paste(tax.clean[i,3], "_X", sep = "")
      tax.clean[i, 4] <- class
      tax.clean[i, 5] <- paste(class,"X", sep = "")
      tax.clean[i, 6] <- paste(class,"XX", sep = "")
      tax.clean[i, 7] <- paste(class,"XX_sp.", sep = "")
      } else 
        if (tax.clean[i,5] == ""){order <- paste(tax.clean[i,4], "_X", sep = "")
        tax.clean[i, 5] <- order
        tax.clean[i, 6] <- paste(order,"X", sep = "")
        tax.clean[i, 7] <- paste(order,"X_sp.", sep = "")
        } else 
          if (tax.clean[i,6] == ""){family <- paste(tax.clean[i,5], "_X", sep = "")
          tax.clean[i, 6] <- family
          tax.clean[i, 7] <- paste(family,"X_sp.", sep = "")
          } else 
            if (tax.clean[i,7] == ""){tax.clean$Species[i] <- paste(tax.clean$Genus[i], "sp.",sep = "_")}
}

tax_table(psf_16S) <- as.matrix(tax.clean)

tax.clean <- psf_18S@tax_table

tax.clean <- data.frame(row.names = row.names(psf_18S@tax_table),
                        Domain = str_replace(psf_18S@tax_table[,1], "D_0__",""),
                        Supergroup = str_replace(psf_18S@tax_table[,2], "D_1__",""),
                        Division = str_replace(psf_18S@tax_table[,3], "D_2__",""),
                        Subdivision = str_replace(psf_18S@tax_table[,4], "D_3__",""),
                        Class = str_replace(psf_18S@tax_table[,5], "D_4__",""),
                        Order = str_replace(psf_18S@tax_table[,6], "D_5__",""),
                        Family = str_replace(psf_18S@tax_table[,7], "D_6__",""),
                        Genus = str_replace(psf_18S@tax_table[,8], "D_7__",""),
                        Species = str_replace(psf_18S@tax_table[,9], "D_8__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""

for (i in 1:9){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fille holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,2] == ""){domain <- paste(tax.clean[i,1],"_X", sep = "")
  tax.clean[i, 2] <- domain
  tax.clean[i, 3] <- paste(domain,"X", sep = "")
  tax.clean[i, 4] <- paste(domain,"XX", sep = "")
  tax.clean[i, 5] <- paste(domain,"XXX", sep = "")
  tax.clean[i, 6] <- paste(domain,"XXXX", sep = "")
  tax.clean[i, 7] <- paste(domain,"XXXXX", sep = "")
  tax.clean[i, 8] <- paste(domain,"XXXXXX", sep = "")
  tax.clean[i, 9] <- paste(domain,"XXXXXX_sp.", sep = "")
  } else 
    if (tax.clean[i,3] == ""){supergroup <- paste(tax.clean[i,2], "_X", sep = "")
    tax.clean[i, 3] <- supergroup
    tax.clean[i, 4] <- paste(supergroup,"X", sep = "")
    tax.clean[i, 5] <- paste(supergroup,"XX", sep = "")
    tax.clean[i, 6] <- paste(supergroup,"XXX", sep = "")
    tax.clean[i, 7] <- paste(supergroup,"XXXX.", sep = "")
    tax.clean[i, 8] <- paste(supergroup,"XXXXX", sep = "")
    tax.clean[i, 9] <- paste(supergroup,"XXXXX_sp.", sep = "")
    } else 
      if (tax.clean[i,4] == ""){division <- paste(tax.clean[i,3], "_X", sep = "")
      tax.clean[i, 4] <- division
      tax.clean[i, 5] <- paste(division,"X", sep = "")
      tax.clean[i, 6] <- paste(division,"XX", sep = "")
      tax.clean[i, 7] <- paste(division,"XXX", sep = "")
      tax.clean[i, 8] <- paste(division,"XXXX", sep = "")
      tax.clean[i, 9] <- paste(division,"XXXX_sp.", sep = "")
      } else 
        if (tax.clean[i,5] == ""){subdivision <- paste(tax.clean[i,4], "_X", sep = "")
        tax.clean[i, 5] <- subdivision
        tax.clean[i, 6] <- paste(subdivision,"X", sep = "")
        tax.clean[i, 7] <- paste(subdivision,"XX", sep = "")
        tax.clean[i, 8] <- paste(subdivision,"XXX", sep = "")
        tax.clean[i, 9] <- paste(subdivision,"XXX_sp.", sep = "")
        } else 
          if (tax.clean[i,6] == ""){class <- paste(tax.clean[i,5], "_X", sep = "")
          tax.clean[i, 6] <- class
          tax.clean[i, 7] <- paste(class,"X", sep = "")
          tax.clean[i, 8] <- paste(class,"XX", sep = "")
          tax.clean[i, 9] <- paste(class,"XX_sp.", sep = "")
          } else 
            if (tax.clean[i,7] == ""){order <- paste(tax.clean[i,6], "_X", sep = "")
            tax.clean[i, 7] <- order
            tax.clean[i, 8] <- paste(order,"X", sep = "")
            tax.clean[i, 9] <- paste(order,"X_sp.", sep = "")
            } else 
              if (tax.clean[i,8] == ""){family <- paste(tax.clean[i,7], "_X", sep = "")
              tax.clean[i, 8] <- family
              tax.clean[i, 9] <- paste(family,"X_sp.", sep = "")
              } else 
                if (tax.clean[i,9] == ""){tax.clean$Species[i] <- paste(tax.clean$Genus[i], "sp.",sep = "_")}
}

tax_table(psf_18S) <- as.matrix(tax.clean)

### Extract OTU and Taxa Tables, then combined for both 16S and 18S interactive tables

psf_16S_combined <- merge_samples(psf_16S, "Combine") 
table_16S <- cbind(data.frame(tax_table(psf_16S_combined)), t(data.frame(otu_table(psf_16S_combined))))
names(table_16S)[names(table_16S) == 'T0_Control'] <- 'Healthy Control'
names(table_16S)[names(table_16S) == 'T0_Treatment'] <- 'Healthy Treatment'
names(table_16S)[names(table_16S) == 'T3_Control'] <- 'Mid-Bleach-1 Control'
names(table_16S)[names(table_16S) == 'T3_Treatment'] <- 'Mid-Bleach-1 Treatment'
names(table_16S)[names(table_16S) == 'T5_Control'] <- 'Mid-Bleach-2 Control'
names(table_16S)[names(table_16S) == 'T5_Treatment'] <- 'Mid-Bleach-2 Treatment'
names(table_16S)[names(table_16S) == 'T7_Control'] <- 'Bleached Control'
names(table_16S)[names(table_16S) == 'T7_Treatment'] <- 'Bleached Treatment'

table_16S <- table_16S %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species, `Healthy Control`, `Mid-Bleach-1 Control`, `Mid-Bleach-2 Control`, `Bleached Control`, `Healthy Treatment`, `Mid-Bleach-1 Treatment`, `Mid-Bleach-2 Treatment`, `Bleached Treatment`)

## This is the searchable table for the 16S/Prokaryome

psf_18S_combined <- merge_samples(psf_18S, "Combine") 
table_18S <- cbind(data.frame(tax_table(psf_18S_combined)), t(data.frame(otu_table(psf_18S_combined)))) ## This is the searchable table for the 18S/Eukaryome
names(table_18S)[names(table_18S) == 'T0_Control'] <- 'Healthy Control'
names(table_18S)[names(table_18S) == 'T0_Treatment'] <- 'Healthy Treatment'
names(table_18S)[names(table_18S) == 'T3_Control'] <- 'Mid-Bleach-1 Control'
names(table_18S)[names(table_18S) == 'T3_Treatment'] <- 'Mid-Bleach-1 Treatment'
names(table_18S)[names(table_18S) == 'T5_Control'] <- 'Mid-Bleach-2 Control'
names(table_18S)[names(table_18S) == 'T5_Treatment'] <- 'Mid-Bleach-2 Treatment'
names(table_18S)[names(table_18S) == 'T7_Control'] <- 'Bleached Control'
names(table_18S)[names(table_18S) == 'T7_Treatment'] <- 'Bleached Treatment'

table_18S <- table_18S %>% select(Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species, `Healthy Control`, `Mid-Bleach-1 Control`, `Mid-Bleach-2 Control`, `Bleached Control`, `Healthy Treatment`, `Mid-Bleach-1 Treatment`, `Mid-Bleach-2 Treatment`, `Bleached Treatment`)

## This is the searchable table for the 18S/Eukaryome

## Buildable Bubble plot
#The idea is to build a bubble plot where users can:
#Plot at the taxonomic level of choice (family, order, genus, etc.)
#Filter the shown plot to their desired median relative abundance
#Download a pdf of that plot.
### First need to create a list so that we can sort the microbial taxa phylogenetically along y-axis

tax.sort <- as.data.frame(psf_16S@tax_table)
tax.sort <- tax.sort %>%
  arrange(Kingdom,	Phylum,	Class,	Order,	Family,	Genus,	Species)
orders_16S <- tax.sort %>%
  pull(Order) %>%
  unique() %>%
  toString()
# Split the string into individual values
orders_16S <- strsplit(orders_16S, ", ")[[1]]
orders_16S <- c(orders_16S, "Below Threshold") # This contains a list of all the Orders sorted phylogenetically

tax.sort <- as.data.frame(psf_16S@tax_table)
tax.sort <- tax.sort %>%
  arrange(Kingdom,	Phylum,	Class,	Order,	Family,	Genus,	Species)
family_16S <- tax.sort %>%
  pull(Family) %>%
  unique() %>%
  toString()
family_16S <- strsplit(family_16S, ", ")[[1]]
family_16S <- c(family_16S, "Below Threshold")
class_16S <- tax.sort %>%
  pull(Class) %>%
  unique() %>%
  toString()
class_16S <- strsplit(class_16S, ", ")[[1]]
class_16S <- c(class_16S, "Below Threshold")
genus_16S <- tax.sort %>%
  pull(Genus) %>%
  unique() %>%
  toString()
genus_16S <- strsplit(genus_16S, ", ")[[1]]
genus_16S <- c(genus_16S, "Below Threshold")
phylum_16S <- tax.sort %>%
  pull(Phylum) %>%
  unique() %>%
  toString()
phylum_16S <- strsplit(phylum_16S, ", ")[[1]]
phylum_16S <- c(phylum_16S, "Below Threshold")
species_16S <- tax.sort %>%
  pull(Species) %>%
  unique() %>%
  toString()
species_16S <- strsplit(species_16S, ", ")[[1]]
species_16S <- c(species_16S, "Below Threshold")

tax.sort <- as.data.frame(psf_18S@tax_table)
tax.sort <- tax.sort %>%
  arrange(Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species)

supergroup_18S <- tax.sort %>%
  pull(Supergroup) %>%
  unique() %>%
  toString()
supergroup_18S <- strsplit(supergroup_18S, ", ")[[1]]
supergroup_18S <- c(supergroup_18S, "Below Threshold")
division_18S <- tax.sort %>%
  pull(Division) %>%
  unique() %>%
  toString()
division_18S <- strsplit(division_18S, ", ")[[1]]
division_18S <- c(division_18S, "Below Threshold")
subdivision_18S <- tax.sort %>%
  pull(Subdivision) %>%
  unique() %>%
  toString()
subdivision_18S <- strsplit(subdivision_18S, ", ")[[1]]
subdivision_18S <- c(subdivision_18S, "Below Threshold")
class_18S <- tax.sort %>%
  pull(Class) %>%
  unique() %>%
  toString()
class_18S <- strsplit(class_18S, ", ")[[1]]
class_18S <- c(class_18S, "Below Threshold")
order_18S <- tax.sort %>%
  pull(Order) %>%
  unique() %>%
  toString()
order_18S <- strsplit(order_18S, ", ")[[1]]
order_18S <- c(order_18S, "Below Threshold")
family_18S <- tax.sort %>%
  pull(Family) %>%
  unique() %>%
  toString()
family_18S <- strsplit(family_18S, ", ")[[1]]
family_18S <- c(family_18S, "Below Threshold")
genus_18S <- tax.sort %>%
  pull(Genus) %>%
  unique() %>%
  toString()
genus_18S <- strsplit(genus_18S, ", ")[[1]]
genus_18S <- c(genus_18S, "Below Threshold")
species_18S <- tax.sort %>%
  pull(Species) %>%
  unique() %>%
  toString()
species_18S <- strsplit(species_18S, ", ")[[1]]
species_18S <- c(species_18S, "Below Threshold")

## Color Palette to color these at phylum level

n <- 100
palette <- distinctColorPalette(n)


##From table part code

names(table_16S) <- gsub("-", " ", names(table_16S))
names(table_18S) <- gsub("-", " ", names(table_18S))
coralstatus <- c("Healthy Control", "Mid-Bleach-1 Control",   "Mid-Bleach-2 Control", "Bleached Control", "Healthy Treatment", "Mid-Bleach-1 Treatment", "Mid-Bleach-2 Treatment", "Bleached Treatment")

##From heatmap part

bmin <- readxl::read_excel("bmin.eggnog.emapper.annotations.filtered.xlsx")

dtre <- readxl::read_excel("dtre.eggnog.emapper.annotations.filtered.nots.xlsx")

ofav <- readxl::read_excel("ofav.anno.edited.xlsx")

bmindf <- data.frame(bmin$query, bmin$Preferred_name, bmin$Description, bmin$PFAMs)
dtredf <- data.frame(dtre$query, dtre$Preferred_name, dtre$Description, dtre$PFAMs)
ofavdf <- data.frame(ofav$GeneID, ofav$Name, ofav$Product, ofav$PFAM)

newcn <- c("Gene ID", "Gene name", "Product description", "PFAMs")
colnames(bmindf) <- newcn
colnames(dtredf) <- newcn
colnames(ofavdf) <- newcn

#bmindf$`Gene ID` <- str_replace(bmindf$`Gene ID`, "\\:[^*]*$", "")
#dtredf$`Gene ID` <- str_replace(dtredf$`Gene ID`, "\\:[^*]*$", "")
#ofavdf$`Gene ID` <- str_replace(ofavdf$`Gene ID`, "\\:[^*]*$", "")
ofavdf$`Gene ID`[ofavdf$`Gene ID`== "" | is.na(ofavdf$`Gene ID`)] <- paste0("Unnamed gene ", seq(1, length(ofavdf$`Gene ID`[is.na(ofavdf$`Gene ID`) | ofavdf$`Gene ID` == ""])))

ofavdf <- ofavdf %>% replace_na(list(`Gene ID` = "-", `Gene name` = "-", `Product description` = "-", PFAMs = "-"))
dtredf <- dtredf %>% replace_na(list(`Gene ID` = "-", `Gene name` = "-", `Product description` = "-", PFAMs = "-"))
bmindf <- bmindf %>% replace_na(list(`Gene ID` = "-", `Gene name` = "-", `Product description` = "-", PFAMs = "-"))


ofav.combined.labeled <- readRDS("ofav.combined.labeled.rds")
my_levels <- c("Unidentified", "Calicoblastic Cells", "Cnidocytes", "Digestive Filaments", "Epidermis", "Gastrodermis (Muscle-like)", "Gastrodermis 1", "Gastrodermis 2", "Gland 1", "Gland 2", "Gland 3", "Immune", "Neurons 1", "Neurons 2", "Neurons 3", "Precursors", "Extracellular Breviolum", "Breviolum-hosting Cells", "Durusdinium-hosting Cells")
ofav.combined.labeled@active.ident <- factor(x = ofav.combined.labeled@active.ident, levels = my_levels)
cluster_colors <- c("azure2", "darkorchid4", "plum1", "red1", "seagreen1", "darkolivegreen4", "darkolivegreen1", "darkolivegreen2", "darkorange1", "darkorange2", "darkorange", "darkkhaki", "deepskyblue3", "deepskyblue4", "deepskyblue2", "bisque3", "darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")
ofav.combined.labeled$cluster_colors <- cluster_colors[Idents(ofav.combined.labeled)]
my_levels <- c("Unidentified", "Calicoblastic Cells", "Cnidocytes", "Digestive Filaments", "Epidermis", "Gastrodermis (Muscle-like)", "Gastrodermis 1", "Gastrodermis 2", "Gland 1", "Gland 2", "Gland 3", "Immune", "Neurons 1", "Neurons 2", "Neurons 3", "Precursors", "Extracellular Breviolum", "Breviolum-hosting Cells", "Durusdinium-hosting Cells")
ofav.combined.labeled@meta.data$orig.ident <- recode(ofav.combined.labeled@meta.data$orig.ident, "Ofav-T0-T" = "Ofav Healthy-T", "Ofav-T3-T" = "Ofav Mid-Bleach 1-T", "Ofav-T5-T" = "Ofav-Mid-Bleach 2-T", "Ofav-T7-T" = "Ofav-Bleached-T")
ofav.combined.labeled$seurat_clusters.Timepoints <- paste(Idents(ofav.combined.labeled), ofav.combined.labeled$orig.ident, sep = " ")
DefaultAssay(ofav.combined.labeled) <- "RNA"
all.genes <- rownames(ofav.combined.labeled)
ofav.combined.labeled.hm<- ScaleData(ofav.combined.labeled, features = all.genes)
levels(ofav.combined.labeled.hm)
my_levels <- levels(ofav.combined.labeled.hm)
colors_hm <- c("azure2", "darkorchid4", "plum1", "red1", "seagreen1", "darkolivegreen4", "darkolivegreen1", "darkolivegreen2", "darkorange1", "darkorange2", "darkorange", "darkkhaki", "deepskyblue3", "deepskyblue4", "deepskyblue2", "bisque3", "darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4", "azure2", "darkorchid4", "plum1", "red1", "seagreen1", "darkolivegreen4", "darkolivegreen1", "darkolivegreen2", "darkorange1", "darkorange2", "darkorange", "darkkhaki", "deepskyblue3", "deepskyblue4", "deepskyblue2", "bisque3", "darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4", "azure2", "darkorchid4", "plum1", "red1", "seagreen1", "darkolivegreen4", "darkolivegreen1", "darkolivegreen2", "darkorange1", "darkorange2", "darkorange", "darkkhaki", "deepskyblue3", "deepskyblue4", "deepskyblue2", "bisque3", "darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4", "azure2", "darkorchid4", "plum1", "red1", "seagreen1", "darkolivegreen4", "darkolivegreen1", "darkolivegreen2", "darkorange1", "darkorange2", "darkorange", "darkkhaki", "deepskyblue3", "deepskyblue4", "deepskyblue2", "bisque3", "darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")
ofav.combined.labeled.hm@active.ident <- factor(x = ofav.combined.labeled.hm@active.ident, levels = my_levels)
ofav.combined.labeled.hm@meta.data[,c("seurat_clusters.Timepoints", "orig.ident"),drop=FALSE] %>% tibble::rownames_to_column() %>% group_by(seurat_clusters.Timepoints) %>% slice_sample(n= 25) -> index

#UI

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Home", tabName = "home", icon = icon("house")),
    menuItem(" ASV Tables", tabName = "tables", icon = icon("table"), startExpanded = FALSE,
             menuSubItem("Prokaryome", tabName = "tabprokaryome"),
             menuSubItem("Eukaryome", tabName = "tabeukaryome")),
    menuItem("Abundance Bubbleplots", icon = icon("circle-dot"), tabName = "bubplots", startExpanded = FALSE ,
             menuSubItem("Prokaryome", tabName = "bubprokaryome"),
             menuSubItem("Eukaryome", tabName = "bubeukaryome")),
    menuItem("scRNA-seq Heatmap", icon = icon("chess-board"), tabName = "heatmap")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "home",
            fluidRow(
              tags$h1("Single-Cell Ecology of a Multi-Partner Symbiosis through Disruption", style = "margin-left: 20px;"),
              tags$h4("Anthony M. Bonacolta, Grace Snyder, Richard Karp, Emily Yeager, Alexandra Wen, Carly Dennison, Jordi Nonell, Nikki Traylor-Knowles, Andrew C. Baker, & Javier del Campo", style = "margin-left: 20px;"),
              tags$h3("Abstract", style = "margin-left: 20px;"),
              tags$hr(style = "margin-top: 20px; margin-bottom: 20px; margin-left: 20px; margin-right: 600px; border-top: 2px solid #333;"),
              tags$h4("The coral–algae symbiosis is disrupted by marine heatwaves in a process popularly known as bleaching, threatening coral reefs worldwide. Bleaching is the manifestation of a response to heat stress at a cellular scale that has a planetary impact. 
                      Combining innovative techniques such as single-cell transcriptomics with metabarcoding we can dissect the microbial ecology of bleaching, one cell at a time. We subjected replicate fragments of a colony of the endangered Caribbean coral", tags$em("Orbicella faveolata"), "which harbors two co-dominant algal symbionts (", tags$em("Durusdinium"), " and  ", tags$em("Breviolum"), "; Family Symbiodiniaceae), to thermal stress. 
                      We then characterized the simultaneous gene expression and shifts in the taxonomic composition of microbes and coral cells.  Our study identifies differential gene expression patterns specific to certain holobiont cell clusters, including a differential suppression of heat stress genes in coral gastrodermal cells associated with different algal symbionts, which may explain the higher thermal tolerance of", tags$em("Durusdinium"), ". 
                      Additionally, we identified key genes and microbial abundance changes linked to alterations in nutrient dynamics during the breakdown of the coral-algal symbiosis, emphasizing nitrogen cycling as crucial to understanding the single-cell ecology of coral bleaching. Our findings underscore the necessity of understanding the responses of the different cellular components within the coral holobiont during bleaching to fully comprehend the mechanisms leading to symbiosis breakdown, coral mortality, and ultimately, reef decline.", style = "margin-left: 20px; margin-right: 600px; text-align: justified;"),
              HTML("<br>"),
              tags$div(style = "width: 45vw; margin-left: 150px;", box(width = NULL, style =  "overflow: hidden;", tags$img(src = "graphical_abstract_v2.png", style = "width: 100%; height: auto; display: block;"))),
            )),
    tabItem(tabName = "tabprokaryome",
            fluidRow(
             tags$h1("Prokaryote ASV abundance table", style = "margin-left: 20px;"),  
             tags$h4("The following table displays all the different prokaryote ASV recovered from the", tags$em("O. faveolata"), " control and treated cells and their relative abundance throughout each timepoint of its bleaching. This data can be filtered, downloaded as .csv, copied or printed.", style = "margin-left: 20px;"),
             HTML("<br>"),
             box(
               DTOutput('asvprok', width = "100%"), width = 12, height = "auto", solidHeader = TRUE, status = "primary"))),
    tabItem(tabName = "tabeukaryome",
            fluidRow(
              tags$h1("Eukaryote ASV abundance table", style = "margin-left: 20px;"),  
              tags$h4("The following table displays all the different eukaryote ASV recovered from the", tags$em("O. faveolata"), " control and treated cells and their relative abundance throughout each timepoint of its bleaching. This data can be filtered, downloaded as .csv, copied or printed.", style = "margin-left: 20px;"),
              HTML("<br>"),
              box(
                DTOutput('asveuk', width = "100%"), width = 12, height = "auto", solidHeader = TRUE, status = "primary"))),
    tabItem(tabName = "bubprokaryome",
            tags$style(HTML("#proktaxonlevel .control-label {
                            display: inline-block;
                            padding-bottom: 20px; 
                            }
                            ")),
            fluidRow(
              tags$h1(tags$em("O. faveolata"), "prokaryome abundance bubbleplot", style = "margin-left: 20px;"),
              tags$h4("Generate a bubbleplot to visualize the median relative abundance of prokaryotic taxa across the different stages of bleaching in control and treated", tags$em("O. faveolata"), "cells.", style = "margin-left: 20px;"),
              HTML("<br>"),
              box(
                fluidRow(
                  tags$h5("Select a taxonomic level to group the ASV and a cutoff to filter these upon their median relative abundance.", style = "margin-left: 15px;"),
                  column(4, radioButtons('proktaxonlevel', 'Taxonomic level of display', list("Phylum" = "Phylum", 
                                                                                              "Class" = "Class", 
                                                                                              "Order" = "Order", 
                                                                                              "Family" = "Family",
                                                                                              "Genus" = "Genus"
                  ), "Family", inline = TRUE)),
                  column(7, sliderInput('cutoffprokaryome', 'Cutoff level', min = 0.00, max = 15, value = 0.25, step = 0.05, ticks = TRUE)),
                  column(1)), height = 230, width = 12, solidHeader = TRUE, status = "primary",
                footer = actionButton("loadbub16S", "Load a bubbleplot"))),
            fluidRow(
              box(title = "Prokaryome ASV bubbleplot", uiOutput("bubplot16S_ui"), width = 12, solidHeader = TRUE, status = "primary", downloadButton("downloadBubplotprok", "Download this bubble plot")))
    ),
    tabItem(tabName = "bubeukaryome",
            tags$style(HTML("#euktaxonlevel .control-label {
                            display: inline-block;
                            padding-bottom: 20px; 
                            }
                            ")),
            fluidRow(
              tags$h1(tags$em("O. faveolata"), "eukaryome abundance bubbleplot", style = "margin-left: 20px;"),
              tags$h4("Generate a bubbleplot to visualize the median relative abundance of prokaryotic taxa across the different stages of bleaching in control and treated", tags$em("O. faveolata"), "cells.", style = "margin-left: 20px;"),
              HTML("<br>"),
              box(
                fluidRow(
                  tags$h5("Select a taxonomic level to group the ASV and a cutoff to filter these upon their median relative abundance.", style = "margin-left: 15px;"),
                  column(4, radioButtons('euktaxonlevel', 'Taxonomic level of display', list("Division" = "Division", 
                                                                                             "Subdivision" = "Subdivision", 
                                                                                             "Class" = "Class", 
                                                                                             "Order" = "Order",
                                                                                             "Family" = "Family",
                                                                                             "Genus" = "Genus"
                  ), "Genus", inline = TRUE)),
                  column(7, sliderInput('cutoffeukaryome', 'Cutoff level', min = 0.00, max = 2.5, value = 0.25, step = 0.05, ticks = TRUE)),
                  column(1)), height = 230, width = 12, solidHeader = TRUE, status = "primary",
                footer = actionButton("loadbub18S", "Load a bubbleplot"))),
            fluidRow(
              box(title = "Eukaryome ASV bubbleplot", uiOutput("bubplot18S_ui"), width = 12, solidHeader = TRUE, status = "primary", downloadButton("downloadBubploteuk", "Download this bubble plot")))
    ),
    tabItem(tabName = "heatmap",
            fluidRow(
              tags$h1("scRNA-seq gene expression heatmap", style = "margin-left: 20px;"),
              tags$h4("Plot a heatmap that shows the expression levels of genes belonging to ", tags$em("O. faveolata"), ",", tags$em("D. trenchii"), "and", tags$em("B. minutum"), "and their evolution throughout the bleaching process in control and treated", tags$em("O. faveolata"), "cells", style = "margin-left: 20px;"),
              tabBox(
                tabPanel(HTML("<em>O. faveolata</em>"),
                         tags$h5("Select genes of your interest by clicking on their row. You can navigate through", tags$em("O. faveolata"), ",", tags$em("D. trenchii"), "and", tags$em("B. minutum"), "genes and load them to the selection table by pressing the 'Add selected genes' button."),
                         DTOutput('ofavtab', width = "100%")),
                tabPanel(HTML("<em>D. trenchii</em>"),
                         tags$h5("Select genes of your interest by clicking on their row. You can navigate through", tags$em("O. faveolata"), ",", tags$em("D. trenchii"), "and", tags$em("B. minutum"), "genes and load them to the selection table by pressing the 'Add selected genes' button."),
                         DTOutput('dtretab', width = "100%")),
                tabPanel(HTML("<em>B. minutum</em>"),      
                         tags$h5("Select genes of your interest by clicking on their row. You can navigate through", tags$em("O. faveolata"), ",", tags$em("D. trenchii"), "and", tags$em("B. minutum"), "genes and load them to the selection table by pressing the 'Add selected genes' button."),
                         DTOutput('bmintab', width = "100%")),
                type = "tabs",
                header = tags$style(type="text/css", "
                                        div.dataTables_wrapper div.dataTables_filter input {
                                          width: 300px;
                                          margin-right: 20px; 
                                        }
                                         div.data"
                ),
                footer = list(actionButton("add", "Add selected genes", style = "margin-left: 10px;"), actionButton("ofavadd", HTML("Add all <em>O. faveolata</em> genes"), style = "margin-left: 10px;"), actionButton("dtreadd", HTML("Add all <em>D. trenchii</em> genes"), style = "margin-left: 10px;"), actionButton("bminadd", HTML("Add all <em>B. minutum</em> genes"), style = "margin-left: 10px;"))
                , width = 12, height = 890)
            ),
            fluidRow(
              box(
                DTOutput('selectedtab', width = "100%"), solidHeader = TRUE, status = "primary",
                header = tags$style(type="text/css", "
                                        div.dataTables_wrapper div.dataTables_filter input {
                                          width: 300px;
                                          margin-right: 20px; 
                                        }
                                         div.data"
                ),
                footer = list(actionButton("remove", "Remove selected genes"), actionButton("ofavremove", HTML("Remove all <em>O. faveolata</em> genes"), style = "margin-left: 10px;"), actionButton("dtreremove", HTML("Remove all <em>D. trenchii</em> genes"), style = "margin-left: 10px;"), actionButton("bminremove", HTML("Remove all <em>B. minutum</em> genes"), style = "margin-left: 10px;"), actionButton("loadhm", "Load selection into heatmap", style = "margin-left: 10px;"))
                ,width = 12, height = 810)
            ),
            fluidRow(
              box(title = "Gene expression heatmap", solidHeader = TRUE, status = "primary", uiOutput("heatmap_ui"), width = 12,downloadButton("downloadHeatmap", "Download this heatmap")))
    )
  )
)



ui <- dashboardPage(
  dashboardHeader(title = shiny::h3("Single-cell ecology of coral bleaching", style = "font-size: 20px; line-height: 7px;"), titleWidth = 500),
  sidebar,
  body
)

server <- function(input, output, session) {
  
  #servertab
  
  table_16S_filter = reactive(table_16S)
  table_18S_filter = reactive(table_18S)
  
  output$asvprok <- DT::renderDataTable({ #Prokaryome ASV table
    table_16S_react = table_16S_filter()
    DT::datatable(table_16S_react, 
                  extensions = c("Scroller", "Buttons"),
                  style = "bootstrap4",
                  options = list(
                    searching = TRUE,
                    autoWidth = TRUE,
                    scrollX = TRUE,
                    scrollY = "750px",
                    dom = 'l<"sep">Bfrtip',
                    fillContainter = TRUE,
                    buttons = list(
                      extend = "csv",
                      text = "<i class='glyphicon glyphicon-download-alt'></i> Download the displayed table as .csv"
                    ),
                    pageLength = 50,
                    lengthMenu = list(c(50,100,200, nrow(table_16S)), c("50", "100", "200", "All")),
                    orderClasses=TRUE,
                    columnDefs = list(list(width = '140px', targets = c(8:15)))
                  )
    )  
  })
  
  output$asveuk <- DT::renderDataTable({ #Eukaryome ASV table
    table_18S_react = table_18S_filter()
    DT::datatable(table_18S_react, 
                  extensions = c("Scroller", "Buttons"),
                  style = "bootstrap4",
                  options = list(
                    searching = TRUE,
                    autoWidth = TRUE,
                    scrollX = TRUE,
                    scrollY = "750px",
                    dom = 'l<"sep">B + <br> + frtip',
                    fillContainer = TRUE,
                    buttons = list(
                      extend = "csv",
                      text = "<i class='glyphicon glyphicon-download-alt'></i> Download the displayed table as .csv"
                    ),
                    pageLength = 50,
                    lengthMenu = list(c(50,100,200, nrow(table_18S)), c("50", "100", "200", "All")),
                    orderClasses=TRUE,
                    columnDefs = list(list(width = '100px', targets = c(10:17)))
                    
                  )
    )  
  })
  
  #serverbub
  
  
  datatable16S <- reactive({ #reactive prokaryome datatable to bubbleplot. Changes according to the taxon and cutoff level of interest
    
    taxonchoice <- input$proktaxonlevel
    
    if(input$proktaxonlevel == "Phylum"){
      sd16S <- psf_16S %>% 
        tax_glom(taxrank = "Phylum") %>%
        transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
        psmelt() %>% group_by(OTU, Combine)  %>% 
        summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
      bub16S <- merge_samples(psf_16S, "Combine")
      type_taxa16S <- bub16S %>% 
        tax_glom(taxrank = "Phylum") %>% 
        transform_sample_counts(function(x) {x/sum(x)*100} )  
      type_taxa16S <- type_taxa16S %>%  
        psmelt() %>% 
        arrange(OTU) 
      dat16S <- data.table(type_taxa16S) 
      dat16S$SD <- sd16S$sd 
      dat16S$Phylum <- as.character(dat16S$Phylum)
      dat16S[, median := median(Abundance, na.rm = FALSE), 
             by = Phylum]
      dat16S[(median <= input$cutoffprokaryome), Phylum := "Below Threshold"] 
    } else {
      if(input$proktaxonlevel == "Class"){
        sd16S <- psf_16S %>% 
          tax_glom(taxrank = "Class") %>%
          transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
          psmelt() %>% group_by(OTU, Combine)  %>% 
          summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
        bub16S <- merge_samples(psf_16S, "Combine")
        type_taxa16S <- bub16S %>% 
          tax_glom(taxrank = "Class") %>% 
          transform_sample_counts(function(x) {x/sum(x)*100} )  
        type_taxa16S <- type_taxa16S %>%  
          psmelt() %>% 
          arrange(OTU) 
        dat16S <- data.table(type_taxa16S) 
        dat16S$SD <- sd16S$sd 
        dat16S$Class <- as.character(dat16S$Class)
        dat16S[, median := median(Abundance, na.rm = FALSE), 
               by = Class]
        dat16S[(median <= input$cutoffprokaryome), Class := "Below Threshold"] 
        dat16S[(Class == "Below Threshold"), Phylum := "Other"]  
      } else {
        if(input$proktaxonlevel == "Order"){
          sd16S <- psf_16S %>% 
            tax_glom(taxrank = "Order") %>%
            transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
            psmelt() %>% group_by(OTU, Combine)  %>% 
            summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
          bub16S <- merge_samples(psf_16S, "Combine")
          type_taxa16S <- bub16S %>% 
            tax_glom(taxrank = "Order") %>% 
            transform_sample_counts(function(x) {x/sum(x)*100} )  
          type_taxa16S <- type_taxa16S %>%  
            psmelt() %>% 
            arrange(OTU) 
          dat16S <- data.table(type_taxa16S) 
          dat16S$SD <- sd16S$sd 
          dat16S$Order <- as.character(dat16S$Order)
          dat16S[, median := median(Abundance, na.rm = FALSE), 
                 by = Order]
          dat16S[(median <= input$cutoffprokaryome), Order := "Below Threshold"] 
          dat16S[(Order == "Below Threshold"), Phylum := "Other"]
        } else {
          if(input$proktaxonlevel == "Family"){
            sd16S <- psf_16S %>% 
              tax_glom(taxrank = "Family") %>%
              transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
              psmelt() %>% group_by(OTU, Combine)  %>% 
              summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
            bub16S <- merge_samples(psf_16S, "Combine")
            type_taxa16S <- bub16S %>% 
              tax_glom(taxrank = "Family") %>% 
              transform_sample_counts(function(x) {x/sum(x)*100} )  
            type_taxa16S <- type_taxa16S %>%  
              psmelt() %>% 
              arrange(OTU) 
            dat16S <- data.table(type_taxa16S) 
            dat16S$SD <- sd16S$sd 
            dat16S$Family <- as.character(dat16S$Family)
            dat16S[, median := median(Abundance, na.rm = FALSE), 
                   by = Family]
            dat16S[(median <= input$cutoffprokaryome), Family := "Below Threshold"] 
            dat16S[(Family == "Below Threshold"), Phylum := "Other"]
          } else {
            if(input$proktaxonlevel == "Genus"){
              sd16S <- psf_16S %>% 
                tax_glom(taxrank = "Genus") %>%
                transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
                psmelt() %>% group_by(OTU, Combine)  %>% 
                summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
              bub16S <- merge_samples(psf_16S, "Combine")
              type_taxa16S <- bub16S %>% 
                tax_glom(taxrank = "Genus") %>% 
                transform_sample_counts(function(x) {x/sum(x)*100} )  
              type_taxa16S <- type_taxa16S %>%  
                psmelt() %>% 
                arrange(OTU) 
              dat16S <- data.table(type_taxa16S) 
              dat16S$SD <- sd16S$sd 
              dat16S$Genus <- as.character(dat16S$Genus)
              dat16S[, median := median(Abundance, na.rm = FALSE), 
                     by = Genus]
              dat16S[(median <= input$cutoffprokaryome), Genus := "Below Threshold"] 
              dat16S[(Genus == "Below Threshold"), Phylum := "Other"]
            } 
          }
        }
      }
    }
    return(dat16S)
  }) %>% bindCache(input$proktaxonlevel, input$cutoffprokaryome)
  
  
  bubplot_16S <- reactive({ #reactive generation of all the ggplot2 bubbleplots. Changes upon selection of taxon level.
    
    dat16S <- datatable16S()
    ggplot_phylum16S <- ggplot(dat16S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Phylum, levels=phylum_16S))) + geom_point(aes(size = Abundance, fill=Phylum, color=Phylum), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Prokaryotic Phylum", size = "Relative Abundance %", fill="Prokaryotic Phylum") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Phyla within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffprokaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Phylum, color=Phylum), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Phylum, color=Phylum), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_class16S <- ggplot(dat16S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Class, levels=class_16S))) + geom_point(aes(size = Abundance, fill=Phylum, color=Phylum), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Prokaryotic Class", size = "Relative Abundance %", fill="Prokaryotic Phylum") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Classes within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffprokaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Phylum, color=Phylum), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Phylum, color=Phylum), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_order16S <- ggplot(dat16S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Order, levels=orders_16S))) + geom_point(aes(size = Abundance, fill=Phylum, color=Phylum), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Prokaryotic Order", size = "Relative Abundance %", fill="Prokaryotic Phylum") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Orders within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffprokaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Phylum, color=Phylum), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Phylum, color=Phylum), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_family16S <- ggplot(dat16S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Family, levels=family_16S))) + geom_point(aes(size = Abundance, fill=Phylum, color=Phylum), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Prokaryotic Family", size = "Relative Abundance %", fill="Prokaryotic Phylum") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Families within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffprokaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Phylum, color=Phylum), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Phylum, color=Phylum), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_genus16S <- ggplot(dat16S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Genus, levels=genus_16S))) + geom_point(aes(size = Abundance, fill=Phylum, color=Phylum), alpha = 0.5, shape = 21) + #Creació del bubble plot
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Prokaryotic Genus", size = "Relative Abundance %", fill="Prokaryotic Phylum") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Genera within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffprokaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Phylum, color=Phylum), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Phylum, color=Phylum), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    if (input$proktaxonlevel == "Phylum") {
      return(ggplot_phylum16S)
    } else {
      if (input$proktaxonlevel == "Class") {
        return(ggplot_class16S)
      } else {
        if (input$proktaxonlevel == "Order") {
          return(ggplot_order16S)
        } else {
          if (input$proktaxonlevel == "Family") {
            return(ggplot_family16S)
          } else {
            if (input$proktaxonlevel == "Genus") {
              return(ggplot_genus16S)
            } }
        }
      }
    }
  }) 
  
  observeEvent(input$proktaxonlevel, {isolate({ #exclusive for the prokaryome. Forces a change in the range of the cutoff sliderInput each time there's a change in taxa level selection
    if (input$proktaxonlevel == "Phylum") {
      updateSliderInput(session, "cutoffprokaryome", max = 17.5)
    } else {
      if (input$proktaxonlevel == "Class") {
        updateSliderInput(session, "cutoffprokaryome", max = 17.5)
      } else {
        if (input$proktaxonlevel == "Order") {
          updateSliderInput(session, "cutoffprokaryome", max = 7.5)
        } else {
          if (input$proktaxonlevel == "Family") {
            updateSliderInput(session, "cutoffprokaryome", max = 7.5)
          } else {
            if (input$proktaxonlevel == "Genus") {
              updateSliderInput(session, "cutoffprokaryome", max = 5)
            }
          }
        }
      }
    }
  })
  }) 
  
  dynamicheight16S <- reactive({ #reactive value for the height of the prokaryome bubbleplot to ensure its correct display
    bubplot_16S_build <- ggplot_build(bubplot_16S())
    bubplot_16S_yticks <- bubplot_16S_build$layout$panel_params[[1]]$y.sec$breaks
    bubplot_16S_num_ticks <- length(bubplot_16S_yticks)
    return(bubplot_16S_num_ticks*22.72 + 500)
  })
  
  output$bubplot16S_ui <- renderUI({ #reactiveUI. When no bubbleplot is loading at the start, displays a message. Once the first bubbleplot is requested it displays the spinner and a box which's height will vary depending on the plot height.
    if(input$loadbub16S > 0) {
      withSpinner(plotOutput("bubplot_16S", height = "auto", width = "auto"), type = 6)
    } else {
      HTML(
        '<div style="display: flex; justify-content: center; align-items: center; height: auto; min-height: 200px;">',
        '<p>Load your prokaryome abundance bubbleplot!</p>',
        '</div>'
      )
    }
  })
  
  output$bubplot_16S <- renderPlot({bubplot_16S()}, height = dynamicheight16S, width = 1800) %>% bindEvent(input$loadbub16S) #output of the prokaryome bubbleplot. The event is bound to the loading button so a change of plot can only be triggered by it
  
  dynamicheightdownload16S <- reactive({ #returns the height of the png image that contains the downloaded bubbleplot.
    bubplot_16S_build <- ggplot_build(bubplot_16S())
    bubplot_16S_yticks <- bubplot_16S_build$layout$panel_params[[1]]$y.sec$breaks
    bubplot_16S_num_ticks <- length(bubplot_16S_yticks)
    return(bubplot_16S_num_ticks + 5)
  })
  
  output$downloadBubplotprok <- downloadHandler( #download button for the prokaryome bubbleplot. 
    filename = function() { "prokaryome_bubbleplot.png"},
    content = function(file) {
      ggsave(file, bubplot_16S(), device = "png", width = 10, height = dynamicheightdownload16S(), limitsize = FALSE)
    }
  )
  
  datatable18S <- reactive({  #reactive eukaryome datatable to bubbleplot. Changes according to the taxon and cutoff level of interest
    
    taxonchoice <- input$euktaxonlevel
    
    if(input$euktaxonlevel == "Division"){
      sd18S <- psf_18S %>% 
        tax_glom(taxrank = "Division") %>%
        transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
        psmelt() %>% group_by(OTU, Combine)  %>% 
        summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
      bub18S <- merge_samples(psf_18S, "Combine")
      type_taxa18S <- bub18S %>% 
        tax_glom(taxrank = "Division") %>% 
        transform_sample_counts(function(x) {x/sum(x)*100} )  
      type_taxa18S <- type_taxa18S %>%  
        psmelt() %>% 
        arrange(OTU) 
      dat18S <- data.table(type_taxa18S) 
      dat18S$SD <- sd18S$sd 
      dat18S$Division <- as.character(dat18S$Division)
      dat18S[, median := median(Abundance, na.rm = FALSE), 
             by = Division]
      dat18S[(median <= input$cutoffeukaryome), Division := "Below Threshold"] 
    } else {
      if(input$euktaxonlevel == "Subdivision"){
        sd18S <- psf_18S %>% 
          tax_glom(taxrank = "Subdivision") %>%
          transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
          psmelt() %>% group_by(OTU, Combine)  %>% 
          summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
        bub18S <- merge_samples(psf_18S, "Combine")
        type_taxa18S <- bub18S %>% 
          tax_glom(taxrank = "Subdivision") %>% 
          transform_sample_counts(function(x) {x/sum(x)*100} )  
        type_taxa18S <- type_taxa18S %>%  
          psmelt() %>% 
          arrange(OTU) 
        dat18S <- data.table(type_taxa18S) 
        dat18S$SD <- sd18S$sd 
        dat18S$Subdivision <- as.character(dat18S$Subdivision)
        dat18S[, median := median(Abundance, na.rm = FALSE), 
               by = Subdivision]
        dat18S[(median <= input$cutoffeukaryome), Subdivision := "Below Threshold"] 
        dat18S[(Subdivision == "Below Threshold"), Division := "Other"]  
      } else {
        if(input$euktaxonlevel == "Class"){
          sd18S <- psf_18S %>% 
            tax_glom(taxrank = "Class") %>%
            transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
            psmelt() %>% group_by(OTU, Combine)  %>% 
            summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
          bub18S <- merge_samples(psf_18S, "Combine")
          type_taxa18S <- bub18S %>% 
            tax_glom(taxrank = "Class") %>% 
            transform_sample_counts(function(x) {x/sum(x)*100} )  
          type_taxa18S <- type_taxa18S %>%  
            psmelt() %>% 
            arrange(OTU) 
          dat18S <- data.table(type_taxa18S) 
          dat18S$SD <- sd18S$sd 
          dat18S$Class <- as.character(dat18S$Class)
          dat18S[, median := median(Abundance, na.rm = FALSE), 
                 by = Class]
          dat18S[(median <= input$cutoffeukaryome), Class := "Below Threshold"] 
          dat18S[(Class == "Below Threshold"), Division := "Other"]  
        } else {
          if(input$euktaxonlevel == "Order"){
            sd18S <- psf_18S %>% 
              tax_glom(taxrank = "Order") %>%
              transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
              psmelt() %>% group_by(OTU, Combine)  %>% 
              summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
            bub18S <- merge_samples(psf_18S, "Combine")
            type_taxa18S <- bub18S %>% 
              tax_glom(taxrank = "Order") %>% 
              transform_sample_counts(function(x) {x/sum(x)*100} )  
            type_taxa18S <- type_taxa18S %>%  
              psmelt() %>% 
              arrange(OTU) 
            dat18S <- data.table(type_taxa18S) 
            dat18S$SD <- sd18S$sd 
            dat18S$Order <- as.character(dat18S$Order)
            dat18S[, median := median(Abundance, na.rm = FALSE), 
                   by = Order]
            dat18S[(median <= input$cutoffeukaryome), Order := "Below Threshold"] 
            dat18S[(Order == "Below Threshold"), Division := "Other"]
          } else {
            if(input$euktaxonlevel == "Family"){
              sd18S <- psf_18S %>% 
                tax_glom(taxrank = "Family") %>%
                transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
                psmelt() %>% group_by(OTU, Combine)  %>% 
                summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
              bub18S <- merge_samples(psf_18S, "Combine")
              type_taxa18S <- bub18S %>% 
                tax_glom(taxrank = "Family") %>% 
                transform_sample_counts(function(x) {x/sum(x)*100} )  
              type_taxa18S <- type_taxa18S %>%  
                psmelt() %>% 
                arrange(OTU) 
              dat18S <- data.table(type_taxa18S) 
              dat18S$SD <- sd18S$sd 
              dat18S$Family <- as.character(dat18S$Family)
              dat18S[, median := median(Abundance, na.rm = FALSE), 
                     by = Family]
              dat18S[(median <= input$cutoffeukaryome), Family := "Below Threshold"] 
              dat18S[(Family == "Below Threshold"), Division := "Other"]
            } else {
              if(input$euktaxonlevel == "Genus"){
                sd18S <- psf_18S %>% 
                  tax_glom(taxrank = "Genus") %>%
                  transform_sample_counts(function(x) {x/sum(x)*100} ) %>% 
                  psmelt() %>% group_by(OTU, Combine)  %>% 
                  summarize(sd = sd(Abundance, na.rm=FALSE)) %>% arrange(OTU) 
                bub18S <- merge_samples(psf_18S, "Combine")
                type_taxa18S <- bub18S %>% 
                  tax_glom(taxrank = "Genus") %>% 
                  transform_sample_counts(function(x) {x/sum(x)*100} )  
                type_taxa18S <- type_taxa18S %>%  
                  psmelt() %>% 
                  arrange(OTU) 
                dat18S <- data.table(type_taxa18S) 
                dat18S$SD <- sd18S$sd 
                dat18S$Genus <- as.character(dat18S$Genus)
                dat18S[, median := median(Abundance, na.rm = FALSE), 
                       by = Genus]
                dat18S[(median <= input$cutoffeukaryome), Genus := "Below Threshold"] 
                dat18S[(Genus == "Below Threshold"), Division := "Other"]
              } 
            }
          }
        }
      }
    }
    return(dat18S)
  }) %>% bindCache(input$euktaxonlevel, input$cutoffeukaryome)
  
  
  bubplot_18S <- reactive({ #reactive ggplot eukaryome bubbleplot. Returns only one plot based on the taxa chosen by the user.
    
    dat18S <- datatable18S()
    
    ggplot_division18S <- ggplot(dat18S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Division, levels=division_18S))) + geom_point(aes(size = Abundance, fill=Division, color=Division), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Eukaryotic Division", size = "Relative Abundance %", fill="Eukaryotic Division") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Divisions within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffeukaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Division, color=Division), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Division, color=Division), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_subdivision18S <- ggplot(dat18S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Subdivision, levels=subdivision_18S))) + geom_point(aes(size = Abundance, fill=Division, color=Division), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Eukaryotic Subdivision", size = "Relative Abundance %", fill="Eukaryotic Division") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Subdivisions within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffeukaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Division, color=Division), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Division, color=Division), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_class18S <- ggplot(dat18S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Class, levels=class_18S))) + geom_point(aes(size = Abundance, fill=Division, color=Division), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Eukaryotic Class", size = "Relative Abundance %", fill="Eukaryotic Division") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Classes within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffeukaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Division, color=Division), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Division, color=Division), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_order18S <- ggplot(dat18S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Order, levels=order_18S))) + geom_point(aes(size = Abundance, fill=Division, color=Division), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Eukaryotic Order", size = "Relative Abundance %", fill="Eukaryotic Division") + theme_bw() +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Orders within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffeukaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Division, color=Division), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Division, color=Division), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_family18S <- ggplot(dat18S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Family, levels=family_18S))) + geom_point(aes(size = Abundance, fill=Division, color=Division), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Eukaryotic Family", size = "Relative Abundance %", fill="Eukaryotic Division") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Families within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffeukaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Division, color=Division), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Division, color=Division), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    ggplot_genus18S <- ggplot(dat18S[Abundance > 0], aes(x = factor(Sample, levels= coralstatus), y = factor(Genus, levels=genus_18S))) + geom_point(aes(size = Abundance, fill=Division, color=Division), alpha = 0.5, shape = 21) + 
      scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(1,10,50,75)) + 
      labs(x="Sample", y = "Eukaryotic Genus", size = "Relative Abundance %", fill="Eukaryotic Division") + theme_bw(base_size = 15) +
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) + 
      ggtitle(label=expression(paste("Microbial Genera within ", italic("O. faveolata"))), subtitle = paste0("With median relative abundance above ", input$cutoffeukaryome, "%")) + geom_point(aes(size=Abundance+SD, fill = Division, color=Division), alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-SD, fill = Division, color=Division), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
    
    if (input$euktaxonlevel == "Division") {
      return(ggplot_division18S)
    } else {
      if (input$euktaxonlevel == "Subdivision") {
        return(ggplot_subdivision18S)
      } else {
        if (input$euktaxonlevel == "Class") {
          return(ggplot_class18S)
        } else {
          if (input$euktaxonlevel == "Order") {
            return(ggplot_order18S)
          } else {
            if (input$euktaxonlevel == "Family") {
              return(ggplot_family18S)
            } else {
              if (input$euktaxonlevel == "Genus") {
                return(ggplot_genus18S)
              } else {
                if (input$euktaxonlevel == "Species") {
                  return(ggplot_species18S)
                }
              }
            }
          }
        }
      }
    }
    
  }) 
  
  
  dynamicheight18S <- reactive({ #height for the bubbleplot output.
    bubplot_18S_build <- ggplot_build(bubplot_16S())
    bubplot_18S_yticks <- bubplot_18S_build$layout$panel_params[[1]]$y.sec$breaks
    bubplot_18S_num_ticks <- length(bubplot_18S_yticks)
    return(bubplot_18S_num_ticks*22.72 + 500)
  })
  
  output$bubplot18S_ui <- renderUI({ #reactive UI.
    if(input$loadbub18S > 0) {
      withSpinner(plotOutput("bubplot_18S", height = "auto", width = "auto"), type = 6)
    } else {
      HTML(
        '<div style="display: flex; justify-content: center; align-items: center; height: auto; min-height: 200px;">',
        '<p>Load your eukaryome abundance bubbleplot!</p>',
        '</div>'
      )
    }
  })
  
  output$bubplot_18S <- renderPlot({bubplot_18S()}, height = dynamicheight18S, width = 1800) %>% bindEvent(input$loadbub18S) #eukaryome bubbleplot output.
  
  dynamicheightdownload18S <- reactive({ #height for the eukaryome bubbleplot png downloadable image.
    bubplot_18S_build <- ggplot_build(bubplot_18S())
    bubplot_18S_yticks <- bubplot_18S_build$layout$panel_params[[1]]$y.sec$breaks
    bubplot_18S_num_ticks <- length(bubplot_18S_yticks)
    return(bubplot_18S_num_ticks + 5)
  })
  
  output$downloadBubploteuk <- downloadHandler( #download button for the eukaryome bubbleplot png.
    filename = function() { "eukaryome_bubbleplot.png"},
    content = function(file) {
      ggsave(file, bubplot_18S(), device = "png", width = 10, height = dynamicheightdownload18S(), limitsize = FALSE)
    })
  
  #Heatmap
  
  ofavproxy <- dataTableProxy('ofavtab') #proxy of the o. faveolata gene datatable that will allow us to manipulate it without changing the whole source table.
  output$ofavtab <- DT::renderDataTable({#base output for the o. faveolata gene datatable.
    DT::datatable(ofavdf, 
                  extensions = c("Scroller", "Buttons"),
                  style = "bootstrap4",
                  options = list(
                    searching = TRUE,
                    autoWidth = TRUE,
                    deferRender = TRUE,
                    scrollX = TRUE,
                    scrollY = "500px",
                    dom = 'l<"sep">B + <br> + frtip',
                    fillContainer = TRUE,
                    buttons = list(
                      extend = "csv",
                      text = "<i class='glyphicon glyphicon-download-alt'></i> Download this table as .csv"
                    ),
                    pageLength = 50,
                    lengthMenu = list(c(50,100,200, nrow(ofavdf)), c("50", "100", "200", "All")),
                    orderClasses=TRUE,
                    columnDefs = list(list(width = '406px', targets = c(1:4)))
                  )
    )  
  })
  
  dtreproxy <- dataTableProxy('dtretab') #proxy of the d. trenchii gene datatable that will allow us to manipulate it without changing the whole source table.
  output$dtretab <- DT::renderDataTable({#base output for the d. trenchii gene datatable.
    DT::datatable(dtredf, 
                  extensions = c("Scroller", "Buttons"),
                  style = "bootstrap4",
                  options = list(
                    searching = TRUE,
                    autoWidth = TRUE,
                    deferRender = TRUE,
                    scrollX = TRUE,
                    scrollY = "500px",
                    dom = 'l<"sep">B + <br> + frtip',
                    fillContainer = TRUE,
                    buttons = list(
                      extend = "csv",
                      text = "<i class='glyphicon glyphicon-download-alt'></i> Download this table as .csv"
                    ),
                    pageLength = 50,
                    lengthMenu = list(c(50,100,200, nrow(dtredf)), c("50", "100", "200", "All")),
                    orderClasses=TRUE,
                    columnDefs = list(list(width = '406px', targets = c(1:4)))
                  )
    )  
  })
  
  bminproxy <- dataTableProxy('bmintab') #proxy of the b. minutum gene datatable that will allow us to manipulate it without changing the whole source table.
  output$bmintab <- DT::renderDataTable({#base output for the b. minutum gene datatable.
    DT::datatable(bmindf, 
                  extensions = c("Scroller", "Buttons"),
                  style = "bootstrap4",
                  options = list(
                    searching = TRUE,
                    autoWidth = TRUE,
                    deferRender = TRUE,
                    scrollX = TRUE,
                    scrollY = "500px",
                    dom = 'l<"sep">B + <br> + frtip',
                    fillContainer = TRUE,
                    buttons = list(
                      extend = "csv",
                      text = "<i class='glyphicon glyphicon-download-alt'></i> Download this table as .csv"
                    ),
                    pageLength = 50,
                    lengthMenu = list(c(50,100,200, nrow(bmindf)), c("50", "100", "200", "All")),
                    orderClasses=TRUE,
                    columnDefs = list(list(width = '406px', targets = c(1:4)))
                  )
    )  
  })
  
  ofavadd <-  eventReactive(input$add, { #event triggered by the "add selected genes" button that adds to an in-between dataframe any o. faveolata selected gene.
    baseofavadd <- matrix(nrow=0, ncol = 4)
    baseofavadd <- data.frame(baseofavadd)
    colnames(baseofavadd) <- c("Gene ID", "Gene name", "Product description", "PFAMs")
    rowsofavadd <- ofavdf[input$ofavtab_rows_selected,]
    colnames(rowsofavadd) <- c("Gene ID", "Gene name", "Product description", "PFAMs")
    rowsofavadd["Organism"] <- rep("O. faveolata", nrow(rowsofavadd))
    if(length(input$ofavtab_rows_selected) == 0){ #if there are no rows selected, it just returns a blank dataframe. If there are 1 or more rows selected, they are added onto the in-between dataframe.
      return(baseofavadd)
    } else if(length(input$ofavtab_rows_selected) != 0) {
      return(rowsofavadd)
    }
  })
  
  dtreadd <-  eventReactive(input$add, {#event triggered by the "add selected genes" button that adds to an in-between dataframe any d. trenchii selected gene.
    basedtreadd <- matrix(nrow=0, ncol = 4)
    basedtreadd <- data.frame(basedtreadd)
    colnames(basedtreadd) <- c("Gene ID", "Gene name", "Product description", "PFAMs")
    rowsdtreadd <- dtredf[input$dtretab_rows_selected,]
    colnames(rowsdtreadd) <- c("Gene ID", "Gene name", "Product description", "PFAMs")
    rowsdtreadd["Organism"] <- rep("D. trenchii", nrow(rowsdtreadd))
    if(length(input$dtretab_rows_selected) == 0){
      return(basedtreadd)
    } else if(length(input$dtretab_rows_selected) != 0) {
      return(rowsdtreadd)
    }
  })
  
  bminadd <-  eventReactive(input$add, {#event triggered by the "add selected genes" button that adds to an in-between dataframe any b. minutum selected gene.
    basebminadd <- matrix(nrow=0, ncol = 4)
    basebminadd <- data.frame(basebminadd)
    colnames(basebminadd) <- c("Gene ID", "Gene name", "Product description", "PFAMs")
    rowsbminadd <- bmindf[input$bmintab_rows_selected,]
    colnames(rowsbminadd) <- c("Gene ID", "Gene name", "Product description", "PFAMs")
    rowsbminadd["Organism"] <- rep("B. minutum", nrow(rowsbminadd))
    if(length(input$bmintab_rows_selected) == 0){
      return(basebminadd)
    } else if(length(input$bmintab_rows_selected) != 0) {
      return(rowsbminadd)
    }
  })
  
  hmreactval <- reactiveValues() #generation of a reactivevalues item that will store in-between objects for the heatmap generation.
  selectedprovamatrix <- matrix(ncol=5, nrow=0) #generation of the selection table. First we create a matrix
  selectedprovadf <- data.frame(selectedprovamatrix) #generation of the selection table. The matrix is turned into a dataframe
  colnames(selectedprovadf) <- c("Gene ID", "Gene name", "Product description", "PFAMs", "Organism") #addition of colnames to the dataframe
  hmreactval$df <- selectedprovadf #we store the selection table dataframe in this pocket.
  ofavfulladd <- ofavdf #preparation of the dataframe that will be added in case someone uses the "Add all O. faveolata genes" button.
  ofavfulladd$Organism <- rep("O. faveolata", nrow(ofavdf))
  hmreactval$ofav <- ofavfulladd
  dtrefulladd <- dtredf #preparation of the dataframe that will be added in case someone uses the "Add all D. trenchii genes" button.
  dtrefulladd$Organism <- rep("D. trenchii", nrow(dtredf))
  hmreactval$dtre <- dtrefulladd
  bminfulladd <- bmindf #preparation of the dataframe that will be added in case someone uses the "Add all B. minutum genes" button.
  bminfulladd$Organism <- rep("B. minutum", nrow(bmindf))
  hmreactval$bmin <- bminfulladd
  
  observeEvent(input$add, { #triggered by the "Add selected genes" button, this event joins the products of the events ofavadd, dtreadd and bmin add, generating the selection table, stored in hmreactval$df.
    hmreactval$df <- rbind(hmreactval$df, ofavadd(), dtreadd(), bminadd())
    hmreactval$df <- unique(hmreactval$df) #any possible duplicate rows are deleted.
    row.names(hmreactval$df) <- 1:nrow(hmreactval$df) #the rownames are reset to be between 1 and the total number of rows in the selection table.
    ofavproxy %>% selectRows(NULL) #the selection made by the user in the gene datatables is reset using their proxies.
    dtreproxy %>% selectRows(NULL)
    bminproxy %>% selectRows(NULL)
  })
  
  observeEvent(input$ofavadd, {#event for "Add all O. faveolata genes" button
    hmreactval$df <- rbind(hmreactval$df, hmreactval$ofav)
    row.names(hmreactval$df) <- 1:nrow(hmreactval$df)
  })
  
  observeEvent(input$dtreadd, {#event for "Add all D. trenchii genes" button
    hmreactval$df <- rbind(hmreactval$df, hmreactval$dtre)
    row.names(hmreactval$df) <- 1:nrow(hmreactval$df)
  })
  
  observeEvent(input$bminadd, {#event for "Add all B. minutum genes" button
    hmreactval$df <- rbind(hmreactval$df, hmreactval$bmin)
    row.names(hmreactval$df) <- 1:nrow(hmreactval$df)
  })
  
  observeEvent(input$remove, {#event for the removal of selected rows in the selection table upon pressing the "Remove selected genes" button.
    if(is.null(input$selectedtab_rows_selected) == FALSE){#the event won't trigger if there are no selected rows in the selection table. 
      hmreactval$df <- hmreactval$df[-input$selectedtab_rows_selected,]#if there are selected rows, they are removed.
      if(nrow(hmreactval$df) != 0) {#in case there are still rows in the selection table once the "Remove all O. faveolata genes" button is pressed, the row names are reset to a 1:nrow logic.
        row.names(hmreactval$df) <- 1:nrow(hmreactval$df) #the condition is set because otherwise it may generate a fatal error.
      }}
  })
  
  observeEvent(input$ofavremove, {#event of removal of o. faveolata genes upon pressing the "Remove all O. faveolata genes" button.
    if(is.null(hmreactval$df) == FALSE){#the event only triggers if there is at least one row in the selection table, otherwise any attempt of row removal would produce an error
      hmreactval$df <- hmreactval$df[!hmreactval$df$Organism == "O. faveolata",] #removal of rows which have "O. faveolata" as value for the row "Organism".
      if (nrow(hmreactval$df) != 0){#in case there are still rows in the selection table once the "Remove selected genes" button is pressed, the row names are reset to a 1:nrow logic.
        row.names(hmreactval$df) <- 1:nrow(hmreactval$df)
      }}
  })
  
  observeEvent(input$dtreremove, {#event of removal of d. trenchii genes upon pressing the "Remove all D. trenchii genes" button.
    if(is.null(hmreactval$df) == FALSE){#the event only triggers if there is at least one row in the selection table, otherwise any attempt of row removal would produce an error
      hmreactval$df <- hmreactval$df[!hmreactval$df$Organism == "D. trenchii",] #removal of rows which have "D. trenchii" as value for the row "Organism".
      if(nrow(hmreactval$df) != 0){#in case there are still rows in the selection table once the "Remove all D. trenchii genes" button is pressed, the row names are reset to a 1:nrow logic.
        row.names(hmreactval$df) <- 1:nrow(hmreactval$df)
      }}
  })
  
  observeEvent(input$bminremove, {#event of removal of b. minutum genes upon pressing the "Remove all B. minutum genes" button.
    if(is.null(hmreactval$df) == FALSE){#the event only triggers if there is at least one row in the selection table, otherwise any attempt of row removal would produce an error
      hmreactval$df <- hmreactval$df[!hmreactval$df$Organism == "B. minutum",] #removal of rows which have "B. minutum" as value for the row "Organism".
      if(nrow(hmreactval$df) != 0){#in case there are still rows in the selection table once the "Remove B. minutum genes" button is pressed, the row names are reset to a 1:nrow logic.
        row.names(hmreactval$df) <- 1:nrow(hmreactval$df)
      }}
  })
  
  
  output$selectedtab <- DT::renderDataTable({#output for the selection table.
    DT::datatable(hmreactval$df, 
                  extensions = c("Scroller", "Buttons"),
                  style = "bootstrap4",
                  options = list(
                    searching = TRUE,
                    autoWidth = TRUE,
                    language = list(zeroRecords = "No records to display - Select rows from the tables above!"),
                    deferRender = TRUE,
                    scrollX = TRUE,
                    scrollY = "500px",
                    dom = 'l<"sep">B + <br> + frtip',
                    fillContainer = TRUE,
                    buttons = list(
                      extend = "csv",
                      text = "<i class='glyphicon glyphicon-download-alt'></i> Download this table as .csv"
                    ),
                    pageLength = 50,
                    lengthMenu = list(c(50,100,200, nrow(bmindf)), c("50", "100", "200", "All")),
                    orderClasses=TRUE,
                    columnDefs = list(list(width = '325px', targets = c(1:5))))) %>% formatStyle("Organism", fontStyle = "italic")
  })
  
  
  genes_to_use <- reactive({as.character(hmreactval$df$'Gene ID') #reactive value that returns a vector that contains the id of all the genes present in the selection table 
  }) 
  
  heatmaps_prep <- reactive({#reactive function made to ease the heatmap generation. It returns a list of 4 heatmaps, each for a certain timepoint, which use the id genes present in the genes_to_use vector
    list(
      hmt1 = DoHeatmap(ofav.combined.labeled[,ofav.combined.labeled.hm@meta.data$orig.ident=="Ofav Healthy-T"], features = genes_to_use(),label=T,cells=index$rowname[index$orig.ident=="Ofav Healthy-T"], draw.lines=T, lines.width=2, group.colors=colors_hm, angle=90, group.bar = T, combine = T) + scale_fill_gradientn(colors = c("blue", "#fffff2", "red")) + NoLegend() + theme(axis.title.x = element_text(size=25)) + xlab("Healthy"),
      hmt3 = DoHeatmap(ofav.combined.labeled[,ofav.combined.labeled.hm@meta.data$orig.ident=="Ofav Mid-Bleach 1-T"], features = genes_to_use(),label=T,cells=index$rowname[index$orig.ident=="Ofav Mid-Bleach 1-T"], draw.lines=T, lines.width=2, group.colors=colors_hm, angle=90, group.bar = T, combine = T) + NoLegend() + scale_fill_gradientn(colors = c("blue", "#fffff2", "red")) + theme(axis.text.y = element_blank(), axis.title.x = element_text(size=25)) + xlab("Mid-bleach 1"),
      hmt5 = DoHeatmap(ofav.combined.labeled[,ofav.combined.labeled.hm@meta.data$orig.ident=="Ofav-Mid-Bleach 2-T"], features = genes_to_use(),label=T,cells=index$rowname[index$orig.ident=="Ofav-Mid-Bleach 2-T"], draw.lines=T, lines.width=2, group.colors=colors_hm, angle=90, group.bar = T, combine = T) + NoLegend() + scale_fill_gradientn(colors = c("blue", "#fffff2", "red")) + theme(axis.text.y = element_blank(), axis.title.x = element_text(size=25)) + xlab("Mid-bleach 2"),
      hmt7 = DoHeatmap(ofav.combined.labeled[,ofav.combined.labeled.hm@meta.data$orig.ident=="Ofav-Bleached-T"], features = genes_to_use(),label=T,cells=index$rowname[index$orig.ident=="Ofav-Bleached-T"], draw.lines=T, lines.width=2, group.colors=colors_hm, angle=90, group.bar = T, combine = T) + scale_fill_gradientn(colors = c("blue", "#fffff2", "red")) + NoLegend() + theme(axis.text.y = element_blank(), axis.title.x = element_text(size=25)) + xlab("Bleached")
    )
  }) %>% bindCache(genes_to_use()) %>% bindEvent(input$loadhm) #cache is bound to reutilize the heatmaps in case they are requested again with the same genes_to_use and is bound to the "Load this heatmap" button to optimize
  
  hmreactval$prehm <- reactiveVal(NULL) #a reactivevalues object is created under hmreactval$prehm to allow the individual use of each heatmap listed in heatmaps_prep
  
  output$heatmap_ui <- renderUI({#reactive UI for the heatmap. Before the "Load this heatmap" button is pressed, a text is shown asking the user to load it. After that, a reactive box is set to contain the heatmap.
    if(input$loadhm > 0) {
      withSpinner(plotOutput("heatmap", height = "auto", width = "auto"), type = 6)
    } else {
      HTML(
        '<div style="display: flex; justify-content: center; align-items: center; height: auto; min-height: 200px;">',
        '<p>Load your scRNA-seq expression heatmap!</p>',
        '</div>'
      )}
  })
  
  observeEvent(input$loadhm, {#loading of the heatmap in the app after the pressing of the "Load this heatmap" button.
    if(nrow(hmreactval$df) != 0){#only triggered if there is at least one row in the selection table, otherwise does nothing.
      hmreactval$prehm(heatmaps_prep()) #heatmaps_prep() list of 4 heatmaps is requested in the variable hmreactval$prehm
      hmreactval$hmcomb <- arrangeGrob(grobs = list(hmreactval$prehm()$hmt1, hmreactval$prehm()$hmt3, hmreactval$prehm()$hmt5, hmreactval$prehm()$hmt7), ncol = 4, widths = c(1.20,1,1,1)) #this will be our heatmap, an arrangeGrob of the 4 heatmaps in the heatmap_pre list.
      hm_build <- ggplot_build(hmreactval$prehm()$hmt1)
      hm_yticks <- hm_build$layout$panel_params[[1]]$y.sec$breaks
      hm_num_ticks <- length(hm_yticks)
      dynamicheighthm <- hm_num_ticks * 5 + 1000 #value calculated from the number of genes used in the heatmap to ensure a good height display of it.
      output$heatmap <- renderPlot({#output of the heatmap.
        hmreactval$hmfinal <- grid.arrange(hmreactval$hmcomb) #grid.arrange of the arrangeGrob hmreactval$hmcomb.
        print(hmreactval$hmfinal)}, width = 1750, height = dynamicheighthm) #printing of the result of the grid.arrange.
    }}) 
  
  
  dynamicheightdownloadhm <- reactive({#the png image that is downloaded when the user wants to download the heatmap needs a variable height that depends on the number of genes used for its proper visualization.
    hmreactval$hmcomb <- arrangeGrob(grobs = list(hmreactval$prehm()$hmt1, hmreactval$prehm()$hmt3, hmreactval$prehm()$hmt5, hmreactval$prehm()$hmt7), ncol = 4, widths = c(1.20,1,1,1))
    hm_build <- ggplot_build(hmreactval$prehm()$hmt1)
    hm_yticks <- hm_build$layout$panel_params[[1]]$y.sec$breaks
    hm_num_ticks <- length(hm_yticks)
    return(hm_num_ticks + 10)
  })
  
  output$downloadHeatmap <- downloadHandler(#download of the heatmap png image after pressing the "Donwload this heatmap" button.
    filename = function() { "scRNA_seq_heatmap.png"},
    content = function(file) {
      ggsave(file, hmreactval$hmfinal, device = "png", width = 30, height = dynamicheightdownloadhm(),  limitsize = FALSE)
    })
}

shinyApp(ui, server)

