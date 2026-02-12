
# load packages
library(readxl)
library(gplots)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(data.table)
library(ape)
library(phytools)
library(metafor)


# set working directory

#####################
##### PREP DATA #####
#####################

### effect size data & moderator variables
# data
df_dataset<-read_excel("Meta-analysis_Dataset.xlsx")

# change "NR" (Not Reported) to NA
df_dataset[df_dataset == "NR"] <- NA

# make sure R reads effect sizes as numeric
df_dataset$EffectSize <- as.numeric(df_dataset$EffectSize)
df_dataset$EffectSizeSE <- as.numeric(df_dataset$EffectSizeSE)
df_dataset$SampleSize <- as.numeric(df_dataset$SampleSize)
df_dataset$`95CI_lwr` <- as.numeric(df_dataset$`95CI_lwr`)
df_dataset$`95CI_upr` <- as.numeric(df_dataset$`95CI_upr`)

#calculate weighted variances from SE
# if SE is missing, calculate variance from confidence intervals
df_dataset$EffectSizeVariance <- ifelse(is.na(df_dataset$EffectSizeSE),
                                        sqrt(df_dataset$SampleSize)*((df_dataset$`95CI_upr`- df_dataset$`95CI_lwr`)/3.92)^2,
                                        df_dataset$EffectSizeSE^2)

# remove columns about literature (e.g., paper name)
df_studies <- df_dataset[,c(2, 5, 16:24,27:56)]

# create column for scientific name
df_studies$sci_name <- as.factor(paste(df_studies$Genus, df_studies$Species, sep = "_"))

# create table of response and precipitation metrics
data_vars <- data.frame(table(df_studies$ResponseCode, df_studies$PrecipCode))
data_vars <- data_vars[data_vars$Freq > 0,]
colnames(data_vars) <- c("ResponseVariable", "PrecipVariable", "N_Studies")
data_vars <- data_vars[data_vars$ResponseVariable != "individual",]

# run contingency test on rainfall and reproduction
con.test <- df_dataset[,c("Relationship","PrecipCode")]
con.test <- con.test[con.test$Relationship != "quad",]
con.test <- con.test[!is.na(con.test$PrecipCode),]
con.test$Relationship <- ifelse(con.test$Relationship == "neg", "Negative response", 
                   ifelse(con.test$Relationship == "none", "No Response","Positive response"))
con.test$PrecipCode <- ifelse(con.test$PrecipCode == "direct", "Precipitation during nesting period", 
                   ifelse(con.test$PrecipCode == "nesting", "Precipitation during breeding season","Precipitation before start of breeding season"))
con.test <- table(con.test)
chisq.test(con.test)
balloonplot(con.test, dotcol = "skyblue1", 
            main = "",
            xlab = "", ylab = "")
options(scipen = 999)


### categorize variables for analyses

# if species are known for being biparental or cooperative care, change to multi-parental
df_studies$ParentalCare[df_studies$ParentalCare %in% c("B", "C", "B/C")] <- "M"

# make sure latitude, longitude, and elevation are numeric
df_studies$Lat <- as.numeric(df_studies$Lat)
df_studies$Long <- as.numeric(df_studies$Long)
df_studies$Elevation <- as.numeric(df_studies$Elevation)

#scale continuous variables
df_studies$AbsLat <- round(abs(df_studies$Lat))
df_AbsLat_z <-data.frame(AbsLat = seq(0,78, by = 1), AbsLat_z = scale(seq(0,78, by = 1))[,1])
df_studies <- left_join(df_studies, df_AbsLat_z, by = "AbsLat")

df_Elevation_z <-data.frame(Elevation = seq(0,2600, by = 0.5), Elevation_z = scale(seq(0,2600, by = 0.5))[,1])
df_studies <- left_join(df_studies, df_Elevation_z, by = "Elevation")

df_AvgBodySize_z <-data.frame(AvgBodySize = seq(1,10000, by = 0.5), AvgBodySize_z = scale(seq(1,10000, by = 0.5))[,1])
df_studies <- left_join(df_studies, df_AvgBodySize_z, by = "AvgBodySize")

df_AvgClutchSize_z <-data.frame(AvgClutchSize = seq(1,13, by = 0.5), AvgClutchSize_z = scale(seq(1,13, by = 0.5))[,1])
df_studies <- left_join(df_studies, df_AvgClutchSize_z, by = "AvgClutchSize")

# if species eat a variety of vertebrates, condense to one "vertebrate diet" variable
df_studies$Diet_Vertebrates <- ifelse(df_studies$Diet_Amphibians == "y" | df_studies$Diet_Birds == "y" | df_studies$Diet_Carrion == "y" | df_studies$Diet_Fish == "y" | df_studies$Diet_Mammals == "y" | df_studies$Diet_Reptiles == "y", "y","n")

# categorize nest types into three bins: minimal, structured, covered
df_studies$NestCode <- ifelse(df_studies$NestType %in% c("scrape", "platform"), "Minimal",
                              ifelse(df_studies$NestType %in% c("cup", "bowl"), "Structured", "Protected"))

# remove few studies with quadratic effect
df_studies <- df_studies[df_studies$Relationship != "quad",]
df_studies$Relationship <- ifelse(df_studies$Relationship == "pos", 1,
                                  ifelse(df_studies$Relationship == "neg", -1, 0))

# change diet to just specialist or generalist
df_studies$Diet <-paste(df_studies$Diet_Inverts, df_studies$Diet_Plant, df_studies$Diet_Fruit, df_studies$Diet_Vertebrates, sep = "")
df_studies$Diet <- ifelse(str_count(df_studies$Diet, "y") <= 1 , "Specialist", "Generalist")


### create table of species included in quantitative analyses (Table S1)
df_table <- df_studies 
df_table$Diet <- paste(df_table$Diet_Inverts, df_table$Diet_Plant, df_table$Diet_Fruit, df_table$Diet_Vertebrates, sep = "")
df_table$Diet <- ifelse(str_count(df_table$Diet, "y") <= 1 , "Specialist", "Generalist")
df_table <- df_table[,c("Family","Genus","Species","Common Name","GroundNest", 
                        "NestCode", "NGDevelopment", "ParentalCare","AvgBodySize",
                        "Diet", "AvgClutchSize")]

df_table <- df_table[!duplicated(df_table),]
df_table$GroundNest <- ifelse(df_table$GroundNest == "y", "Yes", "No")
df_table$NGDevelopment <- ifelse(df_table$NGDevelopment == "P", "Precocial", "Altricial")
df_table$ParentalCare <- ifelse(df_table$ParentalCare == "U", "Uniparental", 
                                ifelse(df_table$ParentalCare == "B", "Biparental",  
                                       ifelse(df_table$ParentalCare == "C", "Cooperative", "Biparental/Cooperative")))

write.csv(df_table, "SpeciesInQuantAnalyses.csv")


# SUMMARY STATISTICS ------------------------------------------------------

### number of studies in qualitative synthesis
df_qual <- df_dataset[,2:5]
df_qual <- df_qual[!duplicated(df_qual),]


### number of studies in vote counting analyses
df_votecounting <- df_dataset[!is.na(df_dataset$Relationship),]
df_votecounting <- df_votecounting[df_votecounting$Relationship != "quad",]
df_votecounting <- df_votecounting[,2:5]
df_votecounting <- df_votecounting[!duplicated(df_votecounting),]
length(df_votecounting$Authors)

### number of studies in phylogenetic synthesis
df_phylosynth <- df_dataset[!is.na(df_dataset$EffectSizeVariance),]
df_phylosynth <- df_phylosynth[,2:5]
df_phylosynth <- df_phylosynth[!duplicated(df_phylosynth),]


### species summaries
df_spp <- df_studies[, c("Family","Common Name", "AvgBodySize", "GroundNest","NestCode", "NGDevelopment", "ParentalCare", "Diet_Inverts", "AvgClutchSize")]
df_spp <- df_spp[!duplicated(df_spp),]

table(df_dataset$PrecipCode)

length(unique(df_spp$Family))
length(unique(df_spp$`Common Name`))
min(df_spp$AvgBodySize)
max(df_spp$AvgBodySize)
quantile(df_spp$AvgBodySize, 0.75)
quantile(df_spp$AvgBodySize, 0.25)
table(df_spp$GroundNest)
table(df_spp$NestCode)
table(df_spp$NGDevelopment)
table(df_spp$ParentalCare)
nrow(df_spp[df_spp$Diet_Inverts == "y",])
mean(df_spp$AvgClutchSize)
range(df_spp$AvgClutchSize)


### site summaries
df_sites <- df_studies[, c("Lat", "Long", "Elevation", "BiomeCode")]
df_sites <- df_sites[!duplicated(df_sites[,c("Lat","Long"),]),]
df_sites$AbsLat <- abs(df_sites$Lat)
range(df_sites$Elevation, na.rm = TRUE)
mean(df_sites$Elevation, na.rm = TRUE)
range(df_sites$AbsLat, na.rm = TRUE)
mean(df_sites$AbsLat, na.rm = TRUE)
table(df_sites$BiomeCode)


### plot responses to precipitation from full dataset (FIGURE 2)

df_temp_plot <- df_dataset %>% drop_na(EffectSize) %>% drop_na(PrecipCode)
df_temp_plot$scaled_b <- scale(df_temp_plot$EffectSize)
df_temp_plot <- df_temp_plot[df_temp_plot$scaled_b < 3.29 & df_temp_plot$scaled_b > -3.29,]

df_bar <- data.frame(table(df_dataset$Relationship, df_dataset$PrecipCode))
colnames(df_bar) <- c("Relationship", "PrecipCode", "Freq")

#reorder factor levels
df_bar$Relationship <- factor(df_bar$Relationship, levels = c("quad", "pos", "none", "neg"))
df_bar$PrecipCode <- factor(df_bar$PrecipCode, levels = c("direct", "nesting", "lagged"))
df_bar <- df_bar[df_bar$Relationship != "quad",]

ggplot(df_bar, aes(fill=Relationship, y=Freq, x=PrecipCode)) + 
  geom_bar(position="fill", stat="identity", alpha = 0.7) +
  scale_fill_manual(values=c("palegreen4","ivory3","indianred4"), labels = c("Positive", "None", "Negative")) +
  scale_x_discrete(labels = c("During \n nesting period", "Throughout \n breeding season", "Before \n breeding season")) +
  xlab("") +
  ylab("Frequency") +
  theme_pubr() +
  theme(legend.position = "right") 

ggsave("Figure_2.tiff", units = "in",
       width = 9, height = 6, 
       dpi = 300, compression = "lzw")

# VOTE COUNTING ANALYSES --------------------------------------------------


#### ALL RESPONSES ####

### correlations between predictors
df_corr <- df_studies[,c("BiomeCode","Elevation","NestCode","GroundNest","NGDevelopment","ParentalCare","AvgBodySize","Diet", "AvgClutchSize")]
df_corr$BiomeCode <- rank(df_corr$BiomeCode)
df_corr$Elevation <- rank(df_corr$Elevation)
df_corr$NestCode <- rank(df_corr$NestCode)
df_corr$GroundNest <- rank(df_corr$GroundNest)
df_corr$NGDevelopment <- rank(df_corr$NGDevelopment)
df_corr$ParentalCare <- rank(df_corr$ParentalCare)
df_corr$AvgBodySize <- rank(df_corr$AvgBodySize)
df_corr$Diet <- rank(df_corr$Diet)
df_corr <- df_corr[!duplicated(df_corr),]

corrplot(cor(df_corr), method = "number", type = "upper")


#adjusting reference levels to elucidate relationships
df_studies$BiomeCode <- relevel(factor(df_studies$BiomeCode), ref = "wooded")
df_studies$PrecipCode <- relevel(factor(df_studies$PrecipCode), ref = "nesting")

#fit model
mod_global <- glm(Relationship ~ PrecipCode + BiomeCode + AbsLat_z * Elevation_z + AvgBodySize_z + GroundNest + NestCode + ParentalCare  + Diet + AvgClutchSize_z,
                  data = df_studies, family = gaussian)

summary(mod_global)

# extract model coefficients and st. errors
betas_vc <- data.frame(summary(mod_global)$coefficients)
betas_vc$lwr <- betas_vc$Estimate - 1.96*betas_vc$Std..Error
betas_vc$upr <- betas_vc$Estimate + 1.96*betas_vc$Std..Error
betas_vc$var <- rownames(betas_vc)
betas_vc <- betas_vc[-1,]


#fix predictor names
betas_vc$Predictor <- c("Precip. timing: direct","Precip. timing: lagged", 
                        "Biome: open", "Biome: aquatic","Biome: urban", 
                        "Latitude", "Elevation", 
                        "Body size", "Ground nest", 
                        "Nest shape: protected", "Nest shape: structured",
                        "Parental care: uniparental",
                        "Diet: specialist", "Clutch size",
                        "Latitude*Elevation")

betas_vc$PredictorType <- ifelse(betas_vc$Predictor %in% c("Precip. timing: lagged", "Precip. timing: direct"),1,
                                 ifelse(betas_vc$Predictor %in% c("Latitude", "Latitude*Elevation", "Elevation",
                                                                     "Biome: urban", "Biome: open", "Biome: aquatic"),2, 3)) 

#reorder variables to match Figure 1 (hypotheses plot) in manuscript
betas_vc$PlotOrder <- c(1,2,7,6,8,4,3,14,12,10,11,13,9,15,5)

# ggplot of beta estimates
plot_vc_betas <- 
ggplot() +
  theme_pubr() +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_vline(xintercept =  0, col = "gray50", linetype = 2) +
  geom_linerange(data = betas_vc, aes(x = Estimate, y = reorder(Predictor,-PlotOrder), xmin = lwr, xmax = upr, color = as.factor(PredictorType)), size = 1, alpha = 0.9) +
  geom_point(data = betas_vc, aes(x = Estimate, y = reorder(Predictor,-PlotOrder), color = as.factor(PredictorType)),  size = 3, alpha = 0.9) +
  annotate("text", x = 2, y = c("Precip. timing: lagged", "Elevation"), label = "**", size = 8) +
  annotate("text", x = 2, y = c("Latitude*Elevation"), label = "*", size = 8) +
  scale_x_continuous(limits = c(-2.6, 2.6))+
  scale_color_manual(values = c("skyblue3", "darkseagreen3","khaki3"))  +   xlab("Parameter estimate") 


#### DAILY RESPONSE ####

df_daily <- df_studies[df_studies$ResponseCode == "daily",c(7:14,16,25:40,43:50)]
df_daily$PrecipCode <- relevel(factor(df_daily$PrecipCode), ref = "nesting")
df_daily$BiomeCode <- relevel(factor(df_daily$BiomeCode), ref = "wooded")


### correlations between predictors
df_corr <- df_daily[,c("BiomeCode","Elevation","NestCode","GroundNest","NGDevelopment","ParentalCare","AvgBodySize","Diet", "AvgClutchSize")]
df_corr$BiomeCode <- rank(df_corr$BiomeCode)
df_corr$Elevation <- rank(df_corr$Elevation)
df_corr$NestCode <- rank(df_corr$NestCode)
df_corr$GroundNest <- rank(df_corr$GroundNest)
df_corr$NGDevelopment <- rank(df_corr$NGDevelopment)
df_corr$ParentalCare <- rank(df_corr$ParentalCare)
df_corr$AvgBodySize <- rank(df_corr$AvgBodySize)
df_corr$Diet <- rank(df_corr$Diet)
df_corr <- df_corr[!duplicated(df_corr),]

corrplot(cor(df_corr), method = "number", type = "upper")
#NG Development & AvgBodySize highly correlated with parental care and nest code (>0.60) -- remove from analyses

#fit model
mod_daily <- glm(Relationship ~ PrecipCode + BiomeCode + AbsLat_z * Elevation_z + GroundNest + NestCode + ParentalCare  + Diet + AvgClutchSize_z,
                 data = df_daily, family = gaussian)

summary(mod_daily)

# extract model coefficients and st. errors
temp <- data.frame(summary(mod_daily)$coefficients)[,1:2]
temp$lwr <- temp$Estimate - 1.96*temp$Std..Error
temp$upr <- temp$Estimate + 1.96*temp$Std..Error
round(temp, digits = 3)



# NEST LEVEL RESPONSE -----------------------------------------------------

df_nest <- df_studies[df_studies$ResponseCode == "nest",c(7:14,16,25:40,43:50)]
df_nest$PrecipCode <- relevel(factor(df_nest$PrecipCode), ref = "nesting")
df_nest$BiomeCode <- relevel(factor(df_nest$BiomeCode), ref = "wooded")


### correlations between predictors
df_corr <- df_nest[,c("BiomeCode","Elevation","NestCode","GroundNest","NGDevelopment","ParentalCare","AvgBodySize","Diet", "AvgClutchSize")]
df_corr$BiomeCode <- rank(df_corr$BiomeCode)
df_corr$Elevation <- rank(df_corr$Elevation)
df_corr$NestCode <- rank(df_corr$NestCode)
df_corr$GroundNest <- rank(df_corr$GroundNest)
df_corr$NGDevelopment <- rank(df_corr$NGDevelopment)
df_corr$ParentalCare <- rank(df_corr$ParentalCare)
df_corr$AvgBodySize <- rank(df_corr$AvgBodySize)
df_corr$Diet <- rank(df_corr$Diet)
df_corr <- df_corr[!duplicated(df_corr),]

corrplot(cor(df_corr), method = "number", type = "upper")
# NG development == correlation of 0.67, removed


mod_nest <- glm(Relationship ~ PrecipCode + BiomeCode + AbsLat_z * Elevation_z +  GroundNest + NestCode + ParentalCare + AvgBodySize_z + Diet + AvgClutchSize_z,
                data = df_nest, family = gaussian)

summary(mod_nest)

# extract model coefficients and st. errors
temp <- data.frame(summary(mod_nest)$coefficients)[,1:2]
temp$lwr <- temp$Estimate - 1.96*temp$Std..Error
temp$upr <- temp$Estimate + 1.96*temp$Std..Error
round(temp, digits = 3)


# PRODUCTIVITY RESPONSE -----------------------------------------------------

df_prod <- df_studies[df_studies$ResponseCode == "productivity",c(7:14,16,25:40,43:50)]
df_prod$PrecipCode <- relevel(factor(df_prod$PrecipCode), ref = "nesting")
df_prod$BiomeCode <- relevel(factor(df_prod$BiomeCode), ref = "wooded")


### correlations between predictors
df_corr <- df_prod[,c("BiomeCode","Elevation","NestCode","GroundNest","NGDevelopment","ParentalCare","AvgBodySize","Diet", "AvgClutchSize")]
df_corr$BiomeCode <- rank(df_corr$BiomeCode)
df_corr$Elevation <- rank(df_corr$Elevation)
df_corr$NestCode <- rank(df_corr$NestCode)
df_corr$GroundNest <- rank(df_corr$GroundNest)
df_corr$NGDevelopment <- rank(df_corr$NGDevelopment)
df_corr$ParentalCare <- rank(df_corr$ParentalCare)
df_corr$AvgBodySize <- rank(df_corr$AvgBodySize)
df_corr$Diet <- rank(df_corr$Diet)
df_corr <- df_corr[!duplicated(df_corr),]

corrplot(cor(df_corr), method = "number", type = "upper")
#nest code highly correlated (-0.80) with avg body size -- body size removed from analyses

mod_prod <- glm(Relationship ~ PrecipCode + BiomeCode + AbsLat_z * Elevation_z + GroundNest + NestCode + NGDevelopment + ParentalCare + Diet + AvgClutchSize_z,
                data = df_prod, family = gaussian)

summary(mod_prod)

# extract model coefficients and st. errors
temp <- data.frame(summary(mod_prod)$coefficients)[,1:2]
temp$lwr <- temp$Estimate - 1.96*temp$Std..Error
temp$upr <- temp$Estimate + 1.96*temp$Std..Error
round(temp, digits = 3)


# META-ANALYSES -----------------------------------------------------------

### only include studies with effect sizes, SEs, and sample sizes
df_meta <- df_studies
df_meta <- df_meta[!duplicated(df_meta),]
df_meta <- df_meta[!is.na(df_meta$EffectSizeSE),]
df_meta <- df_meta[!is.na(df_meta$SampleSize),]
df_meta <- df_meta[df_meta$EffectSizeSE != 0,]
df_meta <- df_meta[df_meta$EffectSizeVariance < 10,] #a few studies with very large SEs

# remove studies with effect sizes and response is "failure" (i.e. a binary response)
df_meta <- df_meta[!df_meta$ResponseVariable %like% "fail", ]
df_meta <- df_meta[!df_meta$ResponseVariable %like% "loss", ]

# create variable for absolute latitude
df_meta$AbsLat <- abs(df_meta$Lat)

# adjust scientific names to be recognized by tree
df_meta$sci_name <- as.character(df_meta$sci_name)
df_meta$sci_name[df_meta$sci_name == "Ammospiza_maritima"] <- "Ammodramus_maritimus"
df_meta$sci_name[df_meta$sci_name == "Anser_caerulescens"] <- "Chen_caerulescens"
df_meta$sci_name[df_meta$sci_name == "Antigone_antigone"] <- "Grus_antigone"
df_meta$sci_name[df_meta$sci_name == "Ardea_alba"] <- "Casmerodius_albus"
df_meta$sci_name[df_meta$sci_name == "Bucorvus_leadbeateri"] <- "Bucorvus_cafer"
df_meta$sci_name[df_meta$sci_name == "Charadrius_nivosus"] <- "Charadrius_alexandrinus"
df_meta$sci_name[df_meta$sci_name == "Dendrocoptes_medius"] <- "Dendrocopos_medius"
df_meta$sci_name[df_meta$sci_name == "Parkesia_motacilla"] <- "Seiurus_motacilla"
df_meta$sci_name[df_meta$sci_name == "Sternula_antillarum"] <- "Sterna_antillarum"
df_meta$sci_name[df_meta$sci_name == "Leuconotopicus_borealis"] <- "Picoides_borealis"
df_meta$sci_name[df_meta$sci_name == "Leiothlypis_celata"] <- "Vermivora_celata"
df_meta$sci_name[df_meta$sci_name == "Iduna_caligata"] <- "Hippolais_caligata"
df_meta$sci_name[df_meta$sci_name == "Poecile_gambeli"] <- "Parus_gambeli"
df_meta$sci_name[df_meta$sci_name == "Artemisiospiza_belli"] <- "Amphispiza_belli"
df_meta$species.id.phy <- df_meta$sci_name

length(unique(df_meta$`Common Name`))
length(unique(df_meta$Authors))

#plot all species in meta-analyses

# # of species represented from each family
df_fams <- df_meta[,c("Family","Genus","Species")]
df_fams <- df_fams[!duplicated(df_fams),]
df_fams <- df_fams[order(df_fams$Family),]

col.pal<-colorRampPalette(c("thistle3", "darkslategray"))

df_summ_fams <- data.frame(table(df_fams$Family))
df_summ_fams <- arrange(df_summ_fams,desc(Freq))

spp_plot <- ggplot(data = df_fams, aes(x = reorder(Family, Family, function(x)+length(x)))) + 
  geom_bar(fill = "darkslategray4", stat="count") +
  xlab("") +
  ylab("Number of species") +
  theme_minimal() + 
  theme(legend.position = "none") +
  coord_flip()

df_fams$sci <- paste(df_fams$Genus, df_fams$Species)

# load trees
tree<-read.nexus("all_species_quant_analyses.nex")

tree_plot <- 
  plot.phylo(tree[[1]], type = "phylogram", cex = 0.6, label.offset = 2, no.margin = TRUE)
  
  
ggarrange(spp_plot, tree_plot, nrow = 1, ncol = 2,
          labels = c("A","B"), align = "v")


# all species -------------------------------------------------------------

all_meta <- df_meta

### correlations between predictors
df_corr <- all_meta[,c("BiomeCode","Elevation_z","AbsLat_z","NestCode","GroundNest","NGDevelopment","ParentalCare","AvgBodySize_z","Diet", "AvgClutchSize_z")]
df_corr$BiomeCode <- rank(df_corr$BiomeCode)
df_corr$NestCode <- rank(df_corr$NestCode)
df_corr$GroundNest <- rank(df_corr$GroundNest)
df_corr$NGDevelopment <- rank(df_corr$NGDevelopment)
df_corr$ParentalCare <- rank(df_corr$ParentalCare)
df_corr$Diet <- rank(df_corr$Diet)
df_corr <- df_corr[!duplicated(df_corr),]

corrplot(cor(df_corr), method = "number", type = "upper")
# nestling development and parental care correlated (0.66). remove ng development


### prep phylogenetic tree data

# load trees
tree<-read.nexus("all_spp.nex")

tree_1 <- compute.brlen(tree$tree_5474)
A <- vcv(tree_1, corr=TRUE)


### fit multilevel phylogenetic meta-analytic model

#calculate variance
all_meta <- escalc(measure = "MN", yi = EffectSize, sei = EffectSizeSE, ni = SampleSize, data = all_meta)

# make forest plot
## forest plot info: https://wviechtb.github.io/metafor/reference/forest.rma.html; https://s4be.cochrane.org/blog/2016/07/11/tutorial-read-forest-plot/
forest(phylo_mod_all, addpred = TRUE, shade = "zebra")

# make funnel plot
par(mar=c(5, 4, 4, 2) + 0.5)
funnel(phylo_mod_all, yaxis = "seinv")

# fit Egger's regression test for funnel plot symmetry
regtest(phylo_mod_all$yi, phylo_mod_all$vi)

### loop through all phylogenetic trees

#save model outputs in a list
all_mods_meta <- list()
betas_meta <- list()

#loop through all 1,000 phylogenetic trees, estimate betas

for(i in data.frame(names(tree))[,1]){
  
  tree_1 <- compute.brlen(tree[[i]]) #save tree
  A <- vcv(tree_1, corr=TRUE) #create vcv matrix for all species
  
  # fit model
  mod<-rma.mv(yi, vi,
              mods = ~PrecipCode + BiomeCode  + AbsLat_z * Elevation_z + NestCode + ParentalCare + GroundNest + AvgBodySize_z + Diet + AvgClutchSize_z,
              random = list(~ 1 | species.id.phy),
              R=list(species.id.phy=A), data=all_meta, slab = paste(all_meta$Common.Name, ", ", all_meta$PrecipCode, paste(" (",all_meta$Authors, " ", all_meta$Year, ")", sep = ""), sep = ""))
  
  #save model outputs
  all_mods_meta[[i]] <- mod
  
  betas_meta <- rbind(betas_meta, 
                       data.frame(var = rownames(mod$beta), 
                                  beta = mod$beta,
                                  se = mod$se,
                                  p = mod$pval))
}

avg_beta_meta <- betas_meta %>%
  group_by(var) %>%
  summarize(beta = mean(beta), se = mean(se), p = mean(p))

avg_beta_meta$lwr <- avg_beta_meta$beta - 1.96*avg_beta_meta$se
avg_beta_meta$upr <- avg_beta_meta$beta + 1.96*avg_beta_meta$se
avg_beta_meta$sig <- ifelse(avg_beta_meta$p <= 0.001, "***",
                            ifelse(avg_beta_meta$p <= 0.01, "**", 
                                   ifelse(avg_beta_meta$p <= 0.05, "*", ".")))

betas_meta <- betas_meta[betas_meta$var != "intrcpt",]
avg_beta_meta <- avg_beta_meta[-16,]

#fix predictor names
avg_beta_meta$Predictor <- c("Latitude", "Latitude*Elevation",
                             "Body size", "Clutch size",
                             "Biome: aquatic","Biome: urban", "Biome: open", 
                             "Diet: specialist", "Elevation", "Ground nest", 
                             "Nest shape: protected", "Nest shape: structured",
                             "Parental care: uniparental", 
                             "Precip. timing: direct", "Precip. timing: lagged")


avg_beta_meta$PredictorType <- ifelse(avg_beta_meta$Predictor %in% c("Precip. timing: lagged", "Precip. timing: direct"),1,
                                 ifelse(avg_beta_meta$Predictor %in% c("Latitude", "Latitude*Elevation", "Elevation",
                                                                  "Biome: urban", "Biome: open", "Biome: aquatic"),2, 3)) 

#reorder variables to match predictions plot in manuscript
avg_beta_meta$PlotOrder <- c(4,5,14,15,6,8,7,9,3,12,10,11,13,1,2)  

# ggplot of beta estimates
plot_meta_betas <- 
  ggplot() +
  theme_pubr() +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_vline(xintercept =  0, col = "gray50", linetype = 2) +
  geom_linerange(data = avg_beta_meta, aes(x = beta, y = reorder(Predictor,-PlotOrder), xmin = lwr, xmax = upr, color = as.factor(PredictorType)), size = 1, alpha = 0.9) +
  geom_point(data = avg_beta_meta, aes(x = beta, y = reorder(Predictor,-PlotOrder), color = as.factor(PredictorType)),  size = 3, alpha = 0.9) +
  annotate("text", x = 2, y = c("Precip. timing: lagged","Precip. timing: direct", "Latitude", "Latitude*Elevation"), 
           label = "***", size = 8) +
  scale_color_manual(values = c("skyblue3", "darkseagreen3","khaki3"))  +   xlab("Parameter estimate") +
    xlab("Parameter estimate") +
  scale_x_continuous(limits = c(-2.6, 2.6))


ggarrange(plot_vc_betas, plot_meta_betas, ncol = 2, nrow = 1, align = "h", widths = c(1.3, 1), labels = c("A", "B"))

ggsave("Figure_3.tiff", units = "in",
       width = 9, height = 6, 
       dpi = 300, compression = "lzw")

# daily response ----------------------------------------------------------

# save data frame for species with daily nest survival analyses
daily_meta <- df_meta[which(df_meta$ResponseCode == "daily"),]

### correlations between predictors
df_corr <- daily_meta[,c("BiomeCode","Elevation","NestCode","GroundNest","NGDevelopment","ParentalCare","AvgBodySize","Diet", "AvgClutchSize")]
df_corr$BiomeCode <- rank(df_corr$BiomeCode)
df_corr$Elevation <- rank(df_corr$Elevation)
df_corr$NestCode <- rank(df_corr$NestCode)
df_corr$GroundNest <- rank(df_corr$GroundNest)
df_corr$NGDevelopment <- rank(df_corr$NGDevelopment)
df_corr$ParentalCare <- rank(df_corr$ParentalCare)
df_corr$AvgBodySize <- rank(df_corr$AvgBodySize)
df_corr$Diet <- rank(df_corr$Diet)
df_corr <- df_corr[!duplicated(df_corr),]

corrplot(cor(df_corr), method = "number", type = "upper")
#removed NGDevelopment and Parental care

### load phylogenetic tree data
# load trees
tree<-read.nexus("daily_1000trees.nex")

# save trees as separate .txt files
for(i in data.frame(names(tree))[,1]){
  ape::write.tree(tree[[i]], file = paste(i,"_daily.txt",sep = ""))}
tree_1 <- compute.brlen(tree$tree_1039)
A <- vcv(tree_1, corr=TRUE)


### fit multilevel phylogenetic meta-analytic model

#calculate variance
daily_meta <- escalc(measure = "MN", yi = EffectSize, sei = EffectSizeSE, ni = SampleSize, data = daily_meta)

# make forest plot
## forest plot info: https://wviechtb.github.io/metafor/reference/forest.rma.html; https://s4be.cochrane.org/blog/2016/07/11/tutorial-read-forest-plot/
forest(phylo_mod_d, addpred = TRUE, shade = "zebra")


### loop through all phylogenetic trees

#save model outputs in a list
all_mods_daily <- list()
betas_daily <- list()

#loop through all 1,000 phylogenetic trees, estimate betas
for(i in data.frame(names(tree))[,1]){
  
  tree_1 <- compute.brlen(tree[[i]]) #save tree
  A <- vcv(tree_1, corr=TRUE) #create vcv matrix for all species
  
  # fit model
  mod<-rma.mv(yi, vi,
              mods = ~PrecipCode + BiomeCode + AbsLat_z + Elevation_z + AbsLat_z * Elevation_z + NestCode  + GroundNest + AvgBodySize_z + Diet + AvgClutchSize_z,
              random = list(~ 1 | species.id.phy),
              R=list(species.id.phy=A), data=daily_meta, slab = paste(daily_meta$Common.Name, ", ", daily_meta$PrecipCode, paste(" (",daily_meta$Authors, " ", daily_meta$Year, ")", sep = ""), sep = ""))
  
  #save model outputs
  all_mods_daily[[i]] <- mod
  
  betas_daily <- rbind(betas_daily, 
                      data.frame(var = rownames(mod$beta), 
                                 beta = mod$beta,
                                 se = mod$se,
                                 p = mod$pval))
}

avg_beta_daily <- betas_daily %>%
  group_by(var) %>%
  summarize(beta = mean(beta), se = mean(se), p = mean(p))

avg_beta_daily$lwr <- avg_beta_daily$beta - 1.96*avg_beta_daily$se
avg_beta_daily$upr <- avg_beta_daily$beta + 1.96*avg_beta_daily$se
avg_beta_daily$sig <- ifelse(avg_beta_daily$p <= 0.001, "***",
                            ifelse(avg_beta_daily$p <= 0.01, "**", 
                                   ifelse(avg_beta_daily$p <= 0.05, "*", ".")))

round(avg_beta_daily[,2:5], digits = 3)


# nest response ----------------------------------------------------------

# save data frame for species with daily nest survival analyses
nest_meta <- df_meta[which(df_meta$ResponseCode == "nest"),]

### correlations between predictors
df_corr <- nest_meta[,c("BiomeCode","Elevation","AbsLat", "NestCode","GroundNest","NGDevelopment","ParentalCare","AvgBodySize","Diet", "AvgClutchSize")]
df_corr$BiomeCode <- rank(df_corr$BiomeCode)
df_corr$NestCode <- rank(df_corr$NestCode)
df_corr$GroundNest <- rank(df_corr$GroundNest)
df_corr$NGDevelopment <- rank(df_corr$NGDevelopment)
df_corr$ParentalCare <- rank(df_corr$ParentalCare)
df_corr$Diet <- rank(df_corr$Diet)
df_corr <- df_corr[!duplicated(df_corr),]

corrplot(cor(df_corr), method = "number", type = "upper")
#remove nestling development and body size b/c correlations, remove parental care because no uniparental care in dataset

length(unique(nest_meta$`Common Name`))
length(unique(nest_meta$Authors))

### tree data
# load trees
tree<-read.nexus("nest_1000trees.nex")

# save trees as separate .txt files
for(i in data.frame(names(tree))[,1]){
  ape::write.tree(tree[[i]], file = paste(i,"_nest.txt",sep = ""))}
tree_1 <- compute.brlen(tree$tree_7048)
A <- vcv(tree_1, corr=TRUE)
setdiff(tree$tree_0000$tip.label, nest_meta$sci_name)


### fit multilevel phylogenetic meta-analytic model

#calculate variance
nest_meta <- escalc(measure = "MN", yi = EffectSize, sei = EffectSizeSE, ni = SampleSize, data = nest_meta)


### loop through all phylogenetic trees

#save model outputs in a list
all_mods_nest <- list()
betas_nest <- list()


#loop through all 1,000 phylogenetic trees, estimate betas

for(i in data.frame(names(tree))[,1]){
  
  tree_1 <- compute.brlen(tree[[i]]) #save tree
  A <- vcv(tree_1, corr=TRUE) #create vcv matrix for all species

  # fit model
  mod<-rma.mv(yi, vi,
              mods = ~PrecipCode + BiomeCode + AbsLat_z + Elevation_z + AbsLat_z * Elevation_z + NestCode + GroundNest + Diet + AvgClutchSize_z,
              random = list(~ 1 | species.id.phy),
              R=list(species.id.phy=A), data=nest_meta, slab = paste(nest_meta$Common.Name, ", ", nest_meta$PrecipCode, paste(" (",nest_meta$Authors, " ", nest_meta$Year, ")", sep = ""), sep = ""))
  
  #save model outputs
  all_mods_nest[[i]] <- mod

  betas_nest <- rbind(betas_nest, 
                      data.frame(var = rownames(mod$beta), 
                                beta = mod$beta,
                                se = mod$se,
                                p = mod$pval))
}

avg_beta_nest <- betas_nest %>%
  group_by(var) %>%
  summarize(beta = mean(beta), se = mean(se), p = mean(p))

avg_beta_nest$lwr <- avg_beta_nest$beta - 1.96*avg_beta_nest$se
avg_beta_nest$upr <- avg_beta_nest$beta + 1.96*avg_beta_nest$se
avg_beta_nest$sig <- ifelse(avg_beta_nest$p <= 0.001, "***",
                             ifelse(avg_beta_nest$p <= 0.01, "**", 
                                    ifelse(avg_beta_nest$p <= 0.05, "*", ".")))



# Lat*Elev plots ----------------------------------------------------------


### make matrices over which to predict responses

df_elev_20 <- expand.grid(PrecipCodelagged = 0,
                          PrecipCodedirect = 1,
                          BiomeCodeopen = 1,
                          BiomeCodeaquatic = 0,
                          BiomeCodedeveloped = 0,
                          AbsLat_z = df_AbsLat_z[df_AbsLat_z$AbsLat == 20,2],
                          Elevation_z = df_Elevation_z[seq(1,4001,by = 100),2],
                          GroundNesty = 1,
                          NestCodeProtected = 0,
                          NestCodeStructured = 1,
                          ParentalCareU = 0,
                          AvgBodySize_z = -1.716678,
                          DietSpecialist = 0,
                          AvgClutchSize_z = -0.8152395)

df_elev_20$'AbsLat_z:Elevation_z' <- df_elev_20$AbsLat_z*df_elev_20$Elevation_z


df_elev_40 <- expand.grid(PrecipCodelagged = 0,
                          PrecipCodedirect = 1,
                          BiomeCodeopen = 1,
                          BiomeCodeaquatic = 0,
                          BiomeCodedeveloped = 0,
                          AbsLat_z = df_AbsLat_z[df_AbsLat_z$AbsLat == 40,2],
                          Elevation_z = df_Elevation_z[seq(1,4001,by = 100),2],
                          GroundNesty = 1,
                          NestCodeProtected = 0,
                          NestCodeStructured = 1,
                          ParentalCareU = 0,
                          AvgBodySize_z = -1.716678,
                          DietSpecialist = 0,
                          AvgClutchSize_z = -0.8152395)

df_elev_40$'AbsLat_z:Elevation_z' <- df_elev_40$AbsLat_z*df_elev_40$Elevation_z


df_elev_60 <- expand.grid(PrecipCodelagged = 0,
                          PrecipCodedirect = 1,
                          BiomeCodeopen = 1,
                          BiomeCodeaquatic = 0,
                          BiomeCodedeveloped = 0,
                          AbsLat_z = df_AbsLat_z[df_AbsLat_z$AbsLat == 60,2],
                          Elevation_z = df_Elevation_z[seq(1,4001,by = 100),2],
                          GroundNesty = 1,
                          NestCodeProtected = 0,
                          NestCodeStructured = 1,
                          ParentalCareU = 0,
                          AvgBodySize_z = -1.716678,
                          DietSpecialist = 0,
                          AvgClutchSize_z = -0.8152395)

df_elev_60$'AbsLat_z:Elevation_z' <- df_elev_60$AbsLat_z*df_elev_60$Elevation_z


### predict responses

preds_elev_a20 <- data.frame()

for(i in data.frame(names(all_mods_meta))[,1]){
  
  # select model
  temp_mod <- all_mods_meta[[i]]
  
  # predict values
  preds_elev_a20 <- rbind(preds_elev_a20, data.frame(var = names(all_mods_meta[i]), elev = df_Elevation_z[seq(1,4001,by = 100),1], predict(temp_mod, as.matrix(df_elev_20))))
}

preds_elev_a40 <- data.frame()

for(i in data.frame(names(all_mods_meta))[,1]){
  
  # select model
  temp_mod <- all_mods_meta[[i]]
  
  # predict values
  preds_elev_a40 <- rbind(preds_elev_a40, data.frame(var = names(all_mods_meta[i]), elev = df_Elevation_z[seq(1,4001,by = 100),1], predict(temp_mod, as.matrix(df_elev_40))))
}

preds_elev_a60 <- data.frame()

for(i in data.frame(names(all_mods_meta))[,1]){
  
  # select model
  temp_mod <- all_mods_meta[[i]]
  
  # predict values
  preds_elev_a60 <- rbind(preds_elev_a60, data.frame(var = names(all_mods_meta[i]), elev = df_Elevation_z[seq(1,4001,by = 100),1], predict(temp_mod, as.matrix(df_elev_60))))
}


avg_preds_elev_a20 <- preds_elev_a20 %>%
  group_by(elev) %>%
  summarize(beta = mean(pred), se = mean(se), AbsLat = 20)

avg_preds_elev_a40 <- preds_elev_a40 %>%
  group_by(elev) %>%
  summarize(beta = mean(pred), se = mean(se), AbsLat = 40)

avg_preds_elev_a60 <- preds_elev_a60 %>%
  group_by(elev) %>%
  summarize(beta = mean(pred), se = mean(se), AbsLat = 60)


# save all predictions together
avg_preds_elev_a <- rbind(avg_preds_elev_a20, avg_preds_elev_a40, avg_preds_elev_a60)


#create color palette for plots
col.pal<-colorRampPalette(c("coral1","deepskyblue4"))

plot_meta_a<-
  ggplot(data = avg_preds_elev_a) +
  theme_pubr() + 
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank(), legend.position = "right",
        plot.margin = unit(c(1, .25, 2, .25), "cm")) + 
  ggtitle("") + labs(y="Predicted response",x="Elevation") + 
  geom_line(aes(x = elev, y = beta, color = as.factor(AbsLat)), size = 1.5) +
  geom_ribbon(aes(x = elev, ymin = beta - (1.96*se), max = beta + (1.96*se), fill = as.factor(AbsLat)), alpha = 0.1) +
  scale_color_manual(values = (col.pal(3)), name = "Latitude") +
  scale_fill_manual(values = (col.pal(3)), name = "Latitude") +
  geom_hline(yintercept =  0, col = "gray50", linetype = 2)  +
  scale_y_continuous(limits = c(-8, 8), breaks = c(-8, -4, 0, 4, 8))



### plot distribution of lat/elev in studies

plot_data <- ggplot() +
  theme_pubr() +
  geom_point(data = all_meta, aes(x = Elevation, y = AbsLat, color = AbsLat), alpha = 0.8, size = 4) +
  xlab("Elevation") +
  ylab("Latitude") +
  scale_color_gradient(low = "coral1", high = "deepskyblue4", name = "Latitude")


# FIGURE 4

ggarrange(plot_data, plot_meta_a, ncol = 2, nrow = 1, align = "hv", 
          labels = c("A","B"), legend = "none")

ggsave("Figure_4.tiff", units = "in",
       width = 9, height = 6, 
       dpi = 300, compression = "lzw")

#
