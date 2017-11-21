#################
# PRELIMINARIES #
#################

rm(list=ls(all=T))

options(stringsAsFactors = F)
# Sys.setlocale("LC_ALL", "de_DE")

# install / load packages
sapply(c("dplyr", "ggplot2", "gridExtra"), function(x) 
  if(!is.element(x, installed.packages())) install.packages(x, dependencies = T))

lapply(c("dplyr", "ggplot2", "gridExtra"), function(x) require(x, character.only = T))


#############
# READ DATA #
#############

# read files
fnhd <- read.csv("ung_fnhd.csv", sep="\t", head=T, quote="",
                 encoding = "UTF-8")
gc <- read.csv("ung_gc01.csv", sep="\t", head=T, quote="",
               encoding="UTF-8", fill = T)
dta <- read.csv("ungbaby5.csv",  sep="\t", head=T, quote="",
                encoding = "UTF-8")


# in one table ------------------------------------------------------
fnhd <- mutate(fnhd, CORPUS="MzENHG")
gc <- mutate(gc, CORPUS="GerManC")
dta <- mutate(dta, CORPUS="DTA")

# add "core years" to DTA
dta_core_years <- c(1625, 1675, 1725, 1775, 1825, 1875)

dta$Core_Year <- NA

for(i in 1:length(levels(factor(dta$Period)))) {
  dta[which(dta$Period==levels(factor(dta$Period))[i]),]$Core_Year <- dta_core_years[i]
}

# add "number" to DTA
dta$Numerus <- NA
dta[grep(".*ungen", dta$Key),]$Numerus <- "PL"
dta[grep(".*ungen", dta$Key, invert = T),]$Numerus <- "SG"

# add PREP to DTA
dta$Prep_V_ung_simple <- ifelse(dta$Prep_V_ung=="", "NEIN", "JA")


# one df for all
ung <- data.frame(LEMMA=c(fnhd$Lemma, gc$Lemma, dta$Lemma),
           DET=c(fnhd$Det, gc$Det, dta$Determiner),
           PREP=c(fnhd$als_praep_Kompl, gc$als_praep_Kompl, dta$Prep_V_ung_simple),
           PREP_VALUE=c(fnhd$Prep, gc$Praeposition, dta$Prep_V_ung),
           MOD=c(fnhd$Adj, gc$Adj, dta$Adj),
           NUM=c(fnhd$Numerus, gc$Numerus, dta$Numerus),
           YEAR=c(fnhd$Jahr, gc$Jahrzehnt, dta$Core_Year),
           PREF=c(fnhd$Praefix, gc$Praefix, dta$Prefix),
           CORPUS=c(fnhd$CORPUS, gc$CORPUS, dta$CORPUS),
           GEN=c(fnhd$Genitiv, gc$Genitiv, dta$Genitive))

ung$DET_SIMPLE <- ifelse(ung$DET %in% c("", "KEIN"), "n", "y")
ung$PREF_SIMPLE <- ifelse(is.na(ung$PREF)|ung$PREF==""|ung$PREF=="komplBV", "NO", "YES")
ung$MOD_SIMPLE <- NA

for(i in 1:nrow(ung)) {
  if(is.na(ung$MOD[i])) {
    ung$MOD_SIMPLE[i] <- "NEIN"
  } else if (ung$MOD[i]=="") {
    ung$MOD_SIMPLE[i] <- "NEIN"
  } else {
    ung$MOD_SIMPLE[i] <- "JA"
  }
  
  # print(i)
  
}

ung$GEN <- ifelse(is.na(ung$GEN), "KEIN", ung$GEN)
ung$GEN_SIMPLE <- ifelse(ung$GEN %in% c("", "KEIN"), "NEIN", "JA")

# unify spelling (as all umlauts have been replaced in the MzENHG data)
ung$MOD <- gsub("ä", "ae", ung$MOD)
ung$MOD <- gsub("ö", "oe", ung$MOD)
ung$MOD <- gsub("ü", "ue", ung$MOD)
ung$MOD <- gsub("ß", "ass", ung$MOD)
ung$MOD <- gsub("Ä", "Ae", ung$MOD)
ung$MOD <- gsub("Ö", "Oe", ung$MOD)
ung$MOD <- gsub("Ü", "Ue", ung$MOD)
ung$LEMMA <- gsub("ä", "ae", ung$LEMMA)
ung$LEMMA <- gsub("ö", "oe", ung$LEMMA)
ung$LEMMA <- gsub("ü", "ue", ung$LEMMA)
ung$LEMMA <- gsub("Ä", "Ae", ung$LEMMA)
ung$LEMMA <- gsub("Ö", "Oe", ung$LEMMA)
ung$LEMMA <- gsub("Ü", "Ue", ung$LEMMA)
ung$LEMMA <- gsub("ß", "ass", ung$LEMMA)

# add century
for(i in 1:nrow(ung)) {
  ung$CENTURY[i] <- as.numeric(as.character(paste(unlist(strsplit(as.character(ung$YEAR[i]), ""))[1:2], collapse="")))+1
}



############################
# MULTIDIMENSIONAL SCALING #
############################

# add Lemma+century
ung$CHRONLEMMA <- sapply(1:nrow(ung), function(i) paste(ung$LEMMA[i], "_", ung$CENTURY[i], sep=""))

# create input for MDS
mds_input <- ung %>% group_by(CHRONLEMMA) %>% summarise(
  prep = sum(PREP=="JA") / n(),
  det = sum(DET_SIMPLE=="y") / n(),
  pl = sum(NUM=="PL") / n(),
  mod = sum(MOD_SIMPLE=="JA") / n(),
  gen = sum(GEN_SIMPLE=="JA") / n()
)

# tibble to dataframe
mds_input <- as.data.frame(mds_input)
rownames(mds_input) <- mds_input$CHRONLEMMA
mds_input <- mds_input[,2:length(mds_input)]

# convert to dist; cmdscale
mds_output <- cmdscale(dist(as.matrix(mds_input)))
mds_output <- as.data.frame(mds_output)
mds_output$Lemma <- gsub("_1.", "", rownames(mds_output))

# add Freq column
mds_output$Freq <- sapply(1:nrow(mds_output), function(i) length(which(ung$LEMMA==mds_output$Lemma[i] & ung$CENTURY==as.numeric(gsub("\\D", "", rownames(mds_output)[i], perl=T)))))

# add RelFreq column
mds_output$RelFreq <- sapply(1:nrow(mds_output), function(i) length(which(ung$LEMMA==mds_output$Lemma[i] & ung$CENTURY==as.numeric(gsub("\\D", "", rownames(mds_output)[i], perl=T)))) /
                            length(which(ung$CENTURY==as.numeric(gsub("\\D", "", rownames(mds_output)[i], perl = T)))))


# create one dataframe for each mds output
for(i in c("_16", "_17", "_18", "_19")) {
  assign(paste("p", i, sep="", collapse=""),
         ggplot(mds_output[grep(i, rownames(mds_output)),], aes(x = V1, y = V2, label = Lemma, size=Freq)) + geom_text() + theme_bw() +
           scale_size_continuous(guide = FALSE) +
           ylab("Dim 2") + xlab("Dim 1") + ggtitle(paste(gsub("_", "", i), ". Jahrhundert", sep="", collapse="")) +
           theme(plot.title = element_text(hjust=0.5, face="bold")))
}

# plot mds output
# png("ung_mds_DE3.png", width=20, height=20, un="in", res=300)
grid.arrange(p_16, p_17, p_18, p_19, ncol=2, nrow=2)
# dev.off()



# MOTION CHART #

# googlevis: complex vs. simplex types
library(googleVis)
mds_output$Century <- as.numeric(gsub(".*_", "", rownames(mds_output)))*100
mds_output$Pref <- NA
for(i in 1:nrow(mds_output)) {
  temp <- filter(ung, LEMMA==gsub("_.*", "", rownames(mds_output)[i]))$PREF[1]
  mds_output$Pref[i] <- ifelse(is.na(temp), "SIMPLE", "PREFIXED")
}

# as factor
mds_output$Pref <- factor(mds_output$Pref)

# plot all items
mc <- gvisMotionChart(data = mds_output, id = "Lemma", xvar = "V1", yvar = "V2",
                timevar = "Century", colorvar = "Pref", sizevar = "RelFreq",
                options = list(state='{"showTrails":false};'))

# plot only items with Freq>5
mc2 <- gvisMotionChart(data = filter(mds_output, Freq>5),
                       id = "Lemma", xvar = "V1", yvar = "V2",
                       timevar = "Century", colorvar = "Pref", sizevar = "Freq")

# output plot
plot(mc)
plot(mc2)
