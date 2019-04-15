## March 12th

## Replace michels and cox IDs with appropriate GEO accessions

library(readxl)
library(tidyr)
library(plyr)
library(dplyr)

# read in sample keys
michels <- read.csv('Z:/Victor/Projects/DNAm-Ethnicity-Predictor/data/Michels/DeIDs_to_GEOIDs.csv') %>% as_tibble()
cox <- read_xlsx('Z:/Victor/Projects/DNAm-Ethnicity-Predictor/data/Cox/FromBCox_Sample to GEO ID.xlsx')

# read in current GEO accession (note that this is deleted after this script is ran
geo <- read_xls('Z:/Victor/Publications/DNAm Ethnicity Predictor/GEO submission/Samples metadata.xls',
                skip = 22, n_max = 509)

# join based on sample ID column

#rename columns to join on
cox <- cox %>% 
  rename(`Sample name` = Sample,
         `characteristics: GSM Accessions` = `GEO ID`) %>% 
  mutate(`Sample name` = gsub('_', '', `Sample name`)) %>%
  select(-Phenotype)
sum( geo$`Sample name` %in% cox$`Sample name`)

michels <- michels %>% 
  rename(`Sample name` = DeID,
         `characteristics: GSM Accessions` = `GSM`) %>%
  select(-NumID, -NewDeID)

mic <- geo %>% slice(grep('DeID', `Sample name`)) %>% pull(`Sample name`)
all(mic %in% michels$`Sample name`)# T

geo <- geo %>% 
  mutate(`characteristics: GSM Accessions` = ifelse(`characteristics: GSM Accessions` == 'NA', NA, 
                                                    `characteristics: GSM Accessions`)) %>%
  left_join(bind_rows(michels, cox), by = 'Sample name') %>%
  mutate(`characteristics: GSM Accessions` =
            coalesce(`characteristics: GSM Accessions.x`, `characteristics: GSM Accessions.y`)) %>%
  select(-`characteristics: GSM Accessions.x`, -`characteristics: GSM Accessions.y`)

write.table(geo, 'Z:/Victor/Publications/DNAm Ethnicity Predictor/GEO submission/Samples metadata2.txt',
            quote = F, sep = '\t', row.names = F)
## After this replace the Samples metadata.xlsx with the updated GSM column. Deleted version 2.


## March 13th
# Two issues with the data:
# Cox IDs were linking to the matched expression data
# [2] The Processed_Matrix.txt file contains the following data columns:
#
# PM139_r2
# PM139
#
# However, the Metadata worksheet lists:
#
#   PM139_vc_r3
#   PM139_vc_r2

#Read processed matrix file
mat <- read.table('Z:/Victor/Publications/DNAm Ethnicity Predictor/GEO submission/Processed_Matrix.txt')
colnames(mat)

geo <- read_xls('Z:/Victor/Publications/DNAm Ethnicity Predictor/GEO submission/Samples metadata.xls',
                skip = 22, n_max = 509)


#rearrange mat to match geo
data.frame(colnames(mat), geo$`Sample name`)
mat2 <- mat[,c(grep('^P[LM]', colnames(mat)), grep('C[OP]X', colnames(mat)), grep('GSM194', colnames(mat)),
       grep('^[PQ][0-9]', colnames(mat)), grep('DeID', colnames(mat)))]
data.frame(colnames(mat2), geo$`Sample name`)

# fix to match completely
colnames(mat2) <- gsub('COX_', 'COX', gsub('CPX', 'COX', colnames(mat2)))
colnames(mat2)[grepl('^P[ML][0-9]+$', colnames(mat2), perl = T)] <- # add vc
  paste0(colnames(mat2)[grepl('^P[ML][0-9]+$', colnames(mat2), perl = T)], '_vc')
colnames(mat2)[colnames(mat2)=='PM139_r2'] <- 'PM139_vc_r2'
data.frame(colnames(mat2), geo$`Sample name`)

# pm139 AND pm139_r2 in the geo submission have same idat file paths, I need to fix this
# load in the rgset to see which one i used 
pDat_train <- readRDS('../../Robjects_final/01_pDat.rds')
pDat_train %>% mutate(Case_ID = gsub('_r.*', '', sampleNames)) %>%
  group_by(Case_ID) %>% filter(n() > 1) %>%
  arrange(Case_ID, desc(meanSScor)) %>%
  select(sampleNames, Plate, Array, Slide, meanSScor)
# I used R2 to train (9977525013_R04C01)

# load in the data from omni analysis to see which one I used
pDat_C6_final <- readRDS('../../Robjects_final/04_final_pData_C6.rds')
pDat_C6_final %>% select(sampleNames, Sample_Plate, Sentrix_ID, Sentrix_Position)

# looks like I also used R2 for the omni section: 9977525013_R04C01
# So basically I have a duplicate sample in the processed data nad the metadata
# I remove one sample from the processed data here, and for the metadata I do that manually in excel
mat2 <- mat2[,-which(colnames(mat2) == 'PM139_vc')]
dim(mat2) # 508

# all should be match now
all(colnames(mat2) %in% geo$`Sample name`[which(geo$`Sample name` != 'PM139_vc_r3')])


# last thing to do is to update the cox gsm accessions again to the methylation GSM IDs
cox_k <- read_xlsx('../../data/Cox/2019-03-13 From Katie - Sample to GEO ID.xlsx')
cox_k <- cox_k %>%
  rename(`Sample name` = Sample,
         `characteristics: GSM Accessions` = `Meth GEO ID`) %>% 
  mutate(`Sample name` = gsub('_', '', `Sample name`)) %>%
  select(-Phenotype, -`GE GEO ID`)

geo <- geo %>% 
  filter(`Sample name` != 'PM139_vc_r3') %>% 
  mutate(`characteristics: GSM Accessions` = ifelse(grepl('COX', `Sample name`), NA, 
                                                    `characteristics: GSM Accessions`)) %>%
  left_join(cox_k, by = 'Sample name') %>% 
  mutate(`characteristics: GSM Accessions` = coalesce(`characteristics: GSM Accessions.x`, 
                                                      `characteristics: GSM Accessions.y`)) %>%
  select(-`characteristics: GSM Accessions.x`, -`characteristics: GSM Accessions.y`) 


# now rearrange processed matrix to match geo
all(geo$`Sample name` %in% colnames(mat2)) #T
all(colnames(mat2) %in% geo$`Sample name`) #T
mat2 <- mat2[,geo$`Sample name`]

dim(mat2);dim(geo)

# save metadata and processed data
write.table(geo, 'Z:/Victor/Publications/DNAm Ethnicity Predictor/GEO submission/Samples metadata2.txt',
            quote = F, sep = '\t', row.names = F)
write.table(mat2, 'Z:/Victor/Publications/DNAm Ethnicity Predictor/GEO submission/Processed_Matrix.txt',
            quote = F, sep = '\t')
