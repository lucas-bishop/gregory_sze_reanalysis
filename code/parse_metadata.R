library(tidyverse)

files <- read_tsv("data/mothur_output/stability.paired.files",
                  col_names = c("filename", "r1_fastq", "r2_fastq"))


shared <- read.delim("data/mothur_output/final.0.03.subsample.shared", 
                     header = T, stringsAsFactors = F) %>% 
  select(-label, -numOtus)

# Build first meta data file
temp_meta <- files %>% filter(!grepl("PCRwater", filename) & !grepl("empty", filename)) %>% select(-r2_fastq) %>% 
  separate(col = r1_fastq,
          into = c("sampleid", "kit", "storage", "samplenum", "x1", "x2", "x3", "x4"), sep = "_") %>% select(-contains("x"))

#Need to parse what each storage and kit abbreviation mean
#need to look at plate map to see if any can be excluded
#Need to look at the samplenum key that Marc gave me to make sure everything matches 

actual_meta <- mutate(temp_meta, 
                      # This grabs the ID for samples that links it to the shared
                      sample_id = shared$Group[!grepl("DM", shared$Group) & !grepl("DS", shared$Group)], 
                      amp_cycles = ifelse(grepl("x", temp_meta$V1) == TRUE, temp_meta$V1, "30x"), 
                      bead_beat = ifelse(grepl("min", temp_meta$V2) == TRUE, temp_meta$V2, 
                                         ifelse(grepl("x", temp_meta$V1) == TRUE & 
                                                  grepl("DA", temp_meta$V5) == FALSE, NA, "10min")), 
                      taq = ifelse(V2 %in% c("PHU", "PL", "ACC", "Q5", "k"), temp_meta$V2, 
                                   ifelse(grepl("PHU", temp_meta$V1) == TRUE, temp_meta$V1, "ACC")), 
                      ext_kit = ifelse(grepl("PMM", temp_meta$V3) == TRUE, "PMM", 
                                       ifelse(grepl("x", temp_meta$V1) == TRUE | 
                                                grepl("Zmock", temp_meta$V1) == TRUE, NA, 
                                              ifelse(grepl("PHU", temp_meta$V1) == TRUE, NA, 
                                                     ifelse(V1 %in% c("PMM", "PMS", "Z"), temp_meta$V1, 
                                                            ifelse(grepl("S", temp_meta$V4) == TRUE, "PMS", "PMM"))))), 
                      collection_kit = ifelse(V3 %in% c("Omni", "Zymo", "Frozen"), temp_meta$V3, 
                                              ifelse(V4 %in% c("Omni", "Zymo", "Frozen"), temp_meta$V4, 
                                                     ifelse(grepl("frozens", temp_meta$V2) == TRUE, "Frozen", 
                                                            ifelse(V4 %in% "FS", "Frozen", NA)))), 
                      sample_type = ifelse(grepl("DA", temp_meta$V2) == TRUE | 
                                             grepl("DA", temp_meta$V3) == TRUE | 
                                             grepl("DA", temp_meta$V4) == TRUE | 
                                             grepl("DA", temp_meta$V5) == TRUE, "Stool", 
                                           ifelse(grepl("ZymoSTD", temp_meta$V2) == TRUE, "artificial_com", "DNA")), 
                      stool_id = ifelse(grepl("DA", temp_meta$V2) == TRUE, temp_meta$V2, 
                                        ifelse(grepl("DA", temp_meta$V3) == TRUE, temp_meta$V3, 
                                               ifelse(grepl("DA", temp_meta$V4) == TRUE, temp_meta$V4, 
                                                      ifelse(grepl("DA", temp_meta$V5) ==TRUE, temp_meta$V5, NA))))) %>% 
  select(sample_id, amp_cycles, bead_beat, taq, ext_kit, collection_kit, 
         sample_type, stool_id)


write.csv(actual_meta, "data/process/final.metadata.csv", row.names = F)
