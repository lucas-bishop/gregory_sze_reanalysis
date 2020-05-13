library(tidyverse)

files <- read_tsv("data/mothur_output/stability.paired.files",
                  col_names = c("filename", "r1_fastq", "r2_fastq"))


shared <- read.delim("data/mothur_output/final.0.03.subsample.shared", 
                     header = T, stringsAsFactors = F) %>% 
  select(-label, -numOtus) %>% filter(!grepl("PCRwater", Group) & !grepl("empty", Group) & !grepl("bad", Group))

# Build first metadata file 
temp_meta <- files %>% filter(!grepl("PCRwater", filename) & !grepl("empty", filename) & !grepl("bad", filename)) %>% separate(col = r1_fastq,
          into = c("kit", "storage", "stool_id", "x1", "x2", "x3", "x4"), sep = "_", remove = FALSE) %>% 
  select(filename, kit, storage, stool_id) %>% rename(Group = filename)
  

#Need to parse what each storage and kit abbreviation mean
# Kits: Zymobiomics (z), PowerMag Microbiome (pm), PowerSoil (ps)
# Storage: Room Temperature (rt), Zymo DNA/RNA Shield (zdr), 
#         frozen only (f), Omnigene Genotek (og)

metadata <- temp_meta %>% mutate(kit = str_replace_all(string = temp_meta$kit, 
                                                       pattern = c("z" = "Zymobiomics","pm" = "PowerMag", "ps" = "PowerSoil")),
                                 storage = str_replace_all(string = temp_meta$storage,
                                                           pattern = c("rt" = "RoomTemp", "zdr" = "Zymo_DNARNA_Shield", "f" = "Frozen", "og" = "Omnigene_Genotek")))
#change all external standards to labeled clearly as controls
metadata$storage[metadata$storage == 'ext'] <- "standard"
metadata$stool_id[metadata$storage == 'standard'] <- "control"

#need to look at plate map to see if any can be excluded
#Need to look at the samplenum key that Marc gave me to make sure everything matches 


write.csv(metadata, "data/final.metadata.csv", row.names = F)
