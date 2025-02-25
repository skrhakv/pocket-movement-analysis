# CryptoSite:
# structures <- c('3CHEA', '2AKAA', '2GFCA', '1ALBA', '1NEPA', '3MN9A', '1ALVA', '2YQCA', '2QFOB', '1RTCA', '1RDWX', '1TQOA', '3PUWE', '3L7UC', '1R1WA', '1G4EB', '1G4EB', '3F74C', '1MY1C', '2WGQB', '1DUBD', '1PZTA', '1XMGB', '1IMFA', '2AX9A', '1EXMA', '3KQAB', '2BF3A', '1CLLA', '3H5RA', '1NI6D', '1QLWB', '1HAGE', '1OK8A', '3DXNA', '1K5HC', '1HKAA', '2BU8A', '2OHGA', '2BLSB', '1RHBA', '1ADEA', '1KS9A', '1H09A', '3BL9B', '1BP5A', '2Q8FA', '2IYTA', '2IYTA', '1EX6A', '1RRGA', '1ECJD', '2CM2A', '1NUWA', '3CJ0A', '3CJ0A', '3CJ0A', '3CJ0A', '1UK2A', '2WGBA', '1HOOB', '1FA9A', '1W50A', '3B7DE', '1PKLB', '1PKLB', '1PKLB', '1K3FB', '1FVRA', '3HQDA', '2ZB1A', '3PEOG', '2BRKA', '2H4EB', '2CGAB', '1XCGB', '1FXXA', '4AKEB', '3NNUA', '1A8IA', '1SWXA', '2QLRC', '2F6VA', '1G24D', '1SU4A', '2AIRH', '1E2XA', '1MY0B', '1ZAHB', '4HB2C', '1BSQA', '2GPOA', '3FDLA', '1B6BA', '1KZ7D', '1JWPA', '3GXDB', '1BNCB', '1Z92A', '1JBUH') 

# PocketMiner:
# apo:
# structures <- c('3PPNA', '4W51A', '4I92A', '6E5DA', '5H9AA', '4IC4A', '1EZMA', '1S2OA', '1TVQA', '2OY4A', '5NIAA', '1Y1AA', '3NX1A', '1KX9B', '5NZMA', '6YPKA', '2CEYA', '1URPA', '3FVJA', '3UGKA', '2ZKUB', '5ZA4A', '4V38A', '3P53A', '2HQ8B', '2W9TA', '4R72A', '3QXWB', '3KJEA', '2LAOA', '4P0IB', '1KMOA', '1J8FC', '5UXAA', '6HB0A', '2FJYA', '5G1MA', '6RVMC')
# holo:
# structures <- c('3PPRB', '4W58A', '4I94B', '6E5FA', '6E5LA', '4INQA', '3DBKA', '1U2SA', '1TW4A', '3DPFA', '5NI6A', '1Y1AB', '3NX2A', '1N8VB', '2OEGA', '5OSZA', '6H76A', '2DRIA', '2B03A', '3UH1A', '3VQSA', '5YYBA', '4V3BA', '6I11A', '2W9SC', '4R74D', '3QXVA', '3KJGA', '1LAHE', '5OTAA', '1KMPA', '6QCNB', '5IGYA', '6HBDA', '2P70A', '5G3RA', '5XDTA')
# REMARK: '2HPSA' from holo set throws an error, therefore it was temporarily removed from the set
# highly rigid structures:
# structures <- c('2FD7A', '4HJKA', '1IGDA', '4TQLA', '2ALPA', '1AMMA', '1OFVA', '1QYSA', '5BVLA')
# REMARK: '4TQLB' changed to '4TQLA'
library(frustratometeR)

ResultsDir <- "/home/vit/Projects/flexibility-analysis/data/features/frustration/rigid-dataset/configurational"
PdbFiles <- "/home/vit/Projects/flexibility-analysis/src/feature-extraction/frustration/modified-pdbs"
# dir.create(file.path(ResultsDir))
PdbFilesList <- list.files(PdbFiles, pattern = "*.pdb", full.names = TRUE)
for (PdbFile in PdbFilesList) {
    pdbId <- substr(basename(PdbFile), 1, 4)
    chainId <- substr(basename(PdbFile), 5, 5)
    tryCatch({
        split_name <- strsplit(basename(PdbFile), "\\.")[[1]][1]
        if (dir.exists(file.path(ResultsDir, paste0(split_name, '_', chainId, '.done')))) {
            message(paste("Skip", split_name, "..."))
            next
        }
        Pdb_conf <- calculate_frustration(PdbFile = PdbFile, Chain = chainId, Mode = "configurational", ResultsDir = ResultsDir, Graphics = FALSE)
    }, error = function(e) {
        message(paste("Error processing", PdbFile, ":", e$message))
    })
}

# for (structure in structures) {
#     pdbId <- substr(structure,1,4)
#     chainId <- substr(structure,5,5)
#     PdbPath <- paste(PdbFiles, tolower(pdbId), chainId, '_clean_h.pdb', sep='')
#     Pdb_conf <- calculate_frustration(PdbFile = "/content/P01116-2.pdb", Chain = chainId, Mode = "configurational", ResultsDir = ResultsDir, Graphics = FALSE)
#     #Pdb_conf <- calculate_frustration(PdbFile = \"/content/P01116-2.pdb\", Chain=\"A\", Mode = 'mutational', ResultsDir = \"/content/Results/\")
# }
# 