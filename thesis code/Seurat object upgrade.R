library("Seurat")
load("original_data/results_seurat_180426_Jane_1.Rdata")
so <- UpdateSeuratObject(cells)
saveRDS(so,"updated_data/D_Static_culture.Robj")

load("original_data/results_seurat_180426_Jane_2.Rdata")
so <- UpdateSeuratObject(cells)
saveRDS(so,"updated_data/C_Static_culture+drug_treatment.Robj")

load("original_data/results_seurat_180426_Jane_3.Rdata")
so <- UpdateSeuratObject(cells)
saveRDS(so,"updated_data/B_Shear_culture.Robj")

load("original_data/results_seurat_180426_Jane_4.Rdata")
so <- UpdateSeuratObject(cells)
saveRDS(so,"updated_data/A_Shear_culture+drug_treatment.Robj")


