
/Information about count matrix and metadata


/24M     ./counts.rds
/8.5K    ./metadata.rds

/dim(counts)

/#[1] 11544   291

/colnames(metadata)


/#[1] "specimenID"         "source"             "cases.submitter_id" "sample_type"        "figo_stage"        

/#[6] "primary_diagnosis"  "HPV_type"           "HPV_clade"


/dim(metadata)

/#[1] 290   8



/table(metadata$source)


/#GTEX TCGA 

/#19  271


/table(metadata$sample_type)

/#Primary Tumor  Solid Tissue Normal 

/#268                  22
