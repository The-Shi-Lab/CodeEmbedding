# Analysis Pipeline to Generate Code Embeddings -- an example using MGI data

## adapted codes from Xu Shi, Xianshi Yu, and Kuan-Han Wu ##
## 05/05/22 ##

# STEP 1 #
A data prepossessing step to convert original EHR data into "cooccur" format
- Step1_EHR2LongNumFormat.R
- Convert EHR data (diagnostics, procedures, lab) to long format with 3 columns (PId, numDays, and CId)
- All numeric variables and sorted by PId then numDays

# STEP 2 #
- Step2_EHR2CoOccurMatrix.py
- Calculate cooccurence matrix from long format EHR table
- ex: python /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step2_EHR2CoOccurMatrix.py -i /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/LongFormat_Num/codeRecord.csv -o /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/Dx_CoOccurMatrix.csv -w 0 1 6 13
- To calculate the cooccurence matrix by chunks (add arguments -c and -tc)
- ex: python /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step2_EHR2CoOccurMatrix.py -i /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/LongFormat_Num/codeRecord.csv -o /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/Dx_CoOccurMatrix.csv -tc 20 -c 2 -w 0 1 6 13

# STEP 3 #
- Step3_Convert2SparseMatrix.R
- Merge all chunks into one matrix
- ex: Rscript /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step3_Convert2SparseMatrix.R
- In the future, add DataSource parameters

# STEP 4 #
- Step4_CreateEmbeddings.R
- Generate embeddings for 500 vectors
- ex: Rscript /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step4_CreateEmbeddings.R
- In the future, add DataSource and iter parameters

# STEP 5 #
- Calculate cosine similarity and use phecode to obtain AUC
- Phecode/ICD crosswalk and code to calculate cosine similary can be found in Step5_CosineSimilarity_short.R
