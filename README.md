# Analysis Pipeline to Generate Code Embeddings -- an example using MGI data

## Adapted codes from Xu Shi, Xianshi Yu, and Kuan-Han Wu ##

# STEP 1 #
A data prepossessing step to convert original EHR data into "cooccur" long number format
- Step1_EHR2LongNumFormat.R
- Convert EHR data (diagnostics, procedures, lab) to long format with 3 columns (PId, numDays, and CId)
- All numeric variables and sorted by PId then numDays

# STEP 2 #
Use long number format EHR data to calculate cooccurence matrix
- Step2_EHR2CoOccurMatrix.py
- ex1: python /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step2_EHR2CoOccurMatrix.py -i /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/LongFormat_Num/codeRecord.csv -o /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/Dx_CoOccurMatrix.csv -w 0 1 6 13
- Alternative: Submit parallel jobs, to calculate the cooccurence matrix by chunks (numerous subsets of long number format EHR data, by adding arguments -c and -tc)
- ex2: python /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step2_EHR2CoOccurMatrix.py -i /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/LongFormat_Num/codeRecord.csv -o /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/data/CoOccurMatrix/Dx_CoOccurMatrix.csv -tc 20 -c 2 -w 0 1 6 13

# STEP 3 #
Merge all chunks result together into one triplet format sparse matrix 
- Step3_Convert2SparseMatrix.R
- ex: Rscript /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step3_Convert2SparseMatrix.R
- In the future, add DataSource parameters

# STEP 4 #
Generate embeddings for desired dimension (here, we set dimension of each code embedding vector be 1x500)
- Step4_CreateEmbeddings.R
- ex: Rscript /nfs/turbo/mgi-shixu/project/AnalysisPipeline/CodeEmbedding_pipeline/code/Step4_CreateEmbeddings.R
- In the future, add DataSource and iter parameters

# STEP 5 #
- Calculate cosine similarity and use phecode to obtain AUC
- Phecode/ICD crosswalk and format data
- Code to calculate cosine similary 
- The condensed version of code can be found in Step5_CosineSimilarity_short.R

# References #
[1] Part of the code in Step2 is originated from the package in [LargeScaleClinicalEmbedding](https://github.com/rusheniii/LargeScaleClinicalEmbedding) and modified by the team later.  
[2] Beam, A. L., Kompa, B., Schmaltz, A., Fried, I., Weber, G., Palmer, N., ... & Kohane, I. S. (2019). [Clinical concept embeddings learned from massive sources of multimodal medical data](https://www.worldscientific.com/doi/epdf/10.1142/9789811215636_0027). In PACIFIC SYMPOSIUM ON BIOCOMPUTING 2020 (pp. 295-306).  
