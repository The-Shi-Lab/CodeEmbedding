# Analysis Pipeline to Generate Code Embeddings -- an example using MGI data

## Adapted codes from Xu Shi, Xianshi Yu, and Kuan-Han Wu ##

# STEP 1 #
As shown in Step1_EHR2LongNumFormat.R. This is a data prepossessing step to convert original EHR data into desired format. In this step we convert the raw data DiagnosisCode ProcedureCode and ResultCode from the dataset diagnosis, procedure, and lab respectively into unique long-form and record them in column code_num in new datasets. The new dataset converted from the raw set now has 3 columns: patient id, the day of the visit recorded as the day since birth, and medical codes as code_num.
- Notice: Please note that this procedure is specified with different data resources, that is, you need to create your personalized converting code for Step1 to get ready for further steps. 
- All the data are numeric variables and sorted by Patient ID and then numDays


# STEP 2 #
As shown in Step2_EHR2CoOccurMatrix.py. This step uses data modified from step1 to calculate the co-occurrence matrix.
- Read in data and set up parameters
  * Create command line arguments,this part can be modified accordingly with the needs to convenient the process of running the codes.
   ```
   parser = argparse.ArgumentParser(description='Create co-occurrence matrix for EHR code within different time windows.')
   parser.add_argument('argument',nargs='associates a different number of command-line arguments with a single action',type=check the argument and converse it into integer,help='customized help message')
   args = parser.parse_args()
   ```
  * Define window in terms of day difference
  * Create matrices to store co-occurence count for each windows
  * Read input file(the data set we created in Step1)
  * Input chunk and total chunks through command line
  * Calculate total number of patients and number of patients per chunk 
   ```
    total patients= max(patient id)--patients' id are recorded as a continuous integer
    number of patients per chunk = round(total patients / total number of chuncks)
   ```
  * Calculate minimum patient id and maximum patient id for a chunk 
- Calculate cooccurence matrix
  * subset data to chunks based on the calculated patient ID range
  * calculate cooccurance matrix
   ```
   FOR i in the range of subset chunk
      READ patient_id, day, code to events_perpt.iloc[i]
    FOR j from i+1 to the end of subset chunk
      READ npatient_id, nday, ncode to events_perpt.iloc[j]
        FOUND tempLoc(non-overlapping window interval) to insert the count
        INSERT count
         IF code < ncode 
           matrices[tempLoc][(code,ncode)]+=1
           ELSE
             matrices[tempLoc][(ncode, code)]+=1
         END IF  
        IF patient_id = npatient_id 
              BREAK
        IF nday - day > the largest inputted window 
              BREAK
    ``` 
 - Format output table
   * format table to data frame
   * 
   
  

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
