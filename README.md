## Step 1. Demultiplexing and tabulating BarSeq results
### Getting environment and files ready
1. Upload files to CHPC. We have two levels of multiplexing so that each well has its own combination of forward and reverse adapters that can be used for demultiplexing and keeping individual samples separate after sequencing.
    - The sequencing facility will demultiplex using the first set of adapters (in our case IT001-IT006). This corresponds the the A-F rows of our plates. Zipped sequence files will be in each of the IT001-IT006 folders; upload these folders.
    - Prepare an idxfile.txt for each of the IT00# directories. This file contains the barcodes that correspond to the second demultiplexing step and are unique for each column of our plate (1-12). For each IT00# folder, change the identifier so that IT001 = A, IT002 = B, etc. The barcodes will not change. Example here:
     - In IT001, upload the following .txt file with tab separation between the two columns:
       
        |   |  |
        |----|----|
        | A01P01 | NGCACTANNN |
        | A02P01 | NNTGTAGCNN |
        | A03P01 | NNNCGGATTN |
        | A04P01 | NNNNACCAGT |
        | A05P01 | NGTGACANNN |
        | A06P01 | NNTAACCGNN |
        | A07P01 | NNNCTAGACN |
        | A08P01 | NNNNAGTTCA |
        | A09P01 | NGACTAGNNN |
        | A10P01 | NNTTCGATNN |
        | A11P01 | NNNCATCGGN |
        | A12P01 | NNNNATGTTC |
      - In IT002, upload the following .txt file:
        
        |   |  |
        |----|----|
        | B01P01 | NGCACTANNN |
        | B02P01 | NNTGTAGCNN |
        | B03P01 | NNNCGGATTN |
        | B04P01 | NNNNACCAGT |
        | B05P01 | NGTGACANNN |
        | B06P01 | NNTAACCGNN |
        | B07P01 | NNNCTAGACN |
        | B08P01 | NNNNAGTTCA |
        | B09P01 | NGACTAGNNN |
        | B10P01 | NNTTCGATNN |
        | B11P01 | NNNCATCGGN |
        | B12P01 | NNNNATGTTC |
      - do this for all of the IT00# folders.
3. Create new conda environment called je_env.
     ```
     conda create --name je_env
     ```
     ```
     conda activate je_env
     ```
     Install [je](https://github.com/gbcs-embl/Je).
     ```
     conda install -c bioconda je-suite
     ```
     Install [cutadapt](https://cutadapt.readthedocs.io/en/stable/).
     ```
     conda install bioconda::cutadapt
     ```
5. Upload the [loop.sh](https://github.com/delaney-beals/RB-TnSeq/blob/main/loop.sh) script and make it executable.
      ```
      chmod +x loop.sh
      ```

### For each "IT00#" directory...
0. Activate je_env if needed. 
    ```
    conda activate je_env
    ```
1. Demultiplex samples and put output into a new directory called demultiplexed_NEW by running the following, making sure that F1 and F2 correspond to the names of the .fq.gz files in the IT00# directory.
    ```
    je demultiplex F1=IT001_CKDL240032762-1A_22GKHFLT4_L3_1.fq.gz F2=IT001_CKDL240032762-1A_22GKHFLT4_L3_2.fq.gz BF=idxfile.txt O=demultiplexed_NEW BPOS=READ_1
    ```
2. Move loop.sh into the demultiplexed_NEW directory that you're currently working on.
    ```
    cd ../
    cp IT002/demultiplexed_NEW/loop.sh IT003/demultiplexed_NEW
    ```
3. Run the executable loop.sh in interactive mode. 
    ```
    ./loop.sh
    ```
4. Go into the newly created demultiplexed_NEW
    ```
    cd demultiplexed_NEW/
    ```
5. Check the number of barcodes found for each of the samples by looking in the output file called jemultiplexer_out_stats.txt
    ```
    cat jemultiplexer_out_stats.txt  
    ```
6. Check the distribution of sequence lengths associated with each barcode detected. The majority of reads should be 20 bp long.
    ```
    cat D01P01/read_length_D01P01.txt
    ```
8. Move to the next IT00# directory and repeat this process.


## Step 2. Data exploration
Perform clustering and preliminary DESeq2 analyses to confirm that technical replicates are similar to one another and that any differences between conditions are more likely to be from biological, rather than technical, factors. Use R script [data_exploration.Rmd](https://github.com/delaney-beals/RB-TnSeq/blob/main/data_exploration.Rmd). 

## Step 3. Gene fitness analysis
Estimate the fitness of each gene from experiments. Use R script [FEBA.R](https://github.com/delaney-beals/RB-TnSeq/blob/main/FEBA.R).
