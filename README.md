SNVoter
=======
  
**Improving SNV detection from low coverage nanopore sequencing data (<30x)**
  
Table of Contents
=================

* **[Installation](https://github.com/vahidAK/NanoMethPhase/blob/master/README.md#installation)**
  * [Using pip](https://github.com/vahidAK/SNVoter#using-pip)
  * [From source](https://github.com/vahidAK/SNVoter#from-source)
* **[SNVoter Modules](https://github.com/vahidAK/SNVoter#snvoter-modules)**
  * [prediction](https://github.com/vahidAK/SNVoter#prediction)
  * [extraction](https://github.com/vahidAK/SNVoter#extraction)
  * [train](https://github.com/vahidAK/SNVoter#train)
* **[Tutorial](https://github.com/vahidAK/SNVoter#tutorial)**
  * [Variant Calling](https://github.com/vahidAK/SNVoter#variant-calling)
  * [Improving SNV calling using SNVoter](https://github.com/vahidAK/SNVoter#improving-snv-calling-using-snvoter)
* **[Example](https://github.com/vahidAK/SNVoter#example)**
  
# Installation
**NOTE:** SNVoter uses several fixed versions of its dependencies in [environment.yaml](https://github.com/vahidAK/SNVoter/blob/master/envs/environment.yaml) file . Users are encouraged to use a conda or similar environment to isolate the packages from their
default python instance. Then activate the environment and install SNVoter using pip or you can clone the git repo and use it from source.  
You can make the conda environment and install all dependencies by downloading the [environment.yaml](https://github.com/vahidAK/SNVoter/blob/master/envs/) file and running these lines of codes:  

```
conda env create -f environment.yaml
conda activate snvoter
```
Now you can install SNVoter using pip or use it from source in the dedicated environment with all dependencies installed.  
## Using pip

```
pip install snvoter
```

## From source

```
git clone https://github.com/vahidAK/SNVoter.git
cd SNVoter
./snvoter.py
```
# SNVoter Modules
## prediction:
To predict dtetedte SNVs are true calls or false positives.
```
usage: snvoter prediction [-h] --input INPUT --bam BAM --reference REFERENCE
                          --output OUTPUT [--model_file MODEL_FILE]
                          [--mappingQuality MAPPINGQUALITY] [--depth DEPTH]
                          [--window_bam WINDOW_BAM] [--nine_mer]
                          [--threads THREADS] [--chunk_size CHUNK_SIZE]

Predict based on a model.

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  --input INPUT, -i INPUT
                        The path to the input vcf or bed file. NOTE. Files
                        must end with .bed or .vcf. vcf files are 1-based and
                        beds are zero-based
  --bam BAM, -b BAM     The path to the alignment bam file
  --reference REFERENCE, -r REFERENCE
                        The path to the reference file. File must be indexed
                        by samtools faidx.
  --output OUTPUT, -o OUTPUT
                        The path to the output directory and prefix for output
                        file.

optional arguments:
  --model_file MODEL_FILE, -mf MODEL_FILE
                        Path to the trained model. Default is
                        NA12878_20FC_model.h5
  --mappingQuality MAPPINGQUALITY, -mq MAPPINGQUALITY
                        Cutt off for filtering out low quality mapped reads
                        from bam. Default is 0
  --depth DEPTH, -d DEPTH
                        Cutt off for filtering out regions with low depth to
                        have frequencies. Default >= 1
  --window_bam WINDOW_BAM, -w WINDOW_BAM
                        if you want to only do for a region or chromosom You
                        must insert region like this chr1 or chr1:1000-100000.
  --nine_mer, -nm       Prediction for 9-mer. Default is five mer
  --threads THREADS, -t THREADS
                        Number of threads. Default is 4.
  --chunk_size CHUNK_SIZE, -cs CHUNK_SIZE
                        Chunk size. Default is 100.
```
## extraction:
Extract features to train a new model.
```
usage: snvoter extraction [-h] --input INPUT --mod_status MOD_STATUS --bam BAM
                          --reference REFERENCE
                          [--mappingQuality MAPPINGQUALITY] [--depth DEPTH]
                          [--window_bam WINDOW_BAM] [--nine_mer]
                          [--threads THREADS] [--chunk_size CHUNK_SIZE]

Extract mutation frequencicies in 5-mer window.

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  --input INPUT, -i INPUT
                        The path to the input vcf or bed file. NOTE. Files
                        must end with .bed or .vcf. vcf files are 1-based and
                        beds are zero-based
  --mod_status MOD_STATUS, -ms MOD_STATUS
                        0 or 1. If you are extracting frequencies to train a
                        model, give the modification status for your bed file
                        either it is modified (1) or unmodified (0) regions.
  --bam BAM, -b BAM     The path to the alignment bam file
  --reference REFERENCE, -r REFERENCE
                        The path to the reference file. File must be indexed
                        by samtools faidx

optional arguments:
  --mappingQuality MAPPINGQUALITY, -mq MAPPINGQUALITY
                        Cutt off for filtering out low quality mapped reads
                        from bam. Default is 0
  --depth DEPTH, -d DEPTH
                        Cutt off for filtering out regions with low depth to
                        have frequencies. Default >=1
  --window_bam WINDOW_BAM, -w WINDOW_BAM
                        if you want to only do for a region or chromosom, you
                        must insert region like this chr1 or chr1:1000-100000.
  --nine_mer, -nm       Extraction for 9-mer. Default is five mer
  --threads THREADS, -t THREADS
                        Number of threads
  --chunk_size CHUNK_SIZE, -cs CHUNK_SIZE
                        Number of sites send to each processes for parrallel
                        processing. Default is 50.
```
## train:
To train a new model using extracted features.
```
usage: snvoter train [-h] --train TRAIN --test TEST --out_dir OUT_DIR
                     [--epochs EPOCHS] [--batch_size BATCH_SIZE] [--plot]
                     [--nine_mer]

train a new model

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  --train TRAIN, -tr TRAIN
                        The path to the shuffled and ready file for training
  --test TEST, -te TEST
                        The path to the shuffled and ready file for testing.
  --out_dir OUT_DIR, -o OUT_DIR
                        Output directory and prefix for saving the model and
                        figures

optional arguments:
  --epochs EPOCHS, -e EPOCHS
                        Number of epochs. Default is 100
  --batch_size BATCH_SIZE, -batch BATCH_SIZE
                        batch size for model training. Default is 400.
  --plot, -plt          Select this option if you wish to output training
                        plots.
  --nine_mer, -nm       Training for 9-mer. Default is five mer
```
# Tutorial

## Variant Calling

You first need to call variants using [Clair](https://github.com/HKU-BAL/Clair)

You can call variants for each chromosome using the following command and the
concatenate all files:

```
for i in chr{1..22} chrX chrY; do callVarBam --chkpnt_fn <path to model file> --ref_fn <reference_genome.fa> --bam_fn <sorted_indexed.bam> --ctgName $i --sampleName <your sample name> --call_fn $i".vcf" --threshold 0.2 --samtools <path to executable samtools software> --pypy <path to executable pypy > --threads <number of threads>
```

For the full tutorial please refer to [Clair](https://github.com/HKU-BAL/Clair)
page on GitHub.

## Improving SNV calling using SNVoter:

```
snvoter prediction -i <SNVs_Clair.vcf> -b <sorted_indexed.bam> -mf <path to model file (model.h5)> -r <reference_genome.fa> -t number_of_threads -o output_prefix
```

It will produce two files.

1- Prediction file that includes each prediction for each 5-mer. The first 10
columns are from vcf file and the last seven columns indicate:
   - **chrom**:            the chromosome name
   - **pos_start**:        0-based position of the 5-mer start
   - **pos_end**:          0-based position of the 5-mer end
   - **pos**:              0-based position of the SNV
   - **5-mer sequence**:   sequence of five-mer
   - **Coverage**:         this might be different from Clair's coverage as
                           SNVoter uses different mapping quality threshold
   - **Prediction**

2- The second file is the ready vcf file with weighted qualities. You can plot
the distribution of weighted quality to obtain optimal threshold for filtering.
The plot usually looks like the following plots:
![Quality distribution of 10x coverage data](docs/images/QualDist10x.png)
![Quality distribution of 18x coverage data](docs/images/QualDist18x.png)
![Quality distribution of 22x coverage data](docs/images/QualDist22x.png)

The optimal threshold is the end of the first peak and start of the valley
(highlighted regions).

By default SNVoter will use the model file trained by us using NA12878 20 flow cells and you do not need to specify path to the model if you want to use our model.  

# Example
We have included an example data in the Example_data folder which you can use for a quick prediction.
