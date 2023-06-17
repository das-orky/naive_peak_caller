# Naive Peak Caller

The first step in any peak caller classifier is to utilize ground truth data set to train the model. Since this is an independent project the ground truth is extracted from publically available data sets. To create the ground truth two approaches were employed:

1. Collecting narrorPeak (.narrowPeak) files from the ENCODE project. DNase-seq or ATAC-Seq experiments are considered.
2. To learn the process in totality a pipeline is created with inputs as fastq (DNase-seq or ATAC-Seq) file and outputs are narrowPeak, WIG file and BedGraph format. Which are created by a peak calling tool like macs3.

Usually in a research paper the authors team up with an experimental scientist to get the data and know the underlying enrichment they are searching for. So peaks are usually human annotated or carefully annotated with tools to highlight the enrichment for which the experiment is done.

The pipeline used in 2nd approach is as follows [Please beware the input parameter to the tools are specific to my setup]:

0. To analyze the quality of the fastq file use the tool multiqc. 

1. bowtie2 to build the reference genome (here hg19 is considered): 

```bash
bowtie2-build --noauto --large-index --threads 16 --dcv 4096 --bmax 100000000 hg19.fa hg19
```
2. Align the downloaded fastq file to the reference genome (for single ended read):

```bash
time bowtie2 --large-index -p 16 -q --local -x /path_to_reference/hg19 -U /path_to_fastq/nanog_rep1.fastq -S /path_to_output/nanog_rep1_unsorted.sam
```

3. Samtool to convert the .sam file to .bam file

```bash
time samtools view -h -S -b -@ 16 -o nanog_rep1_unsorted.bam nanog_rep1_unsorted.sam
```

4. sort and remove duplicate entries, using sambamba:

```bash
sambamba view -h -t 2 -f bam -F "[XS] == null and not unmapped  and not duplicate" nanog_rep1_sorted.bam > nanog_rep1_align.bam
```

5. It is better for an experiment to have a control. If control is available run the above steps with the control files. 

6. To call peaks on the .bam file macs3 is used. If control is available use the -c parameter to identify it. It is a good approach to create a log file, which contains meta data. 

```bash
macs3 callpeak -t /__path__/nanog_rep1_sorted.bam -c /__path__/nanog_control_sorted.bam -n output_name -g hs --bdg -q 0.05 -f BAM 2> ~/__path__/output.log 
```

7. Finally .narrowPeak, WIG file and BedGraph files will produced. narrowPeak can be used as a positive or ground truth for the peak calling classifier.

**Various peak calling algorithms are there and the use of different aligner bowtie2 to BWA-MEM can vary the number of peaks called. But which are false positive is hard to define. Better to have a geneticist in the mix. Since for a particular problem at hand different tools and settings are preferable.** 

The narrowPeak file either from the pipeline or directly from ENCODE portal is used as an input. Few parameters that needs to be checked before running the python code on the file. 

1. Window size: Peak calling algorithm will output ranges of nucleotides where the tool thinks enrichment happened. These ranges can be very narrow and quite broad as well, we need to find a uniform window size since DNN (Deep Neural Network) taks as input a fixed length not variable length input. One can plot the distribution of the ranges and select an appropriate window size consisting of most peaks. In this experiment window size is 150.
2. Score: The quality of peak is marked in the narrowPeak file. We don't want to create out positive or ground truth with low quality peaks. This score is usually used to indicate how dark the peak will be displayed in the browser, but used here to judge the quality of the peak. Plot the scores and find a threshold good enough.

The following function converts a narrowPeak file to a file of dimension ( No_of_peaks X window_size) of nucleotides. Where no_of_peaks depends on the selected window size and quality score selected. 

```bash
python3 narrow_to_sequence.py <path_to_narrowPeak_file> <path_to_ref_genome_folder>
```

The output file is created in the same folder as the narrowPeak file folder. The file name: <code>input_file_name_DNA.txt</code>

This file is used as the positive or ground truth. To create the negative value file, i.e samples which will not produce a peak. There can be two approaches to create such a file. 

1. Find ranges of size window_size that are mutually exclusive to the ranges found in the narrowPeak file. In this implementation the narrowPeak file is used to find these mutually exclusive ranges.
2. Another approach can be to use the control file (if it exists). Find ranges of size window_size that are mutually exclusive to the ranges present in the narrowPeak file.

The program <code> narrow_to_sequence.py</code> take care of it automatically and will create a negative sample file with suffix <code> input_file_name_DNA_neg.txt</code>.

Collect all the positive and negative samples file and segregate into two folders. Then use the <code> peak_call_dnn.ipynb</code> to call peaks. Please note this is just an experimental development. There are hard coded paths for positive and negative samples. Please change them.

The notebook to call peaks will be converted into scripts in near future. 