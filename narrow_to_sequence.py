import pandas as pd
from pathlib import Path
from shutil import which
import subprocess
import argparse
import os

def main():

    parser=argparse.ArgumentParser()
    parser.add_argument(help='Path to the MACS narrowpeak file',dest='narrowpeak',type=str)        
    parser.add_argument(help='Path to the reference Genome',dest='ref_genome',type=str)  
    #parser.add_argument(help='Path to the reference Genome',dest='ref_genome',type=Path)                  
    args = parser.parse_args()
    
    #file_name=sys.argv[1]
    file_name=args.narrowpeak
    filename=Path(file_name)
    ref_file_path=Path(args.ref_genome)
    out_file_bed=file_name+"_.bed"
    out_fastq_file=file_name+".fastq"
    out_dna_file=file_name+"_DNA.txt"
    out_chromo_index=file_name+"_chr_index.txt"


    if(os.path.isfile(filename)):
        print(filename)
    else:
        print("Please give a valid narrowpeaks filename")

    narrow_peak_raw=pd.read_csv(filename, sep='\t', lineterminator='\n',header=None)
    narrow_peak_raw.columns = ['chro','start_pos','end_pos','id','score','strand','signalValue','pvalue','qvalue','peak']

    # convert the MACS narrowpeak file into sequences containing positive data, i.e contains peaks

    narrowpeak_to_bed(narrow_peak_raw,window_size=150, score=50,out_file_bed=out_file_bed,out_chromo_index=out_chromo_index)

    bed_to_fastq(out_file_bed,ref_file_path,out_fastq_file)

    fastq_to_txt(out_fastq_file=out_fastq_file,out_dna_file=out_dna_file)

    out_file_bed=file_name+"_neg.bed"
    out_fastq_file=file_name+"_neg.fastq"
    out_dna_file=file_name+"_DNA_neg.txt"

    # convert the MACS narrowpeak file into genome sequence containing negative data, i.e no peak 

    narrowpeak_to_bed_neg(narrow_peak_raw=narrow_peak_raw,out_chromo_index=out_chromo_index,out_file_bed=out_file_bed, window_size=150)

    bed_to_fastq(out_file_bed,ref_file_path,out_fastq_file)

    fastq_to_txt(out_fastq_file=out_fastq_file,out_dna_file=out_dna_file)


def narrowpeak_to_bed(narrow_peak_raw, window_size, score,out_file_bed,out_chromo_index):
    
    # Cut a window_size from the peaks identified by the MACS tool.
    #
    # Select those peaks whose difference is more than equal to the window_size.
    # i.e: end_position - start-position > 150 
    # After the filteration based on the window_size. Next the peaks are filtered by
    # the value called 'score'. The 5th value in narrowpeaks file. It is calculated as
    # int(-10*log10pvalue) OR for q-value. The peak can be very broad, to cut out a 'window_size'
    # portion from a peak, place the window in such a way that center of the window is the value 'peak'
    # i.e the 10th value in the narrowpeaks file of MACS output.

    

    temp=narrow_peak_raw[ (narrow_peak_raw['end_pos'] - narrow_peak_raw['start_pos']) > window_size]
    temp=temp[temp['score'] > score]

    peak_final=pd.DataFrame(columns=['chro','start_pos','end_pos'])
    peak_final['chro']=temp['chro']
    peak_final['start_pos']= (temp['start_pos'] + temp['peak']) -((window_size//2)-1)
    peak_final['end_pos']=((temp['start_pos'] + temp['peak']) -((window_size//2)-1)) + window_size

    peak_final.to_csv(out_file_bed, header=None, index=None, sep='\t',mode='w')    

    #chromosome_index(dataset=peak_final,out_chromo_index=out_chromo_index)
        

def bed_to_fastq(out_file_bed,ref_file_path,out_fastq_file):

    # Convert the bed file to a fastq file
    
    if which('bedtools') is not None:
        command="bedtools getfasta -fi {} -bed {} -fo {}".format(ref_file_path,out_file_bed,out_fastq_file)
        bedtools_output=subprocess.call(command, shell=True)
    else:
        print("Either Bedtools is not installed or environment variable for Bedtools is not set")

def fastq_to_txt(out_fastq_file, out_dna_file):

    # To extract the sequences from the fastq files.
    # NR%2 menas collect two lines at a time and do modulo operation with 2 to find the even line
    # Even lines contains the DNA sequence in a fastq file

    command="paste -d\"\n\" {} | awk '(NR%2==0)' > {}".format(out_fastq_file,out_dna_file)
    to_DNA=subprocess.call(command,shell=True)

def chromosome_index(dataset,out_chromo_index):

    # It will create a metadata file that contains the chromosomes name, the start and the
    # end row of the dataset. i.e starting row and the end row for a certain chromosome that is
    # present in the dataset. e.g chr2, 567 , 1023 : Means chrosome-2 start from row 567 and end
    # at row number 1023. The first if statement is for the last chromosome in the dataset. If the
    # current chromosome is equal to next chromosome then continue. If current is not equal to the
    # next chromosome write to a file.

    start=0
    end=0
    f=open(out_chromo_index,'w')
    f.close()
    chro_num=dataset.iloc[0,0]

    for index in range(1,len(dataset)):
        chro_next=dataset.iloc[index,0]
        if index==(len(dataset)-1):
            with open(out_chromo_index, 'a') as f:
                f.write(chro_num+" , "+str(start)+" , "+ str(len(dataset)-1)+"\n")
        elif chro_num==chro_next:
            continue
        else:
            #print(chro_num,chro_next)
            end=index
            with open(out_chromo_index, 'a') as f:
                f.write(chro_num+" , "+str(start)+" , "+ str(end-1)+"\n")
            start=index
            chro_num=dataset.iloc[index,0]

def narrowpeak_to_bed_neg(narrow_peak_raw, out_chromo_index,out_file_bed,window_size):
    
    # To create negative dataset i.e dataset for of no peak region. 
    chromosome_index(narrow_peak_raw,out_chromo_index)

    chrom_index=pd.read_csv(out_chromo_index, lineterminator='\n',header=None)
    chrom_index.columns=['chro','start_pos','end_pos']
    negative_index=pd.DataFrame(columns=['chro','start_pos','end_pos'])
    end=[]
    for i in range(0,len(chrom_index)):
        end.append(chrom_index.iloc[i,2])
    #window_size=150
    start=0
    #temp=[]
    for i in end:
        for j in range(start,i):
            t=j+1
            if t == end[-1]:
                gap=narrow_peak_raw.iloc[(j+1),1] - narrow_peak_raw.iloc[j,2]
                if gap > window_size*2:
                    mid=gap//2
                    neg_start=narrow_peak_raw.iloc[j,2] + (mid - window_size//2)
                    neg_end=narrow_peak_raw.iloc[j,2] + (mid + window_size//2)    
                    negative_index.loc[j]=[narrow_peak_raw.iloc[j,0],neg_start,neg_end]
                #temp.append(narrow_peak_raw.iloc[j,2])
                #temp.append(narrow_peak_raw.iloc[(j+1),1])
                break
            elif t == i:
                gap=narrow_peak_raw.iloc[(j+1),1] - narrow_peak_raw.iloc[j,2]
                if gap > window_size*2:
                    mid=gap//2
                    neg_start=narrow_peak_raw.iloc[j,2] + (mid - window_size//2)
                    neg_end=narrow_peak_raw.iloc[j,2] + (mid + window_size//2) 
                    negative_index.loc[j]=[narrow_peak_raw.iloc[j,0],neg_start,neg_end]
                #temp.append(narrow_peak_raw.iloc[j,2])
                #temp.append(narrow_peak_raw.iloc[(j+1),1]) 
                break
            gap=narrow_peak_raw.iloc[(j+1),1] - narrow_peak_raw.iloc[j,2]
            if gap > window_size*2:
                mid=gap//2
                neg_start=narrow_peak_raw.iloc[j,2] + (mid - window_size//2)
                neg_end=narrow_peak_raw.iloc[j,2] + (mid + window_size//2) 
                negative_index.loc[j]=[narrow_peak_raw.iloc[j,0],neg_start,neg_end]
            #temp.append(narrow_peak_raw.iloc[j,2])
            #temp.append(narrow_peak_raw.iloc[(j+1),1])                
        start=i+1   

    negative_index.to_csv(out_file_bed, header=None, index=None, sep='\t',mode='a')


if __name__ == "__main__":
    main()
    print("Probably done with the conversion look for file")