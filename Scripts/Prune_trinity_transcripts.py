# Written by Sofia Medina (Nov 13, 2019). 
# This program is intended to filter out spurious trascript assemblies from Trinity.
# To run: module load python/2.7-anaconda-4.4 && source activate Cori_new &&
# Prune_trinity_transcripts.py <bam> -s 

import argparse, os, pysam, sys
import pandas as pd
import numpy as np
import os
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(description='This script splits up single- and multi-exon Trinity transcripts from the input bam. Single-exon transcripts should be later filtered with TransDecoder. Multi-exon transcripts are shortened if the first and/or last exon is less than 150 bp. For stranded RNAseq single exons are keps it the length is >250bp and for non-stranded if the alignment length is >800bp')

    parser.add_argument('-o', '--output', metavar='STR', help='output file prefix [output]', type=str, default='output')
    parser.add_argument('-s', '--stranded',  help='Transcripts come from stranded RNA-seq library', action='store_true')
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')
    required = parser.add_argument_group('required arguments')
    required.add_argument("bam", help="input bam of Trinity transcripts aligned to genome", type=str)
    args = parser.parse_args()
    return(args)
    
    
def main():
    args = parse_args()
    STAR_alignment = args.bam 
    stranded = args.stranded

    print "\n *** Filtering and pruning Trinity transcripts *** "
    print "Input alignment:\t",   STAR_alignment
    print "\t\t\t(All_Trinity.fa -> STARlong to Genome Index #2 --> Alignment.bam)"
    
    Exon_limit_min_size = 25
    
    if stranded == False:
        print "Non-stranded RNA-seq parameters:"
        Single_exon_min_seq = 600
        print "Non-stranded RNA-seq parameters:"
        print "\tExon_limit_min_size: ",Exon_limit_min_size
        print "\tSingle_exon_min_seq: ",Single_exon_min_seq
    if stranded == True:
        Single_exon_min_seq = 300
        print "Stranded RNA-seq parameters:"
        print "\tExon_limit_min_size: ",Exon_limit_min_size
        print "\tSingle_exon_min_seq: ",Single_exon_min_seq
        
    if 'output' == args.output:
        fasta_out_single = STAR_alignment.replace('Aligned.sortedByCoord.out.bam','_Single_exon.fa') 
        fasta_out_multi  = STAR_alignment.replace('Aligned.sortedByCoord.out.bam','_Multi_exon.fa') 
        fasta_out_badMQ  = STAR_alignment.replace('Aligned.sortedByCoord.out.bam','_LowQual.fa')
        summary_out      = STAR_alignment.replace('Aligned.sortedByCoord.out.bam','_summary.out')
        sample_id        = STAR_alignment.replace('Aligned.sortedByCoord.out.bam','').split('/')[-1]
    else:
        fasta_out_single = ''.join((args.output,'_Single_exon.fa'))
        fasta_out_multi  = ''.join((args.output,'_Multi_exon.fa'))
        fasta_out_badMQ  = ''.join((args.output,'_LowQual.fa'))
        summary_out      = ''.join((args.output,'_summary.tab'))
        sample_id        = args.output.split('/')[-1]
        

    samfile = pysam.AlignmentFile(STAR_alignment, "rb")
    prune_transcripts(samfile, Exon_limit_min_size, Single_exon_min_seq, fasta_out_single, fasta_out_multi, fasta_out_badMQ, summary_out, sample_id)
    samfile.close()
    
    print "Files saved as: \n\t", fasta_out_single,"\n\t", fasta_out_multi , "\n\t",fasta_out_badMQ,"\n\t",summary_out,"\n"
    
    return()
    

def obtain_sequence_coordinates_of_pruned_transcript(cigar_type,cigar_len, list_introns, min_len, seq_length):
    min_length_first_exon = min( cigar_len[list_introns[0]]/float(500) + min_len, 150)
    min_length_last_exon = min( cigar_len[list_introns[-1]]/float(500) + min_len, 150)
    length_first_exon = 0
    last_exon_coordinate = seq_length
    first_exon_coordinate = 0
    for i in range(0,list_introns[0]):
        if (cigar_type[i] == 0) | (cigar_type[i] == 4):
            length_first_exon =  length_first_exon + cigar_len[i]
        if (cigar_type[i] == 1): #Deletion
            length_first_exon = length_first_exon + cigar_len[i]
    
    length_last_exon = 0
    for i in range(list_introns[-1]+1,len(cigar_len)):
        if (cigar_type[i] == 0) | (cigar_type[i] == 4):
            length_last_exon =  length_last_exon + cigar_len[i]
        if (cigar_type[i] == 1): 
            length_last_exon = length_last_exon + cigar_len[i]
    
    if length_first_exon< min_length_first_exon:
        first_exon_coordinate = length_first_exon      
    if length_last_exon<min_length_last_exon:
        last_exon_coordinate  = seq_length - length_last_exon
    return(first_exon_coordinate, last_exon_coordinate)


def Obtain_right_seq_orientation(New_sequence, strand):
    if strand == '-':
        return(str(Seq(New_sequence).reverse_complement()))
    else:
        return(New_sequence)
    

def prune_transcripts(samfile, Exon_limit_min_size, Single_exon_min_seq, fasta_out_single, fasta_out_multi, fasta_out_badMQ, summary_out, sample_id):
    
    Min_Aln_Score = 300
    Exon_limit_min_size = 25
    Multi_exon  = open(fasta_out_multi, 'w') 
    Single_exon = open(fasta_out_single, 'w')
    Low_MQ      = open(fasta_out_badMQ, 'w')
    num_alignments, num_seqs_over_255,num_multiple,num_single, too_short,n_Low_MQ = (0,0,0,0,0,0)
    fasta_lq = ''
    
    for read in samfile.fetch():
        New_sequence = ''
                
        num_alignments = num_alignments+1 
        if read.is_secondary or read.is_supplementary:
            continue
    
        if  read.mapping_quality < 255:
            n_Low_MQ= n_Low_MQ+1
            sequence = read.seq
            if read.is_reverse== False:
                sequence =  str(Seq(read.seq).reverse_complement())
            fasta_lq = ''.join(('>',read.qname,' n=',str(len(read.seq)),'\t',read.reference_name,":",str(read.reference_start),'-',str(read.reference_end), '\t',str(dict(read.get_tags())['AS']),' + ',read.cigarstring,'\n',sequence,'\n'))            
            Low_MQ.write(fasta_lq)
                
        if  read.mapping_quality == 255:
            num_seqs_over_255 = num_seqs_over_255+1
            cigar_string = read.cigarstring
            seq_length = len(read.seq)
            alignment_score = dict(read.get_tags())['AS']
            cigar = read.cigar
            cigar_type = [i[0] for i in cigar]
            cigar_len = [i[1] for i in cigar]
            list_Matches  = [i for i in range(len(cigar_type)) if cigar_type[i] == 0] 
            list_introns  = [i for i in range(len(cigar_type)) if cigar_type[i] == 3] 
            
            if (read.is_reverse== True): strand = '+'
            if (read.is_reverse== False): strand = '-'          
            
            if 'N' in cigar_string:  
            # If there is an intron
    
                start_coord , end_coord = obtain_sequence_coordinates_of_pruned_transcript(cigar_type,cigar_len, list_introns, Exon_limit_min_size, seq_length)
                New_sequence = Obtain_right_seq_orientation(read.seq[start_coord:end_coord], strand)
                
                if (len(New_sequence) > Single_exon_min_seq)  & (alignment_score > Min_Aln_Score):
                    extra = ''.join(('[',str(start_coord) ,':',str(end_coord),']','::',str(len(New_sequence))))
                    seq_ = ''.join(('>',read.qname,' n=',str(len(read.seq)),extra,'\t',read.reference_name,":",str(read.reference_start),'-',str(read.reference_end), '\t',str(alignment_score),' ',strand,read.cigarstring,'\n',New_sequence,'\n'))
                    Multi_exon.write(seq_)
                    num_multiple=num_multiple+1
                else:
                    too_short = too_short + 1
                    
            else:  
            # If sequence has no introns:
                list_matches = sum(np.array(cigar_len)[list_Matches])
                if (list_matches > Single_exon_min_seq) & (alignment_score>Min_Aln_Score):
                    num_single=num_single+1
                    New_sequence = Obtain_right_seq_orientation(read.seq, strand)
                    seq_ = ''.join(('>',read.qname,' n=',str(len(read.seq)),'\t',read.reference_name,":",str(read.reference_start),'-',str(read.reference_end), '\t',str(alignment_score),' ',strand,read.cigarstring,'\n',New_sequence,'\n'))
                    Single_exon.write(seq_)
                    
                else:
                    too_short = too_short + 1
                    
    samfile.close()
    Single_exon.close()
    Multi_exon.close()
    Low_MQ.close()
    
    Summary = pd.DataFrame(index=['Num_alignments','Low_MQ','Good_MQ','Good_MQ_Multi-exon','Good_MQ_Single-exon','Good_MQ_too_short'], columns =[sample_id])
    Summary.loc['Num_alignments',sample_id] = num_alignments
    Summary.loc['Low_MQ',sample_id] = n_Low_MQ
    Summary.loc['Good_MQ',sample_id] = num_seqs_over_255
    Summary.loc['Good_MQ_Multi-exon',sample_id] = num_multiple
    Summary.loc['Good_MQ_Single-exon',sample_id] = num_single
    Summary.loc['Good_MQ_too_short',sample_id] = too_short
    Summary.to_csv(summary_out, sep='\t')
    print Summary

    return()



if __name__=='__main__':
    main()
