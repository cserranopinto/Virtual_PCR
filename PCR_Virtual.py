# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 08:29:36 2022

@author: cserrano
"""

import argparse
from Bio import SeqIO
from Bio import Align
from Bio import Seq

parser = argparse.ArgumentParser(prog= 'PCR_Virtual.py',
                                 description='BIENVENIDO A SU PROGRAMA DE SIMULACION DE PCR FAVORITO',
                                 usage = 'python PCR_VIRTUAL.py -i INPUT -f forward -r reverse -p 0.8')


# ARGUMENTOS OPCIONALES
parser.add_argument('-i', '--Infile',  help = 'Archivo fasta/multifasta')
parser.add_argument('-f', '--Forward', help = 'Secuencia forward', type=str)
parser.add_argument('-r', '--Reverse', help = 'Secuencia reverse', type=str)
parser.add_argument('-rc', '--ReverseComplement', help = 'Secuencia reverse complement', type=str)
parser.add_argument('-p', '--Pident',  help = 'Porcentaje identidad mínimo (Default 0.9)', default=0.9, type=float)
parser.add_argument('-l', '--Length', help = 'Largo amplicon. (Default 1-10000000)',default='1-10000000', type=str)
parser.add_argument('-s', '--Show', help = 'Muestra producto PCR (A) o producto sin partidores (a).',default= 'A', type=str)
# Read arguments from command line
args = parser.parse_args()

if args.Infile:
    in_file = args.Infile
    
if args.Length:
    l_min = int(args.Length.split('-')[0])
    l_max = int(args.Length.split('-')[1])
    
if args.Infile == None:
    parser.print_help()
    exit()
if args.Forward == None:
    parser.print_help()
    exit()    
if args.Reverse == None and args.ReverseComplement == None:
    parser.print_help()
    exit()


if args.ReverseComplement != None:
    args.Reverse = Seq.Seq(args.ReverseComplement).reverse_complement()


def PCR(ref, seq):    
    
    matches = []
    
    # SE DEFINE ALINEAMIENTO Y SCORES
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'

    aligner.match_score = 1 
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1000

    large = len(str(seq)) 

    # EJECUTA ALINEAMIENTO
    alignments = aligner.align(ref.upper(), seq.upper())

    for ali in alignments:

        pident_ali = round(ali.score/large, 2)
        # print(pident_ali)
        
        # Filtro por porcentage identidad
        if pident_ali >= args.Pident:
    
            c_1 = 0
            c_2 = 0
            
            # Captura coordenadas de alineamiento en target
            # for subali_t,subali_s in ali.aligned:
            subali_t = ali.aligned[0][0]
            subali_s = ali.aligned[1][0]
            dif_L = 0 + subali_s[0]
            dif_R = len(seq) - subali_s[-1]
                
            if len(seq) == (subali_s[1]-subali_s[0]):
                
                if subali_t[0] > c_1:
                    c_1 = subali_t[0]
                if subali_t[1] > c_2:
                    c_2 = subali_t[1]
                    
            elif len(seq) != (subali_s[1]-subali_s[0]):
                
                
                if subali_s[0] != 0 and subali_s[-1] == len(seq):
                    if subali_t[0] > c_1:
                        c_1 = subali_t[0] - dif_L
                    if subali_t[1] > c_2:
                        c_2 = subali_t[1]
                        
                elif subali_s[0] != 0 and subali_s[-1] != len(seq):
                    if subali_t[0] > c_1:
                        c_1 = subali_t[0] - dif_L
                    if subali_t[1] > c_2:
                        c_2 = subali_t[1] + dif_R
                        
                     
                elif subali_s[0] == 0 and subali_s[-1] != len(seq):
                    if subali_t[0] > c_1:
                        c_1 = subali_t[0]
                    if subali_t[1] > c_2:
                        c_2 = subali_t[1] + dif_R
                        
            # Secuencia correspondiente en ref     
            target = ali.target[c_1:c_2]
            
            matches.append([c_1, c_2, str(target)])
            
    return matches
        
#### MAIN ####
def main():
    resume = {}
    multifasta = SeqIO.parse(in_file, 'fasta')

    for fasta in multifasta:
        ID = fasta.description
        forward = PCR(fasta.seq, args.Forward)
        reverse = PCR(fasta.seq, args.Reverse)

        
        if type(forward and reverse) == list:
      
            for f in forward:
                
                for r in reverse:
                    
                    if f[1] < r[0]:
    
                        amp = (fasta.seq)[f[1]: r[0]]
                        
                        if len(amp) >= l_min and len(amp) <= l_max:
                            
                            largo_amp = '(%spb)' %(str(len(amp)))
                            print('>'+ID+ '-' + largo_amp)
                            
                            # Muestra solo amplicon
                            if args.Show == 'a':
                                print(amp.upper())
                            # Muestra primer_F + amplicons + primer_R
                            elif args.Show == 'A':
                                print(f[2].lower() + amp.upper() + r[2].lower())
    
                            
                            if ID not in resume:  
                                resume[ID] = 1
                            else:
                                resume[ID] += 1
                                
        
        forward = PCR(Seq.Seq(fasta.seq).reverse_complement(), args.Forward)
        reverse = PCR(Seq.Seq(fasta.seq).reverse_complement(), args.Reverse)
                                
        if type(forward and reverse) == list:
      
            for f in forward:
                
                for r in reverse:
                    
                    if f[1] < r[0]:
    
                        amp = (Seq.Seq(fasta.seq)).reverse_complement()[f[1]: r[0]]
                        
                        if len(amp) >= l_min and len(amp) <= l_max:
                        # if len(amp) >= 66 and len(amp) <= 74:
                            
                            largo_amp = '(%spb)' %(str(len(amp)))
                            # largo_amp = '_C6'
                            print('>RC_'+ID+ '-' + largo_amp)
                            
                            # Muestra solo amplicon
                            if args.Show == 'a':
                                print(amp.upper())
                            # Muestra primer_F + amplicons + primer_R
                            elif args.Show == 'A':
                                print(f[2].lower() + amp.upper() + r[2].lower())
    
                            
                            if ID not in resume:  
                                resume[ID] = 1
                            else:
                                resume[ID] += 1 
            
if __name__=='__main__':
    main()
