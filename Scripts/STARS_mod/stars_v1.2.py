#Performs STARS analysis for genetic perturbation screens
#Input: 1. Input file 
#       2. Chip File
#       3. Direction("P" or "N")
#       4. Threshold 
#       5. Null distribution file for threshold specified
#       6. Use first perturbation value(Y/N)
import pandas as pd
import numpy as np
from scipy import stats
from math import log, isnan, isinf, log10
import csv, argparse, os, sys, re
from datetime import datetime
from statsmodels.distributions.empirical_distribution import ECDF

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file',
        type=str,
        help='File containing sgRNA sequence and score(s) columns with headers')
    parser.add_argument('--chip-file',
        type=str,
        help='Chip file')
    parser.add_argument('--dir',
        type=str,
        help='Direction of scores')
    parser.add_argument('--thr',
        type=int,
        help='Threshold')
    parser.add_argument('--null',
        type=str, 
        help='Null distribution file')
    parser.add_argument('--use-first-pert',
        type=str,
        default='N',
        help='Y/N for value of first perturbation to be used in final STARS score calculation for a gene')
    return parser

'''
Sorts a data frame on the column specified and re-indexes
Argument: Dataframe, column to sort on, direction of sorting
Return Value: Sorted and re-index dataframe
'''
def sort_reindex(df, col, direction):
    if direction == 'P' or direction == 'p':
        df = df.sort_values(by=col, ascending=False)
        df.index = range(0, len(df))
        #df['Rank'] = map(lambda x:x+1, range(len(df)))
    elif direction == 'N' or direction == 'n':
        df = df.sort_values(by=col, ascending=True)
        df.index = range(0, len(df))
    else:
        print 'Please enter a relevant direction; P for positive and N for negative'
        sys.exit(1)
    return df

'''
Changes p-values and FDRs of 0 to <mininimum_value
Argument: List of p-values or FDRs
Return Value: String with appended minimum value
'''
def get_min(val_list):
    temp = val_list
    if len(temp) > 1:
        min_val = min(filter(lambda a: a != 0, temp))
        for i,n in enumerate(val_list):
            if n == 0:
                val_list[i] = '<'+ str(min_val)
    return val_list

'''
Calculates the q-value 
Argument: List of FDRs
Return Value: List of q-values
'''
def get_qval(val_list):
    temp = val_list
    sp_qval = list()
    if len(temp) > 1:
        for i,t in enumerate(temp[:-1]):
            qval = min(t,min(temp[(i+1):]))
            sp_qval.append(qval)
    sp_qval.append(temp[-1])
    return sp_qval

if __name__ == '__main__':
    args = get_parser().parse_args()
    inputfile = args.input_file
    input_df = pd.read_table(inputfile)
    cols = list(input_df.columns)[1:]
    cols = [re.sub('[^a-zA-Z0-9 \n\.]', '_', x) for x in cols]
    cols.insert(0,'Spacer Sequence')
    input_df.columns = cols
    cols_iter = cols[1:]
    ref = pd.read_table(args.chip_file)
    ref_colnames = list(ref.columns)
    ref_colnames[0:2] = ['Spacer Sequence', 'Gene Symbol']
    ref.columns = ref_colnames
    direction = args.dir
    thr = args.thr
    first_pert = args.use_first_pert
    for ci,c in enumerate(cols_iter):
        outputfile = c+'_STARSOutput_'+direction+'_'+str(thr)+'_'+str(datetime.now().strftime("%y-%m-%d-%H-%M"))+'.txt'
        st_in = input_df[['Spacer Sequence',c]]
        st_in = st_in.rename(columns={c:'Score'})
        merged = pd.merge(st_in, ref, on='Spacer Sequence')
        merged = sort_reindex(merged, 'Score', direction)
        st_in = sort_reindex(st_in, 'Score', direction)
        grouped = merged.groupby('Gene Symbol')
        sps = list(st_in['Spacer Sequence'])
        tot_sps = len(sps)

        sp_rank = {}
        for i in range(len(sps)):
            sp_rank[sps[i]] = i+1

        thr_rank = (thr*tot_sps)/100
        ge=list(merged['Gene Symbol'])
        ge=list(set(ge))

        i=1
        with open(outputfile,'w') as o:
            w = csv.writer(o, delimiter='\t', lineterminator='\n')
            w.writerow(('Gene Symbol', 'Number of perturbations', 'Ranks of perturbations', 'Perturbations','Most enriched perturbation', 'STARS Score', 'Average Score'))
            print 'Analyzing ...'
            for g in ge:
                #print i
                i=i+1
                ge_df = grouped.get_group(g)
                ge_df = ge_df.drop_duplicates('Spacer Sequence')
                gene = g
                spacers = list(ge_df['Spacer Sequence'])
                spacer_list = ';'.join(spacers)
                length_sps = len(spacers)
                rank = list()
                for s in spacers:
                    r=sp_rank[s]
                    rank.append(r)
                rank.sort()          
                all_ranks = ';'.join(str(r) for r in rank)
                if first_pert == 'N' or first_pert == 'n':
                    f_s = stats.binom.pmf(1, length_sps, rank[0]/float(tot_sps))  #Calculates probability of first perturbation
                    avg = -log10(f_s)                                             #Include value of first perturbation in Average Score
                    rank=rank[1:]                                                 #Ranks of perturbation #2 onwards
                    j=2
                    index_add=2
                elif first_pert == 'Y' or first_pert == 'y':
                    avg = 0
                    j=1    
                    index_add=1
                else:
                    print >> sys.stderr,'Please specify Y/y for value of first perturbation to be included else specify N/n.'
                    os.remove(outputfile)
                    sys.exit(1)
                indices = [r for r in rank if r <= thr_rank]
                if len(indices) == 0:
                    continue
                else:
                    sc_array = [stats.binom.pmf(indices.index(x)+index_add, length_sps, x/float(tot_sps)) for x in indices]
                    spacer_rank = sc_array.index(min(sc_array))+index_add
                    log_score = -log10(min(sc_array))
                    avg = (avg + sum([-log10(x) for x in sc_array]))/(len(indices)+1)
                    #sp_stars.append(log_score)
                    w.writerow((gene, length_sps, all_ranks, spacer_list,spacer_rank, log_score, avg))

        sp_stars = pd.read_table(outputfile)
        sp_stars = sp_stars.sort_values(by='STARS Score', ascending=False)
        stars_score = list(sp_stars['STARS Score'])
        if len(stars_score) == 0:
            print "No hit genes found for "+c
            #os.remove(outputfile)
            continue 
        else:
            null = [line.strip() for line in open(args.null)]
            null = [float(x) for x in null]
            null.sort(reverse=True)
            #Calculate the test statistic for your null distribution
            null_cdf = ECDF(null)
            min_n_cdf = 1-null_cdf(null[1])
            obs_cdf = ECDF(stars_score)
            sp_pval = list()
            sp_fdr = list()
            min_p = 0
            min_fdr = 0

            for ss in stars_score:
                n_cdf = null_cdf(ss-0.000001)
                o_cdf = obs_cdf(ss-0.000001)
                sp_pval.append((1-n_cdf))
                fdr = (1-n_cdf)/(1-o_cdf)
                if isnan(fdr) == True or isinf(fdr) == True:
                    fdr = 0
                if fdr > 1:
                    fdr = 1
                if fdr == 0:
                    min_p = (1-null_cdf(null[0]-0.000001))
                    min_fdr = min_p/(1-o_cdf)
                sp_fdr.append(fdr)
            #sp_pval_1 = getMin(sp_pval)
            #p_fdr_1 = getMin(sp_fdr)
            sp_stars['p-value'] = sp_pval
            sp_stars['FDR'] = sp_fdr
            qvals = get_qval(sp_fdr)
            sp_stars['q-value'] = qvals
            #min_fdr = min(filter(lambda a: a != 0, sp_fdr))
            with open(outputfile, 'w') as o:
                w = csv.writer(o,delimiter='\t', lineterminator='\n')
#                if min_p != 0:
#                    w.writerow((['p-value=0 means <'+str(min_p)+', FDR=0 means <'+str(min_fdr)]))
                w.writerow((sp_stars.columns))
                for row in sp_stars.iterrows():
                    w.writerow((row[1][0],row[1][1],row[1][2],row[1][3],row[1][4],row[1][5],row[1][6],row[1][7],row[1][8],row[1][9]))
            #sp_stars.to_csv(outputfile, sep='\t', index=None)

