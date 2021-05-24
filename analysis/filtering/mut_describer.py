'''
mut_describer.py - create detailed table of muts
'''

import csv
import re
from glob import glob
import annotation_table
import ness_vcf
import mutation
import vcf # pyvcf, not cyvcf2
import argparse
from tqdm import tqdm
from Bio import SeqIO

def args():
    parser = argparse.ArgumentParser(
        description='create detailed table of muts', 
        usage='python mut_describer.py [options]')

    parser.add_argument('-f', '--fname', required=True,
        type=str, help='Tab-separated file containing input muts')
    parser.add_argument('-v', '--vcf_path', required=True,
        type=str, help='Path to dir containing sample-specific VCF files')
    parser.add_argument('-r', '--ref_fasta', required=True,
        type=str, help='Path to reference genome FASTA')
    parser.add_argument('-t', '--ant_file', required=True,
        type=str, help='Path to annotation table')
    parser.add_argument('-p', '--pop_vcf_file', required=False,
        type=str, help='Path to VCF file for pi calc')
    parser.add_argument('-o', '--outname', required=True,
        type=str, help='File to write to')

    args = parser.parse_args()

    return args.fname, args.vcf_path, args.ref_fasta, args.ant_file, \
        args.pop_vcf_file, args.outname

def describe_mut(d_in, vcf_path, ref_fasta, ant_file, pop_vcf_file) -> dict:
    """
    convert single mutation line (as list) to 'mut describer format'

    mut list should be in format
    [fname, chrom, pos, ref, alt, gt_bases, gt_quals, gt_depths, GQ1, GQ2]

    Parameters
    ----------
    d_in : dict
        mut info from csv.DictReader
    vcf_path : str
        path to dir containing per-sample VCF files
    ref_fasta : str
        path to reference genome FASTA file
    ant_file : str
        path to annotation table file - should be tabixed
    pop_vcf_file : str
        path to VCF file to be used for pi calc [optional]

    Returns
    -------
    mut_dict : dict
        dict with 'mut describer fields' to be written to outfile
    """
    if not isinstance(d_in, dict):
        raise TypeError()
    ref_dict = SeqIO.to_dict(SeqIO.parse(ref_fasta, 'fasta'))
    d = {}

    # create mutation.mutation object
    fname = d_in['fname']
    vcf_file = glob(vcf_path + fname.replace('samples', '') + '*.vcf.gz')[0] # ugh
    tqdm.write(vcf_file)
    chrom = d_in['chrom']
    pos = int(d_in['pos'])
    mut = mutation.mutation(chrom, pos)

    # populate dict
    mut.ref = d_in['ref']
    mut.alt = eval(d_in['alt'])
    
    d['chromosome'] = chrom
    d['position'] = pos
    d['ref'] = d_in['ref']
    d['alt'] = eval(d_in['alt'])
    d['qual'] = d_in['QUAL']
    d['MQ'] = d_in['MQ']

    # mut properties
    try:
        d['mutation'] = '>'.join(mut.mutation(vcf_file, mut_quality=10))
    except:
        print(fname)
        print(chrom, pos)
        print(mut.vcf_record(vcf_file))
        print(mut.vcf_record(vcf_file).samples)
        print(chrom, pos, d_in['fname'])
    d['type'] = mut.mutation_type(vcf_file)
    d['mutant_sample'] = mut.mutant_sample(vcf_file, mut_quality=10)

    # purity
    purity = mut.purity(vcf_file)
    if len(purity) > 0:
        d['min_purity'] = min(purity)
    else:
        d['min_purity'] = None

    # variant stats
    d['het_count'] = mut.het_count(vcf_file)
    d['number_mutants'] = mut.number_mutants(vcf_file)
    d['number_genotypes'] = mut.number_genotypes(vcf_file)

    # genomic stats
    d['gc20'] = mut.gc_content(ref_dict, 10)
    d['gc2000'] = mut.gc_content(ref_dict, 1000)
    if pop_vcf_file:
        d['pi1000'] = mut.pi(pop_vcf_file, 500, ref_dict)

    # nonsynonymous vs synonymous changes
    ant_pos = mut.annotation_position(ant_file)
    if ant_pos.CDS == '1' and d['type'] == 'SNP':
        d['nonsyn_v_syn'], d['anc_codon'], d['mut_codon'] = mut.nonsynonymous(ant_file, vcf_file)
    else:
        d['nonsyn_v_syn'], d['anc_codon'], d['mut_codon'] = None, None, None
    d['nonsense'] = mut.nonsense(vcf_file, ant_file)

    ant_cols = [
    'chromosome', 'position', 'reference_base', 'genic', 'exonic',
    'intronic', 'intergenic', 'utr5', 'utr3', 'fold0', 'fold4',
    'fold2', 'fold3', 'CDS', 'mRNA', 'rRNA', 'tRNA', 'feature_names',
    'feature_types', 'feature_ID', 'cds_position', 'strand', 'frame',
    'codon', 'aa', 'degen', 'FPKM', 'rho', 'FAIRE', 'recombination',
    'mutability', 'all_quebec_alleles']
    ant_vals = ant_pos.output_line().split('\t')
    for i, _ in enumerate(ant_cols):
        d[ant_cols[i]] = ant_vals[i]

    return d


def parse_muts(fname, vcf_path, ref_fasta, ant_file, pop_vcf_file, outname) -> None:
    """
    parse input muts and 'describe' them

    passes all args besides fname and outname to `describe_muts`

    Parameters
    ----------
    fname : str
        path to tsv containing mut info (generated from filter scripts)
    vcf_path : str
        path to dir containing per-sample VCF files
    ref_fasta : str
        path to reference genome FASTA file
    ant_file : str
        path to annotation table file - should be tabixed
    pop_vcf_file : str
        path to VCF file to be used for pi calc [optional]
    outname : str
        file to write to

    Returns
    -------
    mut_dict : dict
        dict with 'mut describer fields' to be written to outfile
    """
    if not vcf_path.endswith('/'):
        vcf_path += '/'

    fieldnames = [
    'chromosome', 'position', 'ref', 'alt', 'qual', 'MQ',
    'mutation', 'type', 'mutant_sample', 'min_purity',
    'het_count', 'number_mutants', 'number_genotypes',
    'gc20', 'gc2000', 'nonsyn_v_syn', 'anc_codon',
    'mut_codon', 'nonsense']

    if pop_vcf_file:
        fieldnames.append('pi1000')

    # may become optional, so keeping separate
    ant_cols = [
    'reference_base', 'genic', 'exonic',
    'intronic', 'intergenic', 'utr5', 'utr3', 'fold0', 'fold4',
    'fold2', 'fold3', 'CDS', 'mRNA', 'rRNA', 'tRNA', 'feature_names',
    'feature_types', 'feature_ID', 'cds_position', 'strand', 'frame',
    'codon', 'aa', 'degen', 'FPKM', 'rho', 'FAIRE', 'recombination',
    'mutability', 'all_quebec_alleles']
    fieldnames.extend(ant_cols)

    total_lines = 0
    with open(fname, 'r') as f:
        for line in f:
            total_lines += 1

    with open(outname, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        with open(fname, 'r', newline='') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            for mut in tqdm(reader, total=total_lines):
                mut_dict = describe_mut(mut, vcf_path, ref_fasta, ant_file, pop_vcf_file)
                if mut_dict:
                    writer.writerow(mut_dict)
        
def main():
    parse_muts(*args())

if __name__ == '__main__':
    main()

        

