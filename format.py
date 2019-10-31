#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import csv
import gzip

print '\t'.join(['Locus', 'dbSNP', 'Ref', 'Alt', 'Type', 'Filter', \
        'Gene', 'Transcript', 'Exon', 'DNA_change', 'protein_change', 'CGI_input', 'MAF', \
        'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score', 'Polyphen2_HDIV_pred', \
        'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred', 'FATHMM_score', 'FATHMM_pred', \
        'gnomAD_exomes_AF', 'gnomAD_exomes_AC', 'gnomAD_exomes_AN', \
        'clinvar_clnsig', 'clinvar_trait', 'Cancer_Gene_Census', \
        'Biocarta_Pathway', 'KEGG_Pathway'])

census_file = '/home/pcusco/lab/pipelines/cancer_gene_census.csv'
pathways_file = '/home/pcusco/data/hg38/dbnsfp/dbNSFP4.0_gene.gz'
exonic_variant_function = '4_ann/' + sys.argv[1] + '.exonic_variant_function'
variant_function = '4_ann/' + sys.argv[1] + '.variant_function'

def parse_variants(line, offset):
    chrom = sline[2 + offset].lstrip('chr')
    start = sline[3 + offset]
    end = sline[4 + offset]
    locus = '%s:%s-%s' % (chrom, start, end)
    dbsnp = sline[9 + offset]
    ref = sline[5 + offset]
    alt = sline[6 + offset]
    var_type = sline[0 + offset]
    filter_type = sline[13 + offset]

    # TODO: remove split when bcftools norm implemented.
    annotation = sline[1 + offset].split(',', 1)[0]
    protein, cgi = ['NA'] * 2
    if '(' in annotation:
        gene, annotation = annotation.split('(')
        annotation = annotation.rstrip(')')
        transcript, exon, dna = annotation.split(':')
    else:
        split_annotation = annotation.split(':')
        if len(split_annotation) == 5:
            gene, transcript, exon, dna, protein = split_annotation
            cgi = '%s:%s' % (gene, protein.lstrip('p.'))
        elif len(split_annotation) == 4:
            gene, transcript, exon, dna = split_annotation
        elif len(split_annotation) == 3 and 'wholegene' in split_annotation:
            gene, transcript, _ = split_annotation
            exon = dna = split_annotation[-1]
        elif len(split_annotation) == 1:
            gene = split_annotation[0]
            transcript, exon, dna = ['NA'] * 3
        else:
            raise ValueError('Malformed annotation: %s' % annotation)

    format_field = sline[16 + offset].split(':')
    maf = format_field[2].split(',', 1)[0] # TODO: remove split when bcftools norm implemented.

    info_field = dict((e.split('=') for e in sline[14 + offset].split(';') if '=' in e))

    sift_score = info_field['SIFT_score'] if 'SIFT_score' in info_field else 'NA'
    soft_pred = info_field['SIFT_pred'] if 'SIFT_pred' in info_field else 'NA'
    hdiv_score = info_field['Polyphen2_HDIV_score'] if 'Polyphen2_HDIV_score' in info_field else 'NA'
    hdiv_pred = info_field['Polyphen2_HDIV_pred'] if 'Polyphen2_HDIV_pred' in info_field else 'NA'
    hvar_score = info_field['Polyphen2_HVAR_score'] if 'Polyphen2_HVAR_score' in info_field else 'NA'
    hvar_pred = info_field['Polyphen2_HVAR_pred'] if 'Polyphen2_HVAR_pred' in info_field else 'NA'
    fathmm_score = info_field['FATHMM_score'] if 'FATHMM_score' in info_field else 'NA'
    fathmm_pred = info_field['FATHMM_pred'] if 'FATHMM_pred' in info_field else 'NA'
    gnomad_af = info_field['gnomAD_exomes_AF'] if 'gnomAD_exomes_AF' in info_field else 'NA'
    gnomad_ac = info_field['gnomAD_exomes_AC'] if 'gnomAD_exomes_AC' in info_field else 'NA'
    gnomad_an = info_field['gnomAD_exomes_AN'] if 'gnomAD_exomes_AN' in info_field else 'NA'
    clinvar_clnsig = info_field['clinvar_clnsig'] if 'clinvar_clnsig' in info_field else 'NA'
    clinvar_trait = info_field['clinvar_trait'] if 'clinvar_trait' in info_field else 'NA'

    census = 'Tier %s' % cgc[gene] if gene in cgc else 'NA'
    biocarta, kegg = pathways[gene] if gene in pathways else ('NA', 'NA')
    if biocarta == '.': biocarta = 'NA'
    if kegg == '.': kegg = 'NA'

    return '\t'.join([locus, dbsnp, ref, alt, var_type, filter_type, \
            gene, transcript, exon, dna, protein, cgi, maf, \
            sift_score, soft_pred, hdiv_score, hdiv_pred, hvar_score, hvar_pred, \
            fathmm_score, fathmm_pred, gnomad_af, gnomad_ac, gnomad_an, \
            clinvar_clnsig, clinvar_trait, census, \
            biocarta, kegg])

cgc = {}
with open(census_file) as fin:
    for i, row in enumerate(csv.reader(fin)):
        if i == 0: continue
        cgc[row[0]] = row[4]

pathways = {}
with gzip.open(pathways_file) as fin:
    _ = fin.readline()
    for line in fin:
        sline = line.split('\t')
        gene = sline[0]
        biocarta = sline[15]
        kegg = sline[19]
        pathways[gene] = (biocarta, kegg)

with open(exonic_variant_function) as fin:
    for line in fin:
        sline = line.split('\t')

        var_type = sline[1]
        #if var_type == 'synonymous SNV':
        if var_type == 'unknown':
            continue
        filter_stat = sline[14]
        if filter_stat == 'PASS' or filter_stat == 'base_quality':
            print parse_variants(sline, 1)

with open(variant_function) as fin:
    for line in fin:
        sline = line.split('\t')

        var_type = sline[0]
        filter_stat = sline[13]
        if var_type == 'splicing' and (filter_stat == 'PASS' or filter_stat == 'base_quality'):
            print parse_variants(sline, 0)
