#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : ngs
# AUTHOR  : codeunsolved@gmail.com
# CREATED : March 10 2018
# VERSION : v0.0.1a7
# UPDATE  : [v0.0.1a1] May 16 2018
# 1. add `sequence_complement()`;
# 2. add :AnnCoordinates:;
# UPDATE  : [v0.0.1a2] July 21 2018
# 1. :AnnCoordinates: add transcript and rank annotation;
# UPDATE  : [v0.0.1a3] August 30 2018
# 1. add :VcfParser:;
# 2. :AnnCoordinates: add gene name map option and multiple genes select policy;
# 3. :AnnCoordinates: amend `self.rank_hits` data structure to dict in list;
# UPDATE  : [v0.0.1a4] September 5 2018
# 1. :AnnCoordinates: add `exon_num`, `exon_total`, `cds_min` and `cds_max`;
# 2. :AnnCoordinates: optimize entire query pipeline;
# UPDATE  : [v0.0.1a5] September 18 2018
# 1. :AnnCoordinates: add `exon_start`, `exon_end`;
# UPDATE  : [v0.0.1a6] November 13 2018
# 1. add logger to :AnnCoordinates:;
# UPDATE  : [v0.0.1a7] November 27 2018
# 1. [BugFix] fix tag content extraction regex in `parse_meta()`;
# 2. [BugFix] fix bug in `extract_format()`;

import os
import re
from collections import defaultdict

import numpy as np

from .base import colour
from .base import color_term
from .connector import MysqlConnector


def sequence_complement(sequence, reverse=True):
    base_comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq_comp = ''.join([base_comp[x.upper()] for x in sequence])
    if reverse:
        seq_comp = seq_comp[::-1]
    return seq_comp


class VcfParser(object):
    def __init__(self, path_vcf):
        self.path_vcf = os.path.abspath(path_vcf)
        self.version = None
        self.header = ''
        self.meta_info = {}
        self.meta_format = {}
        self.th = None
        self.tbody = []
        self.PASS = None

        self.read_vcf()
        self.meta_filter = self.parse_meta('FILTER')
        self.meta_info = self.parse_meta('INFO')
        self.meta_format = self.parse_meta('FORMAT')
        self.meta_contig = self.parse_meta('contig')
        self.extract_info()
        self.extract_format()

        self.set_PASS()

    def read_vcf(self):
        def handle_th(line):
            th = re.sub('^#', '', line).rstrip('\n').split('\t')
            return th

        if not os.path.isfile(self.path_vcf):
            color_term("[VCFPARSER] Error! Invalid VCF file path: {}".format(self.path_vcf), 'ERR')
        else:
            with open(self.path_vcf, 'r') as vcf:
                for line in vcf:
                    if re.match('##fileformat', line):
                        self.version = re.search('=VCF(.+)\n', line).group(1)

                    if re.match('##', line):
                        self.header += line
                    elif re.match('#CHROM', line):
                        self.th = handle_th(line)
                    else:
                        line = line.rstrip('\n').split('\t')
                        self.tbody.append(dict(zip(self.th, line)))

    def parse_meta(self, tag):
        def handle_double_quotes_attr(key, content, line_dict):
            val = re.search('{}="([^"]+)"'.format(key), content).group(1)
            line_dict[key] = val.replace('"', '')
            content = re.sub('{}="[^"]+"'.format(key), '', content)
            content = content.rstrip(',')
            return content

        assert tag in ['FILTER', 'INFO', 'FORMAT', 'contig']

        meta_dict = {}
        for line in self.header.rstrip('\n').split('\n'):
            if re.match("##{}".format(tag), line):
                line_dict = {}
                content = re.search('=<(.+)>', line).group(1)

                # Handle 'Description', 'Source' and 'Version' separately,
                # because they could contain ',' in double-quotes.
                # Refer: [VCFv4.2]
                if re.search('Description="', content):
                    content = handle_double_quotes_attr('Description', content, line_dict)
                if re.search('Source="', content):
                    content = handle_double_quotes_attr('Source', content, line_dict)
                if re.search('Version="', content):
                    content = handle_double_quotes_attr('Version', content, line_dict)

                for attr in content.split(','):
                    k, v = attr.split('=')
                    line_dict[k] = v

                meta_dict[line_dict['ID']] = line_dict
        return meta_dict

    def extract_info(self):
        for row in self.tbody:
            info_dict = {}
            for entry in row['INFO'].split(';'):
                if re.search('=', entry):
                    k, v = entry.split('=')
                else:
                    k, v = entry, True
                info_dict[k] = v

                if k not in self.meta_info:
                    color_term("[VCFPARSER] Couldn't find INFO '{}' in header".format(k), 'WRN')
            row['INFO_dict'] = info_dict

    def extract_format(self):
        for row in self.tbody:
            for sample in self.th[9:]:
                format_keys = row['FORMAT'].split(':')
                format_vals = row[sample].split(':')
                sample_dict = dict(zip(format_keys, format_vals))
                for k, v in sample_dict.items():
                    if k not in self.meta_format:
                        color_term("[VCFPARSER] Couldn't find FORMAT '{}' in header".format(k), 'WRN')
                row[sample] = sample_dict

    def set_PASS(self):
        self.PASS = [row for row in self.tbody if row['FILTER'] == 'PASS']


class AnnCoordinates(object):

    def __init__(self, config=None, mysql_con=None,
                 gene_map={}, transcript_map={},
                 refer='hg19', ann_source='GENCODE',
                 use_cache=True, verbose=False, logger=color_term):
        # Configs
        self.config = config
        self.m_con = mysql_con

        # Options
        self.gene_map = gene_map
        self.transcript_map = transcript_map
        self.refer = refer
        self.ann_source = ann_source
        self.use_cache = use_cache
        self.verbose = verbose
        self.logger = logger

        self.ann_tab = None

        # Coordinates
        self.chr_ = None  # original chr
        self.chr_index = None
        self.chr = None
        self.pos = None

        # Annotations
        # Gene hits contains {key: val}:
        # 'id': Ensembl gene ID,
        # 'name': GENECODE gene name,
        # 'type': GENECODE gene type,
        self.gene_hits = []
        # Rank hits contains {key: val}:
        # 'id': Ensembl trasncript ID
        # 'strand': '+'/'-',
        # 'size': transcript size(end - start),
        # 't_start': transcript start
        # 't_end': transcript end
        # If an exon or part of an exon belongs to UTR,
        # it will annotate as UTR.
        # 'rank': 'Exon{N}', 'Intron{N}', "5'UTR", "3'UTR",
        # 'exon_num': {N},
        # 'exon_total': {N},
        # 'cds_min':  {N},
        # 'cds_max':  {N},
        # 'exon_start': exon start, it will be next exon start if rank is 'Intron{N}',
        # 'exon_end': exon_end, it will be last exon end if rank is 'Intron{N}',
        self.rank_hits = []
        self.gene_info = {}  # Choosen one in self.gene_hits
        self.rank_info = {}  # Choosen one in self.rank_hits
        self.default_transcript = None
        self.cache = defaultdict(dict)

        self.ann_tabs = {
            'hg19': {
                'GENCODE': 'GENCODE_Human_v28lift37_annotation_GFF3',
            }
        }

        # Check params validation
        self.check_refer()

        self.connect_db()
        assert self.m_con is not None

    def log(self, msg, level=None, verbose=None):
        if level is None:
            level = 'INFO'
        if verbose is None:
            verbose = self.verbose
        if verbose:
            self.logger(msg, level)

    def check_refer(self):
        assert self.refer in ['hg19', 'hg38']
        assert self.ann_source in ['GENCODE']

        self.ann_tab = self.ann_tabs[self.refer][self.ann_source]

    def check_chr(self):
        self.chr_index = re.sub('^chr', '', self.chr_.lower()).upper()
        self.chr = 'chr' + self.chr_index

        if not re.match('[12]?\d|[XY]$', self.chr_index):
            self.log("[CHECK_CHR] Unrecognized Chr: {}".format(self.chr_), 'WRN', verbose=True)

    def connect_db(self):
        if self.config:
            self.m_con = MysqlConnector(self.config, self.verbose)

    def init_ann(self):
        self.chr_ = None
        self.chr_index = None
        self.chr = None
        self.pos = None

        self.gene_hits = []
        self.rank_hits = []
        self.gene_info = {}
        self.rank_info = {}
        self.default_transcript = None

    def query(self, chr_, pos):
        self.init_ann()

        self.chr_ = chr_
        self.pos = int(pos)
        self.check_chr()

        key = "{}:{}".format(self.chr, self.pos)
        if self.use_cache and key in self.cache:
            self.log("[QUERY] Use cache for {}".format(key), 'WRN')
            self.gene_hits = self.cache[key]['gene_hits']
            self.rank_hits = self.cache[key]['rank_hits']
            self.gene_info = self.cache[key]['gene_info']
            self.rank_info = self.cache[key]['rank_info']
        else:
            try:
                self.query_gene()
            except Exception as e:
                self.log("[QUERY_GENE] Error! [{}:{}] {}".format(
                    self.chr, self.pos, e), 'ERR', verbose=True)
            else:
                self.log("• Query gene ... {}! {}hits genes: {}".format(
                    colour("OK", 'green'),
                    colour(len(self.gene_hits), 'green'),
                    [x['name'] for x in self.gene_hits]))

            if self.gene_info['id']:
                try:
                    self.query_gene_relative()
                except Exception as e:
                    self.log("[QUERY_GENE_RELATIVE] Error! [{}:{}] {}".format(
                        self.chr, self.pos, e), 'ERR', verbose=True)
                else:
                    self.log("• Query gene relative ... {}! {}hits transcripts: {}".format(
                            colour("OK", 'green'),
                            colour(len(self.rank_hits), 'green'),
                            [x['id'] for x in self.rank_hits]))

            if self.gene_info or self.rank:
                self.cache[key]['gene_hits'] = self.gene_hits
                self.cache[key]['rank_hits'] = self.rank_hits
                self.cache[key]['gene_info'] = self.gene_info
                self.cache[key]['rank_info'] = self.rank_info

    def get_attr_value(self, attr, key):
        s = re.search("{}=([^;]+);".format(key), attr)
        if s is not None:
            return s.group(1)
        else:
            self.log("[GET_ATTR_VALUE] No '{}' in attributes: {}:{} - ({})".format(
                key, self.chr, self.pos, attr), 'WRN')
            return None

    def query_gene(self):
        def get_gene_info(attr):
            gene_info = {'id': self.get_attr_value(attr, 'ID'),
                         'name': self.get_attr_value(attr, 'gene_name'),
                         'type': self.get_attr_value(attr, 'gene_type')}
            return gene_info

        def query_gene_intergenic():
            sql_l = "SELECT * FROM `{tab}` " \
                    "WHERE SeqID = %s AND Feature = 'gene' AND " \
                    "End < %s ORDER BY End DESC LIMIT 1".format(
                        tab=self.ann_tab)
            r_l = self.m_con.query(sql_l, [self.chr, self.pos]).fetchone()
            sql_r = "SELECT * FROM `{tab}` " \
                    "WHERE SeqID = %s AND Feature = 'gene' AND " \
                    "Start > %s ORDER BY Start LIMIT 1".format(
                        tab=self.ann_tab)
            r_r = self.m_con.query(sql_r, [self.chr, self.pos]).fetchone()

            gene_name_l = self.get_attr_value(r_l[-1], 'gene_name') if r_l else '-'
            gene_name_r = self.get_attr_value(r_r[-1], 'gene_name') if r_r else '-'
            name = "intergenic({}, {})".format(gene_name_l, gene_name_r)
            self.gene_info = {'id': None, 'name': name, 'type': 'intergenic'}

        sql = "SELECT * FROM `{tab}` " \
              "WHERE SeqID = %s AND Feature = 'gene' AND " \
              "Start <= %s AND End >= %s".format(
                tab=self.ann_tab)

        c = self.m_con.query(sql, [self.chr, self.pos, self.pos])
        r = c.fetchall()

        if len(r) == 0:
            self.log("[QUERY_GENE] No found for {}:{}, "
                     "try to annotate as intergenic region".format(
                        self.chr, self.pos), 'WRN')
            query_gene_intergenic()
        else:
            self.gene_hits = [get_gene_info(x[-1]) for x in r]
            self.gene_info = self.choose_gene()
            if len(r) > 1:
                self.log("[QUERY_GENE] {} entries({}) found for {}:{}, choose: {}".format(
                    len(r), ', '.join([x['name'] for x in self.gene_hits]),
                    self.chr, self.pos, self.gene_info['name']), 'WRN', verbose=True)
            if self.gene_info['name'] in self.gene_map:
                self.gene_info['name'] = self.gene_map[self.gene_info['name']]

    def query_gene_relative(self):
        def set_rank_hits(r):
            for entry in r:
                feature = entry[3]
                attr = entry[-1]
                if feature == 'transcript':
                    transcript_id = self.get_attr_value(attr, 'ID')
                    start = entry[4]
                    end = entry[5]
                    strand = entry[7]
                    size = end - start
                    if transcript_id is not None:
                        self.rank_hits.append({'id': transcript_id, 'strand': strand, 'size': size,
                                               't_start': start, 't_end': end,
                                               'rank': None, 'exon_num': None, 'exon_total': None,
                                               'cds_min': None, 'cds_max': None,
                                               'exon_start': None, 'exon_end': None})

        def set_rank(r):
            def gen_rank(feature, exon_num, exon=True):
                if feature == 'CDS':
                    if exon:
                        return "Exon{}".format(exon_num)
                    else:
                        return "Intron{}".format(exon_num)
                elif feature == 'five_prime_UTR':
                    return "5'UTR"
                elif feature == 'three_prime_UTR':
                    return "3'UTR"
                else:
                    return feature

            for x in self.rank_hits:
                transcript_id = x['id']

                loc_flag = None
                for entry in r:
                    feature = entry[3]
                    start = entry[4]
                    end = entry[5]
                    strand = entry[7]
                    attr = entry[-1]
                    if not re.search("transcript_id={}".format(transcript_id), attr):
                        continue

                    if feature in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
                        if start <= self.pos and self.pos <= end:
                            exon_num = int(self.get_attr_value(attr, 'exon_number'))
                            rank = gen_rank(feature, exon_num)
                            loc_flag = 0
                            break
                        else:
                            if self.pos < start:
                                if loc_flag == 1:
                                    exon_num = int(self.get_attr_value(attr, 'exon_number'))
                                    if strand == '+':
                                        exon_num -= 1
                                    rank = gen_rank(feature, exon_num, exon=False)
                                    loc_flag = 0
                                    break
                                else:
                                    loc_flag = -1
                            elif end < self.pos:
                                loc_flag = 1
                if loc_flag != 0:
                    if transcript_id == self.default_transcript:
                        self.log("[SET_RANK] [{}:{}] seems {} than the range of transcript: {} {}-{}".format(
                            self.chr, self.pos,
                            'LESS' if loc_flag == -1 else 'LARGE',
                            transcript_id, x['t_start'], x['t_end']), 'WRN')
                    exon_num = 'N/A'
                    rank = 'N/A'
                x['rank'] = rank
                x['exon_num'] = exon_num

        def set_exon_and_cds(r):
            for x in self.rank_hits:
                transcript_id = x['id']

                max_exon_num = 0
                min_cds_num = np.inf
                max_cds_num = 0
                for entry in r:
                    feature = entry[3]
                    start = entry[4]
                    end = entry[5]
                    attr = entry[-1]
                    if not re.search("transcript_id={}".format(transcript_id), attr):
                        continue

                    if feature == 'exon':
                        exon_num = int(self.get_attr_value(attr, 'exon_number'))
                        if exon_num > max_exon_num:
                            max_exon_num = exon_num

                        if re.match('Intron', x['rank']):
                            if exon_num == x['exon_num']:
                                x['exon_end'] = end
                            elif exon_num == x['exon_num']+1:
                                x['exon_start'] = start
                        else:
                            if exon_num == x['exon_num']:
                                x['exon_start'] = start
                                x['exon_end'] = end
                    elif feature == 'CDS':
                        cds_num = int(self.get_attr_value(attr, 'exon_number'))
                        if cds_num < min_cds_num:
                            min_cds_num = cds_num
                        if cds_num > max_cds_num:
                            max_cds_num = cds_num
                x['exon_total'] = max_exon_num
                x['cds_min'] = min_cds_num
                x['cds_max'] = max_cds_num

        sql = "SELECT * FROM `{tab}` " \
              "WHERE SeqID = %s AND Feature IN " \
              "('transcript', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR') " \
              "AND Attribute LIKE %s ORDER BY Start".format(
                tab=self.ann_tab)

        c = self.m_con.query(sql, [self.chr, "%gene_id={}%".format(self.gene_info['id'])])
        r = c.fetchall()
        if len(r):
            set_rank_hits(r)
            self.default_transcript = self.choose_transcript()
            set_rank(r)
            set_exon_and_cds(r)
            self.set_rank_info()
        else:
            self.log("[QUERY_GENE_RELATIVE] Found no relative entry for '{}'({})".format(
                self.gene_info['name'], self.gene_info['id']), 'ERR')

    def choose_gene(self):
        # Prefer (first) 'protein_coding' gene, or first gene
        for x in self.gene_hits:
            gene_type = x['type']
            if gene_type == 'protein_coding':
                return x
        return self.gene_hits[0]

    def choose_transcript(self):
        # Prefer (longest) transcript in `transcript_map`, or longest transcript
        rank_hits_sorted = sorted(self.rank_hits, key=lambda x: x['size'], reverse=True)
        for x in rank_hits_sorted:
            transcript_id = x['id']
            transcript_id_main = re.sub('\.\d+$', '', transcript_id)
            if transcript_id_main in self.transcript_map:
                return x['id']
        self.log("[CHOOSE_TRANSCRIPT] Use longest transcript as default: {} for {}".format(
            rank_hits_sorted[0]['id'], self.gene_info['name']), 'WRN')
        return rank_hits_sorted[0]['id']

    def set_rank_info(self):
        for x in self.rank_hits:
            if x['id'] == self.default_transcript:
                self.rank_info = x

    def get_gene_name(self):
        if 'name' in self.gene_info:
            return self.gene_info['name']
        else:
            return None

    def get_rank_info(self):
        if len(self.rank_info):
            return self.rank_info
        else:
            return None
