#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : ngs
# AUTHOR  : codeunsolved@gmail.com
# CREATED : March 10 2018
# VERSION : v0.0.9
# UPDATE  : [v0.0.1] May 16 2018
# 1. add `sequence_complement()`;
# 2. add :AnnCoordinates:;
# UPDATE  : [v0.0.2] July 21 2018
# 1. :AnnCoordinates: add transcript and rank annotation;
# UPDATE  : [v0.0.3] August 30 2018
# 1. add :VcfParser:;
# 2. :AnnCoordinates: add gene name map option and multiple genes select policy;
# 3. :AnnCoordinates: amend `self.rank_hits` data structure to dict in list;
# UPDATE  : [v0.0.4] September 5 2018
# 1. :AnnCoordinates: add `exon_num`, `exon_total`, `cds_min` and `cds_max`;
# 2. :AnnCoordinates: optimize entire query pipeline;
# UPDATE  : [v0.0.5] September 18 2018
# 1. :AnnCoordinates: add `exon_start`, `exon_end`;
# UPDATE  : [v0.0.6] November 13 2018
# 1. add logger to :AnnCoordinates:;
# UPDATE  : [v0.0.7] November 27 2018
# 1. [BugFix] fix tag content extraction regex in `parse_meta()`;
# 2. [BugFix] fix bug in `extract_format()`;
# UPDATE  : [v0.0.8] April 8 2019
# 1. enhance :VcfParser:;
# UPDATE  : [v0.0.9] May 10 2019
# 1. enhance :AnnCoordinates: by using cache for entire gene region;

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
    class Row(object):
        def __init__(self, th, row):
            self.th = th
            self.row = row
            self.SAMPLES = []

            self.init_attr()

        def __repr__(self):
            content = "[CHROM: {};POS: {}; REF: {};ALT: {}]".format(
                self.CHROM, self.POS, self.REF, self.ALT)
            return content

        def init_attr(self):
            class Obj(object):
                def __init__(self, raw):
                    self.raw = raw

                def __repr__(self):
                    return self.raw

            def handle_info(info):
                obj = Obj(info)
                for entry in info.split(';'):
                    if re.search('=', entry):
                        k, v = entry.split('=')
                    else:
                        k, v = entry, True
                    setattr(obj, k, v)
                return obj

            def handle_format(fmt):
                format_keys = self.FORMAT.split(':')
                format_vals = fmt.split(':')
                obj = Obj(fmt)
                for k, v in zip(format_keys, format_vals):
                    setattr(obj, k, v)
                return obj

            assert isinstance(self.th, list)
            assert isinstance(self.row, list)

            for k in ['CHROM', 'POS', 'REF', 'ALT', 'FORMAT']:
                self.set_attr(k, '')

            for k, v in zip(self.th, self.row):
                if k in ['CHROM', 'POS', 'ID',
                         'REF', 'ALT',
                         'QUAL', 'FILTER',
                         'FORMAT']:
                    self.set_attr(k, v)
                elif k == 'INFO':
                    obj = handle_info(v)
                    self.set_attr(k, obj)
                else:
                    self.SAMPLES.append(k)
                    obj = handle_format(v)
                    self.set_attr(k, obj)

            delattr(self, 'th')
            delattr(self, 'row')

        def set_attr(self, k, v):
            setattr(self, k, v)

    def __init__(self, path_vcf):
        self.path_vcf = os.path.abspath(path_vcf)
        self.version = None
        self.header = ''
        self.meta_filter = {}
        self.meta_info = {}
        self.meta_format = {}
        self.meta_contig = {}
        self.th = None
        self.tbody = []
        self.PASS = []

        self.read_vcf()

        self.meta_filter = self.parse_meta('FILTER')
        self.meta_info = self.parse_meta('INFO')
        self.meta_format = self.parse_meta('FORMAT')
        self.meta_contig = self.parse_meta('contig')

        self.check_meta('INFO', self.meta_info)
        self.check_meta('FORMAT', self.meta_format)

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
                        row = self.Row(self.th, line)
                        self.tbody.append(row)

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

    def check_meta(self, field, meta):
        missing = set()
        for row in self.tbody:
            field_obj = getattr(row, field)
            for k in dir(field_obj):
                if re.match('(?:__.+__|raw)$', k):
                    continue
                if callable(getattr(field_obj, k)):
                    continue
                if k not in meta:
                    missing.add(k)
        for k in missing:
            color_term("[VCFPARSER] Couldn't find {} '{}' in header".format(
                field, k), 'WRN')

    def set_PASS(self):
        self.PASS = [row for row in self.tbody if row.FILTER == 'PASS']


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
        self.chr_raw = None  # original chr
        self.chr_index = None
        self.chr = None
        self.pos = None

        # Annotations
        # Gene hits contains {key: val}:
        # 'id': Ensembl gene ID
        # 'name': GENECODE gene name
        # 'type': GENECODE gene type
        # 'chr': chromosome
        # 'source': entry source
        # 'start': {N}
        # 'end': {N}
        # 'strand': +/-/None
        self.gene_hits = []
        self.gene_info = {}  # Choosen one in self.gene_hits
        # Transcript hits contains {key: val}:
        # 'id': Ensembl trasncript ID
        # 'start': transcript start
        # 'end': transcript end
        # 'size': transcript size(end - start)
        # 'strand': +/-
        # 'components': all exon, CDS, five_prime_UTR, three_prime_UTR info
        # 'rank_info': see below
        self.trans_hits = []
        # Rank info  contains {key: val}:
        # 'trans_id': Ensembl trasncript ID
        # 'trans_start': transcript start
        # 'trans_end': transcript end
        # 'trans_strand': +/-
        # If an exon or part of an exon belongs to UTR,
        # it will annotate as UTR.
        # 'rank': 'Exon{N}', 'Intron{N}', "5'UTR", "3'UTR"
        # 'exon_num': {N}
        # 'exon_total': {N}
        # 'cds_min': {N}
        # 'cds_max': {N}
        # 'exon_start': exon start, it will be next exon start if rank is 'Intron{N}'
        # 'exon_end': exon end, it will be last exon end if rank is 'Intron{N}'
        self.rank_info = {}
        self.default_transcript = None
        self.cache_coord = defaultdict(dict)
        self.cache_region = []

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
        self.chr_index = re.sub('^chr', '', self.chr_raw.lower()).upper()
        self.chr = 'chr' + self.chr_index

        if not re.match(r'[12]?\d|[XY]$', self.chr_index):
            self.log("[CHECK_CHR] Unrecognized Chr: {}".format(self.chr_raw), 'WRN', verbose=True)

    def connect_db(self):
        if self.config:
            self.m_con = MysqlConnector(self.config, self.verbose)

    def init_ann(self):
        self.chr_raw = None
        self.chr_index = None
        self.chr = None
        self.pos = None

        self.gene_hits = []
        self.gene_info = {}
        self.trans_hits = []
        self.rank_info = {}
        self.default_transcript = None

    def query(self, chr_, pos):
        def match_cache_region(chr_, pos):
            for region in self.cache_region:
                if region['chr'] == chr_ and \
                   region['start'] <= pos and pos <= region['end']:
                    return region
            return None

        def _query_gene():
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

        def _query_trans():
            try:
                self.query_trans()
            except Exception as e:
                self.log("[QUERY_TRANS] Error! [{}:{}] {}".format(
                    self.chr, self.pos, e), 'ERR', verbose=True)
            else:
                self.log("• Query transcripts {}! {}hits transcripts: {}".format(
                        colour("OK", 'green'),
                        colour(len(self.trans_hits), 'green'),
                        [x['id'] for x in self.trans_hits]))

        def _cache_coord(coord):
            self.cache_coord[coord]['gene_hits'] = self.gene_hits
            self.cache_coord[coord]['gene_info'] = self.gene_info
            self.cache_coord[coord]['trans_hits'] = self.trans_hits
            self.cache_coord[coord]['rank_info'] = self.rank_info

        def _cache_region():
            region = match_cache_region(self.chr, self.pos)
            if region is None:
                region = {}
                region['chr'] = self.chr
                region['start'] = self.gene_info['start']
                region['end'] = self.gene_info['end']
                region['gene_info'] = self.gene_info
                region['trans_hits'] = self.trans_hits
                self.cache_region.append(region)

        self.init_ann()

        self.chr_raw = chr_
        self.pos = int(pos)
        self.check_chr()

        coord = "{}:{}".format(self.chr, self.pos)
        if self.use_cache and coord in self.cache_coord:
            self.log("[QUERY] Use coord cache for {}".format(coord), 'WRN')
            self.gene_hits = self.cache_coord[coord]['gene_hits']
            self.trans_hits = self.cache_coord[coord]['trans_hits']
            self.gene_info = self.cache_coord[coord]['gene_info']
            self.rank_info = self.cache_coord[coord]['rank_info']
        else:
            region = match_cache_region(self.chr, self.pos)

            if self.use_cache and region is not None:
                self.log("[QUERY] Use region cache for {} hit by {}".format(
                    coord, region['gene_info']['name']), 'WRN')
                self.gene_info = region['gene_info']
                if self.gene_info['id']:
                    self.trans_hits = region['trans_hits']
                    self.set_rank_info()
            else:
                _query_gene()

                if self.gene_info['id']:
                    _query_trans()

            # Cache
            _cache_coord(coord)
            _cache_region()

    def query_gene(self):
        def gen_r():
            sql = "SELECT * FROM `{tab}` " \
                  "WHERE SeqID = %s AND Feature = 'gene' AND " \
                  "Start <= %s AND End >= %s".format(
                    tab=self.ann_tab)

            c = self.m_con.query(sql, [self.chr, self.pos, self.pos])
            r = c.fetchall()
            return r

        def gen_gene_hits(entries):
            gene_hits = []
            for entry in entries:
                attr = entry[-1]
                gene_hits.append({
                    'id': self._get_attr_value(attr, 'ID'),
                    'name': self._get_attr_value(attr, 'gene_name'),
                    'type': self._get_attr_value(attr, 'gene_type'),
                    'chr': entry[1],
                    'source': entry[2],
                    'start': entry[4],
                    'end': entry[5],
                    'strand': entry[7],
                })
            return gene_hits

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

            gene_name_l = '-' if r_l is None else self._get_attr_value(r_l[-1], 'gene_name')
            gene_name_r = '-' if r_r is None else self._get_attr_value(r_r[-1], 'gene_name')
            start = None if r_l is None else r_l[5]  # Left gene end
            end = None if r_r is None else r_r[4]    # Right gene start
            name = "intergenic({}, {})".format(gene_name_l, gene_name_r)
            self.gene_info = {
                'id': None,
                'name': name,
                'type': 'intergenic',
                'chr': self.chr,
                'start': start,
                'end': end,
                'strand': None,
            }

        r = gen_r()

        if len(r):
            self.gene_hits = gen_gene_hits(r)
            self.gene_info = self._choose_gene_info()
            if len(r) > 1:
                self.log("[QUERY_GENE] {} entries({}) found for {}:{}, choose: {}".format(
                    len(r), ', '.join([x['name'] for x in self.gene_hits]),
                    self.chr, self.pos, self.gene_info['name']), 'WRN', verbose=True)
            # Map gene name
            if self.gene_info['name'] in self.gene_map:
                self.gene_info['name'] = self.gene_map[self.gene_info['name']]
        else:
            self.log("[QUERY_GENE] No found for {}:{}, "
                     "try to annotate as intergenic region".format(
                        self.chr, self.pos), 'WRN')
            query_gene_intergenic()

    def query_trans(self):
        def gen_r():
            sql = "SELECT * FROM `{tab}` " \
                  "WHERE SeqID = %s AND Feature IN " \
                  "('transcript', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR') " \
                  "AND Attribute LIKE %s ORDER BY Start".format(
                    tab=self.ann_tab)

            c = self.m_con.query(sql, [self.chr, r"%gene_id={}%".format(self.gene_info['id'])])
            r = c.fetchall()
            return r

        def gen_trans_hits(entries):
            trans_hits = {}
            for entry in entries:
                feature = entry[3]
                attr = entry[-1]
                if feature == 'transcript':
                    transcript_id = self._get_attr_value(attr, 'ID')
                    if transcript_id is not None:
                        trans_hits[transcript_id] = {
                            'id': transcript_id,
                            'start': entry[4],
                            'end': entry[5],
                            'size': entry[5]-entry[4],
                            'strand': entry[7],
                            'components': [],
                            'rank_info': {
                                'trans_id': None,
                                'trans_start': None,
                                'trans_end': None,
                                'trans_strand': None,
                                'rank': None,
                                'exon_num': None,
                                'exon_total': None,
                                'cds_min': None,
                                'cds_max': None,
                                'exon_start': None,
                                'exon_end': None,
                            },
                        }

            for entry in entries:
                attr = entry[-1]
                transcript_id = self._get_attr_value(attr, 'transcript_id')
                transcript_id = re.sub(r'_\d+$', '', transcript_id)

                trans_hits[transcript_id]['components'].append({
                    'feature': entry[3],
                    'start': entry[4],
                    'end': entry[5],
                    'strand': entry[7],
                    'attr': entry[-1],
                })
            return list(trans_hits.values())

        r = gen_r()

        if len(r):
            self.trans_hits = gen_trans_hits(r)
            self.set_rank_info()
        else:
            self.log("[QUERY_TRANS] Found no relative entry for '{}'({})".format(
                self.gene_info['name'], self.gene_info['id']), 'ERR')

    def set_rank_info(self):
        def set_rank():
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

            for tran in self.trans_hits:
                transcript_id = tran['id']
                rank_info = tran['rank_info']

                loc_flag = None
                for entry in tran['components']:
                    feature = entry['feature']
                    start = entry['start']
                    end = entry['end']
                    strand = entry['strand']
                    attr = entry['attr']

                    if feature in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
                        if start <= self.pos and self.pos <= end:
                            exon_num = int(self._get_attr_value(attr, 'exon_number'))
                            rank = gen_rank(feature, exon_num)
                            loc_flag = 0
                            break
                        else:
                            if self.pos < start:
                                if loc_flag == 1:
                                    exon_num = int(self._get_attr_value(attr, 'exon_number'))
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
                            transcript_id, tran['start'], tran['end']), 'WRN')
                    rank = 'N/A'
                    exon_num = 'N/A'
                rank_info['rank'] = rank
                rank_info['exon_num'] = exon_num

        def set_exon_and_cds():
            for tran in self.trans_hits:
                rank_info = tran['rank_info']

                max_exon_num = 0
                min_cds_num = np.inf
                max_cds_num = 0
                for entry in tran['components']:
                    feature = entry['feature']
                    start = entry['start']
                    end = entry['end']
                    attr = entry['attr']

                    if feature == 'exon':
                        exon_num = int(self._get_attr_value(attr, 'exon_number'))
                        if exon_num > max_exon_num:
                            max_exon_num = exon_num

                        if re.match('Intron', rank_info['rank']):
                            if exon_num == rank_info['exon_num']:
                                rank_info['exon_end'] = end
                            elif exon_num == rank_info['exon_num']+1:
                                rank_info['exon_start'] = start
                        else:
                            if exon_num == rank_info['exon_num']:
                                rank_info['exon_start'] = start
                                rank_info['exon_end'] = end
                    elif feature == 'CDS':
                        cds_num = int(self._get_attr_value(attr, 'exon_number'))
                        if cds_num < min_cds_num:
                            min_cds_num = cds_num
                        if cds_num > max_cds_num:
                            max_cds_num = cds_num
                rank_info['exon_total'] = max_exon_num
                rank_info['cds_min'] = min_cds_num
                rank_info['cds_max'] = max_cds_num

        def set_trans_info():
            for tran in self.trans_hits:
                rank_info = tran['rank_info']
                rank_info['trans_id'] = tran['id']
                rank_info['trans_start'] = tran['start']
                rank_info['trans_end'] = tran['end']
                rank_info['trans_strand'] = tran['strand']

        self.default_transcript = self._choose_default_transcript()
        set_rank()
        set_exon_and_cds()
        set_trans_info()
        self.rank_info = self._choose_rank_info()

    def _get_attr_value(self, attr, key):
        s = re.search("{}=([^;]+);".format(key), attr)
        if s is not None:
            return s.group(1)
        else:
            self.log("[GET_ATTR_VALUE] No '{}' in attributes: {}:{} - ({})".format(
                key, self.chr, self.pos, attr), 'WRN')
            return None

    def _choose_gene_info(self):
        # Prefer (first) 'protein_coding' gene, or first one
        for x in self.gene_hits:
            gene_type = x['type']
            if gene_type == 'protein_coding':
                return x
        return self.gene_hits[0]

    def _choose_default_transcript(self):
        # Prefer (longest) transcript in `transcript_map`, or longest transcript
        trans_hits_sorted = sorted(self.trans_hits, key=lambda x: x['size'], reverse=True)
        if len(trans_hits_sorted):
            for x in trans_hits_sorted:
                transcript_id = x['id']
                transcript_id_main = re.sub(r'\.\d+$', '', transcript_id)
                if transcript_id_main in self.transcript_map:
                    return x['id']
            self.log("[CHOOSE_TRANSCRIPT] Use longest transcript as default: {} for {}".format(
                trans_hits_sorted[0]['id'], self.gene_info['name']), 'WRN')
            return trans_hits_sorted[0]['id']

    def _choose_rank_info(self):
        # Prefer default transcript, or first one
        for x in self.trans_hits:
            if x['id'] == self.default_transcript:
                return x['rank_info']
        return self.trans_hits[0]['rank_info']

    def get_gene_name(self):
        if 'name' in self.gene_info:
            return self.gene_info['name']
        else:
            return None

    def get_rank_info(self):
        if len(self.rank_info):
            return self.rank_info
        return None
