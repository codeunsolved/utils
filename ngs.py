#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : ngs
# AUTHOR  : codeunsolved@gmail.com
# CREATED : March 10 2018
# VERSION : v0.0.17
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
# 1. [BugFix] fix tag content extraction regex in :VcfParser: `parse_meta()`;
# 2. [BugFix] fix bug in :VcfParser: `extract_format()`;
# UPDATE  : [v0.0.8] April 8 2019
# 1. enhance :VcfParser:;
# UPDATE  : [v0.0.9] May 10 2019
# 1. enhance :AnnCoordinates: by using cache for entire gene region;
# UPDATE  : [v0.0.10] June 17 2019
# 1. [BugFix] fix cross modification problem of 'rank_info' dict in :AnnCoordinates:;
# UPDATE  : [v0.0.11] June 19 2019
# 1. :AnnCoordinates: optimize gene, transcript, rank selection policy and entire work flow;
# 2. :AnnCoordinates: [BugFix] fix UTR rank annotation which does not handle intron;
# 3. :AnnCoordinates: optimize logger;
# UPDATE  : [v0.0.12] June 20 2019
# 1. :AnnCoordinates: add option to filter non-protein-coding gene;
# UPDATE  : [v0.0.13] September 7 2019
# 1. handle intergenic IGH as gene not intergenic;
# UPDATE  : [v0.0.14] November 14 2019
# 1. optimize attr index in :AnnCoordinates:;
# UPDATE  : [v0.0.15] February 26 2020
# 1. :AnnCoordinates: set rank to 'intragenic' if the coordinate out of transcript;
# UPDATE  : [v0.0.16] March 19 2020
# 1. :AnnCoordinates: adjust the framework to prepare to support a new source (RefSeq);
# UPDATE  : [v0.0.17] March 21 2020
# 1. :AnnCoordinates: add new source: RefSeq;

import os
import re
import logging
from copy import deepcopy as copy
from collections import defaultdict

import numpy as np
import pandas as pd

from .base import colour
from .base import color_term
from .connector import MysqlConnector

__VERSION__ = 'v0.0.17'


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

    def __init__(self,
                 mysql_config=None, mysql_con=None,
                 gene_map={}, transcript_map={},
                 refer='hg19',
                 ann_source='GENCODE',
                 ann_tab=None,
                 ann_gene_type='protein_coding',
                 ann_refseq_files=None,  # GenePred table format
                 use_cache=True,
                 logger_name=None, logger_level=logging.INFO):

        # Configs
        self.mysql_config = mysql_config
        self.mysql_con = mysql_con

        # Options
        self.gene_map = gene_map
        self.transcript_map = transcript_map
        self.refer = refer
        self.ann_source = ann_source
        self.ann_tab = ann_tab
        self.ann_gene_type = ann_gene_type
        self.ann_refseq_files = ann_refseq_files
        self.ann_refseq_df = None
        self.use_cache = use_cache
        self.logger_name = logger_name or self.__class__.__name__
        self.logger_level = logger_level

        # Coordinates
        self.chr_raw = None  # original chr
        self.chr_index = None
        self.chr = None
        self.pos = None

        # Annotations
        # Gene hits contains {key: val}:
        # 'id': Ensembl gene ID
        # 'name': GENCODE gene name
        # 'type': GENCODE gene type
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
        # 'rank': 'Exon{N}', 'Intron{N}', "5'UTR", "3'UTR", 'intragenic'
        # 'exon_num': {N}
        # 'exon_total': {N}
        # 'cds_min': {N}
        # 'cds_max': {N}
        # 'exon_start': exon start, it will be next exon start if rank is 'Intron{N}'
        # 'exon_end': exon end, it will be last exon end if rank is 'Intron{N}'
        self.rank_info = {}
        self.default_transcript = None
        self.coord_cache = defaultdict(dict)
        self.region_cache = []

        self.protein_coding_types = [
            'protein_coding',
            r'IG_\w_gene',
            r'TR_\w_gene',
            'polymorphic_pseudogene',
        ]

        self.ann_tabs = {
            'hg19': {
                'GENCODE': 'GENCODE_Human_v28lift37_annotation_GFF3',
            }
        }

        # Init logger
        self.logger = self.init_logger()

        # Check params validation
        self.check_params()

        # Set annotation table (only used by GENCODE)
        if self.ann_source == 'GENCODE':
            self.set_ann_tab()

        # Init MySQL connection (only used by GENCODE)
        if self.ann_source == 'GENCODE':
            if self.mysql_con is None:
                self.connect_db()
            assert self.mysql_con is not None

        # Set annotation RefSeq DataFrame (only used by RefSeq)
        if self.ann_source == 'RefSeq':
            self.ann_refseq_df = self.gen_refseq_df()

    def init_logger(self):
        logger = logging.getLogger(self.logger_name)
        if logger.parent.name == 'root':
            logging.basicConfig(
                format="[{}] %(message)s".format(
                    self.__class__.__name__))
        logger.setLevel(self.logger_level)
        return logger

    def check_params(self):
        valid_refers = ['hg19', 'hg38']
        assert self.refer in valid_refers, \
            "Unexpected reference genome: {}! Valid choices: {}".format(
                self.refer, valid_refers)

        valid_ann_sources = ['GENCODE', 'RefSeq']
        assert self.ann_source in valid_ann_sources, \
            "Unexpected annotation source: {}! Valid choices: {}".format(
                self.ann_source, valid_ann_sources)

        valid_ann_gene_type = ['protein_coding', 'all']
        assert self.ann_gene_type in valid_ann_gene_type, \
            "Unexpected annotation gene type: {}! Valid choices: {}".format(
                self.ann_gene_type, valid_ann_gene_type)

    def set_ann_tab(self):
        if self.ann_tab is None:
            self.ann_tab = self.ann_tabs[self.refer][self.ann_source]

    def format_chr(self):
        self.chr_index = re.sub('^chr', '', self.chr_raw.lower()).upper()
        self.chr = 'chr' + self.chr_index

        if not re.match(r'[12]?\d|[XY]$', self.chr_index):
            self.logger.warning(
                color_term("[CHECK_CHR] Unrecognized Chr: {}".format(self.chr_raw), 'WRN'))

    def connect_db(self):
        if self.mysql_config:
            self.mysql_con = MysqlConnector(
                self.mysql_config,
                verbose=self.logger_level < logging.INFO)

    def gen_refseq_df(self):
        def read_csv(path_file):
            return pd.read_csv(path_file, header=0, sep='\t')

        assert isinstance(self.ann_refseq_files, (str, list, set)), \
            "Invalid RefSeq files type: {}, str or list or set expected!".format(
                type(self.ann_refseq_files))

        if isinstance(self.ann_refseq_files, str):
            assert os.path.isfile(self.ann_refseq_files), \
                "Invalid RefSeq file: {}".format(self.ann_refseq_files)
            return read_csv(self.ann_refseq_files)

        elif isinstance(self.ann_refseq_files, (list, set)):
            ann_refseq_df = None
            for ann_refseq_file in self.ann_refseq_files:
                assert os.path.isfile(ann_refseq_file), \
                    "Invalid RefSeq file: {}".format(ann_refseq_file)

                ann_refseq_file = os.path.abspath(ann_refseq_file)
                if ann_refseq_df is None:
                    ann_refseq_df = read_csv(ann_refseq_file)
                else:
                    ann_refseq_df_new = read_csv(ann_refseq_file)
                    assert list(ann_refseq_df.columns) == list(ann_refseq_df_new.columns), \
                        "{} has a header that is inconsistent with the others".format(
                            os.path.basename(ann_refseq_file))
                    ann_refseq_df = pd.concat([ann_refseq_df, ann_refseq_df_new])
            return ann_refseq_df

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

    def query(self, *args):
        def handle_bp(args):
            chr_raw = None
            pos = 0

            if len(args) == 1:
                bp = args[0]
                if ':' in bp:
                    chr_raw, pos = bp.split(':')
                else:
                    self.logger.error(colour(
                        "[HANDLE_BP] Invalid Breakpoint: {}".format(
                            bp), 'ERR'))
            elif len(args) == 2:
                chr_raw = args[0]
                pos = args[1]

            self.chr_raw = chr_raw
            if isinstance(pos, str):
                pos = re.sub(',', '', pos)
            self.pos = int(pos)
            self.format_chr()

        def retrieve_region_cache(chr_, pos):
            for region in self.region_cache:
                if region['chr'] == chr_ and \
                   region['start'] <= pos and pos <= region['end']:
                    return region
            return None

        def cache_coord(coord):
            self.coord_cache[coord]['gene_info'] = copy(self.gene_info)
            self.coord_cache[coord]['rank_info'] = copy(self.rank_info)

        def cache_region():
            region = retrieve_region_cache(self.chr, self.pos)
            if region is None:
                region = {}
                region['chr'] = self.chr
                region['start'] = self.gene_info['start']
                region['end'] = self.gene_info['end']
                region['trans_hits'] = copy(self.trans_hits)
                region['gene_info'] = copy(self.gene_info)
                self.region_cache.append(region)

        self.init_ann()

        handle_bp(args)

        if None in (self.chr, self.pos):
            return

        coord = "{}:{}".format(self.chr, self.pos)
        if self.use_cache and coord in self.coord_cache:
            self.logger.info("[QUERY] Use coord cache for {}".format(coord))
            self.gene_info = copy(self.coord_cache[coord]['gene_info'])
            self.rank_info = copy(self.coord_cache[coord]['rank_info'])
        else:
            region = retrieve_region_cache(self.chr, self.pos)

            if self.use_cache and region is not None:
                self.logger.info(
                    "[QUERY] Use region cache for {} hit by {}".format(
                        coord, region['gene_info']['name']))

                self.trans_hits = copy(region['trans_hits'])

                self.gene_info = copy(region['gene_info'])
                self.default_transcript = self._choose_default_transcript()
                self.set_rank_info()
                self.rank_info = self._choose_rank_info()
            else:
                self.query_gene()
                self.query_trans()

                self.gene_info = self._choose_gene_info()
                self.default_transcript = self._choose_default_transcript()
                self.set_rank_info()
                self.rank_info = self._choose_rank_info()

            # Cache
            cache_coord(coord)
            cache_region()

    def query_gene(self):
        try:
            if self.ann_source == 'GENCODE':
                self._query_gene_via_gencode()
            elif self.ann_source == 'RefSeq':
                self._query_gene_and_trans_via_refseq()
        except Exception as e:
            self.logger.error(colour(
                "[QUERY_GENE] Error! {}:{}, {}".format(
                    self.chr, self.pos, e), 'ERR'))
        else:
            if self.ann_source == 'GENCODE':
                self.logger.info(
                    "• Query gene {}! {} hits genes: {}".format(
                        colour("OK", 'green'),
                        colour(len(self.gene_hits), 'green'),
                        [x['name'] for x in self.gene_hits]))
            elif self.ann_source == 'RefSeq':
                self.logger.info(
                    "• Query gene {}! {} hits genes: {}; {} hits transcripts: {}".format(
                        colour("OK", 'green'),
                        colour(len(self.gene_hits), 'green'),
                        [x['name'] for x in self.gene_hits],
                        colour(len(self.trans_hits), 'green'),
                        [x['id'] for x in self.trans_hits]))

    def query_trans(self):
        # No need for RefSeq
        if self.ann_source == 'RefSeq':
            return

        try:
            self._query_trans_via_gencode()
        except Exception as e:
            self.logger.error(colour(
                "[QUERY_TRANS] Error! {}:{}, {}".format(
                    self.chr, self.pos, e), 'ERR'))
        else:
            self.logger.info(
                "• Query transcripts {}! {} hits transcripts: {}".format(
                    colour("OK", 'green'),
                    colour(len(self.trans_hits), 'green'),
                    [x['id'] for x in self.trans_hits]))

    def set_rank_info(self):
        try:
            if self.ann_source == 'GENCODE':
                self._set_rank_info_via_gencode()
            elif self.ann_source == 'RefSeq':
                self._set_rank_info_via_refseq()
        except Exception as e:
            self.logger.error(colour(
                "[SET_RANK_INFO] Error! {}:{}, {}".format(
                    self.chr, self.pos, e), 'ERR'))

    def _query_gene_via_gencode(self):
        def gen_gene_type_regex():
            if self.ann_gene_type == 'protein_coding':
                pc_regex = '|'.join(self.protein_coding_types)
                pc_regex = re.sub(r'\\w', '[[:alpha:]]', pc_regex)
                gene_type_regex = " AND Attribute REGEXP 'gene_type=({});'".format(
                    pc_regex)
            else:
                gene_type_regex = ''
            return gene_type_regex

        def gen_r(gene_type_regex):
            sql = ("SELECT * FROM `{tab}` "
                   "WHERE SeqID = %s AND Feature = 'gene' AND "
                   "Start <= %s AND End >= %s{gene_type_regex}")
            sql = sql.format(
                tab=self.ann_tab,
                gene_type_regex=gene_type_regex)

            c = self.mysql_con.query(sql, [self.chr, self.pos, self.pos])
            r = c.fetchall()
            return r

        def gen_gene_hits(entries):
            gene_hits = []
            for entry in entries:
                attr = entry[9]

                gene_id = self._get_attr_value(attr, 'ID')
                gene_name = self._get_attr_value(attr, 'gene_name')
                gene_type = self._get_attr_value(attr, 'gene_type')

                gene_hits.append({
                    'id': gene_id,
                    'name': gene_name,
                    'type': gene_type,
                    'source': entry[2],
                    'chr': entry[1],
                    'start': int(entry[4]),
                    'end': int(entry[5]),
                    'strand': entry[7],
                    'transcripts': [],
                })
            return gene_hits

        def query_gene_intergenic(gene_type_regex):
            sql_l = ("SELECT * FROM `{tab}` "
                     "WHERE SeqID = %s AND Feature = 'gene' AND "
                     "End < %s{gene_type_regex} ORDER BY End DESC LIMIT 1")
            sql_l = sql_l.format(
                tab=self.ann_tab,
                gene_type_regex=gene_type_regex)

            r_l = self.mysql_con.query(sql_l, [self.chr, self.pos]).fetchone()
            sql_r = ("SELECT * FROM `{tab}` "
                     "WHERE SeqID = %s AND Feature = 'gene' AND "
                     "Start > %s{gene_type_regex} ORDER BY Start LIMIT 1")
            sql_r = sql_r.format(
                tab=self.ann_tab,
                gene_type_regex=gene_type_regex)
            r_r = self.mysql_con.query(sql_r, [self.chr, self.pos]).fetchone()

            gene_name_l = '-' if r_l is None else self._get_attr_value(r_l[9], 'gene_name')
            gene_name_r = '-' if r_r is None else self._get_attr_value(r_r[9], 'gene_name')
            start = 0 if r_l is None else r_l[5]  # Left gene's end
            end = 0 if r_r is None else r_r[4]    # Right gene's start
            name = "intergenic({},{})".format(gene_name_l, gene_name_r)

            gene_hits = [{
                'id': None,
                'name': name,
                'type': 'intergenic',
                'source': None,
                'chr': self.chr,
                'start': int(start),
                'end': int(end),
                'strand': None,
                'transcripts': [],
            }]
            return gene_hits

        gene_type_regex = gen_gene_type_regex()
        r = gen_r(gene_type_regex)

        if len(r):
            self.gene_hits = gen_gene_hits(r)
        else:
            self.logger.warning(colour(
                "[QUERY_GENE] No gene found for {}:{}, "
                "try to annotate as intergenic region".format(
                    self.chr, self.pos), 'WRN'))
            self.gene_hits = query_gene_intergenic(gene_type_regex)
            assert len(self.gene_hits) == 1
            # Special case: intergenic IGH
            self._handle_intergenic_igh(self.gene_hits[0])

    def _query_trans_via_gencode(self):
        def gen_r():
            sql = ("SELECT * FROM `{tab}` "
                   "WHERE SeqID = %s AND Feature IN "
                   "('transcript', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR') "
                   "AND Attribute LIKE %s ORDER BY Start")
            sql = sql.format(tab=self.ann_tab)

            r = []
            for gene in self.gene_hits:
                if gene['id'] is not None:
                    c = self.mysql_con.query(sql, [self.chr, r"%gene_id={}%".format(gene['id'])])
                    r += c.fetchall()
            return r

        def append_trans_id(gene_id, trans_id):
            for gene in self.gene_hits:
                if gene['id'] == gene_id:
                    gene['transcripts'].append(trans_id)

        def gen_trans_hits(entries):
            trans_hits = {}
            for entry in entries:
                feature = entry[3]
                attr = entry[9]
                if feature == 'transcript':
                    gene_id = self._get_attr_value(attr, 'gene_id')
                    if gene_id is not None:
                        gene_id = re.sub(r'_\d+$', '', gene_id)
                        transcript_id = self._get_attr_value(attr, 'ID')
                        if transcript_id is not None:
                            trans_hits[transcript_id] = {
                                'id': transcript_id,
                                'gene_id': gene_id,
                                'start': int(entry[4]),
                                'end': int(entry[5]),
                                'size': int(entry[5] - entry[4]),
                                'strand': entry[7],
                                'components': [],
                                'rank_info': {
                                    'trans_id': transcript_id,
                                    'trans_start': int(entry[4]),
                                    'trans_end': int(entry[5]),
                                    'trans_strand': entry[7],
                                    'rank': None,
                                    'exon_num': 0,
                                    'exon_total': 0,
                                    'cds_min': 0,
                                    'cds_max': 0,
                                    'exon_start': 0,
                                    'exon_end': 0,
                                },
                            }
                            append_trans_id(gene_id, transcript_id)

            for entry in entries:
                attr = entry[9]
                transcript_id = self._get_attr_value(attr, 'transcript_id')
                transcript_id = re.sub(r'_\d+$', '', transcript_id)

                trans_hits[transcript_id]['components'].append({
                    'feature': entry[3],
                    'start': int(entry[4]),
                    'end': int(entry[5]),
                    'strand': entry[7],
                    'attr': entry[9],
                })

            trans_hits = sorted(trans_hits.values(), key=lambda x: x['size'], reverse=True)
            return trans_hits

        r = gen_r()

        if len(r):
            self.trans_hits = gen_trans_hits(r)
        else:
            if len([x for x in self.gene_hits if x['id'] is not None]):
                self.logger.error(colour(
                    "[QUERY_TRANS] No transcripts found for [{}]".format(
                        ', '.join([x['name'] for x in self.gene_hits])), 'ERR'))

    def _set_rank_info_via_gencode(self):
        def set_rank():
            def gen_rank(feature, exon_num, exon=True):
                if feature == 'CDS':
                    if exon:
                        return "Exon{}".format(exon_num)
                    else:
                        return "Intron{}".format(exon_num)
                elif feature == 'five_prime_UTR':
                    if exon:
                        return "5'UTR"
                    else:
                        return "Intron{}".format(exon_num)
                elif feature == 'three_prime_UTR':
                    if exon:
                        return "3'UTR"
                    else:
                        return "Intron{}".format(exon_num)
                else:
                    return feature

            rank = None
            exon_num = 0

            for trans in self.trans_hits:
                transcript_id = trans['id']
                rank_info = trans['rank_info']

                loc_flag = -1
                for entry in sorted(trans['components'], key=lambda x: (x['start'], x['end'])):
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
                    if trans['strand'] == '+':
                        strand_flag = 1
                    elif trans['strand'] == '-':
                        strand_flag = -1
                    else:
                        strand_flag = 0

                    if strand_flag * loc_flag == 1:
                        # rank = "3'UTR"
                        rank = "intragenic"
                        exon_num = -3
                    elif strand_flag * loc_flag == -1:
                        # rank = "5'UTR"
                        rank = "intragenic"
                        exon_num = -5
                    else:
                        rank = ''

                    if transcript_id == self.default_transcript:
                        self.logger.warning(colour(
                            "[SET_RANK] {}:{} seems {} than "
                            "the range of transcript: {}({}) {}-{}, "
                            "use {} instead".format(
                                self.chr, self.pos,
                                'LESS' if loc_flag == -1 else 'LARGE',
                                transcript_id, trans['strand'],
                                trans['start'], trans['end'],
                                rank), 'WRN'))
                rank_info['rank'] = rank
                rank_info['exon_num'] = exon_num

        def set_exon_and_cds():
            for trans in self.trans_hits:
                rank_info = trans['rank_info']

                max_exon_num = 0
                min_cds_num = np.inf
                max_cds_num = 0
                for entry in trans['components']:
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
                            elif exon_num == rank_info['exon_num'] + 1:
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

        set_rank()
        set_exon_and_cds()

    def _query_gene_and_trans_via_refseq(self):
        def get_breakpoint_hits(chrom, pos):
            df_bp_hits = self.ann_refseq_df[
                (self.ann_refseq_df['chrom'] == chrom) &
                (self.ann_refseq_df['txStart'] <= pos) &
                (self.ann_refseq_df['txEnd'] >= pos)]
            return df_bp_hits

        def get_gene_hits(chrom, gene_name):
            df_gene_hits = self.ann_refseq_df[
                (self.ann_refseq_df['chrom'] == chrom) &
                (self.ann_refseq_df['name2'] == gene_name)]
            return df_gene_hits

        def gen_gene_and_trans_hits(df_hits):
            gene_hits_dict = {}
            for idx, row in df_hits.iterrows():
                gene_name = row['name2']
                gene_type = 'protein_coding'  # set to 'protein_coding' as default
                chrom = row['chrom']
                strand = row['strand']

                if gene_name not in gene_hits_dict:
                    gene_hits_dict[gene_name] = {
                        'id': gene_name,
                        'name': gene_name,
                        'type': gene_type,
                        'source': 'RefSeq',
                        'chr': chrom,
                        'start': None,
                        'end': None,
                        'strand': strand,
                        'transcripts': [],
                    }

            gene_hits = []
            trans_hits = {}
            for gene_name, item in sorted(gene_hits_dict.items(), key=lambda x: x[0]):
                df_gene_hits = get_gene_hits(item['chr'], gene_name)
                item['start'] = min(df_gene_hits['txStart'])
                item['end'] = max(df_gene_hits['txEnd'])
                item['transcripts'] = sorted(list(set(df_gene_hits['name'])))
                gene_hits.append(item)

                for idx, row in df_gene_hits.iterrows():
                    transcript_id = row['name']
                    gene_name = row['name2']
                    tx_start = row['txStart']
                    tx_end = row['txEnd']
                    strand = row['strand']

                    trans_hits[transcript_id] = {
                        'id': transcript_id,
                        'gene_id': gene_name,
                        'start': tx_start,
                        'end': tx_end,
                        'size': tx_end - tx_start,
                        'strand': strand,
                        'exonStarts': row['exonStarts'],
                        'exonEnds': row['exonEnds'],
                        'rank_info': {
                            'trans_id': transcript_id,
                            'trans_start': tx_start,
                            'trans_end': tx_end,
                            'trans_strand': strand,
                            'rank': None,
                            'exon_num': 0,
                            'exon_total': row['exonCount'],
                            'cds_start': row['cdsStart'],
                            'cds_end': row['cdsEnd'],
                            'cds_min': 0,
                            'cds_max': 0,
                            'exon_start': 0,
                            'exon_end': 0,
                        },
                    }

            trans_hits = sorted(trans_hits.values(), key=lambda x: x['size'], reverse=True)
            return gene_hits, trans_hits

        def get_intergenic_hits(chrom, pos):
            def get_gene_intergenic_name(hits):
                if len(hits):
                    hit = hits.iloc[0]
                    gene_name = hit['name2']
                    return gene_name
                else:
                    return '-'

            def get_gene_intergenic_range(hits, key):
                if len(hits):
                    hit = hits.iloc[0]
                    return hit[key]
                else:
                    return 0

            hits_left = self.ann_refseq_df[
                (self.ann_refseq_df['chrom'] == chrom) &
                (self.ann_refseq_df['txEnd'] <= pos)].sort_values('txEnd', ascending=False)
            hits_right = self.ann_refseq_df[
                (self.ann_refseq_df['chrom'] == chrom) &
                (self.ann_refseq_df['txStart'] >= pos)].sort_values('txStart', ascending=True)

            gene_left = get_gene_intergenic_name(hits_left)
            gene_right = get_gene_intergenic_name(hits_right)
            gene_name = "intergenic({},{})".format(gene_left, gene_right)

            start = get_gene_intergenic_range(hits_left, 'txEnd')  # left gene's end
            end = get_gene_intergenic_range(hits_right, 'txStart')  # right gene's start

            gene_hits = [{
                'id': None,
                'name': gene_name,
                'type': 'intergenic',
                'source': None,
                'chr': self.chr,
                'start': int(start),
                'end': int(end),
                'strand': None,
                'transcripts': [],
            }]
            return gene_hits

        assert self.ann_refseq_df is not None

        df_bp_hits = get_breakpoint_hits(self.chr, self.pos)
        if len(df_bp_hits):
            self.gene_hits, self.trans_hits = gen_gene_and_trans_hits(df_bp_hits)
        else:
            self.logger.warning(colour(
                "[QUERY_GENE] No gene found for {}:{}, "
                "try to annotate as intergenic region".format(
                    self.chr, self.pos), 'WRN'))
            self.gene_hits = get_intergenic_hits(self.chr, self.pos)
            assert len(self.gene_hits) == 1
            # Special case: intergenic IGH
            self._handle_intergenic_igh(self.gene_hits[0])

    def _set_rank_info_via_refseq(self):
        def format_exons(trans):
            exonStarts = trans['exonStarts'].rstrip(',').split(',')
            exonEnds = trans['exonEnds'].rstrip(',').split(',')
            exons = zip(exonStarts, exonEnds)
            exons = [(int(start), int(end)) for (start, end) in exons]

            return exons

        def get_exon_num(i, exons, strand):
            if strand == '+':
                return i + 1
            elif strand == '-':
                return len(exons) - i

        def set_utr(rank_info, rank, pos):
            if re.match('Exon', rank):
                if pos < rank_info['cds_start']:
                    if rank_info['trans_strand'] == '+':
                        return "5'UTR"
                    elif rank_info['trans_strand'] == '-':
                        return "3'UTR"
                elif pos > rank_info['cds_end']:
                    if rank_info['trans_strand'] == '+':
                        return "3'UTR"
                    elif rank_info['trans_strand'] == '-':
                        return "5'UTR"
            return rank

        for trans in self.trans_hits:
            exons = format_exons(trans)
            transcript_id = trans['id']
            strand = trans['strand']
            rank_info = trans['rank_info']

            exon_num = 0
            exon_start = 0
            exon_end = 0

            loc_flag = -1
            for i, (start, end) in enumerate(exons):
                exon_num_i = get_exon_num(i, exons, strand)
                if start + 1 <= self.pos and self.pos <= end:
                    exon_num = exon_num_i
                    exon_start = start
                    exon_end = end
                    loc_flag = 0
                    rank = "Exon{}".format(exon_num)
                    break
                else:
                    if self.pos < start:
                        if loc_flag == 1:
                            if strand == '+':
                                exon_num = exon_num_i - 1
                            elif strand == '-':
                                exon_num = exon_num_i
                            exon_start = start
                            loc_flag = 0
                            rank = "Intron{}".format(exon_num)
                            break
                        else:
                            loc_flag = -1
                    elif end < self.pos:
                        loc_flag = 1
                        exon_end = end

            if loc_flag != 0:
                if trans['strand'] == '+':
                    strand_flag = 1
                elif trans['strand'] == '-':
                    strand_flag = -1
                else:
                    strand_flag = 0

                if strand_flag * loc_flag == 1:
                    # rank = "3'UTR"
                    rank = "intragenic"
                    exon_num = -3
                elif strand_flag * loc_flag == -1:
                    # rank = "5'UTR"
                    rank = "intragenic"
                    exon_num = -5
                else:
                    rank = ''

                if transcript_id == self.default_transcript:
                    self.logger.warning(colour(
                        "[SET_RANK] {}:{} seems {} than "
                        "the range of transcript: {}({}) {}-{}, "
                        "use {} instead".format(
                            self.chr, self.pos,
                            'LESS' if loc_flag == -1 else 'LARGE',
                            transcript_id, trans['strand'],
                            trans['start'], trans['end'],
                            rank), 'WRN'))
            else:
                rank = set_utr(rank_info, rank, self.pos)

            # Set cds_min & cds_max
            for i, (start, end) in enumerate(exons):
                exon_num_i = get_exon_num(i, exons, strand)
                if start <= rank_info['cds_start'] and \
                   rank_info['cds_start'] <= end:
                    if strand == '+':
                        rank_info['cds_min'] = exon_num_i
                    elif strand == '-':
                        rank_info['cds_max'] = exon_num_i
                if start <= rank_info['cds_end'] and \
                   rank_info['cds_end'] <= end:
                    if strand == '+':
                        rank_info['cds_max'] = exon_num_i
                    elif strand == '-':
                        rank_info['cds_min'] = exon_num_i

            rank_info['rank'] = rank
            rank_info['exon_num'] = exon_num
            rank_info['exon_start'] = exon_start
            rank_info['exon_end'] = exon_end

    def _handle_intergenic_igh(self, gene_hit):
        def extract_igh_from_intergenic(gene):
            if 'intergenic' in gene:
                intergenic_gene_regex = \
                    re.match(r'intergenic\(([^,]+), ?([^,]+)\)', gene)
                if intergenic_gene_regex:
                    gene_a, gene_b = intergenic_gene_regex.groups()
                    if re.match('IGH', gene_a) and re.match('IGH', gene_b):
                        if re.match('IGH[VDJ]', gene_a):
                            return gene_a, True
                        elif re.match('IGH[VDJ]', gene_b):
                            return gene_b, True
                        else:
                            return gene_a, True
            return gene, False

        gene = gene_hit['name']
        gene_igh, is_igh = extract_igh_from_intergenic(gene)
        if is_igh:
            self.logger.warning(
                "[HANDLE_INTERGENIC_IGH] intergenic IGH: {} -> {}".format(
                    gene, gene_igh))
            gene_hit['name'] = gene_igh
            gene_hit['type'] = 'intergenic_igh'
            # Suppose all IGH* gene are minus strand(-)
            gene_hit['strand'] = '-'

    def _get_attr_value(self, attr, key):
        s = re.search("{}=([^;]+);".format(key), attr)
        if s is not None:
            return s.group(1)
        else:
            self.logger.warning(colour(
                "[GET_ATTR_VALUE] No '{}' in attributes: {}:{} - ({})".format(
                    key, self.chr, self.pos, attr), 'WRN'))
            return None

    def _check_if_protein_coding(self, gene_type):
        protein_coding_regex = "{}$".format(
            '|'.join(self.protein_coding_types))

        if re.match(
                protein_coding_regex,
                gene_type):
            return True
        else:
            return False

    def _choose_gene_info(self):
        """
        Selection policy:
           1. Prefer gene which transcript in self.transcript_map
           2. Prefer largest gene with protein coding
              (Protein coding: IGC gene, IGD gene, IG gene,
               IGJ gene, IGLV gene, IGM gene, IGV gene, IGZ gene,
               nonsense mediated decay, nontranslating CDS,
               non stop decay, polymorphic pseudogene,
               TRC gene, TRD gene, TRJ gene.)

        Refer: What do the different biotypes in Ensembl mean?
              (https://asia.ensembl.org/Help/Faq?id=468)
        """

        def sort_gene_hits(x):
            gene_len = x['end'] - x['start']
            if self._check_if_protein_coding(x['type']):
                return 10**8 + gene_len
            else:
                return gene_len

        if len(self.gene_hits) == 0:
            return {}
        elif len(self.gene_hits) == 1:
            return self.gene_hits[0]

        gene_info = None

        # Policy 1
        for x in self.gene_hits:
            for trans_id in x['transcripts']:
                trans_id_main = re.sub(r'\.\d+$', '', trans_id)
                if trans_id_main in self.transcript_map:
                    gene_info = copy(x)
                    self.logger.warning(colour(
                        "[CHOOSE_GENE] Use {}({}) due to {} in 'transcript_map'".format(
                            gene_info['id'], gene_info['name'], trans_id_main), 'WRN'))
                    break
            if gene_info is not None:
                break

        # Policy 2
        if gene_info is None:
            gene_hits_sorted = sorted(self.gene_hits, key=sort_gene_hits, reverse=True)
            gene_info = copy(gene_hits_sorted[0])

            self.logger.warning(colour(
                "[CHOOSE_GENE] Use largest (protein coding) gene "
                "as default: {}({}) {}".format(
                    gene_info['id'], gene_info['name'], gene_info['type']), 'WRN'))

        if gene_info is None:
            gene_info = {}

        # Map gene name
        if gene_info.get('name', None) in self.gene_map:
            self.logger.warning(colour(
                "[CHOOSE_GENE] Gene name changed via 'gene_map': {} -> {}".format(
                    gene_info['name'], self.gene_map[gene_info['name']]), 'WRN'))
            gene_info['name'] = self.gene_map[gene_info['name']]

        if len(self.gene_hits) > 1:
            self.logger.warning(colour(
                "[CHOOSE_GENE] {} entries [{}] found for {}:{}, choose: {}".format(
                    len(self.gene_hits), ', '.join([x['name'] for x in self.gene_hits]),
                    self.chr, self.pos, gene_info['name']), 'WRN'))

        return gene_info

    def _choose_default_transcript(self):
        """
        Selection policy:
           0. Prefer transcripts which in self.gene_info
           1. Prefer transcript which in self.transcript_map
           2. Prefer longest transcript
        """

        default_transcript = None

        # Policy 0
        trans_hits = [x for x in self.trans_hits if x['id'] in self.gene_info['transcripts']]

        # Policy 1
        trans_hits_sorted = sorted(trans_hits, key=lambda x: x['size'], reverse=True)
        for x in trans_hits_sorted:
            transcript_id = x['id']
            transcript_id_main = re.sub(r'\.\d+$', '', transcript_id)
            if transcript_id_main in self.transcript_map:
                default_transcript = x['id']
                break

        # Policy 2
        if default_transcript is None:
            if len(trans_hits_sorted):
                self.logger.warning(colour(
                    "[CHOOSE_TRANSCRIPT] Use longest transcript as default: {} for {}".format(
                        trans_hits_sorted[0]['id'], self.gene_info['name']), 'WRN'))
                default_transcript = trans_hits_sorted[0]['id']

        return default_transcript

    def _choose_rank_info(self):
        """
        Selection policy:
           0. Prefer transcripts which in self.gene_info
           1. Prefer default transcript
           2. Prefer first transcript
        """

        rank_info = None

        # Policy 0
        trans_hits = [x for x in self.trans_hits if x['id'] in self.gene_info['transcripts']]

        # Policy 1
        for x in trans_hits:
            if x['id'] == self.default_transcript:
                rank_info = copy(x['rank_info'])
                break

        # Policy 2
        if rank_info is None:
            if len(trans_hits):
                rank_info = copy(trans_hits[0]['rank_info'])

        if rank_info is None:
            rank_info = {}

        # Special case: intergenic IGH
        if self.gene_info.get('type', '') == 'intergenic_igh':
            rank_info['trans_strand'] = '-'

        return rank_info

    def get_gene_name(self):
        assert isinstance(self.gene_info, dict)
        return self.gene_info.get('name', None)

    def get_rank_info(self):
        assert isinstance(self.rank_info, dict)
        return self.rank_info
