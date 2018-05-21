#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : ngs
# AUTHOR  : codeunsolved@gmail.com
# CREATED : March 10 2018
# VERSION : v0.0.1a1
# UPDATE  : [v0.0.1a1] May 16 2017
# 1. add `sequence_complement()`;
# 2. add :AnnCoordinates:;

import re
from collections import defaultdict

from .base import color_term
from .connector import MysqlConnector


def sequence_complement(sequence, reverse=True):
    base_comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq_comp = ''.join([base_comp[x.upper()] for x in sequence])
    if reverse:
        seq_comp = seq_comp[::-1]
    return seq_comp


class AnnCoordinates(object):

    def __init__(self, config=None, mysql_con=None,
                 ann_gene=True, ann_region=False,
                 refer='hg19', refer_source='GENCODE',
                 cached=True, verbose=False):
        # Configs
        self.config = config
        self.m_con = mysql_con

        # Options
        self.refer = refer
        self.refer_source = refer_source
        self.verbose = verbose
        self.cached = cached
        self.ann_gene = ann_gene
        self.ann_region = ann_region

        # Coordindates
        self.chr = None
        self.pos = None
        self.chr_index = None
        self.refer_tab = None

        # Annotations
        self.gene = None
        # Region values:
        # 'exon': Exon
        # 'intron': Intron
        # '5'UTR': 5'UTR
        # '3'UTR': 3'UTR
        self.region = None
        # Gene hits:
        # None: Init status
        # n: n hits in DB, expect 1, 0 means maybe intergenic
        self.gene_hits = None
        self.cache = defaultdict(dict)

        # Check params validation
        self.check_refer()

        self.connect_db()

    def check_chr(self):
        self.chr_index = re.sub('^chr', '', self.chr_.lower()).upper()
        self.chr = 'chr' + self.chr_index

        if not re.match('[12]?\d|[XY]$', self.chr_index):
            color_term("[CHECK_CHR] Unrecognized Chr: {}".format(self.chr_), 'WRN')

    def check_refer(self):
        assert self.refer in ['hg19', 'hg38']
        assert self.refer_source in ['GENCODE']

        self.refer_tabs = {
            'hg19': {
                'GENCODE': 'GENCODE_Human_v28lift37_annotation_GFF3',
            }
        }

        self.refer_tab = self.refer_tabs[self.refer][self.refer_source]

    def connect_db(self):
        if self.config:
            self.m_con = MysqlConnector(self.config, self.verbose)

        assert self.m_con is not None

    def query(self, chr_, pos):
        self.gene = None
        self.region = None

        self.chr_ = chr_
        self.pos = pos
        self.check_chr()

        key = "{}:{}".format(self.chr, self.pos)
        if key in self.cache:
            if self.verbose:
                color_term("[QUERY] Use cache for {}".format(key), 'WRN')
            self.gene = self.cache[key]['gene']
            self.region = self.cache[key]['region']
        else:
            if self.ann_gene:
                try:
                    self.query_gene()
                except Exception as e:
                    print(e)

            if self.gene or self.region:
                self.cache[key]['gene'] = self.gene
                self.cache[key]['region'] = self.region

    def query_gene(self):
        def get_gene_name(attr):
            s = re.search('gene_name=([^;]+);', attr)
            if s:
                return s.group(1)
            else:
                color_term("[GET_GENE_NAME] No 'gene_name' in attributes({}:{}: {}".format(
                    self.chr, self.pos, attr), 'WRN')
                return None

        def query_gene_intergenic():
            sql_l = "SELECT * FROM `{tab}` " \
                    "WHERE SeqID = %s AND Feature = %s AND End < %s ORDER BY End DESC LIMIT 1".format(
                        tab=self.refer_tab)
            r_l = self.m_con.query(sql_l, [self.chr, 'gene', self.pos]).fetchone()
            sql_r = "SELECT * FROM `{tab}` " \
                    "WHERE SeqID = %s AND Feature = %s AND Start > %s ORDER BY Start LIMIT 1".format(
                        tab=self.refer_tab)
            r_r = self.m_con.query(sql_r, [self.chr, 'gene', self.pos]).fetchone()
            gene_l = get_gene_name(r_l[-1]) if r_l else '-'
            gene_r = get_gene_name(r_r[-1]) if r_r else '-'
            self.gene = "intergenic({}, {})".format(gene_l, gene_r)

        sql = "SELECT * FROM `{tab}` " \
              "WHERE SeqID = %s AND Feature = %s AND Start <= %s AND End >= %s".format(
                tab=self.refer_tab)

        c = self.m_con.query(sql, [self.chr, 'gene', self.pos, self.pos])
        r = c.fetchall()
        self.gene_hits = len(r)

        if self.gene_hits == 0:
            if self.verbose:
                color_term("[QUERY_GENE] No found for {}:{}, try intergenic".format(
                    self.chr, self.pos), 'WRN')
            query_gene_intergenic()
        else:
            if self.verbose and self.gene_hits > 1:
                color_term("[QUERY_GENE] {} entries found for {}:{}, ignore[1:]".format(
                    self.gene_hits, self.chr, self.pos), 'WRN')
            self.gene = get_gene_name(r[0][-1])
