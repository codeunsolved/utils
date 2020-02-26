#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : test_ac
# AUTHOR  : codeunsolved@gmail.com
# CREATED : September 7 2019
# VERSION : v0.0.1

import os
import sys
import random
import unittest
from copy import deepcopy

import numpy as np

PACKAGE_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(PACKAGE_PARENT)

from utils.ngs import AnnCoordinates
from utils.connector import MysqlConnector


class TestAnnCoordinates(unittest.TestCase):

    RANDOM_SEED = 0

    MYSQL_USER = None
    MYSQL_PASS = None
    MYSQL_HOST = None
    MYSQL_PORT = None
    MYSQL_DATABASE = None

    MYSQL_CON = None

    default_args = {
        'config': None,
        'use_cache': True,
        'logger_level': 50,  # logging.CRITICAL
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def set_random_seed(cls, random_seed):
        print("\033[1;32m#\033[0m Random seed: \033[1;35m{}\033[0m".format(
            random_seed))
        np.random.seed(random_seed)

    @classmethod
    def setUpClass(cls):
        # !important! set random seed
        cls.set_random_seed(cls.RANDOM_SEED)

        mysql_config = {
            'user': cls.MYSQL_USER,
            'password': cls.MYSQL_PASS,
            'host': cls.MYSQL_HOST or '127.0.0.1',
            'port': cls.MYSQL_PORT or 3306,
            'database': cls.MYSQL_DATABASE,
            'raise_on_warnings': True
        }

        if not mysql_config['database']:
            cls.skipTest(cls, "Invalid MySQL databse: {}".format(
                mysql_config['database']))

        try:
            mysql_con = MysqlConnector(
                mysql_config)
        except Exception as e:
            # Skip all test
            cls.skipTest(cls, "Invalid MySQL config: {}".format(e))
        else:
            cls.default_args['config'] = mysql_config
            cls.MYSQL_CON = mysql_con

    @classmethod
    def tearDownClass(cls):
        if cls.MYSQL_CON:
            cls.MYSQL_CON.close()

    def test_00_mysql_con(self):
        args = deepcopy(self.default_args)

        # Connect with MYSQL config
        ac = AnnCoordinates(**args)
        self.assertIsNotNone(ac.m_con)

        # Connect with MYSQL connection
        args['config'] = None
        args['mysql_con'] = self.MYSQL_CON
        ac = AnnCoordinates(**args)
        self.assertIsNotNone(ac.m_con)

    def _test_ann_alk_(self, args):
        """
        ALK Gene Structure
        Chromosome: chr2
        Genbank ID: NM_004304 Orientation: -
        Length coding sequence : 4860 nucleotides.
        Region  start       end         region length   phase at end
        UTR     30,143,526  30,144,477
        Exon    30,142,859  30,143,525  667             1
        Exon    29,940,444  29,940,563  120             1
        Exon    29,917,716  29,917,880  165             1
        Exon    29,754,781  29,754,982  202             2
        Exon    29,606,598  29,606,725  128             1
        Exon    29,551,216  29,551,347  132             1
        Exon    29,543,617  29,543,748  132             1
        Exon    29,541,170  29,541,270  101             0
        Exon    29,519,754  29,519,923  170             2
        Exon    29,498,268  29,498,362  95              1
        Exon    29,497,965  29,498,093  129             1
        Exon    29,473,971  29,474,133  163             2
        Exon    29,462,546  29,462,696  151             0
        Exon    29,456,431  29,456,562  132             0
        Exon    29,455,170  29,455,314  145             1
        Exon    29,451,750  29,451,932  183             1
        Exon    29,450,440  29,450,538  99              1
        Exon    29,449,788  29,449,940  153             1
        Exon    29,448,327  29,448,431  105             1
        Exon    29,446,208  29,446,394  187             2
        Exon    29,445,383  29,445,473  91              0
        Exon    29,445,210  29,445,274  65              2
        Exon    29,443,572  29,443,701  130             0
        Exon    29,436,850  29,436,947  98              2
        Exon    29,432,652  29,432,744  93              2
        Exon    29,430,037  29,430,138  102             2
        Exon    29,420,408  29,420,542  135             2
        Exon    29,419,636  29,419,726  91              0
        Exon    29,416,090  29,416,788  699             0
        UTR     29,415,640  29,416,089
        """
        # ALK chr2:29,415,640-30,144,477(GRCh37/hg19) from GeneCards
        # ALK chr2:29,415,640-30,144,432(GRCh37/hg19) from GENCODE/HAVANA
        # ENST00000389048 -> NM_004304

        # Case01: Normal
        ac = AnnCoordinates(**args)
        # 5'UTR
        ac.query('2:30,143,526')
        gene_info = ac.gene_info
        self.assertEqual(gene_info['id'], 'ENSG00000171094.17')
        self.assertEqual(gene_info['name'], 'ALK')
        self.assertEqual(gene_info['type'], 'protein_coding')
        self.assertEqual(gene_info['strand'], '-')
        rank_info = ac.rank_info
        self.assertEqual(rank_info['trans_id'], 'ENST00000389048.7')
        self.assertEqual(rank_info['trans_start'], 29415640)
        self.assertEqual(rank_info['trans_end'], 30144432)
        self.assertEqual(rank_info['trans_strand'], '-')
        self.assertEqual(rank_info['exon_total'], 29)
        self.assertEqual(rank_info['cds_min'], 1)
        self.assertEqual(rank_info['cds_max'], 29)
        self.assertEqual(rank_info['rank'], "5'UTR")
        self.assertEqual(rank_info['exon_num'], 1)
        # Exon20
        bp = np.random.choice(range(29446208, 29446394 + 1))
        ac.query("2:{}".format(bp))
        gene_info = ac.gene_info
        self.assertEqual(gene_info['id'], 'ENSG00000171094.17')
        self.assertEqual(gene_info['name'], 'ALK')
        self.assertEqual(gene_info['type'], 'protein_coding')
        self.assertEqual(gene_info['strand'], '-')
        rank_info = ac.rank_info
        self.assertEqual(rank_info['trans_id'], 'ENST00000389048.7')
        self.assertEqual(rank_info['trans_start'], 29415640)
        self.assertEqual(rank_info['trans_end'], 30144432)
        self.assertEqual(rank_info['trans_strand'], '-')
        self.assertEqual(rank_info['exon_total'], 29)
        self.assertEqual(rank_info['cds_min'], 1)
        self.assertEqual(rank_info['cds_max'], 29)
        self.assertEqual(rank_info['rank'], "Exon20")
        self.assertEqual(rank_info['exon_num'], 20)
        # Intron10
        bp = np.random.choice(range(29498093 + 1, 29498268))
        ac.query("2:{}".format(bp))
        gene_info = ac.gene_info
        self.assertEqual(gene_info['id'], 'ENSG00000171094.17')
        self.assertEqual(gene_info['name'], 'ALK')
        self.assertEqual(gene_info['type'], 'protein_coding')
        self.assertEqual(gene_info['strand'], '-')
        rank_info = ac.rank_info
        self.assertEqual(rank_info['trans_id'], 'ENST00000389048.7')
        self.assertEqual(rank_info['trans_start'], 29415640)
        self.assertEqual(rank_info['trans_end'], 30144432)
        self.assertEqual(rank_info['trans_strand'], '-')
        self.assertEqual(rank_info['exon_total'], 29)
        self.assertEqual(rank_info['cds_min'], 1)
        self.assertEqual(rank_info['cds_max'], 29)
        self.assertEqual(rank_info['rank'], "Intron10")
        self.assertEqual(rank_info['exon_num'], 10)
        # 3'UTR
        ac.query('2:29,415,640')
        gene_info = ac.gene_info
        self.assertEqual(gene_info['id'], 'ENSG00000171094.17')
        self.assertEqual(gene_info['name'], 'ALK')
        self.assertEqual(gene_info['type'], 'protein_coding')
        self.assertEqual(gene_info['strand'], '-')
        rank_info = ac.rank_info
        self.assertEqual(rank_info['trans_id'], 'ENST00000389048.7')
        self.assertEqual(rank_info['trans_start'], 29415640)
        self.assertEqual(rank_info['trans_end'], 30144432)
        self.assertEqual(rank_info['trans_strand'], '-')
        self.assertEqual(rank_info['exon_total'], 29)
        self.assertEqual(rank_info['cds_min'], 1)
        self.assertEqual(rank_info['cds_max'], 29)
        self.assertEqual(rank_info['rank'], "3'UTR")
        self.assertEqual(rank_info['exon_num'], 29)

    def test_01_ann_base(self):
        args = deepcopy(self.default_args)

        self._test_ann_alk_(args)

    def test_02_ann_without_cache(self):
        args = deepcopy(self.default_args)

        args['use_cache'] = False
        self._test_ann_alk_(args)

    def test_03_ann_bp_handling(self):
        args = deepcopy(self.default_args)

        ac = AnnCoordinates(**args)

        ac.query('2:30,143,526')
        self.assertEqual(ac.chr_raw, '2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('2:30143526')
        self.assertEqual(ac.chr_raw, '2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('chr2:30,143,526')
        self.assertEqual(ac.chr_raw, 'chr2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('chr2:30143526')
        self.assertEqual(ac.chr_raw, 'chr2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('2', '30,143,526')
        self.assertEqual(ac.chr_raw, '2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('2', '30143526')
        self.assertEqual(ac.chr_raw, '2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('2', 30143526)
        self.assertEqual(ac.chr_raw, '2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('chr2', '30,143,526')
        self.assertEqual(ac.chr_raw, 'chr2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('chr2', '30143526')
        self.assertEqual(ac.chr_raw, 'chr2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

        ac.query('chr2', 30143526)
        self.assertEqual(ac.chr_raw, 'chr2')
        self.assertEqual(ac.chr, 'chr2')
        self.assertEqual(ac.pos, 30143526)

    def test_04_ann_intergenic_gene(self):
        args = deepcopy(self.default_args)

        ac = AnnCoordinates(**args)

        # intergenic(TMEM99,KRT12)
        ac.query('chr17:39013407')
        self.assertEqual(len(ac.gene_hits), 1)
        gene_hit = ac.gene_hits[0]
        self.assertIsNone(gene_hit['id'])
        self.assertEqual(gene_hit['name'], 'intergenic(TMEM99,KRT12)')
        self.assertEqual(gene_hit['type'], 'intergenic')
        self.assertIsNone(gene_hit['source'])
        self.assertIsNone(gene_hit['strand'])

    def test_05_ann_intergenic_igh_gene(self):
        args = deepcopy(self.default_args)

        ac = AnnCoordinates(**args)
        # intergenic(IGHJ6,IGHJ5)
        ac.query('chr14:106329469')
        self.assertEqual(len(ac.gene_hits), 1)
        gene_hit = ac.gene_hits[0]
        self.assertIsNone(gene_hit['id'])
        self.assertEqual(gene_hit['name'], 'IGHJ6')
        self.assertEqual(gene_hit['type'], 'intergenic_igh')
        self.assertIsNone(gene_hit['source'])
        self.assertEqual(gene_hit['strand'], '-')
        rank_info = ac.rank_info
        self.assertEqual(rank_info['trans_strand'], '-')

    def test_06_ann_beyond_transcript_range(self):
        args = deepcopy(self.default_args)

        # BCL6 without specifying transcript map
        ac = AnnCoordinates(**args)
        ac.query('chr3', 187462599)
        gene_info = ac.gene_info
        self.assertEqual(gene_info['id'], 'ENSG00000113916.17')
        self.assertEqual(gene_info['name'], 'BCL6')
        self.assertEqual(gene_info['type'], 'protein_coding')
        self.assertEqual(gene_info['source'], 'HAVANA')
        self.assertEqual(gene_info['strand'], '-')
        rank_info = ac.rank_info
        self.assertEqual(rank_info['trans_id'], 'ENST00000406870.6')
        self.assertEqual(rank_info['trans_start'], 187439165)
        self.assertEqual(rank_info['trans_end'], 187463515)
        self.assertEqual(rank_info['trans_strand'], '-')
        self.assertEqual(rank_info['exon_total'], 10)
        self.assertEqual(rank_info['cds_min'], 3)
        self.assertEqual(rank_info['cds_max'], 10)
        self.assertEqual(rank_info['rank'], "Intron1")
        self.assertEqual(rank_info['exon_num'], 1)

        # BCL6 with specifying transcript map
        args['transcript_map'] = {
            'ENST00000232014': {'gene': 'BCL6', 'refseq': 'NM_001130845'}
        }
        ac = AnnCoordinates(**args)
        ac.query('chr3', 187462599)
        gene_info = ac.gene_info
        self.assertEqual(gene_info['id'], 'ENSG00000113916.17')
        self.assertEqual(gene_info['name'], 'BCL6')
        self.assertEqual(gene_info['type'], 'protein_coding')
        self.assertEqual(gene_info['source'], 'HAVANA')
        self.assertEqual(gene_info['strand'], '-')
        rank_info = ac.rank_info
        self.assertEqual(rank_info['trans_id'], 'ENST00000232014.8')
        self.assertEqual(rank_info['trans_start'], 187440186)
        self.assertEqual(rank_info['trans_end'], 187454357)
        self.assertEqual(rank_info['trans_strand'], '-')
        self.assertEqual(rank_info['exon_total'], 10)
        self.assertEqual(rank_info['cds_min'], 3)
        self.assertEqual(rank_info['cds_max'], 10)
        self.assertEqual(rank_info['rank'], "intragenic")
        self.assertEqual(rank_info['exon_num'], -5)

    def test_07_ann_special_cases(self):
        args = deepcopy(self.default_args)

        # intergenic at left most side
        ac = AnnCoordinates(**args)
        ac.query('chr2', 100)
        gene_info = ac.gene_info
        self.assertIsNone(gene_info['id'])
        self.assertEqual(gene_info['name'], 'intergenic(-,FAM110C)')
        self.assertEqual(gene_info['type'], 'intergenic')
        self.assertEqual(gene_info['start'], 0)
        self.assertEqual(gene_info['end'], 38814)
        self.assertIsNone(gene_info['source'])
        self.assertIsNone(gene_info['strand'])

        # intergenic at right most side
        ac = AnnCoordinates(**args)
        ac.query('chr2', 243199373 - 100)
        gene_info = ac.gene_info
        self.assertIsNone(gene_info['id'])
        self.assertEqual(gene_info['name'], 'intergenic(AC131097.4,-)')
        self.assertEqual(gene_info['type'], 'intergenic')
        self.assertEqual(gene_info['start'], 242844702)
        self.assertEqual(gene_info['end'], 0)
        self.assertIsNone(gene_info['source'])
        self.assertIsNone(gene_info['strand'])
