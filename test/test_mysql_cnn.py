#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : test_mysql_cnn
# AUTHOR  : codeunsolved@gmail.com
# CREATED : September 7 2019
# VERSION : v0.0.1

import os
import sys
import random
import unittest
from copy import deepcopy

import numpy as np
import mysql.connector

PACKAGE_PARENT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(PACKAGE_PARENT)

from utils.connector import MysqlConnector


class TestMysqlConnector(unittest.TestCase):

    MYSQL_USER = None
    MYSQL_PASS = None
    MYSQL_HOST = None
    MYSQL_PORT = None
    MYSQL_DATABASE = None

    MYSQL_CON = None

    MYSQL_CONFIG = {}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def setUpClass(cls):
        mysql_config = {
            'user': cls.MYSQL_USER,
            'password': cls.MYSQL_PASS,
            'host': cls.MYSQL_HOST or '127.0.0.1',
            'port': cls.MYSQL_PORT or 3306,
            'database': cls.MYSQL_DATABASE,
            'raise_on_warnings': True
        }

        try:
            mysql_con = MysqlConnector(mysql_config)
        except Exception as e:
            # Skip all test
            cls.skipTest(cls, "Invalid MYSQL config: {}".format(e))
        else:
            cls.MYSQL_CONFIG = mysql_config
            cls.MYSQL_CON = mysql_con

    @classmethod
    def tearDownClass(cls):
        if cls.MYSQL_CON:
            cls.MYSQL_CON.close()

    def test_01_retry_conn(self):
        mysql_config = deepcopy(self.MYSQL_CONFIG)

        mysql_config['port'] = 65535 + 3306  # No such port

        with self.assertRaises(mysql.connector.Error) as cm:
            MysqlConnector(mysql_config, retry_n=3, retry_sleep=0.5)

        exception = cm.exception
        self.assertEqual(exception.retry, 3)

    def test_02_retry_exec(self):
        mysql_config = deepcopy(self.MYSQL_CONFIG)

        # Make new connection
        m_con = MysqlConnector(mysql_config, retry_n=3, retry_sleep=0.5)

        # Kill new connection to make '2013 (HY000): Lost connection' error
        self.MYSQL_CON.query("kill {}".format(m_con.conn.connection_id))

        # It will generate a new connection ID after retry
        m_con.query("SELECT CONNECTION_ID()".format(m_con.conn.connection_id))

        self.assertTrue(bool(m_con.stats['exe_retry']))
