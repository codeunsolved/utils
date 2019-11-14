#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PROGRAM : connector
# AUTHOR  : codeunsolved@gmail.com
# CREATED : June 14 2017
# VERSION : v0.0.5
# UPDATE  : [v0.0.1] March 21 2018
# 1. add :PostgresqlConnector: with `get()` and exceptions like Djanngo;
# 2. optimize :MysqlConnector: as :PostgresqlConnector:;
# UPDATE  : [v0.0.2] November 11 2018
# 1. add `reconnect()` to :MysqlConnector: and :PostgresqlConnector:;
# 2. add user-defined logger to :MysqlConnector: and :PostgresqlConnector:;
# 3. optimize error code/number for :MysqlConnector: and :PostgresqlConnector:;
# UPDATE  : [v0.0.3] November 29 2018
# 1. use conditional import;
# 2. change line endings to Unix for minimum requirements;
# UPDATE  : [v0.0.4] December 4 2018
# 1. add :MongoConnector:;
# UPDATE  : [v0.0.5] November 13 2019
# 1. add dictionary support for :MysqlConnector:;


import time


class DoesNotExist(Exception):
    pass


class MultipleObjectsReturned(Exception):
    pass


def stats_performance(func):
    def wrapper(*args, **kw):
        self = args[0]
        START = time.time()
        results = func(*args, **kw)
        COST = time.time() - START
        self.log("[{}] {} rows affected, cost {}s".format(
            self.__class__.__name__,
            self.cursor.rowcount,
            round(COST, 7)))
        return results
    return wrapper


class MysqlConnector(object):
    mysql = __import__('mysql.connector')
    """Connect MySQL and execute SQL on it.

    :param config: config for new class.
        config = {
           'user': 'username',
           'password': 'password',
           'host': '127.0.0.1',
           'port': '',  # default: 3306
           'database': 'database',
           'raise_on_warnings': True
        }
    """

    def __init__(self, config, dictionary=False, verbose=False, logger=print):
        self.config = config
        self.dictionary = dictionary
        self.verbose = verbose
        self.logger = logger
        # Add exceptions
        self.DoesNotExist = DoesNotExist
        self.MultipleObjectsReturned = MultipleObjectsReturned
        # Connect
        self.connect(self.config)
        # Select DB if specified
        if 'database' in self.config:
            self.select_db(self.config['database'])

    def log(self, msg, verbose=False):
        if verbose or self.verbose:
            self.logger(msg)

    def _commit(self):
        self.conn.commit()

    def _rollback(self):
        self.conn.rollback()

    @stats_performance
    def _execute(self, sql, vals, retry=True):
        try:
            self.cursor.execute(sql, vals)
        except Exception as e:
            if e.errno in [2006,   # Error: MySQL server has gone away
                           2013]:  # Error: Lost connection to MySQL server during query
                self.log("[MysqlConnector] {}{}".format(e, ', retry!' if retry else ''), verbose=True)
                if retry:
                    self._execute(sql, vals, retry=False)
                else:
                    raise e
            elif e.errno == 1062:   # Error: Duplicate entry '%s' for key %d
                self.log("[MysqlConnector] {}".format(e), verbose=True)
            else:
                raise e
        else:
            self._commit()

    def connect(self, config):
        try:
            self.conn = self.mysql.connector.connect(**config)
        except self.mysql.connector.Error as e:
            if e.errno == 1045:    # Error: Access denied for user '%s'@'%s' (using password: %s)
                raise Exception("[MysqlConnector] Connect failed! Username or password incorrect.")
            else:
                raise e
        else:
            # cursor 'buffered=True' needed to use .fetch* and .rowcount
            self.cursor = self.conn.cursor(buffered=True, dictionary=self.dictionary)

    def select_db(self, db_name):
        msg = "[MysqlConnector] • Select DB: {} ".format(db_name)
        try:
            self.conn.database = db_name
        except self.mysql.connector.Error as e:
            if e.errno == 1049:    # Error: Unknown database '%s'
                msg += "Failed! {}".format(e.msg)
            else:
                raise e
        else:
            msg += "OK!"
        self.log(msg)

    def create_db(self, db_name):
        msg = "[MysqlConnector] • Create database: {} ".format(db_name)
        try:
            self.cursor.execute("CREATE DATABASE {} DEFAULT CHARACTER SET 'utf8'".format(db_name))
        except self.mysql.connector.Error as e:
            raise e
        else:
            msg += "OK!\n"
        self.log(msg, verbose=True)

    def select(self, tab, where_dict, cols_selected='*', filter_='='):
        cols, vals = zip(*where_dict.items())
        expr_select = ', '.join(cols_selected)
        expr_where = ' AND '.join(["{} {} %s".format(x, filter_) for x in cols])

        sql = "SELECT {} FROM `{}` WHERE {}".format(expr_select, tab, expr_where)
        self._execute(sql, vals)
        return self.cursor

    def insert(self, tab, insert_dict):
        """

        :Example:
        add_employee = ("INSERT INTO employees "
                        "(first_name, last_name, hire_date, gender, birth_date) "
                        "VALUES (%s, %s, %s, %s, %s)")
        data_employee = ('Geert', 'Vanderkelen', tomorrow, 'M', date(1977, 6, 14))
        # Insert new employee
        cursor.execute(add_employee, data_employee)
        emp_no = cursor.lastrowid

        add_salary = ("INSERT INTO salaries "
                      "(emp_no, salary, from_date, to_date) "
                      "VALUES (%(emp_no)s, %(salary)s, %(from_date)s, %(to_date)s)")
        data_salary = {
            'emp_no': emp_no,
            'salary': 50000,
            'from_date': tomorrow,
            'to_date': date(9999, 1, 1),
        }
        # Insert salary information
        cursor.execute(add_salary, data_salary)
        # Make sure data is committed to the database
        cnx.commit()`

        :Refer:
        1. [5.3 Inserting Data Using Connector/Python](http://dev.mysql.com/doc/connector-python/en/
           connector-python-example-cursor-transaction.html)
        """
        cols, vals = zip(*insert_dict.items())
        expr_insert = ', '.join(cols)
        expr_values = ', '.join(['%s'] * len(vals))

        sql = "INSERT INTO `{}` ({}) VALUES ({})".format(tab, expr_insert, expr_values)
        self._execute(sql, vals)

    def update(self, tab, where_dict, update_dict):
        cols_w, vals_w = zip(*where_dict.items())
        cols_u, vals_u = zip(*update_dict.items())
        expr_set = ', '.join(["{} = %s".format(x) for x in cols_u])
        expr_where = ' AND '.join(["{} = %s".format(x) for x in cols_w])

        sql = "UPDATE `{}` SET {} WHERE {}".format(tab, expr_set, expr_where)
        self._execute(sql, vals_u + vals_w)

    def query(self, sql, vals=None):
        self._execute(sql, vals)
        return self.cursor

    def get(self, tab, query_dict, cols_selected='*', filter_='='):
        self.cursor = self.select(tab, query_dict, cols_selected, filter_)
        if self.cursor.rowcount == 1:
            return self.cursor.fetchone()
        elif self.cursor.rowcount == 0:
            raise self.DoesNotExist("No record found!")
        else:
            raise self.MultipleObjectsReturned("Get {} records! Only one excepted.".format(self.cursor.rowcount))

    def close(self):
        self.cursor.close()
        self.conn.close()

    def reconnect(self):
        self.connect(self.config)
        # Select DB if specified
        if 'database' in self.config:
            self.select_db(self.config['database'])


class PostgresqlConnector(object):
    psycopg2 = __import__('psycopg2')
    """Connect PostgreSQL and execute SQL on it.

    :param config: config for new class.
        config = {
           'user': 'username',
           'password': 'password',
           'host': '127.0.0.1',
           'port': '',  # default: 5432
           'database': 'database',
        }
    """

    def __init__(self, config, verbose=False, logger=print):
        self.config = config
        self.verbose = verbose
        self.logger = logger
        # Add exceptions
        self.DoesNotExist = DoesNotExist
        self.MultipleObjectsReturned = MultipleObjectsReturned
        # Connect
        self.connect()

    def log(self, msg, verbose=False):
        if verbose or self.verbose:
            self.logger(msg)

    def _commit(self):
        self.conn.commit()

    def _rollback(self):
        self.conn.rollback()

    @stats_performance
    def _execute(self, sql, vals):
        try:
            self.cursor.execute(sql, vals)
        except self.psycopg2.InternalError as e:
            # pgcode: '25P02' (in_failed_sql_transaction)
            # pgerror: 'ERROR:  current transaction is aborted, commands ignored until end of transaction block'
            if e.pgcode == '25P02':
                self._rollback()
                self._execute(sql, vals)
            else:
                raise e
        else:
            self._commit()

    def connect(self):
        self.conn = self.psycopg2.connect(**self.config)
        self.cursor = self.conn.cursor()

    def select(self, tab, where_dict, cols_selected='*', filter_='='):
        cols, vals = zip(*where_dict.items())
        expr_select = ', '.join(cols_selected)
        expr_where = ' AND '.join(["{} {} %s".format(x, filter_) for x in cols])

        sql = "SELECT {} FROM \"{}\" WHERE {}".format(expr_select, tab, expr_where)
        self._execute(sql, vals)
        return self.cursor

    def insert(self, tab, insert_dict):
        cols, vals = zip(*insert_dict.items())
        expr_insert = ', '.join(cols)
        expr_values = ', '.join(['%s'] * len(vals))

        sql = "INSERT INTO \"{}\" ({}) VALUES ({})".format(tab, expr_insert, expr_values)
        self._execute(sql, vals)

    def update(self, tab, where_dict, update_dict):
        cols_w, vals_w = zip(*where_dict.items())
        cols_u, vals_u = zip(*update_dict.items())
        expr_set = ', '.join(["{} = %s".format(x) for x in cols_u])
        expr_where = ' AND '.join(["{} = %s".format(x) for x in cols_w])

        sql = "UPDATE \"{}\" SET {} WHERE {}".format(tab, expr_set, expr_where)
        self._execute(sql, vals_u + vals_w)

    def query(self, sql, vals=None):
        self._execute(sql, vals)
        return self.cursor

    def get(self, tab, query_dict, cols_selected='*', filter_='='):
        self.cursor = self.select(tab, query_dict, cols_selected, filter_)
        if self.cursor.rowcount == 1:
            return self.cursor.fetchone()
        elif self.cursor.rowcount == 0:
            raise self.DoesNotExist("No record found!")
        else:
            raise self.MultipleObjectsReturned("Get {} records! Only one excepted.".format(self.cursor.rowcount))

    def close(self):
        self.cursor.close()
        self.conn.close()

    def reconnect(self):
        self.connect()


class MongoConnector(object):
    pymongo = __import__('pymongo')
    """Connect MongoDB and execute operations on it.

    :param config: config for new class.
        config = {
           'username': '',
           'password': '',
           'host': '127.0.0.1',
           'port': '',  # default: 27017
           'database': 'database',
        }
    """

    def __init__(self, config, verbose=False, logger=print):
        self.config = config
        self.verbose = verbose
        self.logger = logger
        # Add exceptions
        self.DoesNotExist = DoesNotExist
        self.MultipleObjectsReturned = MultipleObjectsReturned
        # Connect
        self.uri = None
        self.connect(self.config)

    def log(self, msg, verbose=False):
        if verbose or self.verbose:
            self.logger(msg)

    def gen_uri(self, config):
        if 'username' in config and 'password' in config and 'authsource' in config:
            if config['username'] and config['password'] and config['authsource']:
                cred = "{u}:{p}@".format(u=config['username'], p=config['password'])
                auth = "?authSource={a}".format(a=config['authsource'])
            else:
                cred = ''
                auth = ''
        else:
            cred = ''
            auth = ''

        if 'host' not in config or not config['host']:
            host = 'localhost'
        else:
            host = config['host']

        if 'port' not in config or not config['port']:
            port = 27017
        else:
            port = config['port']

        if 'database' not in config or not config['database']:
            database = '/'
        else:
            database = "/{}".format(config['database'])

        uri = "mongodb://{cred}{host}:{port}{database}{auth}".format(
            cred=cred, host=host, port=port, database=database, auth=auth)
        return uri

    def connect(self, config, serverSelectionTimeoutMS=1000):
        try:
            self.uri = self.gen_uri(self.config)
            self.conn = self.pymongo.MongoClient(
                self.uri,
                serverSelectionTimeoutMS=serverSelectionTimeoutMS)
            # A cheap way to find out that the connection is established
            self.conn.admin.command('ismaster')
        except Exception as e:
            raise e
        else:
            pass

    def select_db(self, db_name):
        return self.conn[db_name]
