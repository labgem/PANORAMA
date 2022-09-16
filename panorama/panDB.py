#!/usr/bin/env python3
# coding:utf-8
import concurrent.futures
import logging
# default libraries
from threading import Lock
from typing import Union, List
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

# installed librairies
from neo4j import GraphDatabase, Session

try:
    from neo4j._async.driver import AsyncGraphDatabase
except ModuleNotFoundError:
    print('Error! You should be running neo4j python driver version 5 to use async features')

# local librairies
from panorama.pangenomes import Pangenome


class Neo4jDB:

    def __init__(self, uri, user, pwd):
        self.__uri = uri
        self.__user = user
        self.__pwd = pwd
        self.__driver = None
        try:
            self.__driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
            self.__async_driver = AsyncGraphDatabase.driver(self.__uri, auth=(self.__user, self.__pwd))
        except Exception as e:
            raise Exception("Failed to create the driver:", e)
        self._lock = Lock()

    def close(self):
        assert self.__driver is not None, "Driver not initialized!"
        try:
            self.__driver.close()
        except Exception:
            raise Exception("Close driver failed")

    async def close_async(self):
        assert self.__async_driver is not None, "Async driver not initialized!"
        try:
            await self.__async_driver.close()
        except Exception:
            raise Exception("Close async driver failed")

    def run(self, query: str, parameters: dict = None, db: str = None):
        assert self.__driver is not None, "Driver not initialized!"
        session = None
        response = None
        if len(query) > 0:
            try:
                session = self.__driver.session(database=db) if db is not None else self.__driver.session()
            except Exception as e:
                raise Exception("Init session failed:", e)
            finally:
                try:
                    response = session.run(query, parameters)
                except Exception as e:
                    raise Exception("Request failed", e)
                finally:
                    if session is not None:
                        session.close()
        else:
            raise ValueError("Query is empty")
        return response

    def clean(self, max_rows: int = 1000):
        self.run(query=f"MATCH (n) CALL {{WITH n DETACH DELETE n}} IN TRANSACTIONS OF {max_rows} ROWS;")


if __name__ == "__main__":
    neo4j_db = Neo4jDB(uri="bolt://localhost:7687",
                       user="neo4j",
                       pwd="PANORAMA2022")
