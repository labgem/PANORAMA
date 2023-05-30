from multiprocessing import Lock

db_loading_lock: Lock = None


def init_db_lock(lock):
    global db_loading_lock
    if db_loading_lock is None:
        db_loading_lock = lock