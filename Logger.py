from threading import Lock
import datetime

class Logger:
    """
    A thread safe singleton logger

    Attributes:
        __instance    -- The current instance of logger
        __initialized -- True if logger has been initialized, otherwise false
        __lock        -- Ensures only one process has access to logger at a given time

    Methods:
        __new__     -- Acquires lock and initializes an instannce
        __init__    -- Acquires lock and creates a file with initial signature
        __del__     -- Closes file upon logger destruction
        Log         -- Writes to file and flushes
    """
    __instance = None
    __initialized = False
    __lock = Lock()

    def __new__(cls):
        if not cls.__instance:
            with cls.__lock:
                if not cls.__instance:
                    cls.__instance = super(Logger, cls).__new__(cls)
        return cls.__instance

    def __init__(self):
        if not self.__initialized:
            with self.__lock:
                if not self.__initialized:
                    self.file = open(f"Log_{datetime.datetime.now().strftime(r'%Y-%m-%d')}.log", 'w+')
                    self.file.write(f"--- Logged by thread safe singleton logger on {datetime.datetime.now().strftime(f'%Y-%m-%d')} ---\n")
                    self.__initialized = True
        return

    def __del__(self):
        with self.__lock:
            self.file.close()
        return

    def log(self, msg: str):
        with self.__lock:
            self.file.write(f"{datetime.datetime.now().strftime(r'%Y-%m-%d %H:%M:%S')}: {str(msg)}\n")
            self.file.flush()
        return
