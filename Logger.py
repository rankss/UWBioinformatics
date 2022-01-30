import threading
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
    __lock = threading.Lock()

    def __new__(self):
        if not self.__instance:
            with self.__lock:
                if not self.__instance:
                    self.__instance = super(Logger, self).__new__(self)
        return self.__instance

    def __init__(self):
        if not self.__initialized:
            with self.__lock:
                if not self.__initialized:
                    self.file = open(f"Log_{datetime.datetime.now().strftime(f'%Y-%m-%d')}.txt", 'w+')
                    self.file.write(f"--- Logged by thread safe singleton logger on {datetime.datetime.now().strftime(f'%Y-%m-%d')} ---\n")
                    self.__initialized = True
        return

    def __del__(self):
        with self.__lock:
            self.file.close()
        return

    def Log(self, msg: str):
        with self.__lock:
            self.file.write(f"{datetime.datetime.now().strftime(f'%Y-%m-%d %H:%M:%S')}: {str(msg)}\n")
            self.file.flush()
        return
