import threading
import datetime

class Logger:
    """
    A thread safe singleton logger

    Attributes:
        instance    -- The current instance of logger
        initialized -- True if logger has been initialized, otherwise false
        lock        -- Ensures only one process has access to logger at a given time

    Methods:
        __new__     -- Acquires lock and initializes an instannce
        __init__    -- Acquires lock and creates a file with initial signature
        __del__     -- Closes file upon logger destruction
        Log         -- Writes to file and flushes
    """

    instance = None
    initialized = False
    lock = threading.Lock()

    def __new__(self):
        if not self.instance:
            with self.lock:
                if not self.instance:
                    self.instance = super(CLogger, self).__new__(self)
        return self.instance

    def __init__(self):
        if not self.initialized:
            with self.lock:
                if not self.initialized:
                    self.file = open(f"Log_{datetime.datetime.now().strftime(f'%Y-%m-%d')}.txt", 'w+')
                    self.file.write(f"--- Logged by thread safe singleton logger on {datetime.datetime.now().strftime(f'%Y-%m-%d')} ---\n")
                    self.initialized = True
        return

    def __del__(self):
        with self.lock:
            self.file.close()
        return

    def Log(self, msg: str):
        with self.lock:
            self.file.write(F"{datetime.datetime.now().strftime(f'%Y-%m-%d %H:%M:%S')}: {str(msg)}\n")
            self.file.flush()
        return
