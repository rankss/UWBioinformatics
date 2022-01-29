import threading
import datetime

class CLogger:
    """
    A thread safe singleton logger.
    """

    instance = None
    initialized = False
    lock = threading.Lock()
    date = datetime.datetime.now()

    def __new__(self):
        if not self.instance:
            with self.lock:
                if not self.instance:
                    print("New Instance Created!")
                    self.instance = super(CLogger, self).__new__(self)
        return self.instance

    def __init__(self):
        if not self.initialized:
            print("Constructor Called")
            with self.lock:
                if not self.initialized:
                    self.file = open(F"Log_{self.date.strftime('%Y-%m-%d')}.txt", 'w+')
                    self.file.write(F"--- Logged by thread safe singleton logger on {self.date.strftime('%Y-%m-%d')} ---\n")
                    self.initialized = True
        return

    def __del__(self):
        print("Destructor Called")
        with self.lock:
            self.file.close()
        return

    def log(self, msg):
        print("Log Called")
        with self.lock:
            self.file.write(F"{self.date.strftime('%Y-%m-%d %H:%M:%S')}: {str(msg)}\n")
            self.file.flush()
        return
