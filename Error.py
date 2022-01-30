class Error(Exception):
    """
    Base class for user-defined exceptions

    """

    def __init__(self):
        raise NotImplementedError("Implemented by subclasses")

    def __str__(self):
        raise NotImplementedError("Implemented by subclasses")

class SequenceInvalidError(Error):
    """
    Exception raised for invalid sequence

    Attributes:
        seq -- Sequence which caused the exception
        msg -- Explanation of the error
    """

    def __init__(self, seq, msg="Invalid Sequence"):
        self.seq = seq
        if len(self.seq) > 10:
            self.seq = self.seq[:11] + "..."
        self.msg = msg
    
    def __str__(self):
        return f"{self.msg}: {self.seq}"
    