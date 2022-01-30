from Logger import Logger

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

    def __init__(self, sequence, message="Invalid Sequence"):
        self.sequence = sequence
        if len(self.seq) > 10:
            self.sequence = self.sequence[:11] + "..."
        self.message = sequence
    
    def __str__(self):
        return f"{self.message}: {self.sequence}"
    