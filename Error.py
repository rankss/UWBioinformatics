class Error(Exception):
    """[summary]

    Args:
        Exception ([type]): [description]
    """
    def __init__(self, message: str):
        super().__init__()
        self.message = message

    def __str__(self):
        raise NotImplementedError("Implemented by subclasses")

class InvalidSequenceError(Error):
    """[summary]

    Args:
        Error ([type]): [description]
    """
    def __init__(self, message="Invalid Sequence"):
        super().__init__(message)

    def __str__(self):
        return f"{self.message}"
    
class InvalidSequenceTypeError(Error):
    """[summary]

    Args:
        Error ([type]): [description]
    """
    def __init__(self, message="Invalid Sequence Type"):
        super().__init__(message)

    def __str__(self):
        return f"{self.message}"

class InvalidAlignmentTypeError(Error):
    """_summary_

    Args:
        Error (_type_): _description_
    """
    def __init__(self, message="Invalid Alignment Type"):
        super().__init__(message)

    def __str__(self):
        return f"{self.message}"

class InvalidMatrixError(Error):
    """[summary]

    Args:
        Error ([type]): [description]
    """
    def __init__(self, message="Invalid Matrix"):
        super().__init__(message)

    def __str__(self):
        return f"{self.message}"
    