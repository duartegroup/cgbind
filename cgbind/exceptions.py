class ArchitectureNotFound(Exception):

    def __init__(self, message):
        super().__init__(message)


class NoXYZs(Exception):
    """Exception for a molecule that had no 3D coordinates (xyzs)"""


class CgbindCritical(Exception):
    """Non recoverable error"""
    def __init__(self, message=''):
        super().__init__(message)


class CannotBuildCSComplex(Exception):
    pass


class RequiresAutodE(Exception):
    pass


class FileMalformatted(Exception):
    """Exception for a file in the wrong format"""
