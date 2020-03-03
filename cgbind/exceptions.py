class ArchitectureNotFound(Exception):

    def __init__(self, message):
        super().__init__(message)


class NoXYZs(Exception):
    """Exception for a molecule that had no 3D coordinates (xyzs)"""
