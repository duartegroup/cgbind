class ArchitectureNotFound(Exception):
    """Class for an architecture not found"""


class NoXYZs(Exception):
    """Exception for a molecule that had no 3D coordinates (xyzs)"""


class CgbindCritical(Exception):
    """Non recoverable error"""


class CannotBuildCage(Exception):
    """Exception for a cage-substrate complex failing to build"""


class CannotBuildCSComplex(Exception):
    """Exception for a cage-substrate complex failing to build"""


class RequiresAutodE(Exception):
    """Exception for autode being required"""


class FileMalformatted(Exception):
    """Exception for a file in the wrong format"""


class RequiresOpenBabel(Exception):
    """Exception for Open Babel being required"""
