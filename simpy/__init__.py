import logging as log
import sys

class tcol:
    """
    Terminal unicode colors
    """
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'

    WARNING = '\033[93m'
    FAIL = '\033[91m'

    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    ENDC = '\033[0m'


class simpyFormatter(log.Formatter):
    """
    Custom formatter for logging in simpy module
    """

    error_fmt = tcol.BOLD + tcol.FAIL + "%(levelname)s" + tcol.ENDC + \
        " <%(name)s>: " + tcol.BOLD + "%(message)s" + tcol.ENDC

    warning_fmt = tcol.WARNING + "%(levelname)s" + tcol.ENDC + \
        " <%(name)s>: %(message)s"

    debug_fmt = tcol.BOLD + "%(levelname)s" + tcol.ENDC + \
        " <%(name)s>: %(message)s"

    info_fmt = "%(levelname)s <%(name)s>: %(message)s"


    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        log.Formatter.__init__(self, fmt)


    def format(self, record):
        format_orig = self._fmt

        if record.levelno == log.DEBUG:
            self._fmt = self.debug_fmt

        elif record.levelno == log.INFO:
            self._fmt = self.info_fmt

        elif record.levelno == log.WARNING:
            self._fmt = self.warning_fmt

        elif record.levelno == log.ERROR:
            self._fmt = self.error_fmt

        result = log.Formatter.format(self, record)
        self._fmt = format_orig

        return result


fmt = simpyFormatter()
hdlr = log.StreamHandler(sys.stdout)

hdlr.setFormatter(fmt)
log.root.addHandler(hdlr)
log.root.setLevel(log.INFO)
#log.basicConfig(format="%(levelname)s <%(name)s>: %(message)s", level=10)
