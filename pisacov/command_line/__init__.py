"""
This is PISACov, a program designed to infer quaternary structure
of proteins from evolutionary covariance.
"""

from pisacov import __prog__, __description__, __version__
from pisacov import __author__, __date__, __copyright__

import logging
import sys
import os
import time
import datetime

def unwraptime(delta):
    """
    Convert a :obj:`~datetime.timedelta` into a human readable string.

    :param delta: Time interval.
    :type delta: :obj:`~datetime.timedelta`
    :return: Human-readable description of time interval.
    :rtype: str

    """

    if delta.days > 0:
        wrapped = str(delta.days) + 'days, '
    else:
        wrapped = ''
    seconds = delta.seconds%60
    uwminutes = (delta.seconds - seconds)/60
    minutes = uwminutes%60
    hours = (uwminutes - minutes)/60
    if int(hours) > 0 or wrapped != '':
        wrapped += str(int(hours)) + ' hours, '
    if int(minutes) > 0 or wrapped != '':
        wrapped += str(int(minutes)) + ' minutes, '
    wrapped += str(seconds + delta.microseconds/1000000) + ' seconds.'

    return wrapped


def welcome(command=None):
    """
    Produce execution message.

    :param command: Script name, defaults to None.
    :type command: str, optional
    :return: Execution message.
    :rtype: str

    """
    begin = datetime.datetime.now()
    tzbegin = begin.astimezone().strftime("%d %B %Y, %H:%M:%S %Z")
    if command is not None:
        l = len(command)
        title = '*'*7 + ' '*3 + command + ' '*3 + '*'*(74-l)
    else:
        l = len('PISACov')
        title = '*'*7 + ' '*3 + 'PISACov' + ' '*3 + '*'*(74-l)

    msg = os.linesep + '*'*87 + os.linesep
    msg += title + os.linesep
    msg = '*'*87 + os.linesep
    msg += "** Running "+__prog__+" v."+__version__+' **' + os.linesep
    msg += "** "+__description__ + os.linesep
    msg += "** Developed by "+__author__+". Copyright: "+__copyright__+'.' + os.linesep
    msg += '*'*87 + os.linesep
    msg += '** Time of execution: ' + tzbegin + os.linesep + os.linesep

    return msg, begin

def ok(initime, command=None):
    """
    Produce completion message.

    :param initime: Time of execution.
    :type initime: :class:`~datetime.timedelta`
    :param command: Script name, defaults to None.
    :type command: str, optional
    :return: Completion message.
    :rtype: str

    """
    if command is not None:
        l = len(command)
        title = '*'*7 + ' '*3 + command + ' '*3 + '*'*(67-l)
    else:
        l = len('PISACov')
        title = '*'*7 + ' '*3 + 'PISACov' + ' '*3 + '*'*(67-l)

    msg = os.linesep + '*'*87 + os.linesep
    msg += title + os.linesep
    msg += "** " +__prog__+ " finished **"  + os.linesep
    endtime = datetime.datetime.now()
    tzend = endtime.astimezone().strftime(
        "%d %B %Y, %H:%M:%S %Z")
    msg += '** Time of completion: ' + tzend + os.linesep
    delta = endtime - initime
    seconds = str(delta.days*3600*24 + delta.seconds +
                  delta.microseconds/1000000)
    msg += '** Execution time: ' + seconds + ' seconds or' + os.linesep

    msg += '                   ' + unwraptime(delta) + os.linesep + os.linesep

    return msg

def running(subprocess, done=None):
    """
    Subprocess execution message.

    :param subprocess: Subprocess name.
    :type subprocess: str
    :param done: If concluded, indicate the time of conclusion, defaults to None
    :type done: :class:`~datetime.timedelta`, optional
    :return: Subprocess execution message.
    :rtype: str

    """
    l = len(subprocess)

    if done is None:
        title = '*'*3 + ' '*3 + ' Executing subprocess ' + subprocess + ' '*3 + '*'*(19-l)
        msg = os.linesep + '*'*50 + os.linesep
        msg += title + os.linesep + os.linesep
    else:
        delta = datetime.datetime.now() - done
        title = '*'*3 + ' '*3 + ' Subprocess ' + subprocess + ' completed' + ' '*3 + '*'*(19-l)
        msg = (os.linesep + title + os.linesep +
               '** Subprocess execution time: ' + unwraptime(delta) +
               os.linesep + os.linesep)

    return msg


def pisacov_logger(level="info"):
    """Logger setup.
    :param level: Console logging level, defaults to "info".
                      Options: [ notset | info | debug | warning | error | critical ]
    :type level: str, optional
    :param logfile: The path to a full file log, defaults to None.
    :type logfile: str, optional
    :return: Configured logger.
    :rtype: :obj:`~logging.Logger`
    """
    class PisacovFormatter(logging.Formatter):

        def __init__(self, fmt="%(levelname)s: %(message)s"):
            #logging.Formatter.__init__(self, fmt=fmt)
            super().__init__(fmt="%(levelname)s: %(message)s", datefmt=None, style='%')

        def format(self, record):
            # Save the original format configured by the user
            # when the logger formatter was instantiated
            format_orig = self._style._fmt

            # Replace the original format with one customized by logging level
            if record.levelno == logging.DEBUG:
                self._style._fmt ="%(levelname)s: [%(module)s] (%(lineno)d) %(message)s"
            elif record.levelno == logging.INFO:
                self._style._fmt = "%(message)s"
            else:
                self._style._fmt = "%(levelname)s: %(message)s"

            # Call the original formatter class to do the grunt work
            myformat = logging.Formatter.format(self, record)

            # Restore the original format configured by the user
            self._fmt = format_orig

            return myformat

    logging_levels = {
        "notset": logging.NOTSET,
        "info": logging.INFO,
        "debug": logging.DEBUG,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL
    }

    # Create logger and default settings
    logging.getLogger().setLevel(logging.DEBUG)

    # create console handler with a higher log level
    custom_fmt=PisacovFormatter()
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging_levels.get(level, logging.INFO))
    ch.setFormatter(custom_fmt)
    logging.getLogger().addHandler(ch)

    logging.getLogger().debug("Console logger level: %s", logging_levels.get(level, logging.INFO))

    return logging.getLogger()
