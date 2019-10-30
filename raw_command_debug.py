import CRABAPI
import logging
from CRABClient.ClientUtilities import initLoggers, flushMemoryLogger, removeLoggerHandlers, setConsoleLogLevelVar

def crabCommandDebug(command,*args,**kwargs):
    arguments=[]
    for key, val in kwargs.iteritems():
        if isinstance(val, bool):
            if val:
                arguments.append('--'+str(key))
        else:
            arguments.append('--'+str(key))
            arguments.append(val)
    arguments.extend(list(args))

    return execRawDebug(command, arguments)


def execRawDebug(command,args): #Takes actual command rather than a string
    setConsoleLogLevelVar(logging.DEBUG)
    tblogger,logger,memhandler=initLoggers()
    try:
        cmdobj=command(logger,args)
        res=cmdobj()
    except SystemExit as se:
        if se.code == 2:
            raise CRABAPI.BadArgumentException
    finally:
        flushMemoryLogger(tblogger,memhandler,logger.logfile)
        removeLoggerHandlers(tblogger)
        removeLoggerHandlers(logger)
    return res
