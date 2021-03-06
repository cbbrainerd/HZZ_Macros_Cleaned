from __future__ import print_function
from __future__ import division

import CRABClient.Emulator
from CRABClient import __version__
from CRABClient.ClientUtilities import validateJobids,colors
from CRABClient.UserUtilities import getFileFromURL
from CRABClient.Commands.getcommand import getcommand
from CRABClient.ClientExceptions import RESTCommunicationException, MissingOptionException

from ServerUtilities import getProxiedWebDir
from httplib import HTTPException

class getlog(getcommand):
    """
    Retrieve the log files of a number of jobs specified by the -q/--quantity option.
    -q logfiles per exit code are returned if transferLogs = False; otherwise all the log files
    collected by the LogCollect job are returned. The task is identified by the -d/--dir option.
    """
    name = 'getlog'
    shortnames = ['log']
    visible = True #overwrite getcommand

    def __call__(self):
        if self.options.short:
            #Check if splitting is automatic
            try:
                splitting=self.cachedinfo['OriginalConfig'].Data.splitting
            except AttributeError: #Default setting is 'Automatic'
                splitting='Automatic'
            except KeyError: #crab remade task does not have 'OriginalConfig' key, need to fetch from DB
                splitting='Unknown'
            taskname = self.cachedinfo['RequestName']
            inputlist = {'subresource': 'webdir', 'workflow': taskname}
            serverFactory = CRABClient.Emulator.getEmulator('rest')
            server = serverFactory(self.serverurl, self.proxyfilename, self.proxyfilename, version=__version__)
            uri = self.getUrl(self.instance, resource = 'task')
            webdir=None
            if splitting!='Unknown':
                webdir = getProxiedWebDir(taskname, self.serverurl, uri, self.proxyfilename, self.logger.debug)
            if not webdir:
                dictresult, status, reason =  server.get(uri, data = inputlist)
                if status != 200:
                    msg = "Problem retrieving information from the server:\ninput:%s\noutput:%s\nreason:%s" % (str(inputlist), str(dictresult), str(reason))
                    raise RESTCommunicationException(msg)
                if splitting=='Unknown':
                    splitting=getColumn(dictresult,'tm_split_algo')
                webdir = dictresult['result'][0]
                self.logger.info('Server result: %s' % webdir)
            self.setDestination()
            self.logger.info("Setting the destination to %s " % self.dest)
            #check the format of jobids
            self.options.jobids = validateJobids(self.options.jobids,splitting!='Automatic')
            failed, success = self.retrieveShortLogs(webdir, self.proxyfilename)
            if failed:
                msg = "%sError%s: Failed to retrieve the following files: %s" % (colors.RED,colors.NORMAL,failed)
                self.logger.info(msg)
            else:
                self.logger.info("%sSuccess%s: All files successfully retrieved." % (colors.GREEN,colors.NORMAL))
            returndict = {'success': success, 'failed': failed}
        else:
            # Different from the old getlog code: set 'logs2' as subresource so that 'getcommand' uses the new logic.
            returndict = getcommand.__call__(self, subresource = 'logs2')
            if ('success' in returndict and not returndict['success']) or \
               ('failed'  in returndict and returndict['failed']):
                msg = "You can use the --short option to retrieve a short version of the log files from the Grid scheduler."
                self.logger.info(msg)

        return returndict

    def setOptions(self):
        """
        __setOptions__

        This allows to set specific command options
        """
        self.parser.add_option( '--quantity',
                                dest = 'quantity',
                                help = 'The number of logs you want to retrieve (or "all"). Ignored if --jobids is used.' )
        self.parser.add_option( '--parallel',
                                dest = 'nparallel',
                                help = 'Number of parallel download, default is 10 parallel download.',)
        self.parser.add_option( '--wait',
                                dest = 'waittime',
                                help = 'Increase the sendreceive-timeout in second.',)
        self.parser.add_option( '--short',
                                dest = 'short',
                                default = False,
                                action = 'store_true',
                                help = 'Get the short version of the log file. Use with --dir and --jobids.',)
        getcommand.setOptions(self)


    def validateOptions(self):
        getcommand.validateOptions(self)
        if self.options.short:
            if self.options.jobids is None:
                msg  = "%sError%s: Please specify the job ids for which to retrieve the logs." % (colors.GREEN, colors.NORMAL)
                msg += " Use the --jobids option."
                ex = MissingOptionException(msg)
                ex.missingOption = "jobids"
                raise ex


    def retrieveShortLogs(self, webdir, proxyfilename):
        self.logger.info("Retrieving...")
        success = []
        failed = []
        for _, jobid in self.options.jobids:
            # We don't know a priori how many retries the job had. So we start with retry 0
            # and increase it by 1 until we are unable to retrieve a log file (interpreting
            # this as the fact that we reached the highest retry already).
            retry = 0
            succeded = True
            while succeded:
                filename = 'job_out.%s.%s.txt' % (jobid, retry)
                url = webdir + '/' + filename
                try:
                    getFileFromURL(url, self.dest + '/' + filename, proxyfilename)
                    self.logger.info('Retrieved %s' % (filename))
                    success.append(filename)
                    retry += 1 #To retrieve retried job log, if there is any.
                except HTTPException as hte:
                    succeded = False
                    # Ignore the exception if the HTTP status code is 404. Status 404 means file
                    # not found (see http://www.w3.org/Protocols/rfc2616/rfc2616-sec10.html). File
                    # not found error is expected, since we try all the job retries.
                    if hasattr(hte.args[0], 'status') and hte.args[0].status != 404:
                        self.logger.debug(str(hte))
                        failed.append(filename)

        return failed, success

class nil:
    def __getattr__(self,attr): return self
    def __call__(self,*args,**kwargs): return None
import sys
import trace

from CRABClient.ClientUtilities import initLoggers, flushMemoryLogger, removeLoggerHandlers

tblogger,logger,memhandler=initLoggers()
try:
    trace=trace.Trace()
    #trace.runfunc(getlog(logger,sys.argv[1:]))
    getlog(logger,sys.argv[1:])()
finally:
    flushMemoryLogger(tblogger, memhandler, logger.logfile)
    removeLoggerHandlers(tblogger)
    removeLoggerHandlers(logger)
