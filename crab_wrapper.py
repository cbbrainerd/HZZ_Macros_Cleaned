#!/usr/bin/env python
#Magic debug stuff. Change shebang to #!/usr/bin/env python to remove
''':'
env
ls
touch output.root output_bnn.txt output_bnn.root output_txt.txt output_txt_vbf.txt 
exec python "$0" "$@"
'''

import subprocess
import PSet

fileList=list(PSet.process.source.fileNames)
with open('fileList','w') as f:
    for fn in fileList:
        f.write('%s\n' % fn)

with open('FrameworkJobReport.xml','w') as f:
    f.write(
'''<FrameworkJobReport>
<ReadBranches>
</ReadBranches>
<PerformanceReport>
  <PerformanceSummary Metric="StorageStatistics">
    <Metric Name="Parameter-untracked-bool-enabled" Value="true"/>
    <Metric Name="Parameter-untracked-bool-stats" Value="true"/>
    <Metric Name="Parameter-untracked-string-cacheHint" Value="application-only"/>
    <Metric Name="Parameter-untracked-string-readHint" Value="auto-detect"/>
    <Metric Name="ROOT-tfile-read-totalMegabytes" Value="0"/>
    <Metric Name="ROOT-tfile-write-totalMegabytes" Value="0"/>
  </PerformanceSummary>
</PerformanceReport>
<GeneratorInfo>
</GeneratorInfo>
</FrameworkJobReport>''')

import sys
import os

import traceback

class dotDict(dict):
    def __getattr__(self,attr):
        return self[attr]

def parse_command_line(argv):
    args=dotDict(tuple(x.split('=',1)) for x in argv[2:])
    args['job_number']=argv[1]
    return args

def main(argv=None):
    if argv==None:
        argv=sys.argv
    arguments=parse_command_line(argv)
    #print arguments
    #SITE="$1"
    #YEAR="$2"
    #MC="$3"
    #SAMPLE="$4"
    #n="$5"
    
    print("Environment:",os.environ)
    print(os.listdir("%s/src/ZZMatrixElement/MELA/data/%s" % (os.environ['CMSSW_BASE'],os.environ['SCRAM_ARCH'])))
    os.execl('/bin/bash','bash','submit.sh',arguments.site,arguments.year,arguments.mc,arguments.sample,arguments.n,arguments.dataset)
    #print subprocess.check_output(['ls'])
    #try:
    #if True:
        #subprocess.call(['bash','compilereference.sh',arguments.sample],stderr=subprocess.STDERR,stdout=subprocess.STDOUT)
        #subprocess.call(['bash','submit.sh',arguments.site,arguments.year,arguments.mc,arguments.sample,arguments.n,arguments.dataset],stderr=subprocess.STDERR,stdout=subprocess.STDOUT)
    #    os.execlp('/bin/bash','bash','submit.sh',arguments.sample)
    #except subprocess.CalledProcessError:
    #    traceback.print_exc()
    #    sys.exit(1)

if __name__=='__main__':
    main()
