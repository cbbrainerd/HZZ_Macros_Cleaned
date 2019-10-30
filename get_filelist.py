try:
    from CRABClient.Commands.checkwrite import checkwrite
    from WMCore.Services.PhEDEx.PhEDEx import PhEDEx
except ImportError:
    print 'Warning: CRAB environment is not setup. Unable to search for files to run over.'
    checkwrite=None
import subprocess
import json
import traceback
import httplib

if checkwrite is not None:
    class getfilenames(checkwrite):
        name='checkwrite'
        shortnames=['chk']
        def setOptions(self):
            self.parser.add_option('--site',dest='sitename',default=None)
            self.parser.add_option('--lfn',dest='userlfn',default=None)
        def validateOptions(self): #fixme
            pass
        def __call__(self):
            self.lfnsadds=self.options.userlfn.split(',')
            #pfndict=self.phedex.getPFN(nodes=[self.options.sitename],lfns=[lfnsadd])
            self.rootfiles=[]
            self.failed_dirs=dict()
            while len(self.lfnsadds)>0:
                try:
                    lfnsadd=self.lfnsadds.pop(0)
                    #pfn=pfndict[(self.options.sitename,lfnsadd)]
                    pfn=self.phedex.getPFN(nodes=[self.options.sitename],lfns=[lfnsadd])[(self.options.sitename,lfnsadd)]
                    ls_cmd=['/usr/bin/gfal-ls','-l',pfn]
                    #if __name__ == "__main__":
                    #    print lfnsadd
                    #    print pfn
                    #    print ls_cmd
                    lsprocess=subprocess.Popen(ls_cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    lsout,lserr=lsprocess.communicate()
                    returncode=lsprocess.returncode
                    if lserr:
                        print lserr
                    if returncode!=0:
                        self.failed_dirs[pfn]={ 'pfn' : pfn , 'lfnsadd' : lfnsadd , 'returncode' : returncode }
                        print 'Failed with status code %i' % returncode
                        continue
                    lines=lsout.split('\n')
                    dirs=['%s/%s' % (lfnsadd,line.split()[8]) for line in lines if line and line[0]=='d' and line.split()[8] not in ['log','failed']]
                    #pfndict.update(self.phedex.getPFN(nodes=[self.options.sitename],lfns=dirs))
                    self.lfnsadds+=dirs
                    print '%i files added!' % len([self.rootfiles.append('%s/%s' % (lfnsadd,line.split()[8])) for line in lines if line and line.split()[8][-5:]=='.root'])
                except httplib.HTTPException:
                    print "HTTPException %s" % lfnsadd
                    self.lfnsadds.append(lfnsadd)
            return ({'root_files' : self.rootfiles , 'failed_dirs' : self.failed_dirs})
else:
    def getfilenames(*args,**kwargs):
        raise ImportError('Must have crab environment setup!')
            
class nil:
    def __getattr__(self,attr):
        return self
    def __call__(self,*args,**kwargs):
        return None
#print getfilenames(nil(),['--site','T2_IT_Bari','--lfn','/store/user/cbrainer/DoubleMuon'])()

import atexit
import shutil
class dataset_filenames(object):
    def __call__(self,*args,**kwargs):
        return self.get_filenames(*args,**kwargs)
    def __init__(self):
        self.changed=False
        atexit.register(self.save)
        try:
            with open('dataset_filenames.json') as f:
                self.files=json.load(f)
        except IOError:
            self.changed=True
            self.files={}
    #root_dir can be a list
    def get_filenames(self,primary_dataset,site=None,root_dirs=None):
        try:
            entry=self.files[primary_dataset]
            if root_dirs is not None and root_dirs!=entry['root_dirs']:
                raise KeyError(primary_dataset)
            if site is not None and site!=entry['site']:
                raise KeyError(primary_dataset)
            return entry['files']
        except KeyError:
            if root_dirs is None or site is None:
                raise KeyError(primary_dataset)
            files=[]
            failed_dirs=[]
            for root_dir in root_dirs:
                getfilenames_retval=getfilenames(nil(),['--site',site,'--lfn',root_dir])()
                files+=getfilenames_retval['root_files']
                failed_dirs+=getfilenames_retval['failed_dirs']
                if failed_dirs:
                    print primary_dataset
                    print failed_dirs
            self.files[primary_dataset]={
                'root_dirs':root_dirs,
                'site':site,
                'files':files,
                'failed_dirs':failed_dirs,
            }
            self.changed=True
            return self.files[primary_dataset]['files']
    def save(self):
        if not self.changed:
            return
        try:
            with open('.dataset_filenames.json.swp','w') as f:
                #May as well sort the list
                #Sorts by year (2017 before 2018) then by run (2017B before 2017C) and finally by job number
                for v in self.files.values():
                    try:
                        v['files'].sort(key=lambda x: (x.split('_Run')[1][0:5],int(x.split('_')[-1].split('.root')[0])))
                    except IndexError:
                        pass
                json.dump(self.files,f)
            shutil.move('.dataset_filenames.json.swp','dataset_filenames.json')
            print("Saved list of filenames.")
        except IOError:
            traceback.print_exc()
            print("Failed to save list of filenames.")

get_filelist=dataset_filenames()

def taskname_to_lfn(taskname,primary_dataset=None):
    (number,task)=taskname.split(':',1)
    (username,taskdir)=task.split('_',1)
    if primary_dataset is None:
        (crab,analysis,primary_dataset)=taskdir.split('_')[0:3]
        assert(crab=='crab')
    return '/store/user/%s/%s/%s/%s' % (username,primary_dataset,taskdir,number)

if __name__=="__main__":
    import sys
    dataset=sys.argv[1]
    from CRAB_Datasets import datasets as crab_datasets
    from NONCRAB_Datasets import datasets as noncrab_datasets
    from Datasets_Filelist import datasets as datasets_filelist
    try:
        tasks=crab_datasets[dataset]
        primary_dataset=dataset
        if '/' in primary_dataset:
            primary_dataset=primary_dataset.split('/')[1]
        filelist=get_filelist(dataset,'T3_US_FNALLPC',[taskname_to_lfn(task,primary_dataset) for task in tasks])
    except KeyError:
        lfns=noncrab_datasets[dataset]
        filelist=get_filelist(dataset,'T3_US_FNALLPC',lfns)
    except KeyError:
        filelist=datasets_filelist[dataset]
    print filelist
