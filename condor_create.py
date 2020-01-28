#!/usr/bin/env python

from CRAB_Datasets import datasets as crab_datasets
from NONCRAB_Datasets import datasets as noncrab_datasets
from Datasets_Filelist import datasets as datasets_filelist

condor_cfg="""#!/bin/bash

''':'
[ -z "$CMSSW_BASE" ] || exec env -i HOME="$HOME" bash -l "$0" "$@"
exec /usr/bin/python "$0" "$@"
exit 1
'''

import htcondor, classad
import pickle, os, imp, errno, sys
schedd_ad=imp.load_source('schedd','%s/bin/condor_schedd.py' % os.environ['HOME']).get_schedd()

def getconfig(stage):
    sub=dict()
    sub['universe'] = "vanilla"
    sub['use_x509userproxy'] = "True"
    sub['Submission_Stage'] = str(stage)
    sub['InputFilename'] = "inputs_{jobname}_$(Submission_Stage).txt"
    sub['Unique_Id'] = "$(Submission_Stage).$(Cluster).$(Process)"
    sub['Executable'] = "../../{executable}"
    sub['Arguments'] = "{year} {is_mc} {sample} $(process) {dataset} $(InputFilename) {files_per_job} $(Cluster) {jobname}"
    sub['Should_Transfer_Files'] = "YES"
    sub['WhenToTransferOutput'] = "ON_EXIT"
    sub['Requirements'] = 'TARGET.OpSys == "LINUX"&& (TARGET.Arch != "DUMMY" )'
    sub['Transfer_Output_Files'] = "fallback_$(Cluster).$(Process).tar.gz"
    sub['Transfer_Input_Files'] = "../../compilereference.sh, ../../HZZ4LeptonsAnalysis_{sample}.C, ../../HZZ4LeptonsAnalysis_all.h, ../../Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root, ../../ScaleFactors_mu_Moriond2018_final.root, ../../egammaEffi_txt_EGM2D_Moriond2018v1.root, ../../egammaEffi_txt_EGM2D_Moriond2018v1_gap.root, ../../egammaEffi_txt_EGM2D_runBCDEF_passingRECO_lowEt.root, ../../egammaEffi_txt_EGM2D_runBCDEF_passingRECO.root, ../../PU_Reweight_2017.root, ../../HISTOShapes2HDM_READ_ext.root, ../../HISTOShapesZpB_READ.root, $(InputFilename), ../../compilereference_all.C, ../../ZZMatrixElement.tar.gz, ../../pu_weights_2018.root, ../../DataPileupHistogram2018_69200_100bins.root, ../../pileup_corrector.cpp, ../../pileup_corrector.h, ../../NewNtuple.h, ../../NewNtuple.C" 
#../../MELA_libs/libcollier.so, ../MELA_libs/libjhugenmela.so, ../MELA_libs/libmcfm_705.so, 
    sub['Output'] = "$(Unique_Id).stdout"
    sub['Error'] = "$(Unique_Id).stderr"
    sub['Log'] = "$(Unique_Id).log"
    sub['notify_user'] = "{logname}@{domain}"
    sub['max_retries'] = "5"
    sub['coresize'] = "0"
    return sub
    #Add number of jobs to first submit

with open('.condor_history.pkl') as f:
    condor_history=pickle.load(f)
 
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('stage',nargs='?',type=int)

options=parser.parse_args()
last_stage=condor_history['stage']

if last_stage!=0:
    #Do some status stuff, parse and figure out failed jobs
    #Unfinished
    last_schedd_ad,last_cluster=condor_history['stages'][-1]
    last_job=condor_history['clusters'][(last_schedd,last_cluster)]
    last_classad=last_job['classAd']
    last_schedd=htcondor.Schedd(last_schedd_ad)
    status=schedd.xquery(requirements='ClusterId == %d' % last_cluster)
else: #First submit, no need to do any status stuff
    schedd=htcondor.Schedd(schedd_ad)
    sub_dict=getconfig(1)
    sub=htcondor.Submit(sub_dict)
    number_of_jobs=condor_history['number_of_jobs']
    with schedd.transaction() as txn:
        cluster_number=sub.queue(txn,number_of_jobs)
    condor_history['stage']=1
    condor_history['stages'].append((schedd_ad,cluster_number))
    condor_history['clusters'][(schedd_ad,cluster_number)]={{ 'stage' : 1 , 'status' : 'submitted' , 'classAd' : sub_dict , 'numberOfJobs' : number_of_jobs }}

import shutil
with open('.condor_history.pkl.swp','w') as f:
    pickle.dump(condor_history,f)

shutil.move('.condor_history.pkl.swp','.condor_history.pkl')
"""

import argparse
parser=argparse.ArgumentParser(description='Create condor jobs for submission.')
parser.add_argument('dataset',action='store',help='Dataset to run over.')
filesPerJobDefault=1
parser.add_argument('-j','--filesperjob',action='store',help='How many jobs to run over. Default %i.' % filesPerJobDefault,default=filesPerJobDefault)
parser.add_argument('sample',action='store',help='Sample to run over. One of 4e/4mu/2e2mu')
parser.add_argument('jobname',action='store',help='Name of condor job. Defaults to job_<dataset>_<sample>_condor.',nargs='?')
parser.add_argument('-r','--read_site',action='store',help='What site to search to find ntuples.',default='T3_US_FNALLPC')

def taskname_to_lfn(taskname,primary_dataset=None):
    (number,task)=taskname.split(':',1)
    (username,taskdir)=task.split('_',1)
    if primary_dataset is None:
        (crab,analysis,primary_dataset)=taskdir.split('_')[0:3]
        assert(crab=='crab')
    return '/store/user/%s/%s/%s/%s' % (username,primary_dataset,taskdir,number)

import sys
import os
import math
import pickle

def generate_condor(number_of_files,sample,jobname,site,year,is_mc,dataset,logname=os.environ['LOGNAME'],executable='submit_condor.sh',domain='fnal.gov' if 'lpc' in os.environ['HOSTNAME'] else 'cern.ch',files_per_job=10,**kwargs):
    arguments=locals()
    number_of_files=arguments.pop('number_of_files')
    if 'test' not in arguments:
        arguments['test']=''
    if 'number_of_jobs' not in arguments:
        arguments['number_of_jobs']=int(math.ceil(number_of_files/float(files_per_job)))
    return (condor_cfg.format(**arguments),arguments)

args=parser.parse_args()
dataset=args.dataset
jobname=args.jobname
search_site=args.read_site
files_per_job=args.filesperjob
sample=args.sample

import string

if jobname is None:
    jobname='job_%s_%s_condor' % (dataset,sample)

jobname=jobname.translate(string.maketrans('/','.'))
#from get_filelist import get_filelist
#try:
#    tasks=crab_datasets[dataset]
#    primary_dataset=dataset
#    if '/' in primary_dataset:
#        primary_dataset=primary_dataset.split('/')[1]
#        print primary_dataset
#    filelist=get_filelist(dataset,search_site,[taskname_to_lfn(task,primary_dataset) for task in tasks])
#except KeyError:
#    try:
#        lfns=noncrab_datasets[dataset]
#        filelist=get_filelist(dataset,search_site,lfns)
#    except KeyError:
#        try:
#            filelist=datasets_filelist[dataset]
#        except KeyError:
#            print 'Warning: file list for dataset %s could not be found. Continuing with empty file list.' % dataset
#            filelist=['']
import json
with open('HZZ_dataset_filenames.json') as f:
    filelist=json.load(f)[dataset]['files']
filelist=['root://cms-xrd-global.cern.ch/%s' % f for f in filelist]

import errno
try:
    os.mkdir('condor_jobs2/%s' % jobname)
except OSError as e:
    if e.errno == errno.EEXIST:
        print('A job named "{}" already exists. Use a different name or delete the directory first.'.format(jobname))
    else:
        raise
os.chdir('condor_jobs2/%s' % jobname)


test_args= { 
    'site' : 'FNAL',
    'year' : '2018',
    'is_mc'   : 'data',
    'dataset' : dataset,
    'files_per_job'  : files_per_job,
    'sample' : sample,
}

with open('inputs_{}_1.txt'.format(jobname),'w') as f:
    for fn in filelist:
        f.write('%s\n' % fn)
condor_cfg_f,condor_args=generate_condor(number_of_files=len(filelist),jobname=jobname,**test_args)
with open('condor_{}_cfg.py'.format(jobname),'w') as f:
    f.write(condor_cfg_f)

number_of_jobs=condor_args['number_of_jobs']

with open('.condor_history.pkl','w') as f:
    pickle.dump(dict((
        ('number_of_jobs', number_of_jobs), 
        ('stage' , 0), 
        ('stages' , []), 
        ('clusters' , dict()))),f)

#import shutil
#shutil.copyfile('.condor_history.pkl','test/')

os.mkdir('test')
os.chdir('test')
test_args['files_per_job']=1
test_args['number_of_jobs']=1
with open('test_condor_{}.cfg'.format(jobname),'w') as f:
    f.write(generate_condor(number_of_files=1,jobname=jobname,**test_args)[0])
with open('inputs_{}.txt'.format(jobname),'w') as f:
    f.write('%s\n' % filelist[0])

os.chdir('../../..')
