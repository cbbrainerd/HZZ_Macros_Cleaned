import FWCore.ParameterSet.Config as cms
process=cms.Process("Analysis")
process.source=cms.Source("PoolSource",fileNames=cms.untracked.vstring(''))

