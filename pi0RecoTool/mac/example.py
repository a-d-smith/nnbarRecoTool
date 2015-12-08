import sys
from ROOT import gSystem
gSystem.Load("libnnbarRecoTool_pi0RecoTool")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing pi0RecoTool..."

sys.exit(0)

