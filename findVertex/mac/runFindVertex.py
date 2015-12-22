import sys
import ROOT

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from seltool import ertool
from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
curPath = sys.argv[1];
#for x in range(1,41):
#	curPath = curPath + "/" + `x`
#	my_proc.add_input_file(curPath+"/larlite_mcinfo.root")

#for x in xrange(len(sys.argv)-1):
#	my_proc.add_input_file(sys.argv[x+1])
for x in xrange(100):
	filename=curPath+'/larlite_mcinfo_{}.root'.format(x)
	file = ROOT.TFile(filename)
	if not file.IsZombie():
		my_proc.add_input_file(filename)

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("Default.root")

# Create ERTool algorithm (empty base class for this example)
my_algoempart  = ertool.AlgoEMPart()
my_algorecopi0 = ertool.ERAlgorecoPi0()

# Create ERTool analysis (empty base class for this example)
my_ana = ertool.ERAlgofindVertex()

# Create larlite interfce analysis unit for ERTool
my_anaunit = fmwk.ExampleERSelection()

# Set Producers
# First Argument: True = MC, False = Reco
# Second Argument: producer module label

my_anaunit.SetShowerProducer(True,"mcreco");
my_anaunit.SetTrackProducer(True,"mcreco");

# my_anaunit.SetVtxProducer(True,"generator");

# Implement manager
my_anaunit._mgr.AddAlgo(my_algoempart)
my_anaunit._mgr.AddAlgo(my_algorecopi0)
my_anaunit._mgr.AddAlgo(my_ana)
my_anaunit._mgr._mc_for_ana = True
# my_anaunit._mgr.AddAna(my_ana)
# my_ana._mode =True # True = Select. False = Fill mode
my_proc.add_process(my_anaunit)

# run!
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

#my_algo.StoreParams()
sys.exit(0)

