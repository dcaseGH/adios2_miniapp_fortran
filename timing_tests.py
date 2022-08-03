import reframe as rfm
import reframe.utility.sanity as sn
import os

@rfm.simple_test
class TimingBPTest(rfm.RegressionTest):
    def __init__(self):
        super().__init__()
        self.descr = ('Simple build test')
        self.valid_systems = ['*']
        self.valid_prog_environs = ['*']

        self.executable = './adios2_miniapp.exe'
        self.num_tasks = 1
        self.executable_opts = ['1', '1']

        self.build_system = 'Make'
        self.build_system.flags_from_environ = False
        self.build_system.max_concurrency = 1

        @sn.deferrable
        def extract_walltime(self):
            return sn.extractsingle(r'Wallclock time \(s\):\s+(.+)\s', self.stderr, 1, float)

        self.sanity_patterns = sn.all([
                                       sn.assert_found('Step', self.stderr)
            ])

        @sn.deferrable
        def analyse_bp_json(self, required_attribute):
            '''  Example of reading somethig from a json bp file '''
            with open(self.stagedir + '/initial_dat.bp/profiling.json') as json_file:
                json_dict = eval(json_file.read())
            return int(json_dict[0][required_attribute])

        self.perf_patterns = {'walltime':  extract_walltime(self),
                              'bytes':     analyse_bp_json(self, 'bytes')
                              }

        self.tags = {'timing', 'BP'}

@rfm.simple_test
class TimingHDF5Test(TimingBPTest):
    def __init__(self):
        super().__init__()
        self.tags = {'timing', 'HDF5'}
        self.build_system.cppflags = ['-DHDF5']

        #Why not inherited?
        @sn.deferrable
        def extract_walltime(self):
            return sn.extractsingle(r'Wallclock time \(s\):\s+(.+)\s', self.stderr, 1, float)


        self.perf_patterns = {'walltime':  extract_walltime(self)}