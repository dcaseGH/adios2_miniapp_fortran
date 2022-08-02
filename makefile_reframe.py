import reframe as rfm
import reframe.utility.sanity as sn
import os

@rfm.simple_test
class BuildBPTest(rfm.RegressionTest):
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
        self.build_system.max_concurrency = 1 #much safer to be 1
#        self.prebuild_cmds = ['pwd',
#                              'make clean']

        @sn.deferrable
        def check_BP_files(self):
            subdirs = next(os.walk(self.stagedir))[1]
            return (len(subdirs) == 1 and subdirs[0].endswith('.bp'))

        self.sanity_patterns = sn.all([
                                       sn.assert_found('Step', self.stderr),
                                       check_BP_files(self)
        ])

#        self.cleanup(remove_files=True)
        self.maintainers = ['Jigglypuff']
        self.tags = {'build', 'BP'}


@rfm.simple_test
class BuildHDF5Test(BuildBPTest):
    def __init__(self):
        super().__init__()
        self.tags = {'build', 'HDF5'}
        self.build_system.cppflags = ['-DHDF5']

#        @run_before('compile')
#        def set_preproc(self):
#            self.build_system.cppflags = ['-DHDF5']

        @sn.deferrable
        def check_HDF5_files(self):
            output_files = [file for file in os.listdir(self.stagedir)
                                     if file.endswith('.hdf5')]
            return (len(output_files) == 1)

        self.sanity_patterns = sn.all([
                                       sn.assert_found('Step', self.stderr),
                                       check_HDF5_files(self)
            ])

