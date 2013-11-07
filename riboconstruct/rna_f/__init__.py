import ctypes
import os
import os.path

from riboconstruct import __file__ as riboconstruct_rel_path


# load the library
riboconstruct_path = os.path.dirname(os.path.abspath(riboconstruct_rel_path))
rna_f_path = os.path.join(riboconstruct_path, "rna_f", "libRNAf_wrap.so")
rna_f_lib = ctypes.cdll.LoadLibrary(rna_f_path)

# set return types
rna_f_lib.evaluate.restype = ctypes.c_void_p
rna_f_lib.release.argtypes = ctypes.c_void_p,

save_modi = set(('c', 'd', 'r', 'i', 's'))
hd_modi = set(('s', 'n', 'e'))


class RNAf(object):
    def __init__(self):
        # keep reference to make sure it does exist at least as long as
        # as an instance of this object exists
        self.lib = rna_f_lib
        # set standard modi
        self._hd_mode = 'n'
        self._save_mode = 's'

    def evaluate(self, path, pattern_file):
        path_to_pattern_file = os.path.join(path, pattern_file)
        if not os.path.isfile(path_to_pattern_file):
            raise IOError("Pattern '%s' does not exist." % path_to_pattern_file)
        if path[-1] != '/':
            path += '/'

        p_pattern_str = self.lib.evaluate(ctypes.c_char_p(self._hd_mode),
                                          ctypes.c_char_p(self._save_mode),
                                          ctypes.c_char_p(path),
                                          ctypes.c_char_p(pattern_file))
        # make a copy
        pattern_str = (
            ctypes.cast(p_pattern_str, ctypes.c_char_p).value.strip()[0:-1])
        self.lib.release(p_pattern_str)
        return pattern_str

    def save_mode():
        def fset(self, save_mode):
            for c in save_mode:
                if c not in save_modi:
                    return
            self._save_mode = save_mode

        def fget(self):
            return self._save_mode

        return locals()
    save_mode = property(**save_mode())

    def hd_mode():
        def fset(self, hd_mode):
            if hd_mode not in hd_modi:
                return
            self._hd_mode = hd_mode

        def fget(self):
            return self._hd_mode

        return locals()
    hd_mode = property(**hd_mode())
