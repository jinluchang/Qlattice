#!/usr/bin/env python3

import sys

if not (sys.version_info.major == 3 and sys.version_info.minor >= 6):
    print(sys.argv)
    print("You are using not supported Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

import subprocess as p
import time
import os

########################

class LogCmd:

    def __init__(self):
        log_path = __file__ + ".log"
        time_sec = time.time()
        self.log_prefix = f"{time_sec:30.10f}: "
        self.log_file = open(log_path, 'a')
        self.log("Start.")
        self.log(f"pwd='{os.getcwd()}'")

    def __del__(self):
        self.log("End.")
        self.log_file.write(f"\n");
        self.log_file.close()

    def log(self, msg):
        msg = str(msg)
        msg = msg.replace("\n", f"\n{self.log_prefix}")
        self.log_file.write(f"{self.log_prefix}{msg}\n");

log = LogCmd()

########################

def process_remove_arg(argv):
    name = "--wrapper-remove-arg="
    name_len = len(name)
    args_to_remove = []
    for arg in argv:
        if arg[:name_len] == name:
            args_to_remove.append(arg)
            args_to_remove.append(arg[name_len:])
    argv_new = []
    for arg in argv:
        if arg not in args_to_remove:
            argv_new.append(arg)
    return argv_new

def quote_comma(arg):
    return arg.replace(",", "\\,").replace("$", "\\\\$")

class NvccCmdLine:

    def __init__(self, argv):
        self.argv = argv.copy()
        self.name = None
        self.ccbin = None
        self.cc_only_flags = []
        self.is_compile = None
        self.is_link = None
        self.output = None
        self.common_flags = []
        self.nv_flags = []
        self.cc_flags = []
        self.omit_flags = []
        self.wl_group_flags = []
        self.parse_all()

    def __str__(self):
        lines = [
                "NvccCmdLine={",
                f"argv={self.argv}",
                f"name={self.name}",
                f"ccbin={self.ccbin}",
                f"cc_only_flags={self.cc_only_flags}",
                f"is_compile={self.is_compile}",
                f"is_link={self.is_link}",
                f"output={self.output}",
                f"common_flags={self.common_flags}",
                f"nv_flags={self.nv_flags}",
                f"cc_flags={self.cc_flags}",
                f"omit_flags={self.omit_flags}",
                f"wl_group_flags={self.wl_group_flags}",
                "}",
                ]
        return "\n".join(lines)

    def parse_name(self):
        """
        ``self.argv`` should take the form [ "PATH/NVCC.py", "-ccbin", "CXX.sh", ... ]
        """
        self.name = self.argv[0]
        assert self.name.endswith("NVCC.py")
        self.argv = self.argv[1:]
        ccbin = "-ccbin"
        assert self.argv[0] == "-ccbin"
        assert len(self.argv) >= 2
        self.ccbin = self.argv[1]
        self.argv = self.argv[2:]

    def parse_cc_only_flags(self):
        """
        If there is such flags, call cc instead of nvcc
        """
        opt_pool = set([
            "--version",
            "-Wl,--version",
            "--print-search-dirs",
            "-dM",
            "-E",
            "-P",
            ])
        argv_new = []
        for arg in self.argv:
            if arg in opt_pool:
                self.cc_only_flags.append(arg)
            else:
                argv_new.append(arg)
        self.argv = argv_new

    def parse_compile_link(self):
        self.is_compile = False
        self.is_link = False
        if "--NVCC-compile" in self.argv:
            self.argv = list(filter(lambda x: x != "--NVCC-compile", self.argv))
            for arg in self.argv:
                if arg.endswith(".cpp"):
                    self.is_compile = True
                    break
        if "--NVCC-link" in self.argv:
            self.is_link = True
            self.argv = list(filter(lambda x: x != "--NVCC-link", self.argv))

    def parse_output(self):
        argv_new = []
        is_output = False
        for arg in self.argv:
            if arg == "-o":
                is_output = True
            elif is_output:
                assert self.output is None
                self.output = arg
                is_output = False
            else:
                argv_new.append(arg)
        self.argv = argv_new

    def parse_common_flags(self):
        opt_pool = set([
            "-w",
            "-shared",
            ])
        argv_new = []
        for arg in self.argv:
            if arg in opt_pool:
                self.common_flags.append(arg)
            elif arg.startswith("-std="):
                self.common_flags.append(arg)
            elif arg.startswith("-I"):
                self.common_flags.append(arg)
            elif arg.startswith("-isystem"):
                self.common_flags.append("-I" + arg.removeprefix("-isystem"))
            else:
                argv_new.append(arg)
        self.argv = argv_new

    def parse_nv_flags(self):
        opt_pool = set([
            # "-w",
            "--expt-extended-lambda",
            "--expt-relaxed-constexpr",
            ])
        opt1_pool = set([
            "-Xcudafe",
            ])
        argv_new = []
        n_arg = 0
        for arg in self.argv:
            if n_arg > 0:
                self.nv_flags.append(arg)
                n_arg -= 1
            elif arg in opt_pool:
                self.nv_flags.append(arg)
            elif arg.startswith("-arch=sm_"):
                self.nv_flags.append(arg)
            elif arg in opt1_pool:
                self.nv_flags.append(arg)
                n_arg = 1
            else:
                argv_new.append(arg)
        self.argv = argv_new

    def parse_cc_flags(self):
        opt_pool = set([
            "-fopenmp",
            "-fPIC",
            "-fno-strict-aliasing",
            "-fdiagnostics-color=always",
            "-MD",
            "-fpermissive",
            "-D_FILE_OFFSET_BITS=64",
            "-Wl,--allow-shlib-undefined",
            "-Wl,--as-needed",
            "-Wl,--no-undefined",
            "-Wl,-O1",
            "-fvisibility=hidden",
            "-fvisibility-inlines-hidden",
            ])
        opt1_pool = set([
            "-MQ",
            "-MF",
            ])
        argv_new = []
        n_arg = 0
        for arg in self.argv:
            if n_arg > 0:
                self.cc_flags.append(arg)
                n_arg -= 1
            elif arg in opt_pool:
                self.cc_flags.append(arg)
            elif arg in opt1_pool:
                self.cc_flags.append(arg)
                n_arg = 1
            elif arg.startswith("-Wl,-rpath"):
                self.cc_flags.append("-Wl,-rpath,XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
                self.cc_flags.append(arg)
            else:
                argv_new.append(arg)
        self.argv = argv_new

    def parse_omit_flags(self):
        opt_pool = set([
            # "-w",
            "-Wall",
            "-Winvalid-pch",
            "-Wextra",
            "-Wpedantic",
            "-Xcompiler",
            "-MD",
            ])
        opt1_pool = set([
            "-ccbin",
            "-MQ",
            "-MF",
            ])
        argv_new = []
        n_arg = 0
        for arg in self.argv:
            if n_arg > 0:
                self.omit_flags.append(arg)
                n_arg -= 1
            elif arg in opt_pool:
                self.omit_flags.append(arg)
            elif arg in opt1_pool:
                self.omit_flags.append(arg)
                n_arg = 1
            elif arg.startswith("-Wl,-soname,"):
                self.omit_flags.append(arg)
            else:
                argv_new.append(arg)
        self.argv = argv_new

    def parse_wl_group_flags(self):
        opt_pool = set([
            "-Wall",
            "-Winvalid-pch",
            "-Wextra",
            "-Wpedantic",
            ])
        argv_new = []
        is_wl_group = False
        for arg in self.argv:
            if arg == "-Wl,--start-group":
                is_wl_group = True
            elif arg == "-Wl,--end-group":
                is_wl_group = False
            elif is_wl_group:
                self.wl_group_flags.append(arg)
            else:
                argv_new.append(arg)
        self.argv = argv_new

    def parse_all(self):
        self.parse_name()
        self.parse_cc_only_flags()
        self.parse_compile_link()
        self.parse_output()
        self.parse_omit_flags()
        self.parse_common_flags()
        self.parse_nv_flags()
        self.parse_cc_flags()
        self.parse_wl_group_flags()

    def prepare_cc_flags(self):
        argv_new = []
        for arg in self.cc_flags:
            argv_new += [ "-Xcompiler", quote_comma(arg), ]
        return argv_new

    def prepare_wl_group_flags(self):
        argv_new = []
        for arg in self.wl_group_flags:
            if arg.startswith("-l") or arg.endswith(".a"):
                argv_new.append(arg)
            elif (not arg.startswith("-")) and arg.endswith(".so"):
                dirname = os.path.dirname(arg)
                libname = os.path.basename(arg)
                assert libname.startswith("lib")
                libname = libname.removeprefix("lib").removesuffix(".so")
                argv_new.append(f'-L{dirname}')
                argv_new.append(f'-l{libname}')
            else:
                argv_new.append('-Xcompiler')
                argv_new.append(quote_comma(arg))
        return argv_new

    def make_cc_argv(self):
        if not self.cc_only_flags:
            return None
        argv_new = [ self.ccbin, ] + self.cc_only_flags + self.common_flags + self.cc_flags + self.argv
        return argv_new

    def make_argv(self):
        argv_new = []
        assert not self.cc_only_flags
        argv_new += [ "nvcc", "-ccbin", self.ccbin, ]
        if self.is_compile:
            argv_new += [ "-x", "cu", ]
        if self.is_link:
            argv_new += [ "-link", ]
        argv_new += self.nv_flags
        if self.is_compile:
            # argv_new += [ "-dc", ]
            pass
        argv_new += self.common_flags
        argv_new += self.prepare_cc_flags()
        argv_new += self.argv
        argv_new += self.prepare_wl_group_flags()
        if self.is_link:
            argv_new += [ "-lcudart", "-lcufft", ]
        if self.output is not None:
            argv_new += [ "-o", self.output, ]
        return argv_new

    def call(self, argv):
        status = p.call(argv)
        log.log(f"{' '.join([ repr(arg) for arg in argv ])}")
        log.log(f"status={status}")
        if status != 0:
            print(f"pwd={os.getcwd()}")
            print(f"sys.argv={sys.argv}")
            print(f"argv={argv}")
            print(f"status={status}")
        return status

    def ldd(self):
        out = self.output
        if out is None:
            return
        if out[-2:] == ".o":
            return
        log.log(f"ldd '{out}'")
        try:
            ldd_out = p.check_output([ "ldd", out, ], stderr=p.STDOUT).decode("utf-8")
            log.log(ldd_out)
        except Exception:
            log.log(f"ldd failed.")

    def readelf(self):
        out = self.output
        if out is None:
            return
        if out[-2:] == ".o":
            return
        log.log(f"readelf -d '{out}'")
        try:
            ldd_out = p.check_output([ "readelf", "-d", out, ], stderr=p.STDOUT).decode("utf-8")
            log.log(ldd_out)
        except Exception:
            log.log(f"readelf failed.")

    def call_and_exit(self, argv):
        status = self.call(argv)
        self.ldd()
        self.readelf()
        for k, v in os.environ.items():
            log.log(f"os.environ['{k}'] = '{v}'")
        sys.exit(status)

    def run(self):
        if self.cc_only_flags:
            argv = self.make_cc_argv()
            self.call_and_exit(argv)
        else:
            argv = self.make_argv()
            self.call_and_exit(argv)

##################################################

if __name__ == "__main__":
    argv = sys.argv.copy()
    log.log(f"{argv}")
    argv = process_remove_arg(argv)
    cmd = NvccCmdLine(argv)
    log.log(f"{cmd}")
    cmd.run()
