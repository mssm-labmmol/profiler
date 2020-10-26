#
# This file is part of the profilerTools suite (see
# https://github.com/mssm-labmmol/profiler).
#
# Copyright (c) 2020 mssm-labmmol
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import re
import os
import sys
from shutil import which
from io import StringIO
from subprocess import check_output

def createStreamAfterPreprocessing (fn):
    cpp_path = which('cpp')
    if cpp_path is None:
        answer = 'n'
        while not (answer == 'y'):
            print("Error: The stpParser uses the C preprocessor (cpp) to translate")
            print("directives such as #include and #define. However, no path for cpp")
            print("was found. Please, preprocess the file %s manually. If you have already" % fn)
            print("done this, or if your file does not need preprocessing, answer with")
            print("'y'. If you want to quit, send the kill signal Ctrl+C.")
            answer = input("")
        fp = open(fn, 'r')
    else:
        processed_string = check_output([cpp_path, '-P', '-traditional', fn]).decode('utf-8')
        fp = StringIO(processed_string)
    return fp

def check_if_path_is_fullpath (path):
    if path.startswith('/'):
        return True
    else:
        return False

def extract_filename_from_include_considering_path (string, parent_fn):
    if not (check_if_path_is_fullpath(parent_fn)):
        raise Exception("error in function: %s must be an absolute path" % parent_fn)
    m = re.match(r'^#include\s+"(.*)"', string)
    if m:
        if (check_if_path_is_fullpath (m.group(1))):
            return m.group(1)
        else:
            return os.path.abspath(os.path.dirname(parent_fn) + '/' + m.group(1))
    else:
        return ""

def put_includes_into_stream (filename, stream):
    fp = open (filename, "r")
    for line in fp:
        # extra space for blocks
        if (line.rstrip().startswith('[')):
            stream.write("\n")
        nested_include = extract_filename_from_include_considering_path (line, filename)
        if (nested_include == ""):
            stream.write(line)
        else:
            put_includes_into_stream (nested_include, stream)
    fp.close()

def extract_define_directive (string):
    m = re.match(r'^#define\s+(\S+)\s+(.*)', string)
    if m:
        return (m.group(1), m.group(2))
    else:
        return ""

def expand_includes_in_file (input_fn, output_fn):
    input_fn = os.path.abspath(input_fn)
    output_fn = os.path.abspath(output_fn)
    fp = open (input_fn, "r")
    fp_output = open(output_fn, "w")
    for line in fp:
        # extra space for blocks
        if (line.rstrip().startswith('[')):
            fp_output.write("\n")
        nested_include = extract_filename_from_include_considering_path (line, input_fn)
        if (nested_include == ""):
            fp_output.write(line)
        else:
            put_includes_into_stream (nested_include, fp_output)
    fp_output.close()
    fp.close()

def extract_defines_in_file (input_fn):
    hash_of_defines = {}
    fp = open(input_fn, "r")
    for line in fp:
        this_define = extract_define_directive (line)
        if not (this_define == ""):
            hash_of_defines[this_define[0]] = this_define[1]
    return hash_of_defines

def substitute_defines_in_file (input_fn, output_fn, hash_of_defines):
    fp = open(input_fn, "r")
    fp_output = open(output_fn, "w")
    for line in fp:
        this_define = extract_define_directive (line)
        comment = ""
        if (this_define == ""):
            for d in hash_of_defines:
                old_line = line
                if (d in line.split()):
                    line = line.replace(d, hash_of_defines[d])
                if (old_line != line):
                    comment += " " + d
            if (comment != ""):
                fp_output.write(line.rstrip() + " ;" + comment + "\n")
            else:
                fp_output.write(line)
        else:
            fp_output.write("; COMMENTED-OUT " + line)
    fp_output.close()
    fp.close()

def pre_process_file (input_fn, output_fn, method='cpp'):
    if (method == 'cpp'):
        fpo = open(output_fn, 'w')
        fpi = createStreamAfterPreprocessing(input_fn)
        for line in fpi:
            # add padding before blocks
            if line[0] == '[':
                fpo.write('\n')
            fpo.write(line)
        fpi.close()
        fpo.close()
    elif (method == 'in-house'):
        # WARNING: This does not support #ifdef
        middle_output = output_fn + "_tmp"
        expand_includes_in_file (input_fn, middle_output)
        hash_of_defines = extract_defines_in_file (middle_output)
        substitute_defines_in_file (middle_output, output_fn, hash_of_defines)
        os.remove(middle_output)
