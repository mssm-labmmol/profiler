#!/usr/bin/env python3
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

from modules.pre_process import pre_process_file
from modules.process import itp_to_stp
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    args_g = parser.add_argument_group('arguments')
    args_g.add_argument('-f', dest='input', metavar="<input.itp>", type=str, required=True, help='input *.itp file')
    args_g.add_argument('-o', dest='output', metavar="<output.stp>", type=str, required=False, help='output *.stp file')
    args_g.add_argument('--no-inter', dest='interface', required=False, action='store_false', help='suppress interface to write profilerTools blocks [ refdihedral ] and [ opt* ]')
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    use_interface = args.interface

    if not (input_file.endswith('.itp')):
        raise Exception("error: input file must have .itp extension")

    if output_file == None:
        output_file = os.path.basename(input_file)[:-4] + '.stp'
    else:
        if not (output_file.endswith('.stp')):
            raise Exception("error: output file must have .stp extension")

    temporary_file = input_file + ".tmp.itp"
    pre_process_file (input_file, temporary_file)
    itp_to_stp (temporary_file, output_file, use_interface)
    os.remove(temporary_file)
