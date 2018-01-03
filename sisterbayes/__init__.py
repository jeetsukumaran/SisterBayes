#! /usr/bin/env python

##############################################################################
## Copyright (c) 2017 Jeet Sukumaran.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
## IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
## THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
## AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import subprocess
import os
import sys
import site

###############################################################################
## PACKAGE METADATA
import collections
__project__ = "SisterBayes"
__version__ = "2018.01.01"
__author__ = "Jeet Sukumaran"
__copyright__ = "Copyright 2017 Jeet Sukumaran"

# Return the git revision as a string
# Modified from: NumPy
def git_revision_sha1():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.PIPE,).strip().decode("ascii")
    except:
        return None

def git_revision_date():
    try:
        return subprocess.check_output('git log -1 --pretty=format:%ci'.split(), stderr=subprocess.PIPE,).strip().decode("ascii")
    except:
        return None

def revision_description():
    cwd = os.path.dirname(os.path.abspath(__file__))
    try:
        desc = subprocess.check_output(
                ["git", "log", "-1", "--pretty=format:%h, %ci"],
                stderr=subprocess.PIPE,
                cwd=cwd).strip().decode("ascii")[:-6]
    except:
        return ""
    if desc:
        git_date = git_revision_date()
        revision_text = " (revision: {})".format(desc)
    else:
        revision_text = ""
    return revision_text

def package_id():
    return "{} {}{}".format(__project__, __version__, revision_description())

def homedir():
    import os
    try:
        try:
            __homedir__ = __path__[0]
        except AttributeError:
            __homedir__ = os.path.dirname(os.path.abspath(__file__))
        except IndexError:
            __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except OSError:
        __homedir__ = None
    except:
        __homedir__ = None
    return __homedir__

def description(dest=None):
    if dest is None:
        dest = sys.stdout
    fields = collections.OrderedDict()
    fields["SisterBayes version"] = package_id()
    fields["SisterBayes location"] = homedir()
    fields["Python version"] = sys.version.replace("\n", "")
    fields["Python executable"] = sys.executable
    try:
        fields["Python site packages"] = site.getsitepackages()
    except:
        pass
    max_fieldname_len = max(len(fieldname) for fieldname in fields)
    for fieldname, fieldvalue in fields.items():
        dest.write("{fieldname:{fieldnamewidth}}: {fieldvalue}\n".format(
            fieldname=fieldname,
            fieldnamewidth=max_fieldname_len + 2,
            fieldvalue=fieldvalue))

if __name__ == "__main__":
    description()
