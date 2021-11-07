slow_math_warning_string = """
************************************************************************

                                   WARNING!

It looks like the program could not import the ``geometrypy``
module. This might significantly affect performance!

This probably means you either skipped the ``build`` command when
installing the program, or it failed without your notice.

************************************************************************
"""

def slow_math_warning():
    print(slow_math_warning_string)

