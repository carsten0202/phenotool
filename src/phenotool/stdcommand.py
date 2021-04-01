
##################################################
#
# --%%  CLASS: StdCommand  %%--

import click

class StdCommand(click.Command):
        def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)
                self.params.insert(0, click.Argument(['files'], nargs=-1, type=click.File()))
#               self.epilog = EPILOG.legal

