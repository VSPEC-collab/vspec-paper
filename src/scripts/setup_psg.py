"""
When running the GHA it is important that the PSG API key is set.
"""

import os

import pypsg

import paths

outfile = paths.output / 'psg_status.txt'

def setup_psg():
    try:
        key = os.environ['PSG_API_KEY']
        pypsg.settings.save_settings(api_key=key)
        s = 'The PSG API key was set from an environment variable.'

    except KeyError:
        s = 'No environment variable was found to set the API key.'
    print(s)
    with open(outfile,'wt', encoding='utf-8') as f:
        f.write(s)
