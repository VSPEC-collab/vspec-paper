"""
When running the GHA it is important that the PSG API key is set.
"""

import os

import pypsg

def setup_psg():
    try:
        key = os.environ['PSG_API_KEY']
        pypsg.settings.save_settings(api_key=key)
        print('API Key was successfully set.')

    except KeyError:
        print('API Key was not set.')
        assert pypsg.docker.is_psg_installed()
        print('PSG is installed.')
