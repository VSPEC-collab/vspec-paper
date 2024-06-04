"""
When running the GHA it is important that the PSG API key is set.
"""

import os

import pypsg

try:
    key = os.environ['PSG_API_KEY']
    pypsg.settings.save_settings(api_key=key)

except KeyError:
    assert pypsg.docker.is_psg_installed()